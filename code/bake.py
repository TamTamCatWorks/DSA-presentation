#!/usr/bin/env python3
#!/usr/bin/env python3
"""bake.py - trace a Python function and bake it into ONE play/pause widget.

Examples (run from the folder holding the algorithm files):

  # greedy superstring from the teammate's trace, same input
  python bake.py ssp:greedy --args '[["AAA","AAB","ABB","BBB","BBA"]]' \
      --view strings --var arr --live num --hl idx_i idx_j --cmp i j -o widget.html

  # brute force on the same input (permutation view)
  python bake.py ssp:brute --args '[["AAA","AAB","ABB"]]' --view strings --var perm -o brute.html

  # any other algorithm, e.g. bubble sort as bars
  python bake.py sorts.py:bubble --args '[[5,2,4,1,3]]' --view bars --var a --cmp j -o bubble.html

Then in the .qmd, on ONE slide:

  ```{=html}
  {{< include widget.html >}}
  ```

Views: strings (suffix->prefix overlap graph), bars (numeric array), none (code + variables only).
"""
import argparse, copy, html, os, os, importlib.machinery, importlib.util, inspect, json, sys, uuid


def load(spec):
    path, fn = spec.rsplit(":", 1)
    if not os.path.exists(path) and os.path.exists(path + ".py"):
        path += ".py"  # allow `ssp:greedy` for ssp.py
    loader = importlib.machinery.SourceFileLoader("baked_mod", path)
    mod = importlib.util.module_from_spec(importlib.util.spec_from_loader("baked_mod", loader))
    loader.exec_module(mod)
    return getattr(mod, fn)


def snap(loc):
    out = {}
    for k, v in loc.items():
        if callable(v) or inspect.ismodule(v):
            continue
        out[k] = json.loads(json.dumps(v, default=repr))
    return out


def trace(fn, args):
    steps, code = [], fn.__code__

    def tracer(frame, event, arg):
        if frame.f_code is not code:
            return None
        steps.append({"l": frame.f_lineno, "v": snap(frame.f_locals)})  # 'call'

        def local(frame, event, arg):
            if event in ("line", "return"):
                v = snap(frame.f_locals)
                if event == "return":
                    v["return"] = json.loads(json.dumps(arg, default=repr))
                steps.append({"l": frame.f_lineno, "v": v})
            return local
        return local

    sys.settrace(tracer)
    try:
        fn(*copy.deepcopy(args))
    finally:
        sys.settrace(None)
    return steps


def source_lines(fn):
    lines, start = inspect.getsourcelines(fn)
    out = []
    for i, t in enumerate(lines):
        s = t.rstrip("\n")
        if not s.strip() or s.strip().startswith("#"):
            continue  # same compaction as the slide version: no blanks / comment-only lines
        out.append({"n": start + i, "t": s})
    return out


TEMPLATE = r"""
<div class="vz" id="__ID__" tabindex="0">
<style>
#__ID__{font-family:system-ui,sans-serif;width:1200px;outline:none;color:#1b2a41}
#__ID__ .vz-main{display:flex;gap:24px;align-items:flex-start}
#__ID__ .vz-code{flex:0 0 540px;height:400px;overflow:auto;background:#f8f8f8;border:1px solid #e3e8ee;border-radius:8px;padding:8px 0;font:15px/1.55 Menlo,Consolas,monospace}
#__ID__ .vz-ln{white-space:pre;padding-right:10px}
#__ID__ .vz-ln.cur{background:#fff3b0}
#__ID__ .vz-no{display:inline-block;width:40px;text-align:right;color:#9aa5b1;margin-right:6px}
#__ID__ .vz-ar{display:inline-block;width:18px;color:#e0a800}
#__ID__ .k{color:#008000;font-weight:bold}#__ID__ .s{color:#BA2121}#__ID__ .c{color:#3D7B7B;font-style:italic}#__ID__ .m{color:#666}
#__ID__ .vz-right{flex:1;min-width:0}
#__ID__ svg{width:100%;height:300px;display:block}
#__ID__ .vz-vars{font:14px/1.4 Menlo,Consolas,monospace;border-top:1px solid #e3e8ee;padding-top:6px;height:96px;overflow:auto}
#__ID__ .vz-vars b{display:inline-block;width:15ch;color:#7a4fb5}
#__ID__ .vz-vars .chg{background:#e8f5e9}
#__ID__ .vz-bar{display:flex;gap:10px;align-items:center;margin-top:10px;font-size:16px}
#__ID__ button,#__ID__ select{font-size:18px;padding:4px 10px;border:1px solid #c9d3de;border-radius:6px;background:#fff;cursor:pointer}
#__ID__ input[type=range]{flex:1}
#__ID__ .vz-pp{min-width:64px;background:#d6e9ff}
</style>
<div class="vz-main"><div class="vz-code"></div>
<div class="vz-right"><svg viewBox="0 0 640 300"></svg><div class="vz-vars"></div></div></div>
<div class="vz-bar">
<button data-a="first" title="First (Home)">⏮</button><button data-a="prev" title="Back (←)">◀</button>
<button data-a="pp" class="vz-pp" title="Play / pause (Space)">▶ Play</button>
<button data-a="next" title="Step (→)">▶|</button><button data-a="last" title="Last (End)">⏭</button>
<input type="range" min="0" value="0"><span class="vz-n"></span>
<select title="Speed"><option value="0.5">0.5×</option><option value="1" selected>1×</option><option value="2">2×</option><option value="4">4×</option><option value="10">10×</option></select>
</div>
<script>
(function(){
const D=__DATA__, root=document.getElementById("__ID__"), $=s=>root.querySelector(s);
const code=$(".vz-code"), svg=$("svg"), vars=$(".vz-vars"), rng=$("input"), cnt=$(".vz-n"), pp=$(".vz-pp"), spd=$("select");
const KW=/^(def|return|if|else|elif|for|while|in|not|and|or|is|None|True|False|import|from|break|continue|lambda)$/;
const esc=s=>s.replace(/&/g,"&amp;").replace(/</g,"&lt;");
function hl(t){const ci=t.indexOf("#");let c="";if(ci>=0){c='<span class="c">'+esc(t.slice(ci))+"</span>";t=t.slice(0,ci);}
 return t.split(/(\"[^\"]*\"|'[^']*'|\b\w+\b)/).map(p=>/^["']/.test(p)?'<span class="s">'+esc(p)+"</span>":KW.test(p)?'<span class="k">'+p+"</span>":/^\d+$/.test(p)?'<span class="m">'+p+"</span>":esc(p)).join("")+c;}
code.innerHTML=D.lines.map(l=>`<div class="vz-ln" data-n="${l.n}"><span class="vz-no">${l.n}</span><span class="vz-ar"></span>${hl(l.t)}</div>`).join("");
const rows={};code.querySelectorAll(".vz-ln").forEach(e=>rows[e.dataset.n]=e);
const ov=(a,b)=>{for(let k=Math.min(a.length,b.length);k>0;k--)if(a.slice(-k)===b.slice(0,k))return k;return 0;};
const num=x=>typeof x==="number"?x:null;
function drawStrings(v){
 const V=D.view;let arr=v[V.var];if(!Array.isArray(arr)){svg.innerHTML="";return;}
 const live=V.live&&num(v[V.live])!=null?v[V.live]:arr.length;
 const n=live,cx=320,cy=150,R=Math.min(115,40+n*14),pos=[];
 for(let i=0;i<n;i++){const a=-Math.PI/2+2*Math.PI*i/Math.max(n,1);pos.push([cx+R*1.35*Math.cos(a),cy+R*Math.sin(a)]);}
 const best=(V.hl||[]).map(k=>num(v[k])).filter(x=>x!=null), cmp=(V.cmp||[]).map(k=>num(v[k])).filter(x=>x!=null);
 const w=i=>12*String(arr[i]).length+22;let s='<defs><marker id="'+root.id+'a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10z" fill="#9db7d1"/></marker></defs>';
 if(n<=14)for(let i=0;i<n;i++)for(let j=0;j<n;j++){if(i===j)continue;const o=ov(String(arr[i]),String(arr[j]));if(!o)continue;
  const[x1,y1]=pos[i],[x2,y2]=pos[j],dx=x2-x1,dy=y2-y1,L=Math.hypot(dx,dy)||1,nx=-dy/L*9,ny=dx/L*9;
  const on=best.includes(i)&&best.includes(j)&&best.indexOf(i)<best.indexOf(j);
  const sx=x1+dx*.2+nx,sy=y1+dy*.2+ny,ex=x1+dx*.8+nx,ey=y1+dy*.8+ny;
  s+=`<line x1="${sx}" y1="${sy}" x2="${ex}" y2="${ey}" stroke="${on?"#e07b00":"#9db7d1"}" stroke-width="${on?3.5:1.6}" marker-end="url(#${root.id}a)"/><text x="${(sx+ex)/2+nx*.6}" y="${(sy+ey)/2+ny*.6+4}" font-size="12" text-anchor="middle" fill="${on?"#e07b00":"#7b8794"}" font-weight="700">${o}</text>`;}
 for(let i=0;i<n;i++){const[x,y]=pos[i],W=w(i),isb=best.includes(i),isc=cmp.includes(i);
  s+=`<rect x="${x-W/2}" y="${y-16}" width="${W}" height="32" rx="9" fill="${isb?"#ffe3b8":isc?"#d9f2d9":"#d6e9ff"}" stroke="${isb?"#e07b00":isc?"#4caf50":"#5b9bd5"}" stroke-width="1.5"/><text x="${x}" y="${y+6}" text-anchor="middle" font-size="16" font-family="Menlo,Consolas,monospace">${esc(String(arr[i]))}</text><text x="${x-W/2+3}" y="${y-19}" font-size="11" fill="#9aa5b1">[${i}]</text>`;}
 svg.innerHTML=s;}
function drawBars(v){
 const V=D.view;let a=v[V.var];if(!Array.isArray(a)){svg.innerHTML="";return;}
 const mx=Math.max(1,...a.map(Math.abs)),W=600/a.length,cmp=(V.cmp||[]).concat(V.hl||[]).map(k=>num(v[k]));let s="";
 a.forEach((x,i)=>{const h=Math.abs(x)/mx*240,on=cmp.includes(i);s+=`<rect x="${20+i*W+2}" y="${270-h}" width="${W-4}" height="${h}" fill="${on?"#ffb35c":"#8ab8e6"}"/><text x="${20+i*W+W/2}" y="288" text-anchor="middle" font-size="12">${x}</text>`;});
 svg.innerHTML=s;}
let i=0,timer=null,prev={};
function show(k){
 i=Math.max(0,Math.min(D.steps.length-1,k));const st=D.steps[i];
 if(rows.cur)rows.cur.classList.remove("cur"),rows.cur.querySelector(".vz-ar").textContent="";
 const r=rows[st.l];if(r){r.classList.add("cur");r.querySelector(".vz-ar").textContent="▶";rows.cur=r;
  const t=r.offsetTop-code.clientHeight/2;code.scrollTop=Math.max(0,t);}
 const pv=i?D.steps[i-1].v:{};
 vars.innerHTML=Object.entries(st.v).filter(([k])=>!(D.hide||[]).includes(k)).map(([k,x])=>{const t=JSON.stringify(x),ch=JSON.stringify(pv[k])!==t;return `<div class="${ch&&i?"chg":""}"><b>${esc(k)}</b> = ${esc(t)}</div>`;}).join("");
 if(D.view.kind==="strings")drawStrings(st.v);else if(D.view.kind==="bars")drawBars(st.v);else svg.innerHTML="";
 rng.value=i;cnt.textContent=(i+1)+" / "+D.steps.length;}
function stop(){clearInterval(timer);timer=null;pp.textContent="▶ Play";}
function play(){if(i>=D.steps.length-1)show(0);pp.textContent="⏸ Pause";timer=setInterval(()=>{if(i>=D.steps.length-1)stop();else show(i+1);},700/spd.value);}
const toggle=()=>timer?stop():play();
rng.max=D.steps.length-1;rng.oninput=()=>{stop();show(+rng.value);};
spd.onchange=()=>{if(timer){stop();play();}};
root.querySelectorAll("button").forEach(b=>b.onclick=()=>{const a=b.dataset.a;if(a==="pp")return toggle();stop();show(a==="first"?0:a==="last"?1e9:i+(a==="next"?1:-1));root.focus();});
root.addEventListener("keydown",e=>{const m={" ":toggle,ArrowRight:()=>{stop();show(i+1)},ArrowLeft:()=>{stop();show(i-1)},Home:()=>{stop();show(0)},End:()=>{stop();show(1e9)}}[e.key];
 if(m){m();e.preventDefault();e.stopPropagation();}},true);
root.addEventListener("click",()=>root.focus());
if(window.Reveal&&Reveal.on)Reveal.on("slidechanged",stop);
show(0);
})();
</script></div>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", help="file.py:function  (file need not end in .py, e.g. ssp:greedy)")
    ap.add_argument("--args", default="[]", help="JSON list of positional args")
    ap.add_argument("--view", choices=["strings", "bars", "none"], default="none")
    ap.add_argument("--var", help="variable holding the list to draw")
    ap.add_argument("--live", help="variable holding how many leading items of --var are alive")
    ap.add_argument("--hl", nargs="*", default=[], help="index variables to highlight (orange; strings view draws the chosen i->j edge)")
    ap.add_argument("--cmp", nargs="*", default=[], help="loop-index variables being compared (green)")
    ap.add_argument("--hide", nargs="*", default=[], help="variables to hide from the panel")
    ap.add_argument("--standalone", action="store_true", help="wrap in a full HTML page for quick testing")
    ap.add_argument("-o", "--out", default="widget.html")
    a = ap.parse_args()

    fn = load(a.spec)
    steps = trace(fn, json.loads(a.args))
    data = {"lines": source_lines(fn), "steps": steps, "hide": a.hide,
            "view": {"kind": a.view, "var": a.var, "live": a.live, "hl": a.hl, "cmp": a.cmp}}
    out = TEMPLATE.replace("__ID__", "vz" + uuid.uuid4().hex[:8]).replace("__DATA__", json.dumps(data))
    if a.standalone:
        out = f"<!doctype html><meta charset=utf-8><title>{html.escape(a.spec)}</title><body style='margin:20px'>{out}"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"{a.spec}: {len(steps)} steps -> {a.out}")


if __name__ == "__main__":
    main()