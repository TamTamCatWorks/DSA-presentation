"""bake_graph.py - bake the overlap graph of a whole FASTQ into ONE interactive graph-view widget.

Every read is a dot; a line joins two reads whose ends overlap by >= --show-min letters. Greedy's
picks light up in orange as you play, and edges that can no longer be used disappear. Force-directed
layout (d3-force, inlined, works offline) with zoom, pan, drag and hover, in the spirit of Obsidian's
graph view.

  python bake_graph.py ads1_week4_reads.fq --show-min 50 -o graph.html

Needs vendor/d3-force.bundle.min.js next to this script (d3-dispatch + quadtree + timer + force, ISC).
Then in the .qmd, on ONE slide:

  ```{=html}
  {{< include graph.html >}}
  ```
"""
import argparse, html, json, os, random, sys, uuid
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_greedy import read_fastq, greedy_events, overlaps

W, H = 1200, 560

def spectral_init(n, edges, rows=5, seed=7):
    """Order reads along the 1-D structure of the overlap graph (Fiedler vector), then lay that
    order out as a serpentine so the force layout starts untangled and fills the canvas."""
    A = np.zeros((n, n))
    for i, j, o in edges:
        A[i, j] = A[j, i] = max(A[i, j], o / 100.0)
    Lap = np.diag(A.sum(1)) - A
    try:
        vec = np.linalg.eigh(Lap)[1][:, 1]
    except np.linalg.LinAlgError:
        vec = np.random.default_rng(seed).standard_normal(n)
    rank = np.argsort(np.argsort(vec)) / max(1, n - 1)
    rnd = random.Random(seed); pts = []
    for t in rank:
        r = min(rows - 1, int(t * rows)); u = t * rows - r
        x = 1 - u if r % 2 else u
        pts.append([round(70 + x * (W - 140) + rnd.uniform(-2, 2), 1), round(60 + r * (H - 120) / max(1, rows - 1) + rnd.uniform(-2, 2), 1)])
    return pts

TEMPLATE = r"""
<div class="vz" id="__ID__" tabindex="0">
<script>__D3__</script>
<style>
#__ID__{font-family:system-ui,sans-serif;width:1200px;outline:none;color:#1b2a41}
#__ID__ *{box-sizing:border-box}
#__ID__ .gr-wrap{position:relative;width:1200px;height:560px;border:1px solid #e3e8ee;border-radius:10px;background:#f8fafc;overflow:hidden}
#__ID__ canvas{display:block;width:1200px;height:560px;cursor:grab;touch-action:none}
#__ID__ canvas.drag{cursor:grabbing}
#__ID__ .gr-ov{position:absolute;pointer-events:none;background:rgba(248,250,252,.88);border-radius:8px;padding:8px 12px}
#__ID__ .gr-top{left:12px;top:12px;max-width:620px}
#__ID__ .gr-stats{display:flex;gap:22px;align-items:baseline;font-size:17px;color:#52606d}
#__ID__ .gr-stats b{font-size:24px;color:#1b2a41;font-variant-numeric:tabular-nums}
#__ID__ .gr-cap{font-size:16px;line-height:1.35;margin-top:4px}
#__ID__ .gr-leg{left:12px;bottom:12px;font-size:14px;color:#3e4c59;line-height:1.7}
#__ID__ .gr-leg i{display:inline-block;width:26px;height:0;border-top:2px solid;vertical-align:middle;margin-right:8px}
#__ID__ .gr-leg u{display:inline-block;width:10px;height:10px;border-radius:50%;vertical-align:middle;margin:0 12px 0 8px}
#__ID__ .gr-tip{position:absolute;display:none;pointer-events:none;background:#1b2a41;color:#fff;font-size:14px;line-height:1.4;padding:6px 10px;border-radius:6px;white-space:nowrap}
#__ID__ .vz-bar{display:flex;gap:10px;align-items:center;margin-top:12px;font-size:16px}
#__ID__ button,#__ID__ select{font-size:18px;padding:4px 10px;border:1px solid #c9d3de;border-radius:6px;background:#fff;cursor:pointer}
#__ID__ button:focus-visible,#__ID__ select:focus-visible{outline:2px solid #e07b00;outline-offset:2px}
#__ID__ input[type=range]{flex:1}
#__ID__ .vz-pp{min-width:78px;background:#d6e9ff}
#__ID__ .vz-n{min-width:130px;font-variant-numeric:tabular-nums}
#__ID__ .sm{font-size:15px}
</style>
<div class="gr-wrap"><canvas></canvas>
 <div class="gr-ov gr-top"><div class="gr-stats"><span><b class="s-merge"></b> merges</span><span><b class="s-left"></b> pieces</span><span><b class="s-long"></b> longest</span></div><div class="gr-cap"></div></div>
 <div class="gr-ov gr-leg"><div><i style="border-color:#7896be"></i>reads overlap by at least __SHOWMIN__ letters</div><div><i style="border-color:#e07b00"></i>edge picked by greedy</div><div><u style="background:#5b9bd5"></u>single read<u style="background:hsl(150,60%,52%)"></u>reads already joined (one colour per piece)</div></div>
 <div class="gr-tip"></div></div>
<div class="vz-bar">
<button data-a="first" title="First (Home)">⏮</button><button data-a="prev" title="Back (←)">◀</button>
<button data-a="pp" class="vz-pp" title="Play / pause (Space)">▶ Play</button>
<button data-a="next" title="Step (→)">▶|</button><button data-a="last" title="Last (End)">⏭</button>
<input type="range" min="0" value="0"><span class="vz-n"></span>
<select title="Speed (merges per second)"><option value="5">5 / s</option><option value="20">20 / s</option><option value="100" selected>100 / s</option><option value="400">400 / s</option></select>
<button class="sm" data-a="fit" title="Fit graph to view (F)">Fit</button><button class="sm" data-a="shake" title="Re-run the layout (S)">Shake</button>
</div>
<script>
(function(){
const D=__DATA__, root=document.getElementById("__ID__"), $=s=>root.querySelector(s);
const N=D.n,EI=D.ei,EJ=D.ej,EO=D.eo,EV=D.ev,M=EI.length,E=EV.length,W=1200,H=560;
const cv=$("canvas"),ctx=cv.getContext("2d"),tip=$(".gr-tip"),cap=$(".gr-cap"),rng=$("input"),cnt=$(".vz-n"),pp=$(".vz-pp"),spd=$("select");
const dpr=Math.min(2,window.devicePixelRatio||1);cv.width=W*dpr;cv.height=H*dpr;
const fmt=n=>n.toLocaleString("en-US");
// ---- layout (d3-force, deterministic)
const nodes=D.init.map((p,i)=>({index:i,x:p[0],y:p[1]}));
const links=[];for(let q=0;q<M;q++)links.push({source:EI[q],target:EJ[q],o:EO[q]});
function lcg(){let s=7;return()=>(s=(s*1664525+1013904223)%4294967296)/4294967296;}
const sim=d3.forceSimulation(nodes).randomSource(lcg()).alpha(0.35).alphaDecay(0.025)
 .force("link",d3.forceLink(links).distance(l=>Math.max(3,(100-l.o)*0.45)).strength(0.35))
 .force("charge",d3.forceManyBody().strength(-7).theta(0.9)).stop();
const outL=Array.from({length:N},()=>[]),inL=Array.from({length:N},()=>[]);
for(let q=0;q<M;q++){outL[EI[q]].push([EJ[q],EO[q],q]);inL[EJ[q]].push([EI[q],EO[q],q]);}
// ---- view
let K=1,TX=0,TY=0;
function fit(){let a=1e9,b=-1e9,c=1e9,d=-1e9;for(const n of nodes){if(n.x<a)a=n.x;if(n.x>b)b=n.x;if(n.y<c)c=n.y;if(n.y>d)d=n.y;}
 K=Math.min((W-80)/Math.max(1,b-a),(H-120)/Math.max(1,d-c));TX=W/2-K*(a+b)/2;TY=H/2+14-K*(c+d)/2;}
fit();
// ---- greedy replay
let succ,pred,par,pcnt,applied=0,rootOf=new Int32Array(N),state=new Uint8Array(M);
function reset(){succ=new Int32Array(N).fill(-1);pred=new Int32Array(N).fill(-1);par=Int32Array.from({length:N},(_,i)=>i);pcnt=new Int32Array(N).fill(1);applied=0;}
function find(x){while(par[x]!==x){par[x]=par[par[x]];x=par[x];}return x;}
function goto(k){if(k<applied)reset();while(applied<k){const e=EV[applied];succ[e[0]]=e[1];pred[e[1]]=e[0];par[e[5]]=e[4];pcnt[e[4]]+=pcnt[e[5]];applied++;}
 for(let i=0;i<N;i++)rootOf[i]=find(i);
 for(let q=0;q<M;q++){const a=EI[q],b=EJ[q];state[q]=succ[a]===b?2:(succ[a]<0&&pred[b]<0&&rootOf[a]!==rootOf[b])?1:0;}}
reset();
// ---- drawing
let hover=-1,dirty=true,raf=0,drag=null,pan=null,playing=false;
const sx=n=>n.x*K+TX,sy=n=>n.y*K+TY;
const nodeR=()=>Math.max(2,Math.min(7,1.6+0.55*Math.sqrt(K)));
function draw(){
 ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,W,H);
 const dim=hover>=0;ctx.lineWidth=0.9;
 for(let c=0;c<2;c++){ctx.beginPath();for(let q=0;q<M;q++){if(state[q]!==1||((EO[q]>=75)!==(c===1)))continue;const a=nodes[EI[q]],b=nodes[EJ[q]];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}
  ctx.strokeStyle=c?"rgba(80,118,172,"+(dim?0.14:0.38)+")":"rgba(120,150,190,"+(dim?0.07:0.2)+")";ctx.stroke();}
 ctx.beginPath();for(let k=0;k<applied;k++){const e=EV[k],a=nodes[e[0]],b=nodes[e[1]];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}
 ctx.lineWidth=2.2;ctx.strokeStyle="rgba(224,123,0,"+(dim?0.4:0.95)+")";ctx.stroke();
 const nb=new Set();if(hover>=0){nb.add(hover);for(const[j]of outL[hover])nb.add(j);for(const[j]of inL[hover])nb.add(j);}
 const r=nodeR();
 for(let i=0;i<N;i++){const n=nodes[i],rt=rootOf[i];ctx.globalAlpha=dim&&!nb.has(i)?0.3:1;
  ctx.fillStyle=pcnt[rt]>1?"hsl("+(((rt*137.5)%360)|0)+",62%,50%)":"#5b9bd5";ctx.beginPath();ctx.arc(sx(n),sy(n),r,0,6.2832);ctx.fill();}
 ctx.globalAlpha=1;
 if(applied){const e=EV[applied-1];for(const i of [e[0],e[1]]){const n=nodes[i];ctx.beginPath();ctx.arc(sx(n),sy(n),r+4,0,6.2832);ctx.lineWidth=2.5;ctx.strokeStyle="#e07b00";ctx.stroke();}}
 if(hover>=0){const h=nodes[hover];ctx.lineWidth=2;ctx.strokeStyle="#1b2a41";ctx.fillStyle="#1b2a41";
  const inc=[...outL[hover],...inL[hover]].filter(t=>state[t[2]]!==0||true);
  ctx.beginPath();for(const t of inc){if(state[t[2]]===0)continue;const o=nodes[t[0]];ctx.moveTo(sx(h),sy(h));ctx.lineTo(sx(o),sy(o));}ctx.lineWidth=1.8;ctx.strokeStyle="rgba(27,42,65,.8)";ctx.stroke();
  if(K>1.2&&inc.length<=16){ctx.font="700 12px system-ui,sans-serif";ctx.textAlign="center";for(const t of inc){if(state[t[2]]===0)continue;const o=nodes[t[0]],mx=(sx(h)+sx(o))/2,my=(sy(h)+sy(o))/2;ctx.lineWidth=3;ctx.strokeStyle="#f8fafc";ctx.strokeText(t[1],mx,my-3);ctx.fillText(t[1],mx,my-3);}}
  ctx.beginPath();ctx.arc(sx(h),sy(h),r+2.5,0,6.2832);ctx.lineWidth=2.5;ctx.strokeStyle="#1b2a41";ctx.stroke();}
}
function loop(){raf=0;let act=false;if(sim.alpha()>sim.alphaMin()){sim.tick();sim.tick();act=true;dirty=true;}
 if(dirty){draw();dirty=false;}if(act||playing||drag)raf=requestAnimationFrame(loop);}
function kick(){dirty=true;if(!raf)raf=requestAnimationFrame(loop);}
// ---- interaction
const pt=e=>{const r=cv.getBoundingClientRect();return[(e.clientX-r.left)*W/r.width,(e.clientY-r.top)*H/r.height];};
function pick(mx,my){let b=-1,bd=(nodeR()+5)**2;for(let i=0;i<N;i++){const n=nodes[i],d=(sx(n)-mx)**2+(sy(n)-my)**2;if(d<bd){bd=d;b=i;}}return b;}
function setTip(i,mx,my){if(i<0){tip.style.display="none";return;}
 const o=outL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]),n=inL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]);
 let s="read #"+i;if(o.length)s+="<br>→ read #"+o[0][0]+" ("+o[0][1]+" letters)";if(n.length)s+="<br>← read #"+n[0][0]+" ("+n[0][1]+" letters)";s+="<br>"+(o.length+n.length)+" usable overlaps";
 tip.innerHTML=s;tip.style.display="block";tip.style.left=Math.min(W-170,mx+14)+"px";tip.style.top=Math.min(H-90,my+14)+"px";}
cv.addEventListener("wheel",e=>{e.preventDefault();const[mx,my]=pt(e),nk=Math.max(0.4,Math.min(60,K*Math.exp(-e.deltaY*0.0015)));TX=mx-(mx-TX)*nk/K;TY=my-(my-TY)*nk/K;K=nk;kick();},{passive:false});
cv.addEventListener("mousedown",e=>{const[mx,my]=pt(e),i=pick(mx,my);if(i>=0){drag=i;nodes[i].fx=nodes[i].x;nodes[i].fy=nodes[i].y;sim.alphaTarget(0.25);if(sim.alpha()<0.25)sim.alpha(0.25);}else{pan=[mx,my,TX,TY];cv.classList.add("drag");}kick();});
window.addEventListener("mousemove",e=>{const[mx,my]=pt(e);
 if(drag!==null){const n=nodes[drag];n.fx=(mx-TX)/K;n.fy=(my-TY)/K;kick();}
 else if(pan){TX=pan[2]+mx-pan[0];TY=pan[3]+my-pan[1];kick();}
 else if(e.target===cv){const i=pick(mx,my);if(i!==hover){hover=i;setTip(i,mx,my);kick();}else if(i>=0)setTip(i,mx,my);}});
window.addEventListener("mouseup",()=>{if(drag!==null){nodes[drag].fx=null;nodes[drag].fy=null;drag=null;sim.alphaTarget(0);}pan=null;cv.classList.remove("drag");kick();});
cv.addEventListener("mouseleave",()=>{if(hover>=0){hover=-1;tip.style.display="none";kick();}});
// ---- stepping
let i0=0,timer=null,acc=0;
function show(k){i0=Math.max(0,Math.min(E,Math.round(k)));goto(i0);
 $(".s-merge").textContent=fmt(i0)+" / "+fmt(E);$(".s-left").textContent=fmt(N-i0);$(".s-long").textContent=fmt(D.maxlen[i0])+" letters";
 cap.textContent=i0?"Merge "+fmt(i0)+": read #"+EV[i0-1][0]+" → read #"+EV[i0-1][1]+", "+EV[i0-1][2]+" letters in common. Edges that can no longer be used disappear.":"Each dot is a read. A line joins two reads whose ends overlap. Greedy takes the longest overlaps first.";
 rng.value=i0;cnt.textContent="step "+fmt(i0)+" / "+fmt(E);kick();}
function stop(){clearInterval(timer);timer=null;playing=false;pp.textContent="▶ Play";}
function play(){if(i0>=E)show(0);pp.textContent="⏸ Pause";acc=0;playing=true;timer=setInterval(()=>{acc+=spd.value/20;const st=Math.floor(acc);acc-=st;if(i0>=E){stop();return;}show(i0+Math.max(st,0));if(i0>=E)stop();},50);}
const toggle=()=>timer?stop():play();
function shake(){sim.alpha(0.6);kick();}
function refit(){fit();kick();}
rng.max=E;rng.oninput=()=>{stop();show(+rng.value);};
spd.onchange=()=>{if(timer){stop();play();}};
root.querySelectorAll("button").forEach(b=>b.onclick=()=>{const a=b.dataset.a;if(a==="pp")return toggle();if(a==="fit")return refit();if(a==="shake")return shake();stop();show(a==="first"?0:a==="last"?E:i0+(a==="next"?1:-1));root.focus();});
root.addEventListener("keydown",e=>{const m={" ":toggle,ArrowRight:()=>{stop();show(i0+1)},ArrowLeft:()=>{stop();show(i0-1)},PageDown:()=>{stop();show(i0+50)},PageUp:()=>{stop();show(i0-50)},Home:()=>{stop();show(0)},End:()=>{stop();show(E)},f:refit,F:refit,s:shake,S:shake}[e.key];
 if(m){m();e.preventDefault();e.stopPropagation();}},true);
if(window.Reveal&&Reveal.on)Reveal.on("slidechanged",stop);
show(0);
})();
</script></div>
"""

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq")
    ap.add_argument("--min-overlap", type=int, default=20, help="minimum overlap greedy may use")
    ap.add_argument("--show-min", type=int, default=50, help="draw only edges with at least this overlap")
    ap.add_argument("--rows", type=int, default=5, help="serpentine rows used to seed the layout")
    ap.add_argument("--standalone", action="store_true")
    ap.add_argument("-o", "--out", default="graph.html")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    d3 = open(os.path.join(here, "vendor", "d3-force.bundle.min.js"), encoding="utf-8").read().replace("</script", "<\\/script")
    reads = list(dict.fromkeys(read_fastq(a.fastq)))
    ev, contigs, _ = greedy_events(reads, a.min_overlap)
    edges = [(i, j, o) for o, i, j in overlaps(reads, min(a.show_min, a.min_overlap) if False else a.show_min)]
    have = {(i, j) for i, j, _ in edges}
    edges += [(e["a"], e["b"], e["ov"]) for e in ev if (e["a"], e["b"]) not in have]   # picks are always drawn
    cur = max(len(r) for r in reads); maxlen = [cur]
    for e in ev:
        cur = max(cur, e["len"]); maxlen.append(cur)
    data = {"n": len(reads), "ei": [e[0] for e in edges], "ej": [e[1] for e in edges], "eo": [e[2] for e in edges],
            "ev": [[e["a"], e["b"], e["ov"], e["len"], e["root"], e["rb"]] for e in ev], "maxlen": maxlen,
            "init": spectral_init(len(reads), edges, a.rows)}
    out = (TEMPLATE.replace("__ID__", "vz" + uuid.uuid4().hex[:8]).replace("__SHOWMIN__", str(a.show_min))
           .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__D3__", d3))
    if a.standalone:
        out = f"<!doctype html><meta charset=utf-8><title>{html.escape(a.fastq)}</title><body style='margin:20px'>{out}"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"{len(reads)} reads, {len(edges)} edges shown, {len(ev)} merges -> {a.out} ({len(out)//1024} KB)")

if __name__ == "__main__":
    main()
