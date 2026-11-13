"""bake_graph.py - bake the overlap graph of a whole FASTQ into ONE interactive graph-view widget.

Reads are rectangles; a line joins two reads whose ends overlap by >= --show-min letters.
Greedy's picks turn into stiff rope links (reads of one contig stay as staggered rectangles,
like a read pile-up), edges that can no longer be used disappear, and the swarm starts as a
force-directed cloud (d3-force, inlined) that straightens into chains as merges accumulate.
A lens at the bottom right shows the nucleotides of the two reads being merged.

  python bake_graph.py reads.fq --show-min 50 -o graph.html

Needs vendor/d3-force.bundle.min.js next to this script.
Then in the .qmd, on ONE slide:

  ```{=html}
  {{< include graph.html >}}
  ```
"""
import argparse, html, json, os, sys, uuid
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_greedy import read_fastq, greedy_events, overlaps

TEMPLATE = r"""
<div class="vz" id="__ID__" tabindex="0">
<script>__D3__</script>
<style>
#__ID__{font-family:system-ui,sans-serif;width:1200px;outline:none;color:#1b2a41}
#__ID__ *{box-sizing:border-box}
#__ID__ .gr-wrap{position:relative;width:1200px;height:560px;border:1px solid #e3e8ee;border-radius:10px;background:#f8fafc;overflow:hidden}
#__ID__ canvas.gr-main{display:block;width:1200px;height:560px;cursor:grab;touch-action:none}
#__ID__ canvas.gr-main.drag{cursor:grabbing}
#__ID__ .gr-ov{position:absolute;pointer-events:none;background:rgba(248,250,252,.88);border-radius:8px;padding:8px 12px}
#__ID__ .gr-top{left:12px;top:12px}
#__ID__ .gr-stats{display:flex;gap:26px}
#__ID__ .gr-stats div{display:flex;flex-direction:column}
#__ID__ .gr-stats b{font-size:26px;line-height:1.1;font-variant-numeric:tabular-nums}
#__ID__ .gr-stats span{font-size:12px;color:#7b8794;letter-spacing:.04em}
#__ID__ .gr-sp{right:12px;top:12px;padding:6px 8px}
#__ID__ .gr-sp canvas{display:block;width:260px;height:64px}
#__ID__ .gr-lc{right:12px;bottom:12px;padding:0;opacity:0;transition:opacity .15s;background:rgba(248,250,252,.96);box-shadow:0 0 0 1px #e3e8ee}
#__ID__ .gr-lc canvas{display:block;width:560px;height:106px}
#__ID__ .gr-tip{position:absolute;display:none;pointer-events:none;background:#1b2a41;color:#fff;font-size:14px;line-height:1.4;padding:6px 10px;border-radius:6px;white-space:nowrap;font-variant-numeric:tabular-nums}
#__ID__ .vz-bar{display:flex;gap:10px;align-items:center;margin-top:12px;font-size:16px}
#__ID__ button,#__ID__ select{font-size:18px;padding:4px 10px;border:1px solid #c9d3de;border-radius:6px;background:#fff;cursor:pointer}
#__ID__ button:focus-visible,#__ID__ select:focus-visible{outline:2px solid #e07b00;outline-offset:2px}
#__ID__ input[type=range]{flex:1}
#__ID__ .vz-pp{min-width:78px;background:#d6e9ff}
#__ID__ .vz-n{min-width:130px;font-variant-numeric:tabular-nums}
#__ID__ .sm{font-size:15px}
#__ID__ .sm.on{background:#ffe3b8;border-color:#e07b00}
</style>
<div class="gr-wrap"><canvas class="gr-main"></canvas>
 <div class="gr-ov gr-top"><div class="gr-stats">
  <div><b class="s-merge"></b><span>merges</span></div><div><b class="s-left" style="color:#5b9bd5"></b><span>pieces</span></div>
  <div><b class="s-long" style="color:#e07b00"></b><span>longest nt</span></div><div><b class="s-ov" style="color:#7a4fb5"></b><span>overlap nt</span></div></div></div>
 <div class="gr-ov gr-sp"><canvas></canvas></div>
 <div class="gr-ov gr-lc"><canvas></canvas></div>
 <div class="gr-tip"></div></div>
<div class="vz-bar">
<button data-a="first" title="First (Home)">⏮</button><button data-a="prev" title="Back (←)">◀</button>
<button data-a="pp" class="vz-pp" title="Play / pause (Space)">▶ Play</button>
<button data-a="next" title="Step (→)">▶|</button><button data-a="last" title="Last (End)">⏭</button>
<input type="range" min="0" value="0"><span class="vz-n"></span>
<select title="Merges per second"><option value="5">5 / s</option><option value="20">20 / s</option><option value="100" selected>100 / s</option><option value="400">400 / s</option></select>
<button class="sm" data-a="fit" title="Follow / fit (F)">Fit</button><button class="sm" data-a="unfold" title="Unfold (U)">Unfold</button><button class="sm" data-a="shake" title="Shake (S)">Shake</button>
</div>
<script>
(function(){
const D=__DATA__, root=document.getElementById("__ID__"), $=s=>root.querySelector(s);
const N=D.n,EI=D.ei,EJ=D.ej,EO=D.eo,EV=D.ev,M=EI.length,E=EV.length,W=1200,H=560;
const reads=D.reads,LEN=D.lengths,FO=D.fo,U=0.25,RH=5,LANE=6,MAXL=D.maxlen;
const cv=$("canvas.gr-main"),ctx=cv.getContext("2d"),tip=$(".gr-tip"),rng=$("input"),cnt=$(".vz-n"),pp=$(".vz-pp"),spd=$("select");
const lensBox=$(".gr-lc"),lens=lensBox.querySelector("canvas"),lc=lens.getContext("2d"),spk=$(".gr-sp canvas"),sc=spk.getContext("2d");
const dpr=Math.min(2,window.devicePixelRatio||1);
function sizeCv(c,w,h){c.width=w*dpr;c.height=h*dpr;}
sizeCv(cv,W,H);const LW=560,LH=106,SW=260,SH=64;sizeCv(lens,LW,LH);sizeCv(spk,SW,SH);
const fmt=n=>n.toLocaleString("en-US");
function lcg(s0){let s=s0;return()=>(s=(s*1664525+1013904223)%4294967296)/4294967296;}
// ---- state replayed from greedy's merge log
const outL=Array.from({length:N},()=>[]),inL=Array.from({length:N},()=>[]);
for(let q=0;q<M;q++){outL[EI[q]].push([EJ[q],EO[q],q]);inL[EJ[q]].push([EI[q],EO[q],q]);}
let succ,pred,par,pcnt,plen,applied=0,rootOf=new Int32Array(N),state=new Uint8Array(M),usable=new Int32Array(M),usableN=0;
function reset(){succ=new Int32Array(N).fill(-1);pred=new Int32Array(N).fill(-1);par=Int32Array.from({length:N},(_,i)=>i);pcnt=new Int32Array(N).fill(1);plen=Int32Array.from(LEN);applied=0;}
function find(x){while(par[x]!==x){par[x]=par[par[x]];x=par[x];}return x;}
function goto(k){if(k<applied)reset();while(applied<k){const e=EV[applied];succ[e[0]]=e[1];pred[e[1]]=e[0];par[e[5]]=e[4];pcnt[e[4]]+=pcnt[e[5]];plen[e[4]]=e[3];applied++;}
 for(let i=0;i<N;i++)rootOf[i]=find(i);usableN=0;
 for(let q=0;q<M;q++){const a=EI[q],b=EJ[q];state[q]=succ[a]===b?2:(succ[a]<0&&pred[b]<0&&rootOf[a]!==rootOf[b])?1:0;if(state[q]===1)usable[usableN++]=q;}}
reset();
// ---- unfold target: serpentine, arc length = position in the final contig
const LR=1400,RR=70,CYC=LR+Math.PI*RR;
function pathPt(s){const k=Math.floor(s/CYC),t=s-k*CYC,dir=k%2===0?1:-1,y0=k*2*RR;
 if(t<=LR)return[dir>0?t:LR-t,y0];const a=(t-LR)/RR;return[dir>0?LR+RR*Math.sin(a):-RR*Math.sin(a),y0+RR-RR*Math.cos(a)];}
const TGT=[];{let tot=0;for(let i=0;i<N;i++)tot=Math.max(tot,FO[i]+LEN[i]);const rows=Math.ceil(tot*U/CYC);
 for(let i=0;i<N;i++){const p=pathPt((FO[i]+LEN[i]/2)*U);TGT.push([p[0]-LR/2,p[1]-(rows-1)*RR]);}}
// ---- physics: cloud (charge, centre, pull of usable edges) + rope (picked edges, bending) + unfold
const rnd=lcg(11),R0=16*Math.sqrt(N);
const nodes=Array.from({length:N},(_,i)=>{const a=rnd()*6.2832,r=R0*Math.sqrt(rnd());return{index:i,x:r*Math.cos(a),y:r*Math.sin(a)};});
const ang0=Array.from({length:N},(_,i)=>rnd()*Math.PI);
let unf=0,manV=0,manT=0,unfMove=1,still=0;
function pull(alpha){for(let u=0;u<usableN;u++){const q=usable[u],a=nodes[EI[q]],b=nodes[EJ[q]];
 let dx=b.x+b.vx-a.x-a.vx,dy=b.y+b.vy-a.y-a.vy;const d=Math.hypot(dx,dy)||1e-6,r0=34+(100-EO[q])*0.15,f=(d-r0)/d*alpha*0.3*EO[q]/100;dx*=f;dy*=f;
 b.vx-=dx*.5;b.vy-=dy*.5;a.vx+=dx*.5;a.vy+=dy*.5;}}
function rope(){for(let it=0;it<4;it++){
 for(let k=0;k<applied;k++){const e=EV[k],a=nodes[e[0]],b=nodes[e[1]];let dx=b.x-a.x,dy=b.y-a.y;const d=Math.hypot(dx,dy)||1e-6,c=(d-(LEN[e[0]]-e[2])*U)/d*0.45;
  a.x+=dx*c;a.y+=dy*c;b.x-=dx*c;b.y-=dy*c;}
 for(let i=0;i<N;i++){const p=pred[i],s=succ[i];if(p<0||s<0)continue;const a=nodes[p],b=nodes[i],c=nodes[s],mx=((a.x+c.x)/2-b.x)*0.08,my=((a.y+c.y)/2-b.y)*0.08;
  b.x+=mx;b.y+=my;a.x-=mx/2;a.y-=my/2;c.x-=mx/2;c.y-=my/2;}}}
function unfold(){unfMove=0;if(unf<0.0005)return;for(let i=0;i<N;i++){const n=nodes[i],dx=(TGT[i][0]-n.x)*unf,dy=(TGT[i][1]-n.y)*unf;n.x+=dx;n.y+=dy;unfMove+=Math.abs(dx)+Math.abs(dy);}unfMove/=N;}
const sim=d3.forceSimulation(nodes).randomSource(lcg(7)).alphaDecay(0.02).velocityDecay(0.35).stop()
 .force("charge",d3.forceManyBody().strength(-26).distanceMax(150).theta(0.9))
 .force("cx",d3.forceX(0).strength(0.012)).force("cy",d3.forceY(0).strength(0.012))
 .force("pull",pull).force("rope",rope).force("unfold",unfold);
goto(0);sim.alpha(1);for(let t=0;t<150;t++)sim.tick();sim.alpha(0.05);
// ---- camera
let K=1,TX=0,TY=0,follow=true;
function bbox(){let a=1e9,b=-1e9,c=1e9,d=-1e9;for(const n of nodes){if(n.x<a)a=n.x;if(n.x>b)b=n.x;if(n.y<c)c=n.y;if(n.y>d)d=n.y;}return[a,b,c,d];}
function camTarget(){const[a,b,c,d]=bbox(),k=Math.min(3,(W-80)/Math.max(1,b-a),(H-130)/Math.max(1,d-c));return[k,W/2-k*(a+b)/2,(H+90)/2-k*(c+d)/2];}
function snapCam(){const t=camTarget();K=t[0];TX=t[1];TY=t[2];}
snapCam();
// ---- drawing
let hover=-1,dirty=true,raf=0,drag=null,pan=null,playing=false,i0=0,prevI=0,lensT0=0,lensDur=0;
const sx=n=>n.x*K+TX,sy=n=>n.y*K+TY;
const NW=i=>Math.max(5,LEN[i]*U*K),NH=()=>Math.max(2.4,RH*K),LN=()=>Math.max(NH()+0.8,LANE*K);
function pose(i){const n=nodes[i],p=pred[i],s=succ[i];let th;
 if(p>=0||s>=0){const a=s>=0?nodes[s]:n,b=p>=0?nodes[p]:n;th=Math.atan2(a.y-b.y,a.x-b.x);}else th=ang0[i];
 let x=sx(n),y=sy(n);if(p>=0||s>=0){const lane=(FO_RANK[i]%3)-1,o=lane*LN();x+=-Math.sin(th)*o;y+=Math.cos(th)*o;}return[x,y,th];}
const FO_RANK=D.fr;
function rrect(w,h,r){const q=Math.min(r,w/2,h/2);ctx.beginPath();ctx.moveTo(-w/2+q,-h/2);ctx.arcTo(w/2,-h/2,w/2,h/2,q);ctx.arcTo(w/2,h/2,-w/2,h/2,q);ctx.arcTo(-w/2,h/2,-w/2,-h/2,q);ctx.arcTo(-w/2,-h/2,w/2,-h/2,q);ctx.closePath();}
function place(x,y,th){const c=Math.cos(th)*dpr,s=Math.sin(th)*dpr;ctx.setTransform(c,s,-s,c,x*dpr,y*dpr);}
function draw(){
 ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,W,H);
 const dim=hover>=0;ctx.lineWidth=0.9;
 for(let c=0;c<2;c++){ctx.beginPath();for(let u=0;u<usableN;u++){const q=usable[u];if((EO[q]>=75)!==(c===1))continue;const a=nodes[EI[q]],b=nodes[EJ[q]];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}
  ctx.strokeStyle=c?"rgba(80,118,172,"+(dim?0.14:0.38)+")":"rgba(120,150,190,"+(dim?0.07:0.2)+")";ctx.stroke();}
 ctx.beginPath();for(let k=0;k<applied;k++){const e=EV[k],a=nodes[e[0]],b=nodes[e[1]];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}
 ctx.lineWidth=1.4;ctx.strokeStyle="rgba(224,123,0,"+(dim?0.35:0.8)+")";ctx.stroke();
 const nb=new Set();if(hover>=0){nb.add(hover);for(const[j]of outL[hover])nb.add(j);for(const[j]of inL[hover])nb.add(j);}
 const h=NH(),P=new Array(N);
 for(let i=0;i<N;i++){const p=pose(i);P[i]=p;const rt=rootOf[i];ctx.globalAlpha=dim&&!nb.has(i)?0.3:1;
  ctx.fillStyle=pcnt[rt]>1?"hsl("+(((rt*137.5)%360)|0)+",62%,50%)":"#5b9bd5";place(p[0],p[1],p[2]);rrect(NW(i),h,Math.min(3,h/2));ctx.fill();}
 ctx.globalAlpha=1;ctx.setTransform(dpr,0,0,dpr,0,0);
 if(applied){const e=EV[applied-1];ctx.lineWidth=2.5;ctx.strokeStyle="#e07b00";for(const i of [e[0],e[1]]){const p=P[i];place(p[0],p[1],p[2]);rrect(NW(i)+6,h+6,Math.min(5,(h+6)/2));ctx.stroke();}ctx.setTransform(dpr,0,0,dpr,0,0);}
 if(hover>=0){const n=nodes[hover];ctx.beginPath();for(const t of [...outL[hover],...inL[hover]]){if(state[t[2]]===0)continue;const o=nodes[t[0]];ctx.moveTo(sx(n),sy(n));ctx.lineTo(sx(o),sy(o));}
  ctx.lineWidth=1.8;ctx.strokeStyle="rgba(27,42,65,.8)";ctx.stroke();const p=P[hover];place(p[0],p[1],p[2]);rrect(NW(hover)+4,h+4,Math.min(5,(h+4)/2));ctx.lineWidth=2.5;ctx.strokeStyle="#1b2a41";ctx.stroke();ctx.setTransform(dpr,0,0,dpr,0,0);}
}
// ---- lens: the two reads of the latest merge, letter by letter
const NC={A:"#6bbf59",C:"#5b9bd5",G:"#f2b134",T:"#e0605e"};
function strip(s,x,y,w,h,cw){for(let k=0;k<s.length;k++){lc.fillStyle=NC[s[k]]||"#9aa5b1";lc.fillRect(x+k*cw,y,cw+0.5,h);}}
function drawLens(now){
 lc.setTransform(dpr,0,0,dpr,0,0);lc.clearRect(0,0,LW,LH);
 if(!applied){lensBox.style.opacity=0;return false;}lensBox.style.opacity=1;
 const e=EV[applied-1],A=reads[e[0]],B=reads[e[1]],ov=e[2],x0=116,cw=(LW-x0-12)/(A.length+B.length-ov);
 let p=lensDur>0?Math.min(1,(now-lensT0)/lensDur):1;const q=1-Math.pow(1-p,3);
 const xa=x0,xbF=x0+(A.length-ov)*cw,xb=xbF+(1-q)*B.length*cw*0.6,ya=12,yb=38,ym=72,hh=16;
 lc.font="600 12px ui-monospace,Menlo,monospace";lc.textBaseline="middle";lc.fillStyle="#52606d";
 lc.fillText("#"+e[0]+" · "+fmt(e[6]),8,ya+hh/2);lc.fillText("#"+e[1]+" · "+fmt(e[7]),8,yb+hh/2);lc.fillStyle="#1b2a41";lc.fillText("→ "+fmt(e[3]),8,ym+hh/2);
 lc.save();lc.beginPath();lc.rect(x0,0,LW-x0,LH);lc.clip();
 strip(A,xa,ya,0,hh,cw);strip(B,xb,yb,0,hh,cw);
 const bx=xb,bw=xa+A.length*cw-xb;
 if(bw>0){lc.fillStyle="rgba(224,123,0,.14)";lc.fillRect(bx,ya-4,bw,yb+hh-ya+8);lc.strokeStyle="#e07b00";lc.lineWidth=1.5;lc.strokeRect(bx,ya-4,bw,yb+hh-ya+8);}
 if(p>0.55){lc.globalAlpha=(p-0.55)/0.45;strip(A+B.slice(ov),xa,ym,0,hh,cw);lc.globalAlpha=1;}
 lc.restore();
 if(bw>0){lc.font="700 13px system-ui,sans-serif";lc.textAlign="center";const tx=bx+bw/2,ty=(ya+hh+yb)/2;lc.lineWidth=4;lc.strokeStyle="#f8fafc";lc.strokeText(ov,tx,ty);lc.fillStyle="#e07b00";lc.fillText(ov,tx,ty);lc.textAlign="left";}
 return p<1;}
// ---- sparklines: pieces, longest, overlap
const OVS=EV.map(e=>e[2]),OVMAX=Math.max(...OVS,1),OVMIN=Math.min(...OVS,OVMAX);
const SER=[["#5b9bd5",t=>E?(E-t)/E:0],["#e07b00",t=>E?(MAXL[t]-MAXL[0])/Math.max(1,MAXL[E]-MAXL[0]):0],["#7a4fb5",t=>t?(OVS[t-1]-OVMIN)/Math.max(1,OVMAX-OVMIN):1]];
const CURVE=SER.map(([c,f])=>{const a=[];for(let x=0;x<SW;x++)a.push(f(Math.round(x*E/(SW-1))));return a;});
function drawSpark(){sc.setTransform(dpr,0,0,dpr,0,0);sc.clearRect(0,0,SW,SH);const yy=v=>SH-4-v*(SH-8);
 SER.forEach(([c],s)=>{sc.beginPath();CURVE[s].forEach((v,x)=>x?sc.lineTo(x,yy(v)):sc.moveTo(x,yy(v)));sc.strokeStyle=c;sc.lineWidth=1.6;sc.stroke();});
 const px=E?i0/E*(SW-1):0;sc.strokeStyle="rgba(27,42,65,.35)";sc.lineWidth=1;sc.beginPath();sc.moveTo(px,0);sc.lineTo(px,SH);sc.stroke();
 SER.forEach(([c,f],s)=>{sc.fillStyle=c;sc.beginPath();sc.arc(px,yy(f(i0)),3.2,0,6.2832);sc.fill();});}
// ---- main loop
function loop(now){raf=0;let act=false;
 manV+=(manT-manV)*0.06;if(Math.abs(manT-manV)<0.002)manV=manT;
 const run=sim.alpha()>sim.alphaMin()||(unf>0.0005&&still<25);
 if(run){sim.tick();sim.tick();act=true;dirty=true;still=unf>0.0005&&unfMove<0.01?still+1:0;}
 if(follow){const t=camTarget();if(Math.abs(t[0]-K)>1e-4||Math.abs(t[1]-TX)>0.05||Math.abs(t[2]-TY)>0.05){K+=(t[0]-K)*0.12;TX+=(t[1]-TX)*0.12;TY+=(t[2]-TY)*0.12;act=true;dirty=true;}}
 if(dirty){draw();dirty=false;}
 if(drawLens(now))act=true;
 if(act||playing||drag!==null||manV!==manT)raf=requestAnimationFrame(loop);}
function kick(){dirty=true;if(!raf)raf=requestAnimationFrame(loop);}
function reheat(a){if(sim.alpha()<a)sim.alpha(a);still=0;}
// ---- interaction
const pt=e=>{const r=cv.getBoundingClientRect();return[(e.clientX-r.left)*W/r.width,(e.clientY-r.top)*H/r.height];};
function pick(mx,my){let b=-1,bd=Infinity;for(let i=0;i<N;i++){const p=pose(i),d=Math.hypot(p[0]-mx,p[1]-my)-NW(i)/2;if(d<bd){bd=d;b=i;}}return bd<=7?b:-1;}
function setTip(i,mx,my){if(i<0){tip.style.display="none";return;}
 const o=outL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]),n=inL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]);
 let s="#"+i;if(o.length)s+="<br>→ #"+o[0][0]+" · "+o[0][1];if(n.length)s+="<br>← #"+n[0][0]+" · "+n[0][1];
 const r=rootOf[i];if(pcnt[r]>1)s+="<br>"+fmt(pcnt[r])+" reads · "+fmt(plen[r])+" nt";
 tip.innerHTML=s;tip.style.display="block";tip.style.left=Math.min(W-170,mx+14)+"px";tip.style.top=Math.min(H-90,my+14)+"px";}
cv.addEventListener("wheel",e=>{e.preventDefault();follow=false;fitBtn.classList.remove("on");const[mx,my]=pt(e),nk=Math.max(0.05,Math.min(60,K*Math.exp(-e.deltaY*0.0015)));TX=mx-(mx-TX)*nk/K;TY=my-(my-TY)*nk/K;K=nk;kick();},{passive:false});
cv.addEventListener("mousedown",e=>{const[mx,my]=pt(e),i=pick(mx,my);if(i>=0){drag=i;nodes[i].fx=nodes[i].x;nodes[i].fy=nodes[i].y;sim.alphaTarget(0.2);reheat(0.2);}else{pan=[mx,my,TX,TY];follow=false;fitBtn.classList.remove("on");cv.classList.add("drag");}kick();});
window.addEventListener("mousemove",e=>{const[mx,my]=pt(e);
 if(drag!==null){const n=nodes[drag];n.fx=(mx-TX)/K;n.fy=(my-TY)/K;kick();}
 else if(pan){TX=pan[2]+mx-pan[0];TY=pan[3]+my-pan[1];kick();}
 else if(e.target===cv){const i=pick(mx,my);if(i!==hover){hover=i;setTip(i,mx,my);kick();}else if(i>=0)setTip(i,mx,my);}});
window.addEventListener("mouseup",()=>{if(drag!==null){nodes[drag].fx=null;nodes[drag].fy=null;drag=null;sim.alphaTarget(0);}pan=null;cv.classList.remove("drag");kick();});
cv.addEventListener("mouseleave",()=>{if(hover>=0){hover=-1;tip.style.display="none";kick();}});
// ---- stepping
let timer=null,acc=0;
function show(k){const before=i0;i0=Math.max(0,Math.min(E,Math.round(k)));goto(i0);
 const f=E?MAXL[i0]/MAXL[E]:1,x=Math.min(1,Math.max(0,(f-0.55)/0.4));unf=0.08*Math.max(x*x*(3-2*x),manV);
 if(i0!==before){reheat(i0>before?0.3:0.15);const fast=playing&&spd.value>=20;lensT0=performance.now();lensDur=fast?0:(playing?Math.min(420,800/spd.value):420);}
 $(".s-merge").textContent=fmt(i0)+" / "+fmt(E);$(".s-left").textContent=fmt(N-i0);$(".s-long").textContent=fmt(MAXL[i0]);$(".s-ov").textContent=i0?EV[i0-1][2]:"–";
 rng.value=i0;cnt.textContent="step "+fmt(i0)+" / "+fmt(E);drawSpark();kick();}
function stop(){clearInterval(timer);timer=null;playing=false;pp.textContent="▶ Play";}
function play(){if(i0>=E)show(0);pp.textContent="⏸ Pause";acc=0;playing=true;timer=setInterval(()=>{acc+=spd.value/20;const st=Math.floor(acc);acc-=st;if(i0>=E){stop();return;}show(i0+Math.max(st,0));if(i0>=E)stop();},50);}
const toggle=()=>timer?stop():play();
const fitBtn=$("button[data-a=fit]"),unfBtn=$("button[data-a=unfold]");
function refit(){follow=true;fitBtn.classList.add("on");kick();}
function shake(){reheat(0.8);kick();}
function unfoldBtn(){manT=manT?0:1;unfBtn.classList.toggle("on",!!manT);show(i0);}
rng.max=E;rng.oninput=()=>{stop();show(+rng.value);};
spd.onchange=()=>{if(timer){stop();play();}};
root.querySelectorAll("button").forEach(b=>b.onclick=()=>{const a=b.dataset.a;if(a==="pp")return toggle();if(a==="fit")return refit();if(a==="shake")return shake();if(a==="unfold")return unfoldBtn();stop();show(a==="first"?0:a==="last"?E:i0+(a==="next"?1:-1));root.focus();});
root.addEventListener("keydown",e=>{const m={" ":toggle,ArrowRight:()=>{stop();show(i0+1)},ArrowLeft:()=>{stop();show(i0-1)},PageDown:()=>{stop();show(i0+50)},PageUp:()=>{stop();show(i0-50)},Home:()=>{stop();show(0)},End:()=>{stop();show(E)},f:refit,F:refit,s:shake,S:shake,u:unfoldBtn,U:unfoldBtn}[e.key];
 if(m){m();e.preventDefault();e.stopPropagation();}},true);
if(window.Reveal&&Reveal.on)Reveal.on("slidechanged",stop);
fitBtn.classList.add("on");
window.__vz={show,nodes,sim,get applied(){return applied;},get unf(){return unf;},get K(){return K;}};
show(0);
})();
</script></div>
"""

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq")
    ap.add_argument("--min-overlap", type=int, default=20, help="minimum overlap greedy may use")
    ap.add_argument("--show-min", type=int, default=50, help="draw only edges with at least this overlap")
    ap.add_argument("--standalone", action="store_true")
    ap.add_argument("-o", "--out", default="graph.html")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    d3 = open(os.path.join(here, "vendor", "d3-force.bundle.min.js"), encoding="utf-8").read().replace("</script", "<\\/script")
    reads = list(dict.fromkeys(read_fastq(a.fastq)))
    n = len(reads)
    ev, contigs, _ = greedy_events(reads, a.min_overlap)
    edges = [(i, j, o) for o, i, j in overlaps(reads, a.show_min)]
    have = {(i, j) for i, j, _ in edges}
    edges += [(e["a"], e["b"], e["ov"]) for e in ev if (e["a"], e["b"]) not in have]   # picks are always drawn
    # chain lengths before each merge (for the lens), running longest piece
    L = {i: len(r) for i, r in enumerate(reads)}
    cur = max(len(r) for r in reads); maxlen = [cur]; evo = []
    for e in ev:
        la, lb = L[e["root"]], L[e["rb"]]
        L[e["root"]] = e["len"]
        cur = max(cur, e["len"]); maxlen.append(cur)
        evo.append([e["a"], e["b"], e["ov"], e["len"], e["root"], e["rb"], la, lb])
    # final position of every read: offset along its final contig (contigs laid end to end), rank along the path
    succ = {e["a"]: (e["b"], e["ov"]) for e in ev}
    has_pred = {e["b"] for e in ev}
    fo, fr, base = [0] * n, [0] * n, 0
    for h in range(n):
        if h in has_pred:
            continue
        c, off, k = h, base, 0
        while True:
            fo[c], fr[c] = off, k
            if c not in succ:
                break
            nx, o = succ[c]; off += len(reads[c]) - o; c = nx; k += 1
        base = off + len(reads[c]) + 200
    data = {"n": n, "reads": reads, "lengths": [len(r) for r in reads],
            "ei": [e[0] for e in edges], "ej": [e[1] for e in edges], "eo": [e[2] for e in edges],
            "ev": evo, "maxlen": maxlen, "fo": fo, "fr": fr}
    out = (TEMPLATE.replace("__ID__", "vz" + uuid.uuid4().hex[:8])
           .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__D3__", d3))
    if a.standalone:
        out = f"<!doctype html><meta charset=utf-8><title>{html.escape(a.fastq)}</title><body style='margin:20px'>{out}"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"{n} reads, {len(edges)} edges shown, {len(ev)} merges -> {a.out} ({len(out)//1024} KB)")

if __name__ == "__main__":
    main()
