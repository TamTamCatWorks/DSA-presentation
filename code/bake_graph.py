"""bake_graph.py - bake the overlap graph of a whole FASTQ into ONE interactive graph-view widget.

Reads are rectangles; a line joins two reads whose ends overlap by >= --show-min letters.
Each contig is a snake: head = first read, tail = last read. Heads hunt the best free tail and
bodies follow the leader (d3-force cloud + follow-the-leader constraints), edges that can no longer
be used disappear. Level of detail: snake strokes -> read rectangles -> nucleotide cells -> letters,
rendered in place and only for what is on screen. Focus (Z) zooms to the latest merge.

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
#__ID__ .gr-zm{right:12px;bottom:12px;padding:3px 9px;font:600 14px ui-monospace,Menlo,monospace;color:#52606d;font-variant-numeric:tabular-nums}
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
 <div class="gr-ov gr-zm"></div>
 <div class="gr-tip"></div></div>
<div class="vz-bar">
<button data-a="first" title="First (Home)">⏮</button><button data-a="prev" title="Back (←)">◀</button>
<button data-a="pp" class="vz-pp" title="Play / pause (Space)">▶ Play</button>
<button data-a="next" title="Step (→)">▶|</button><button data-a="last" title="Last (End)">⏭</button>
<input type="range" min="0" value="0"><span class="vz-n"></span>
<select title="Merges per second"><option value="1">1 / s</option><option value="5">5 / s</option><option value="20">20 / s</option><option value="100" selected>100 / s</option><option value="400">400 / s</option></select>
<button class="sm" data-a="fit" title="Fit (F)">Fit</button><button class="sm" data-a="focus" title="Follow the latest merge (Z)">Focus</button><button class="sm" data-a="unfold" title="Unfold (U)">Unfold</button><button class="sm" data-a="shake" title="Shake (S)">Shake</button>
</div>
<script>
(function(){
const D=__DATA__, root=document.getElementById("__ID__"), $=s=>root.querySelector(s);
const N=D.n,EI=D.ei,EJ=D.ej,EO=D.eo,EV=D.ev,M=EI.length,E=EV.length,W=1200,H=560;
const reads=D.reads,LEN=D.lengths,FO=D.fo,LNE=D.ln,U=0.25,RH=5,MAXL=D.maxlen,STR=0.12;
const cv=$("canvas.gr-main"),ctx=cv.getContext("2d"),tip=$(".gr-tip"),rng=$("input"),cnt=$(".vz-n"),pp=$(".vz-pp"),spd=$("select"),zm=$(".gr-zm");
const spk=$(".gr-sp canvas"),sc=spk.getContext("2d");
const dpr=Math.min(2,window.devicePixelRatio||1);
cv.width=W*dpr;cv.height=H*dpr;const SW=260,SH=64;spk.width=SW*dpr;spk.height=SH*dpr;
const fmt=n=>n.toLocaleString("en-US");
function lcg(s0){let s=s0;return()=>(s=(s*1664525+1013904223)%4294967296)/4294967296;}
// ---- state replayed from greedy's merge log. A contig is a snake: head = first read (free predecessor slot), tail = last read (free successor slot).
const outL=Array.from({length:N},()=>[]),inL=Array.from({length:N},()=>[]);
for(let q=0;q<M;q++){outL[EI[q]].push([EJ[q],EO[q],q]);inL[EJ[q]].push([EI[q],EO[q],q]);}
let succ,pred,par,pcnt,plen,ovOf,applied=0,rootOf=new Int32Array(N),state=new Uint8Array(M),usable=new Int32Array(M),usableN=0;
const glide=new Array(N).fill(null);let glideN=0,GDUR=100;
const bestIn=new Int32Array(N),headOf=new Int32Array(N),ord=new Int32Array(N),heads=new Int32Array(N),tails=new Int32Array(N);let ordN=0,headsN=0;
function reset(){succ=new Int32Array(N).fill(-1);pred=new Int32Array(N).fill(-1);par=Int32Array.from({length:N},(_,i)=>i);pcnt=new Int32Array(N).fill(1);plen=Int32Array.from(LEN);ovOf=new Int32Array(N);glide.fill(null);applied=0;}
function find(x){while(par[x]!==x){par[x]=par[par[x]];x=par[x];}return x;}
reset();let INIT=null;
function goto(k){const wasReset=k<applied;if(wasReset)reset();const from=applied;
 while(applied<k){const e=EV[applied];succ[e[0]]=e[1];pred[e[1]]=e[0];ovOf[e[1]]=e[2];par[e[5]]=e[4];pcnt[e[4]]+=pcnt[e[5]];plen[e[4]]=e[3];applied++;}
 for(let i=0;i<N;i++)rootOf[i]=find(i);usableN=0;bestIn.fill(-1);
 for(let q=0;q<M;q++){const a=EI[q],b=EJ[q];state[q]=succ[a]===b?2:(succ[a]<0&&pred[b]<0&&rootOf[a]!==rootOf[b])?1:0;
  if(state[q]===1){usable[usableN++]=q;if(bestIn[b]<0||EO[q]>EO[bestIn[b]])bestIn[b]=q;}}
 ordN=0;headsN=0;for(let h=0;h<N;h++)headOf[h]=h;
 for(let h=0;h<N;h++){if(pred[h]>=0||succ[h]<0)continue;heads[headsN]=h;let c=h;while(succ[c]>=0){c=succ[c];ord[ordN++]=c;headOf[c]=h;}tails[headsN++]=c;}
 if(!INIT)return;
 if(wasReset)restore();
 if(applied-(wasReset?0:from)>40)settle();           // big jump: lay snakes out as calm rings instead of animating hundreds of merges
 else for(let j=wasReset?0:from;j<applied;j++){const e=EV[j],nb=nodes[e[1]];glide[e[1]]={a:e[0],x0:nb.x,y0:nb.y,s:0,d:GDUR};glideN++;}}
function restore(){for(let i=0;i<N;i++){const n=nodes[i];n.x=INIT[i][0];n.y=INIT[i][1];n.vx=n.vy=0;}glide.fill(null);}
function settle(){restore();for(let s=0;s<headsN;s++){const h=heads[s];let L=0,c=h;while(succ[c]>=0){const nx=succ[c];L+=(LEN[c]-ovOf[nx])*U;c=nx;}
  const H0=INIT[h],d=Math.hypot(H0[0],H0[1])||1,hs=((h*2654435761)>>>0)/4294967296,psi=Math.atan2(d>1?-H0[1]/d:0,d>1?-H0[0]/d:1)+(hs-0.5)*1.6,k=(hs<0.5?-1:1)*0.9/Math.max(L,1e-3);
  let x=H0[0],y=H0[1],sArc=0;c=h;for(;;){nodes[c].x=x;nodes[c].y=y;const nx=succ[c];if(nx<0)break;const ds=(LEN[c]-ovOf[nx])*U,a=psi+k*(sArc+ds/2);x+=Math.cos(a)*ds;y+=Math.sin(a)*ds;sArc+=ds;c=nx;}}}   // gentle arcs
// ---- unfold target: serpentine, arc length = position in the final contig
const LR=1400,RR=70,CYC=LR+Math.PI*RR;
function pathPt(s){const k=Math.floor(s/CYC),t=s-k*CYC,dir=k%2===0?1:-1,y0=k*2*RR;
 if(t<=LR)return[dir>0?t:LR-t,y0];const a=(t-LR)/RR;return[dir>0?LR+RR*Math.sin(a):-RR*Math.sin(a),y0+RR-RR*Math.cos(a)];}
const TGT=[];{let tot=0;for(let i=0;i<N;i++)tot=Math.max(tot,FO[i]+LEN[i]);const rows=Math.ceil(tot*U/CYC);
 for(let i=0;i<N;i++){const p=pathPt((FO[i]+LEN[i]/2)*U);TGT.push([p[0]-LR/2,p[1]-(rows-1)*RR]);}}
// ---- physics. The cloud is laid out once and then stays calm. A merge animates only the joining snake: its head glides on a curve to the other tail, its body follows the leader.
const rnd=lcg(11),R0=16*Math.sqrt(N);
const nodes=Array.from({length:N},(_,i)=>{const a=rnd()*6.2832,r=R0*Math.sqrt(rnd());return{index:i,x:r*Math.cos(a),y:r*Math.sin(a)};});
const ang0=Array.from({length:N},()=>rnd()*Math.PI);
let unf=0,manV=0,manT=0,unfMove=1,still=0;
function glideTarget(g,b){const a=g.a,A=nodes[a],p=pred[a];let dx,dy;
 if(p>=0){dx=A.x-nodes[p].x;dy=A.y-nodes[p].y;}else{dx=g.x0-A.x;dy=g.y0-A.y;}
 const d=Math.hypot(dx,dy)||1;dx/=d;dy/=d;const r0=(LEN[a]-ovOf[b])*U;return[A.x+dx*r0,A.y+dy*r0,dx,dy];}
function glidePos(g,b,e){const T=glideTarget(g,b),dist=Math.hypot(T[0]-g.x0,T[1]-g.y0),c=Math.max(30,dist*0.5),cx=T[0]+T[2]*c,cy=T[1]+T[3]*c,u=1-e;
 return[u*u*g.x0+2*u*e*cx+e*e*T[0],u*u*g.y0+2*u*e*cy+e*e*T[1]];}
function rope(){glideN=0;
 for(let k=0;k<ordN;k++){const v=ord[k],u=pred[v],nv=nodes[v],g=glide[v];
  if(g){g.s=Math.min(1,g.s+1/g.d);const p=glidePos(g,v,g.s*g.s*(3-2*g.s));nv.x=p[0];nv.y=p[1];nv.vx=nv.vy=0;if(g.s>=1)glide[v]=null;else glideN++;continue;}
  const nu=nodes[u];let dx=nv.x-nu.x,dy=nv.y-nu.y,d=Math.hypot(dx,dy);if(d<1e-6){dx=1;dy=0;d=1;}dx/=d;dy/=d;
  const p=pred[u];if(p>=0){const np=nodes[p];let tx=nu.x-np.x,ty=nu.y-np.y;const td=Math.hypot(tx,ty)||1;dx=dx*(1-STR)+tx/td*STR;dy=dy*(1-STR)+ty/td*STR;const m=Math.hypot(dx,dy)||1;dx/=m;dy/=m;}
  const r0=(LEN[u]-ovOf[v])*U;nv.x=nu.x+dx*r0;nv.y=nu.y+dy*r0;nv.vx=nu.vx;nv.vy=nu.vy;}
 unfMove=0;if(unf>0.0005){const w=unf*1.5;for(let i=0;i<N;i++){const n=nodes[i],dx=(TGT[i][0]-n.x)*w,dy=(TGT[i][1]-n.y)*w;n.x+=dx;n.y+=dy;unfMove+=Math.abs(dx)+Math.abs(dy);}unfMove/=N;}}
function spring(alpha){for(let u=0;u<usableN;u++){const q=usable[u],a=nodes[EI[q]],b=nodes[EJ[q]];
 let dx=b.x+b.vx-a.x-a.vx,dy=b.y+b.vy-a.y-a.vy;const d=Math.hypot(dx,dy)||1e-6,r0=34+(100-EO[q])*0.15,f=(d-r0)/d*alpha*0.3*EO[q]/100;dx*=f;dy*=f;
 b.vx-=dx*.5;b.vy-=dy*.5;a.vx+=dx*.5;a.vy+=dy*.5;}}
const sim=d3.forceSimulation(nodes).randomSource(lcg(7)).alphaDecay(0.02).velocityDecay(0.35).stop()
 .force("charge",d3.forceManyBody().strength(-26).distanceMax(150).theta(0.9))
 .force("cx",d3.forceX(0).strength(0.012)).force("cy",d3.forceY(0).strength(0.012)).force("spring",spring);
goto(0);sim.alpha(1);for(let t=0;t<150;t++)sim.tick();sim.alpha(0);INIT=nodes.map(n=>[n.x,n.y]);
// ---- camera: 0 free, 1 fit, 2 focus on the latest merge
let K=1,TX=0,TY=0,mode=1;
function fitTarget(){let a=1e9,b=-1e9,c=1e9,d=-1e9;for(const n of nodes){if(n.x<a)a=n.x;if(n.x>b)b=n.x;if(n.y<c)c=n.y;if(n.y>d)d=n.y;}
 const k=Math.min(3,(W-80)/Math.max(1,b-a),(H-130)/Math.max(1,d-c));return[k,W/2-k*(a+b)/2,(H+90)/2-k*(c+d)/2];}
function focusTarget(){if(!applied)return fitTarget();const e=EV[applied-1],a=nodes[e[0]],b=nodes[e[1]],k=Math.min(60,1000/((LEN[e[0]]+LEN[e[1]]-e[2])*U)),cx=(a.x+b.x)/2,cy=(a.y+b.y)/2;return[k,W/2-k*cx,H/2-k*cy];}
function camTarget(){return mode===2?focusTarget():fitTarget();}
{const t=fitTarget();K=t[0];TX=t[1];TY=t[2];}
// ---- drawing with level of detail: ribbons -> pile-up of reads -> nucleotide cells -> letters; everything off screen is skipped
let hover=-1,dirty=true,raf=0,drag=null,dragInfo=null,pan=null,playing=false,i0=0;
const sx=n=>n.x*K+TX,sy=n=>n.y*K+TY,hue=(r,s,l)=>"hsl("+(((r*137.5)%360)|0)+","+s+"%,"+l+"%)";
const NW=i=>Math.max(5,LEN[i]*U*K);
const NH=()=>K<1.6?Math.max(2.4,RH*K):K<6.4?Math.min(10,Math.max(4,RH*K*0.5)):Math.max(14,U*K*1.6);
const LCEN=(Math.max.apply(null,LNE)+1)/2;
const NUC={A:"#4caf50",C:"#3f7fc1",G:"#e0a21b",T:"#d9534f"};
const ribW=()=>Math.max(3,Math.min(16,RH*K*1.2+1.5));
function pose(i){const n=nodes[i],p=pred[i],s=succ[i],chain=p>=0||s>=0;let th;
 if(chain){const a=s>=0?nodes[s]:n,b=p>=0?nodes[p]:n;th=Math.atan2(a.y-b.y,a.x-b.x);}else th=ang0[i];
 let x=sx(n),y=sy(n);if(chain&&K>=1.6){const o=(LNE[i]-LCEN)*(NH()*1.2+0.5);x+=-Math.sin(th)*o;y+=Math.cos(th)*o;}return[x,y,th];}
function rrect(w,h,r){const q=Math.min(r,w/2,h/2);ctx.beginPath();ctx.moveTo(-w/2+q,-h/2);ctx.arcTo(w/2,-h/2,w/2,h/2,q);ctx.arcTo(w/2,h/2,-w/2,h/2,q);ctx.arcTo(-w/2,h/2,-w/2,-h/2,q);ctx.arcTo(-w/2,-h/2,w/2,-h/2,q);ctx.closePath();}
function place(x,y,th){const c=Math.cos(th)*dpr,s=Math.sin(th)*dpr;ctx.setTransform(c,s,-s,c,x*dpr,y*dpr);}
const flat=()=>ctx.setTransform(dpr,0,0,dpr,0,0);
const inView=(x,y,m)=>x>-m&&x<W+m&&y>-m&&y<H+m;
function chainPath(h){const p=new Path2D();let c=h;const xs=[sx(nodes[c])],ys=[sy(nodes[c])];
 while(succ[c]>=0&&!glide[succ[c]]){c=succ[c];const x=sx(nodes[c]),y=sy(nodes[c]),dx=x-xs[xs.length-1],dy=y-ys[ys.length-1];if(succ[c]<0||glide[succ[c]]||dx*dx+dy*dy>=9){xs.push(x);ys.push(y);}}
 const n=xs.length;p.moveTo(xs[0],ys[0]);if(n===2)p.lineTo(xs[1],ys[1]);else{for(let i=1;i<n-1;i++)p.quadraticCurveTo(xs[i],ys[i],(xs[i]+xs[i+1])/2,(ys[i]+ys[i+1])/2);p.lineTo(xs[n-1],ys[n-1]);}return [p,c];}
function draw(){
 flat();ctx.clearRect(0,0,W,H);
 const dim=hover>=0,h=NH(),cw=U*K,PILE=K>=1.6,CELLS=cw>=1.6,LETTERS=cw>=6.5,BODY=K<6.4,wR=ribW();
 ctx.lineWidth=0.9;
 for(let c=0;c<2;c++){ctx.beginPath();for(let u=0;u<usableN;u++){const q=usable[u];if((EO[q]>=75)!==(c===1))continue;const a=nodes[EI[q]],b=nodes[EJ[q]];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}
  ctx.strokeStyle=c?"rgba(80,118,172,"+(dim?0.14:0.38)+")":"rgba(120,150,190,"+(dim?0.07:0.2)+")";ctx.stroke();}
 // candidate arrows: free tail -> head, labelled with the overlap
 if(K>=0.7){ctx.fillStyle="rgba(60,100,160,.75)";ctx.font="700 11px system-ui,sans-serif";ctx.textAlign="center";
  for(let i=0;i<N;i++){const q=bestIn[i];if(q<0)continue;const a=nodes[EI[q]],b=nodes[i],ax=sx(a),ay=sy(a),bx=sx(b),by=sy(b);if(!inView(bx,by,40)&&!inView(ax,ay,40))continue;
   const dx=bx-ax,dy=by-ay,d=Math.hypot(dx,dy);if(d<4)continue;const ux=dx/d,uy=dy/d,tx=bx-ux*Math.min(7,d*0.4),ty=by-uy*Math.min(7,d*0.4);
   ctx.beginPath();ctx.moveTo(tx+ux*5,ty+uy*5);ctx.lineTo(tx-ux*3-uy*3.5,ty-uy*3+ux*3.5);ctx.lineTo(tx-ux*3+uy*3.5,ty-uy*3-ux*3.5);ctx.fill();
   if(K>=1.5&&d>26)ctx.fillText(EO[q],(ax+bx)/2-uy*8,(ay+by)/2+ux*8+4);}}
 // contig ribbons; a snake whose head is still gliding to its new tail is drawn as its own ribbon until it arrives
 const segs=[];
 if(BODY){const starts=Array.from(heads.subarray(0,headsN));for(let v=0;v<N;v++)if(glide[v])starts.push(v);
  ctx.lineCap="round";ctx.lineJoin="round";ctx.globalAlpha=dim?0.4:(PILE?0.35:1);
  for(const st of starts){const rt=rootOf[st],r=chainPath(st),P2=r[0];segs.push([st,r[1]]);ctx.lineWidth=wR+2;ctx.strokeStyle=hue(rt,50,30);ctx.stroke(P2);ctx.lineWidth=wR;ctx.strokeStyle=hue(rt,62,56);ctx.stroke(P2);
   if(wR>5){ctx.lineWidth=wR*0.3;ctx.strokeStyle="rgba(255,255,255,.4)";ctx.stroke(P2);}}
  ctx.globalAlpha=1;ctx.lineCap="butt";
  ctx.setLineDash([4,3]);ctx.lineWidth=1.4;ctx.strokeStyle="#e07b00";ctx.beginPath();for(let v=0;v<N;v++){const g=glide[v];if(!g)continue;const a=nodes[g.a],b=nodes[v];ctx.moveTo(sx(a),sy(a));ctx.lineTo(sx(b),sy(b));}ctx.stroke();ctx.setLineDash([]);}
 const nb=new Set();if(hover>=0){nb.add(hover);for(const[j]of outL[hover])nb.add(j);for(const[j]of inL[hover])nb.add(j);}
 const P=new Array(N);
 for(let i=0;i<N;i++){const chain=pred[i]>=0||succ[i]>=0;if(chain&&!PILE)continue;const p=pose(i);P[i]=p;const w=NW(i);if(!inView(p[0],p[1],w))continue;
  ctx.globalAlpha=dim&&!nb.has(i)?0.3:(CELLS&&chain?0.6:1);const col=chain?hue(rootOf[i],62,50):"#5b9bd5";place(p[0],p[1],p[2]);
  if(CELLS){const S=reads[i],x0=-w/2;for(let k=0;k<S.length;k++){ctx.fillStyle=NUC[S[k]]||"#9aa5b1";ctx.fillRect(x0+k*cw,-h/2,cw+0.4,h);}
   ctx.lineWidth=1.5;ctx.strokeStyle=col;ctx.strokeRect(x0,-h/2,w,h);
   if(LETTERS){ctx.fillStyle="#fff";ctx.font="700 "+Math.min(h*0.75,cw*0.95)+"px ui-monospace,Menlo,monospace";ctx.textAlign="center";ctx.textBaseline="middle";for(let k=0;k<S.length;k++)ctx.fillText(S[k],x0+(k+0.5)*cw,0.5);}}
  else{ctx.fillStyle=col;rrect(w,h,Math.min(3,h/2));ctx.fill();}}
 ctx.globalAlpha=1;flat();
 // the latest merge, in place: its two reads on top, the shared letters boxed with the overlap length
 if(CELLS&&applied){const LN=h*1.2+0.5,e=EV[applied-1],a=e[0],b=e[1],pa=P[a]||pose(a),pb=P[b]||pose(b),ov=e[2],bw=ov*cw,wa=NW(a),yo=(LNE[b]-LNE[a])*LN,y0=Math.min(0,yo)-h/2-3,y1=Math.max(0,yo)+h/2+3;
  for(const i of [a,b]){const p=i===a?pa:pb,w=NW(i);if(!inView(p[0],p[1],w))continue;place(p[0],p[1],p[2]);ctx.lineWidth=3;ctx.strokeStyle="#e07b00";ctx.strokeRect(-w/2,-h/2,w,h);}
  place(pa[0],pa[1],pa[2]);ctx.fillStyle="rgba(224,123,0,.22)";ctx.fillRect(wa/2-bw,y0,bw,y1-y0);ctx.lineWidth=2.5;ctx.strokeStyle="#e07b00";ctx.strokeRect(wa/2-bw,y0,bw,y1-y0);flat();
  ctx.textAlign="center";ctx.textBaseline="middle";ctx.font="700 16px system-ui,sans-serif";const lx=wa/2-bw/2,ly=y0-12,c=Math.cos(pa[2]),s=Math.sin(pa[2]),tx=pa[0]+c*lx-s*ly,ty=pa[1]+s*lx+c*ly;ctx.lineWidth=5;ctx.strokeStyle="#f8fafc";ctx.strokeText(ov,tx,ty);ctx.fillStyle="#e07b00";ctx.fillText(ov,tx,ty);}
 // head (filled) and tail (open) of each snake
 if(BODY){const rh=wR*0.7+1.6,rt=wR*0.5+1;
  for(const[st,en] of segs){if(K<1.6&&plen[rootOf[st]]*U*K<16)continue;const hd=nodes[st],tl=nodes[en],hx=sx(hd),hy=sy(hd),tx=sx(tl),ty=sy(tl);
   if(inView(hx,hy,20)){ctx.beginPath();ctx.arc(hx,hy,rh,0,6.2832);ctx.fillStyle="#1b2a41";ctx.fill();ctx.beginPath();ctx.arc(hx,hy,rh*0.38,0,6.2832);ctx.fillStyle="#fff";ctx.fill();}
   if(en!==st&&inView(tx,ty,20)){ctx.beginPath();ctx.arc(tx,ty,rt,0,6.2832);ctx.lineWidth=1.8;ctx.strokeStyle="#1b2a41";ctx.fillStyle="#f8fafc";ctx.fill();ctx.stroke();}}}
 if(applied&&!CELLS){const e=EV[applied-1],a=nodes[e[0]],b=nodes[e[1]];ctx.lineWidth=2.5;ctx.strokeStyle="#e07b00";
  if(PILE){for(const i of [e[0],e[1]]){const p=P[i]||pose(i);place(p[0],p[1],p[2]);rrect(Math.max(NW(i),12)+6,h+6,Math.min(5,(h+6)/2));ctx.stroke();}flat();}
  else{ctx.beginPath();ctx.arc((sx(a)+sx(b))/2,(sy(a)+sy(b))/2,Math.max(7,wR*1.8),0,6.2832);ctx.stroke();}}
 if(hover>=0){const n=nodes[hover];ctx.beginPath();for(const t of [...outL[hover],...inL[hover]]){if(state[t[2]]===0)continue;const o=nodes[t[0]];ctx.moveTo(sx(n),sy(n));ctx.lineTo(sx(o),sy(o));}
  ctx.lineWidth=1.8;ctx.strokeStyle="rgba(27,42,65,.8)";ctx.stroke();ctx.beginPath();ctx.arc(sx(n),sy(n),Math.max(5,wR),0,6.2832);ctx.lineWidth=2.5;ctx.strokeStyle="#1b2a41";ctx.stroke();}
 zm.textContent=(K*U<10?(K*U).toFixed(2):(K*U).toFixed(1))+" px/nt";
}
// ---- sparklines: pieces, longest, overlap
const OVS=EV.map(e=>e[2]),OVMAX=Math.max(...OVS,1),OVMIN=Math.min(...OVS,OVMAX);
const SER=[["#5b9bd5",t=>E?(E-t)/E:0],["#e07b00",t=>E?(MAXL[t]-MAXL[0])/Math.max(1,MAXL[E]-MAXL[0]):0],["#7a4fb5",t=>t?(OVS[t-1]-OVMIN)/Math.max(1,OVMAX-OVMIN):1]];
const CURVE=SER.map(([c,f])=>{const a=[];for(let x=0;x<SW;x++)a.push(f(Math.round(x*E/(SW-1))));return a;});
function drawSpark(){sc.setTransform(dpr,0,0,dpr,0,0);sc.clearRect(0,0,SW,SH);const yy=v=>SH-4-v*(SH-8);
 SER.forEach(([c],s)=>{sc.beginPath();CURVE[s].forEach((v,x)=>x?sc.lineTo(x,yy(v)):sc.moveTo(x,yy(v)));sc.strokeStyle=c;sc.lineWidth=1.6;sc.stroke();});
 const px=E?i0/E*(SW-1):0;sc.strokeStyle="rgba(27,42,65,.35)";sc.lineWidth=1;sc.beginPath();sc.moveTo(px,0);sc.lineTo(px,SH);sc.stroke();
 SER.forEach(([c,f])=>{sc.fillStyle=c;sc.beginPath();sc.arc(px,yy(f(i0)),3.2,0,6.2832);sc.fill();});}
// ---- main loop
function loop(){raf=0;let act=false;
 manV+=(manT-manV)*0.06;if(Math.abs(manT-manV)<0.002)manV=manT;
 if(manV!==manT||manV>0){const f=E?MAXL[i0]/MAXL[E]:1,x=Math.min(1,Math.max(0,(f-0.55)/0.4));unf=0.08*Math.max(x*x*(3-2*x),manV);}
 const simOn=sim.alpha()>sim.alphaMin(),unfOn=unf>0.0005&&still<25;
 if(simOn||unfOn||glideN>0||drag!==null){for(let r=0;r<2;r++){if(simOn)sim.tick();rope();}act=true;dirty=true;still=unf>0.0005&&unfMove<0.01?still+1:0;}
 if(mode){const t=camTarget(),lk=mode===2?0.14:0.12;if(Math.abs(t[0]/K-1)>1e-3||Math.abs(t[1]-TX)>0.05||Math.abs(t[2]-TY)>0.05){K*=Math.pow(t[0]/K,lk);TX+=(t[1]-TX)*lk;TY+=(t[2]-TY)*lk;act=true;dirty=true;}}
 if(dirty){draw();dirty=false;}
 if(act||playing||drag!==null||manV!==manT)raf=requestAnimationFrame(loop);}
function kick(){dirty=true;if(!raf)raf=requestAnimationFrame(loop);}
function reheat(a){if(sim.alpha()<a)sim.alpha(a);still=0;}
// ---- interaction
const fitBtn=$("button[data-a=fit]"),focBtn=$("button[data-a=focus]"),unfBtn=$("button[data-a=unfold]");
function setMode(m){mode=m;fitBtn.classList.toggle("on",m===1);focBtn.classList.toggle("on",m===2);kick();}
const pt=e=>{const r=cv.getBoundingClientRect();return[(e.clientX-r.left)*W/r.width,(e.clientY-r.top)*H/r.height];};
function pick(mx,my){let b=-1,bd=Infinity;const tol=K>=1.6?7:Math.max(7,ribW()/2+4);for(let i=0;i<N;i++){const p=pose(i),chain=pred[i]>=0||succ[i]>=0,d=Math.hypot(p[0]-mx,p[1]-my)-(chain&&K<1.6?0:NW(i)/2);if(d<bd){bd=d;b=i;}}return bd<=tol?b:-1;}
function setTip(i,mx,my){if(i<0){tip.style.display="none";return;}
 const o=outL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]),n=inL[i].filter(t=>state[t[2]]!==0).sort((a,b)=>b[1]-a[1]);
 let s="#"+i;if(o.length)s+="<br>→ #"+o[0][0]+" · "+o[0][1];if(n.length)s+="<br>← #"+n[0][0]+" · "+n[0][1];
 const r=rootOf[i];if(pcnt[r]>1)s+="<br>"+fmt(pcnt[r])+" reads · "+fmt(plen[r])+" nt";
 tip.innerHTML=s;tip.style.display="block";tip.style.left=Math.min(W-170,mx+14)+"px";tip.style.top=Math.min(H-90,my+14)+"px";}
function landGlides(h){let c=h;for(;;){const g=glide[c];if(g){const p=glidePos(g,c,1);nodes[c].x=p[0];nodes[c].y=p[1];glide[c]=null;}if(succ[c]<0)break;c=succ[c];}}
cv.addEventListener("wheel",e=>{e.preventDefault();setMode(0);const[mx,my]=pt(e),nk=Math.max(0.05,Math.min(120,K*Math.exp(-e.deltaY*0.0015)));TX=mx-(mx-TX)*nk/K;TY=my-(my-TY)*nk/K;K=nk;kick();},{passive:false});
cv.addEventListener("mousedown",e=>{const[mx,my]=pt(e),i=pick(mx,my);
 if(i>=0){const chain=pred[i]>=0||succ[i]>=0;if(chain)landGlides(headOf[i]);drag=i;dragInfo={type:!chain?"node":(pred[i]<0?"head":"chain"),mx,my};}   // head: the body follows; body: the whole contig moves
 else{pan=[mx,my,TX,TY];setMode(0);cv.classList.add("drag");}kick();});
window.addEventListener("mousemove",e=>{const[mx,my]=pt(e);
 if(drag!==null){const dx=(mx-dragInfo.mx)/K,dy=(my-dragInfo.my)/K;dragInfo.mx=mx;dragInfo.my=my;
  if(dragInfo.type==="chain"){let c=headOf[drag];for(;;){nodes[c].x+=dx;nodes[c].y+=dy;if(succ[c]<0)break;c=succ[c];}}else{nodes[drag].x+=dx;nodes[drag].y+=dy;nodes[drag].vx=nodes[drag].vy=0;}kick();}
 else if(pan){TX=pan[2]+mx-pan[0];TY=pan[3]+my-pan[1];kick();}
 else if(e.target===cv){const i=pick(mx,my);if(i!==hover){hover=i;setTip(i,mx,my);kick();}else if(i>=0)setTip(i,mx,my);}});
window.addEventListener("mouseup",()=>{drag=null;dragInfo=null;pan=null;cv.classList.remove("drag");kick();});
cv.addEventListener("mouseleave",()=>{if(hover>=0){hover=-1;tip.style.display="none";kick();}});
// ---- stepping
let timer=null,acc=0;
function show(k){const before=i0;i0=Math.max(0,Math.min(E,Math.round(k)));GDUR=playing?Math.max(24,Math.min(100,2500/spd.value)):100;goto(i0);
 const f=E?MAXL[i0]/MAXL[E]:1,x=Math.min(1,Math.max(0,(f-0.55)/0.4));unf=0.08*Math.max(x*x*(3-2*x),manV);still=0;
 $(".s-merge").textContent=fmt(i0)+" / "+fmt(E);$(".s-left").textContent=fmt(N-i0);$(".s-long").textContent=fmt(MAXL[i0]);$(".s-ov").textContent=i0?EV[i0-1][2]:"–";
 rng.value=i0;cnt.textContent="step "+fmt(i0)+" / "+fmt(E);drawSpark();kick();}
function stop(){clearInterval(timer);timer=null;playing=false;pp.textContent="▶ Play";}
function play(){if(i0>=E)show(0);pp.textContent="⏸ Pause";acc=0;playing=true;timer=setInterval(()=>{acc+=spd.value/20;const st=Math.floor(acc);acc-=st;if(i0>=E){stop();return;}show(i0+Math.max(st,0));if(i0>=E)stop();},50);}
const toggle=()=>timer?stop():play();
const refit=()=>setMode(1),refocus=()=>setMode(mode===2?1:2);
function shake(){reheat(0.8);kick();}
function unfoldBtn(){manT=manT?0:1;unfBtn.classList.toggle("on",!!manT);show(i0);}
rng.max=E;rng.oninput=()=>{stop();show(+rng.value);};
spd.onchange=()=>{if(timer){stop();play();}};
root.querySelectorAll("button").forEach(b=>b.onclick=()=>{const a=b.dataset.a;if(a==="pp")return toggle();if(a==="fit")return refit();if(a==="focus")return refocus();if(a==="shake")return shake();if(a==="unfold")return unfoldBtn();stop();show(a==="first"?0:a==="last"?E:i0+(a==="next"?1:-1));root.focus();});
root.addEventListener("keydown",e=>{const m={" ":toggle,ArrowRight:()=>{stop();show(i0+1)},ArrowLeft:()=>{stop();show(i0-1)},PageDown:()=>{stop();show(i0+50)},PageUp:()=>{stop();show(i0-50)},Home:()=>{stop();show(0)},End:()=>{stop();show(E)},f:refit,F:refit,z:refocus,Z:refocus,s:shake,S:shake,u:unfoldBtn,U:unfoldBtn}[e.key];
 if(m){m();e.preventDefault();e.stopPropagation();}},true);
if(window.Reveal&&Reveal.on)Reveal.on("slidechanged",stop);
fitBtn.classList.add("on");
window.__vz={show,nodes,sim,setMode,cam(k,x,y){mode=0;K=k;TX=W/2-k*x;TY=H/2-k*y;kick();},get K(){return K;},get applied(){return applied;},get heads(){return headsN;},get glides(){return glideN;},pick,pose,chainNodes(h){const r=[h];let c=h;while(succ[c]>=0){c=succ[c];r.push(c);}return r;},chains(){const r=[];for(let s=0;s<headsN;s++){let c=heads[s],n=1;while(succ[c]>=0){c=succ[c];n++;}r.push([heads[s],n]);}return r;},get headsArr(){return Array.from(heads.slice(0,headsN));}};
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
    fo, fr, ln, base = [0] * n, [0] * n, [0] * n, 0
    for h in range(n):
        if h in has_pred:
            continue
        c, off, k, lane_end = h, base, 0, []
        while True:
            fo[c], fr[c] = off, k
            lane = next((j for j, e_ in enumerate(lane_end) if e_ <= off), len(lane_end))   # pile-up lane: first free row
            if lane == len(lane_end): lane_end.append(0)
            lane_end[lane] = off + len(reads[c]); ln[c] = lane
            if c not in succ:
                break
            nx, o = succ[c]; off += len(reads[c]) - o; c = nx; k += 1
        base = off + len(reads[c]) + 200
    data = {"n": n, "reads": reads, "lengths": [len(r) for r in reads],
            "ei": [e[0] for e in edges], "ej": [e[1] for e in edges], "eo": [e[2] for e in edges],
            "ev": evo, "maxlen": maxlen, "fo": fo, "fr": fr, "ln": ln}
    out = (TEMPLATE.replace("__ID__", "vz" + uuid.uuid4().hex[:8])
           .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__D3__", d3))
    if a.standalone:
        out = f"<!doctype html><meta charset=utf-8><title>{html.escape(a.fastq)}</title><body style='margin:20px'>{out}"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"{n} reads, {len(edges)} edges shown, {len(ev)} merges -> {a.out} ({len(out)//1024} KB)")

if __name__ == "__main__":
    main()
