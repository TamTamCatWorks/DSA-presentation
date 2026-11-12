"""Efficient greedy overlap-layout: index prefixes, collect all suffix->prefix overlaps >= MIN,
then take edges longest-first (classic greedy), skipping edges that would give a read two
successors/predecessors or close a cycle. Emits an event log of merges."""
import collections, json, sys

def read_fastq(path):
    return [l.strip() for i,l in enumerate(open(path)) if i%4==1]

def overlaps(rs, m):
    idx=collections.defaultdict(list)
    for j,r in enumerate(rs): idx[r[:m]].append(j)
    out=[]
    for i,a in enumerate(rs):
        for s in range(1,len(a)-m+1):
            for j in idx.get(a[s:s+m],()):
                if j!=i and rs[j].startswith(a[s:]): out.append((len(a)-s,i,j))
    return out

def greedy_events(rs, m=20):
    E=overlaps(rs,m); E.sort(key=lambda e:(-e[0],e[1],e[2]))
    n=len(rs); succ=[None]*n; pred=[None]*n
    parent=list(range(n))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    length={i:len(rs[i]) for i in range(n)}   # keyed by union-find root
    events=[]
    for o,a,b in E:
        if succ[a] is not None or pred[b] is not None: continue
        ra,rb=find(a),find(b)
        if ra==rb: continue
        succ[a]=(b,o); pred[b]=a
        parent[rb]=ra
        length[ra]=length[ra]+length[rb]-o
        events.append({"a":a,"b":b,"ov":o,"len":length[ra],"root":ra,"rb":rb})
    # build contigs
    heads=[i for i in range(n) if pred[i] is None]
    contigs=[]
    for h in heads:
        s=rs[h]; cur=h; k=1
        while succ[cur] is not None:
            nxt,o=succ[cur]; s+=rs[nxt][o:]; cur=nxt; k+=1
        contigs.append((s,k,h))
    return events, contigs, len(E)

if __name__=="__main__":
    m=int(sys.argv[1]) if len(sys.argv)>1 else 20
    rs=list(dict.fromkeys(read_fastq('/mnt/user-data/uploads/ads1_week4_reads.fq')))
    ev,ct,ne=greedy_events(rs,m)
    L=sorted((len(c[0]) for c in ct),reverse=True)
    tot=sum(L); acc=0
    for x in L:
        acc+=x
        if acc>=tot/2: n50=x;break
    print(f"min_overlap={m} edges={ne} merges={len(ev)} contigs={len(ct)} lengths={L[:10]} total={tot} N50={n50}")
    ov=[e['ov'] for e in ev]
    print("overlap of last 10 merges:",ov[-10:], " first 5:",ov[:5])
