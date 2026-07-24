// Checkpointed decision solver: does the full N=11 graph contain a clique of size > THR?
// Root-splitting: any clique has a minimum-index vertex i (in fixed order);
// subproblem i searches within {v_i} ∪ (N(v_i) ∩ {v_j : j>i}).  Batches saved to ck_state.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
static int N,V,W,THR;
static uint64_t *adj; static int *vmask; static int found=0;
static double tlim; static clock_t t0; static int to=0;
static inline void setbit(uint64_t*s,int i){s[i>>6]|=1ULL<<(i&63);}
static inline void clrbit(uint64_t*s,int i){s[i>>6]&=~(1ULL<<(i&63));}
static int ap(int m){ if(!m)return 0; if(__builtin_popcount(m)<=2)return 1;
  int p[32],k=0; for(int b=0;b<N;b++) if(m>>b&1)p[k++]=b;
  int d=p[1]-p[0]; for(int i=2;i<k;i++) if(p[i]-p[i-1]!=d) return 0; return 1; }
static int curlen;
static void expand(uint64_t*P,int pc,int rmask){
  if(found||to) return;
  if(((double)(clock()-t0))/CLOCKS_PER_SEC>tlim){to=1;return;}
  // star pruning: if (rmask & AND of all candidate masks) != 0, every completion
  // has a common element, hence size <= t_star(N) <= THR: prune.
  { int andp=(1<<N)-1;
    for(int w=0;w<W;w++){ uint64_t bits=P[w];
      while(bits){ int u=(w<<6)+__builtin_ctzll(bits); bits&=bits-1; andp&=vmask[u];
        if(!(rmask&andp)) goto nostar; } }
    if(rmask&andp) return;
    nostar: ; }
  int *ord=malloc(pc*4),*col=malloc(pc*4);
  uint64_t *un=malloc(W*8),*q=malloc(W*8); memcpy(un,P,W*8);
  int cnt=0,color=0;
  while(1){ int any=0; for(int w=0;w<W;w++) if(un[w]){any=1;break;} if(!any)break;
    color++; memcpy(q,un,W*8);
    while(1){ int v=-1; for(int w=0;w<W;w++) if(q[w]){v=(w<<6)+__builtin_ctzll(q[w]);break;}
      if(v<0)break; clrbit(q,v); clrbit(un,v); ord[cnt]=v; col[cnt]=color; cnt++;
      uint64_t*av=adj+(size_t)v*W; for(int w=0;w<W;w++) q[w]&=~av[w]; } }
  free(q); free(un);
  uint64_t *NP=malloc(W*8);
  for(int i=cnt-1;i>=0&&!found&&!to;i--){
    if(curlen+col[i]<=THR) break;
    int v=ord[i]; curlen++;
    uint64_t*av=adj+(size_t)v*W; int e=1,npc=0;
    for(int w=0;w<W;w++){NP[w]=P[w]&av[w]; if(NP[w])e=0;}
    if(e){ if(curlen>THR){found=curlen;} }
    else { for(int w=0;w<W;w++) npc+=__builtin_popcountll(NP[w]); expand(NP,npc,rmask&vmask[v]); }
    curlen--; clrbit(P,v);
  }
  free(NP); free(ord); free(col);
}
int main(int argc,char**argv){
  N=atoi(argv[1]); THR=atoi(argv[2]); int i0=atoi(argv[3]), i1=atoi(argv[4]); tlim=atof(argv[5]);
  int M=1<<N; V=M-1; W=(V+63)>>6;
  unsigned char *ia=malloc(M); for(int m=1;m<M;m++) ia[m]=ap(m);
  vmask=malloc(V*4); int*deg=calloc(V,4);
  for(int m=1;m<M;m++) vmask[m-1]=m;
  adj=calloc((size_t)V*W,8);
  for(int i=0;i<V;i++)for(int j=i+1;j<V;j++){int mm=vmask[i]&vmask[j];
    if(mm&&ia[mm]){setbit(adj+(size_t)i*W,j);setbit(adj+(size_t)j*W,i);deg[i]++;deg[j]++;}}
  // 58-core peeling: clique of size THR+1 needs internal degree >= THR
  int changed=1; unsigned char*alive=malloc(V); memset(alive,1,V);
  while(changed){ changed=0;
    for(int i=0;i<V;i++) if(alive[i]&&deg[i]<THR){ alive[i]=0; changed=1;
      uint64_t*ai=adj+(size_t)i*W;
      for(int j=0;j<V;j++) if((ai[j>>6]>>(j&63))&1){ deg[j]--; clrbit(adj+(size_t)j*W,i);} 
      memset(ai,0,W*8); deg[i]=0; } }
  int na=0; for(int i=0;i<V;i++) na+=alive[i];
  fprintf(stderr,"alive after %d-core peel: %d of %d\n",THR,na,V);
  // order alive vertices by degree ascending (small subproblems first)
  int *idx=malloc(V*4),nn=0;
  for(int i=0;i<V;i++) if(alive[i]) idx[nn++]=i;
  for(int a=0;a<nn;a++)for(int b=a+1;b<nn;b++) if(deg[idx[b]]<deg[idx[a]]){int t=idx[a];idx[a]=idx[b];idx[b]=t;}
  t0=clock();
  uint64_t *P=malloc(W*8);
  int done_to=i0;
  for(int a=i0;a<i1&&a<nn;a++){
    int v=idx[a];
    // subproblem: cliques with min-order-index a: P = N(v) ∩ {idx[b]: b>a}
    memset(P,0,W*8);
    uint64_t*av=adj+(size_t)v*W;
    for(int b=a+1;b<nn;b++){int u=idx[b]; if((av[u>>6]>>(u&63))&1) setbit(P,u);}
    // subproblem peeling: iterate removing u in P with |N(u)∩P| < THR-1
    int ch=1;
    while(ch){ ch=0;
      for(int w=0;w<W;w++){ uint64_t bits=P[w];
        while(bits){ int u=(w<<6)+__builtin_ctzll(bits); bits&=bits-1;
          uint64_t*au=adj+(size_t)u*W; int du=0;
          for(int ww=0;ww<W;ww++) du+=__builtin_popcountll(au[ww]&P[ww]);
          if(du<THR-1){ clrbit(P,u); ch=1; } } } }
    curlen=1; int pc=0; for(int w=0;w<W;w++) pc+=__builtin_popcountll(P[w]);
    if(pc>=THR) expand(P,pc,vmask[v]);
    if(found){ printf("FOUND clique of size %d (subproblem %d)\n",found,a); return 0; }
    if(to){ printf("TIMEOUT at subproblem %d (completed [%d,%d))\n",a,i0,a); return 2; }
    done_to=a+1;
  }
  printf("NO clique > %d with min-index in [%d,%d)  (n_alive=%d)\n",THR,i0,done_to,nn);
  return 0;
}
