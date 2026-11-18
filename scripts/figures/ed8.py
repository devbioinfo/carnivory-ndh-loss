import csv,json,numpy as np,matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage,leaves_list
from figstyle import *
setup()
G=['ndhA','ndhB','ndhC','ndhD','ndhE','ndhF','ndhG','ndhH','ndhI','ndhJ','ndhK']
P='/home/claude/work/pan/'
tv={r['tip']:r for r in csv.DictReader(open(P+'phaseC/phaseC_traits_cov311.csv'))}
ev=json.load(open(P+'phaseC/loss_clades.json'))
sc={'F':0,'P':1,'-':2}
rows=[]
for nut,k,p,tl in ev:
    m=np.nanmean([[sc.get(c,np.nan) for c in tv[x]['ndh_calls']] for x in tl],axis=0)
    t0=tl[0].replace('_',' ').split(); name=t0[0] if k>1 else f'{t0[0][0]}. {t0[1]}'
    rows.append((nut,name+(f' ({k})' if k>1 else ''),m,p))
order=['carnivorous','parasitic','mycoheterotrophic','CAM','aquatic','autotrophic']
rows.sort(key=lambda r:(order.index(r[0]),-r[2].mean()))
M=np.array([r[2] for r in rows]); lab=np.array([r[0] for r in rows])
fig=plt.figure(figsize=(W2,235*MM))
gs=fig.add_gridspec(2,2,width_ratios=[1.05,1],height_ratios=[1.7,1],wspace=0.45,hspace=0.3)
# a heatmap
ax=fig.add_subplot(gs[0,0]); cmap=LinearSegmentedColormap.from_list('s',[GS['F'],GS['P'],GS['-']])
ax.imshow(M,cmap=cmap,vmin=0,vmax=2,aspect='auto',interpolation='nearest')
for i,(nut,nm,m,p) in enumerate(rows):
    ax.add_patch(plt.Rectangle((-1.3,i-0.5),0.6,1,color=C.get(nut,INK2) if nut in C else '#bdbcb7',clip_on=False))
ax.set_yticks(range(len(rows)));ax.set_yticklabels([r[1] for r in rows],fontsize=5.3,style='italic')
ax.tick_params(axis='y',pad=1,length=0);ax.set_xlim(-1.35,10.5)
ax.set_xticks(range(11));ax.set_xticklabels(G,rotation=90,fontsize=6,style='italic');ax.xaxis.tick_top()
ax.tick_params(axis='x',length=0,pad=1)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xticks(np.arange(-.5,11),minor=True);ax.set_yticks(np.arange(-.5,len(rows)),minor=True)
ax.grid(which='minor',color='white',lw=0.6);ax.tick_params(which='minor',length=0)
for b in np.where(lab[1:]!=lab[:-1])[0]: ax.axhline(b+0.5,color=INK,lw=0.7)
for v,t in [(0,'Intact'),(1,'Pseudogene'),(2,'Absent')]:
    ax.add_patch(plt.Rectangle((-1.3+v*3.6,len(rows)+0.8),0.7,0.8,color=cmap(v/2),clip_on=False,ec=INK3,lw=0.3))
    ax.text(-0.4+v*3.6,len(rows)+1.2,t,fontsize=5.5,va='center')
for j,(g,t) in enumerate([('carnivorous','Carnivorous'),('parasitic','Parasitic'),('mycoheterotrophic','Mycoheterotr.'),('CAM','Other: CAM, aquatic, autotrophic')]):
    xo=[-1.3,2.3,5.9,-1.3][j]; yo=len(rows)+(2.3 if j<3 else 3.8); ax.add_patch(plt.Rectangle((xo,yo),0.6,0.8,color=C[g],clip_on=False)); ax.text(xo+0.9,yo+0.4,t,fontsize=5.3,va='center')
ax.text(-1.3,len(rows)+6.3,'Cell: mean state across tips of each loss clade\n(posterior ≥ 0.5); n tips in brackets',fontsize=5.2,color=INK2)
ax.set_ylim(len(rows)+8.2,-0.5)
panel(ax,'a',x=-0.42,y=1.08)
# b carn vs para
ax=fig.add_subplot(gs[0,1])
c=M[lab=='carnivorous'].mean(0);p=M[lab=='parasitic'].mean(0)
ax.plot([0.3,1.8],[0.3,1.8],color=GRID,lw=0.8,zorder=0)
ax.scatter(c,p,s=22,color='white',edgecolor=INK,lw=0.7,zorder=3)
off={'ndhK':(4,-6),'ndhC':(4,-6),'ndhB':(4,2),'ndhI':(-14,5),'ndhH':(6,-2),'ndhG':(-4,7),'ndhA':(-17,-2),'ndhD':(-17,3),'ndhF':(4,0),'ndhE':(4,0),'ndhJ':(4,0)}
for g,x,y in zip(G,c,p): ax.annotate(g,(x,y),xytext=off.get(g,(3,3)),textcoords='offset points',fontsize=5.8,style='italic')
r,_=spearmanr(c,p);rng=np.random.default_rng(1)
pp=np.mean([spearmanr(c,rng.permutation(p))[0]>=r for _ in range(10000)])
ax.text(0.04,0.95,f'Spearman $\\rho$ = {r:.2f}\npermutation $P$ = {pp:.3f}\n(10,000 permutations)',transform=ax.transAxes,va='top',fontsize=6)
ax.set_xlabel(f'Mean degradation, carnivorous events (n = {sum(lab=="carnivorous")})')
ax.set_ylabel(f'Mean degradation, parasitic events (n = {sum(lab=="parasitic")})')
ax.set_xlim(0.3,1.8);ax.set_ylim(0.3,1.8)
ax.text(0.97,0.04,'0 = intact, 2 = absent\n$ndhB$ (inverted repeat) degrades last',transform=ax.transAxes,ha='right',fontsize=5.6,color=INK2)
panel(ax,'b')
# c co-loss correlation
ax=fig.add_subplot(gs[1,0])
R=np.corrcoef(M.T);o=leaves_list(linkage(M.T,'average',metric='correlation'))
Ro=R[np.ix_(o,o)];div=LinearSegmentedColormap.from_list('d',['#3e5373','#f3f2ee','#a8714a'])
im=ax.imshow(np.where(np.eye(11,dtype=bool),np.nan,Ro),cmap=div,vmin=-1,vmax=1)
ax.set_xticks(range(11));ax.set_xticklabels([G[i] for i in o],rotation=90,fontsize=5.8,style='italic')
ax.set_yticks(range(11));ax.set_yticklabels([G[i] for i in o],fontsize=5.8,style='italic')
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.04,pad=0.03);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=5.5,length=1.5)
cb.set_label('Pearson r (co-degradation across 27 events)',fontsize=5.8)
panel(ax,'c',x=-0.3)
print('mean offdiag r',np.nanmean(R[~np.eye(11,dtype=bool)]))
# d length vs nfun
ax=fig.add_subplot(gs[1,1])
d=list(csv.DictReader(open(P+'phaseC/len311.csv')))
def grp(r):
    if r['carn']=='1':return 'carnivorous'
    if r['para']=='1':return 'parasitic'
    if r['myco']=='1':return 'mycoheterotrophic'
    return 'autotrophic'
rng=np.random.default_rng(2)
for g in ['autotrophic','mycoheterotrophic','parasitic','carnivorous']:
    s=[r for r in d if grp(r)==g and r['nfun']]
    x=np.array([float(r['nfun']) for r in s])+rng.uniform(-.22,.22,len(s));y=[float(r['len']) for r in s]
    ax.scatter(x,y,s=9 if g!='autotrophic' else 6,color=C[g],edgecolor='white',lw=0.3,alpha=0.9 if g!='autotrophic' else 0.6,label=f'{g.capitalize()} (n = {len(s)})',zorder=3 if g!='autotrophic' else 2)
xx=np.array([0,11]);ax.plot(xx,145.76+1.122*xx,color=INK,lw=1)
ax.text(0.97,0.04,'PGLS ($\\lambda$ = 1.00, n = 310), parasitism and\nmycoheterotrophy as covariates:\n+1.12 kb per intact $ndh$ gene ($P$ = 1.4 × 10$^{-5}$)\ncarnivory × $ndh$ interaction $P$ = 0.30',transform=ax.transAxes,fontsize=5.6,va='bottom',ha='right')
ax.set_xlabel('Intact $ndh$ genes (of 11)');ax.set_ylabel('Plastome length (kb)');ax.set_xticks([0,3,6,9,11])
ax.legend(loc='upper left',bbox_to_anchor=(0.2,1.02),fontsize=5.6,handletextpad=0.2,borderaxespad=0,markerscale=1.3)
ax.set_ylim(0,ax.get_ylim()[1]);panel(ax,'d')
save(fig,'Fig6_loss_order')
