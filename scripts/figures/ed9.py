import csv,re,numpy as np,pandas as pd,matplotlib.pyplot as plt
from figstyle import *
setup();P='/home/claude/work/pan/'
W=pd.read_csv(P+'results/dating_final/results/node_ages_wide.csv').set_index('node')
def parse(s):
    m=re.match(r'([\d.]+) \(([\d.]+)-([\d.]+)\)',s);return tuple(map(float,m.groups()))
# literature (stem/crown ages as published; equivalent nodes only)
lit={'stem Droseraceae':(84.8,'F18'),'crown Droseraceae (X1)':(53.4,'F18'),'Dionaea+Aldrovanda':(48.0,'F18'),
     'crown Nepenthes':(20,'N19'),'stem Triphyophyllum':(6.9,'F18'),'stem Sarraceniaceae':(47,'E12'),
     'crown Sarraceniaceae':(35,'E12'),'stem Roridula (C4)':(39,'G15'),'stem Byblis':(44.5,'F18'),
     'stem Lentibulariaceae':(43.4,'F18')}
nice={'stem Droseraceae':'Droseraceae (stem)','crown Droseraceae (X1)':'Droseraceae (crown)','Dionaea+Aldrovanda':'Dionaea–Aldrovanda',
 'crown Nepenthes':'Nepenthes (crown)','stem Drosophyllum':'Drosophyllum (stem)','stem Triphyophyllum':'Triphyophyllum (stem)',
 'stem Sarraceniaceae':'Sarraceniaceae (stem)','crown Sarraceniaceae':'Sarraceniaceae (crown)','stem Roridula (C4)':'Roridula (stem; C4)',
 'stem Cephalotus':'Cephalotus (stem)','stem Byblis':'Byblis (stem)','stem Lentibulariaceae':'Lentibulariaceae (stem)',
 'crown Lentibulariaceae':'Lentibulariaceae (crown)','crown Utricularia':'Utricularia (crown)','stem Triantha':'Triantha (stem)'}
fig=plt.figure(figsize=(W2,W2*0.74))
gs=fig.add_gridspec(2,2,width_ratios=[1.15,1],wspace=0.5,hspace=0.45)
ax=fig.add_subplot(gs[:,0])
def bands(ax,top):
    geobands(ax,110,labels=False,alpha_f=0.85)
    for x,t in [(66,'K–Pg'),(33.9,'EOT')]: ax.text(x,top,t,ha='center',va='bottom',fontsize=5.6,color=INK2)
nodes=list(nice);n=len(nodes)
for i,k in enumerate(nodes):
    g=parse(W.loc[k,'sens_GBM_r1']);l=parse(W.loc[k,'post_ML_r1'])
    ax.plot([g[1],g[2]],[i-0.13]*2,color=C['carnivorous'],lw=1.6,solid_capstyle='round');ax.plot(g[0],i-0.13,'o',ms=3.6,color=C['carnivorous'],mec='white',mew=0.4)
    ax.plot([l[1],l[2]],[i+0.17]*2,color=INK3,lw=0.9);ax.plot(l[0],i+0.17,'o',ms=3,mfc='white',mec=INK2,mew=0.6)
    if k in lit: ax.plot(lit[k][0],i,'D',ms=3.3,color=INK,zorder=4);ax.text(lit[k][0],i-0.36,lit[k][1],fontsize=4.6,ha='center',color=INK2)
ax.set_yticks(range(n));ax.set_yticklabels([nice[k] for k in nodes],fontsize=6);ax.set_ylim(n+2.0,-1.2)
for t in ax.get_yticklabels(): t.set_style('italic') if not t.get_text().startswith(('Dro','Sarr','Lent')) or 'phyllum' in t.get_text() else None
ax.set_xlim(105,0);ax.set_xlabel('Age (Ma)');bands(ax,-1.15)
ax.plot([],[],'o-',color=C['carnivorous'],ms=3.6,lw=1.6,label='This study, GBM clock (95% HPD)')
ax.plot([],[],'o-',color=INK3,mfc='white',mec=INK2,ms=3,lw=0.9,label='This study, ILN clock')
ax.plot([],[],'D',color=INK,ms=3.3,ls='none',label='Published estimate')
ax.legend(loc='lower left',fontsize=5.6,borderaxespad=0.2)
panel(ax,'a',x=-0.42,y=1.0)
# b ILN vs GBM all nodes
ax=fig.add_subplot(gs[:,1])
xs=[];ys=[]
for k in W.index:
    if k.startswith('root') or 'Brocchinia reducta+hechtioides' in k and k.startswith('crown'):continue
    xs.append(parse(W.loc[k,'post_ML_r1'])[0]);ys.append(parse(W.loc[k,'sens_GBM_r1'])[0])
xs,ys=np.array(xs),np.array(ys)
ax.plot([0,140],[0,140],color=GRID,lw=0.8,zorder=0)
cal=[k for k in W.index if re.search(r'\(C\d\)',k) and not k.startswith('root')]
ax.scatter(xs,ys,s=12,color='white',edgecolor=INK,lw=0.6,zorder=3)
for k in cal:
    x,y=parse(W.loc[k,'post_ML_r1'])[0],parse(W.loc[k,'sens_GBM_r1'])[0];ax.scatter(x,y,s=14,color=INK,zorder=4)
    ax.annotate(re.search(r'C\d',k).group(),(x,y),xytext=(4,-6),textcoords='offset points',fontsize=5.5)
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(2,160);ax.set_ylim(2,160)
ax.set_xlabel('Node age, ILN clock (Ma)');ax.set_ylabel('Node age, GBM clock (Ma)')
ax.text(0.97,0.30,'GBM preferred: log BF = 48.8 ± 0.8\n(stepping-stone, 8 β points)\nYoung uncalibrated nodes older under\nGBM; calibrated nodes stable',transform=ax.transAxes,va='top',ha='right',fontsize=5.6)
ax.text(0.97,0.04,f'n = {len(xs)} nodes; filled = calibrations',transform=ax.transAxes,ha='right',fontsize=5.4,color=INK2)
panel(ax,'b',x=-0.3)
save(fig,'ed/ED9_timescale')
print('ok')
