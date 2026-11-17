import sys,subprocess,collections; sys.path.insert(0,'/home/claude/work/pan')
from genbank_utils import read_genbank, extract
import pangene as P
from Bio.Seq import Seq
fu=next(read_genbank('/home/claude/work/pan/genbank_retired/NC_069066.1.gb')); ch=next(read_genbank('/home/claude/work/pan/genbank/PP273159.1.gb'))
rc=lambda s:str(Seq(s).reverse_complement())
def genes(r):
    d={}
    for f in r['features']:
        if f['type'] in ('CDS','gene'):
            g=P.gname(f['qualifiers'])
            if g and (g not in d or (d[g][3]=='gene' and f['type']=='CDS')): d[g]=(extract(r['seq'],f['location']).upper(),f['location'],'pseudo' in f['qualifiers'],f['type'])
    return d
gf,gc=genes(fu),genes(ch)
def ident(a,b):
    open('/tmp/a.fa','w').write('>a\n'+a+'\n'); open('/tmp/b.fa','w').write('>b\n'+b+'\n')
    out=subprocess.run(['mafft','--quiet','--adjustdirection','--auto','/tmp/ab.fa'],capture_output=True,text=True) if False else None
    open('/tmp/ab.fa','w').write(f'>a\n{a}\n>b\n{b}\n')
    o=subprocess.run(['mafft','--quiet','--adjustdirection','--globalpair','/tmp/ab.fa'],capture_output=True,text=True).stdout
    seqs=[''.join(x.split('\n')[1:]).upper() for x in o.split('>')[1:]]
    A,B=seqs; sub=ind=same=0; gap_open=0; prevgap=False
    for x,y in zip(A,B):
        if x=='-' or y=='-':
            ind+=1; gap_open+= (not prevgap); prevgap=True
        else:
            prevgap=False
            if x==y: same+=1
            else: sub+=1
    return same/(same+sub), sub, gap_open, ind, '_R_' in o
SSCg=['ndhF','rpl32','ccsA','ndhD','psaC','ndhE','ndhG','ndhI','ndhA','ndhH','rps15','ycf1']
LSCg=['rbcL','matK','psbA','atpA','rpoB','rpoC2','accD','psaA']
print('gene\tregion\tFu_len\tChang_len\tFu_pseudo\tChang_pseudo\tidentity\tsubstitutions\tindel_events\tindel_bp\tChang_revcomp_needed')
for reg,gl in (('SSC',SSCg),('LSC',LSCg)):
    for g in gl:
        if g in gf and g in gc:
            i,s,go,ib,flip=ident(gf[g][0],gc[g][0])
            print(f'{g}\t{reg}\t{len(gf[g][0])}\t{len(gc[g][0])}\t{gf[g][2]}\t{gc[g][2]}\t{i*100:.2f}\t{s}\t{go}\t{ib}\t{flip}')
        else: print(g,reg,'missing in', 'Fu' if g not in gf else 'Chang')

fssc=fu['seq'].upper()[110594:128682]; open('/tmp/q.fa','w').write('>q\n'+fssc+'\n'); open('/tmp/s.fa','w').write('>s\n'+ch['seq'].upper()+'\n')
for strand in ('plus','minus'):
    o=subprocess.run(['blastn','-query','/tmp/q.fa','-subject','/tmp/s.fa','-strand',strand,'-outfmt','6 qstart qend sstart send pident length mismatch gaps','-evalue','1e-20'],capture_output=True,text=True).stdout
    hs=sorted([l.split() for l in o.splitlines()],key=lambda h:-int(h[5]))
    cov=set(); mm=gp=al=0
    for h in hs:
        new=set(range(int(h[0]),int(h[1])+1))-cov
        if len(new)<0.5*int(h[5]): continue
        cov|=new; mm+=int(h[6]); gp+=int(h[7]); al+=int(h[5])
    print(f'Fu SSC (18,088 bp) vs Chang genome, strand {strand}: query covered {len(cov)/len(fssc)*100:.1f}%, identity of aligned bp {(al-mm-gp)/max(al,1)*100:.2f}%')
