import sys,subprocess,json; sys.path.insert(0,'/home/claude/work/pan')
from genbank_utils import read_genbank, extract
import pangene as P
from Bio.Seq import Seq
fu=next(read_genbank('/home/claude/work/pan/genbank_retired/NC_069066.1.gb')); chs=open('chang.fa').read().split('\n',1)[1].replace('\n','')
act=[l for l in open('/home/claude/work/pan/taxa_selected.csv') if 'Actinidia' in l][0].split(',')[3]
ar=next(read_genbank(f'/home/claude/work/pan/genbank/{act}.gb'))
loc={'ndhF':(128583,130798,'-'),'ndhD':(124472,125993,'+'),'ndhG':(122806,123336,'+')}
res={}
for g,(a,b,strand) in loc.items():
    fcds=[extract(fu['seq'],f['location']).upper() for f in fu['features'] if f['type']=='CDS' and P.gname(f['qualifiers'])==g][0]
    acds=[extract(ar['seq'],f['location']).upper() for f in ar['features'] if f['type']=='CDS' and P.gname(f['qualifiers'])==g][0]
    pad=60; s=chs[a-1-pad:b+pad]
    if strand=='-': s=str(Seq(s).reverse_complement())
    # codon-aware: align Chang region to Fu CDS (nucleotide, global ends free) with mafft
    open('/tmp/x.fa','w').write(f'>fu\n{fcds}\n>act\n{acds}\n>chang\n{s}\n')
    o=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','100','/tmp/x.fa'],capture_output=True,text=True).stdout
    al={x.split('\n')[0]:''.join(x.split('\n')[1:]).upper() for x in o.split('>')[1:]}
    F,A,C=al['fu'],al['act'],al['chang']
    cols=[i for i,ch in enumerate(F) if ch!='-']
    st,en=cols[0],cols[-1]
    Csub=C[st:en+1]; Fsub=F[st:en+1]
    # indels relative to Fu within CDS span
    indels=[]; i=0; fpos=0
    while i<len(Fsub):
        if Fsub[i]=='-' or Csub[i]=='-':
            j=i; kind='ins' if Fsub[i]=='-' else 'del'
            while j<len(Fsub) and ((Fsub[j]=='-')==(kind=='ins')) and (Fsub[j]=='-' or Csub[j]=='-'): j+=1
            indels.append((fpos,kind,j-i)); 
            if kind=='del': fpos+=j-i
            i=j
        else: fpos+=1; i+=1
    cseq=Csub.replace('-','')
    prot=str(Seq(cseq[:len(cseq)//3*3]).translate(table=11))
    fprot=str(Seq(fcds).translate(table=11))
    first_stop=prot.find('*'); 
    fs=[x for x in indels if x[2]%3]
    res[g]=dict(fu_cds_len=len(fcds),fu_internal_stops=fprot[:-1].count('*'),chang_homolog_len=len(cseq),
        chang_frameshift_indels=len(fs),chang_inframe_indels=len(indels)-len(fs),
        chang_first_stop_codon_aa=first_stop, chang_orf_fraction_of_fu=round(first_stop/(len(fcds)/3),2) if first_stop>=0 else None,
        frameshifts=[(p,k,l) for p,k,l in fs][:6], strand=strand, chang_coords=f'{a}-{b}')
    print(g,json.dumps(res[g]))
json.dump(res,open('ndh_orf.json','w'),indent=1)
