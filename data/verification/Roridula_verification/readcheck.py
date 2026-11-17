import mappy as mp, sys, re, json, collections
fu=open('fu.fa').read().split('\n',1)[1].replace('\n',''); ch=open('chang.fa').read().split('\n',1)[1].replace('\n','')
REF={'fu':fu,'chang':ch}
AL={k:mp.Aligner(seq=v,preset='sr') for k,v in REF.items()}
LOCI={'fu':{'ndhF':(110610,112835),'ndhD':(115870,117390),'ndhG':(118507,119037),'SSC':(110595,128682),'LSC':(1,84725)},
      'chang':{'ndhF':(128583,130798),'ndhD':(124472,125993),'ndhG':(122806,123336),'SSC_IRexp':(109485,131423),'LSC':(1,84735)}}
JUNC={'fu':{'LSC/IRb':84725,'IRb/SSC':110594,'SSC/IRa':128682},'chang':{'LSC/IRa':84735,'IRa/SSC':128586,'SSC/IRb':131423}}
L={k:len(v) for k,v in REF.items()}
depth={k:[0]*(L[k]+2) for k in REF}; err={k:[0]*(L[k]+2) for k in REF}
jspan={k:collections.Counter() for k in REF}; ident_hist={k:collections.defaultdict(list) for k in REF}
cs_re=re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)')
n=0
for name,seq,q in mp.fastx_read(sys.argv[1]):
    n+=1
    if n>int(sys.argv[3]): break
    for k,a in AL.items():
        for h in a.map(seq,cs=True):
            idt=1-h.NM/max(h.blen,1)
            if idt<0.90 or h.mlen<80: continue
            p=h.r_st
            for op in cs_re.findall(h.cs):
                if op[0]==':':
                    l=int(op[1:])
                    for i in range(p,p+l): depth[k][i]+=1
                    p+=l
                elif op[0]=='*': depth[k][p]+=1; err[k][p]+=1; p+=1
                elif op[0]=='+': err[k][p]+=1
                else:
                    l=len(op)-1
                    for i in range(p,p+l): depth[k][i]+=1; err[k][i]+=1
                    p+=l
            for lab,(a0,b0) in LOCI[k].items():
                if h.r_st<b0 and h.r_en>a0-1: ident_hist[k][lab].append(idt)
            for lab,j in JUNC[k].items():
                if h.r_st<=j-40 and h.r_en>=j+40 and idt>=0.99: jspan[k][lab]+=1
out={'reads':n}
for k in REF:
    for lab,(a0,b0) in LOCI[k].items():
        d=depth[k][a0-1:b0]; e=err[k][a0-1:b0]
        cov=[x for x in d if x>0]
        disc=sum(1 for x,y in zip(d,e) if x>=5 and y/x>0.5)
        ih=ident_hist[k][lab]
        out[f'{k}:{lab}']=dict(mean_depth=round(sum(d)/len(d),1),frac_bp_covered=round(len(cov)/len(d),3),
            bp_where_majority_of_reads_disagree=disc, pct_bp_disagree=round(100*disc/len(d),2),
            reads=len(ih), frac_reads_ge99=round(sum(x>=0.99 for x in ih)/max(len(ih),1),3))
    out[f'{k}:junction_spanning_reads(>=99%,±40bp)']=dict(jspan[k])
print(json.dumps(out,indent=1)); json.dump(out,open(sys.argv[2],'w'),indent=1)
