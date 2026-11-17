import csv,json,time,urllib.request,urllib.parse
E_="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(4):
        try: return urllib.request.urlopen(u,timeout=60).read().decode()
        except Exception: time.sleep(2)
    return ''
out=open('ath.fa','w'); log=open('ath_check.tsv','w')
for cls,g,loc in csv.reader(open('genes.tsv'),delimiter='\t'):
    q=urllib.parse.quote(f'{loc}[All Fields] AND txid3702[orgn] AND alive[prop]')
    ids_=json.loads(get(E_+f'esearch.fcgi?db=gene&retmode=json&retmax=10&term={q}'))['esearchresult']['idlist']
    gid=None
    if ids_:
        R=json.loads(get(E_+'esummary.fcgi?db=gene&retmode=json&id='+','.join(ids_)))['result']
        for u in ids_:
            if loc in (R[u].get('otheraliases','')+' '+R[u]['name']).replace(' ','').split(','): gid=[u]; s=R[u]; break
    if not gid: log.write(f'{g}\t{loc}\tNO_EXACT_GENE\n'); continue
    lk=json.loads(get(E_+f'elink.fcgi?dbfrom=gene&db=protein&retmode=json&linkname=gene_protein_refseq&id={gid[0]}'))
    pids=[l for ls in lk['linksets'][0].get('linksetdbs',[]) for l in ls['links']]
    seqs=[]
    for p in pids[:6]:
        fa=get(E_+f'efetch.fcgi?db=protein&rettype=fasta&id={p}')
        sq=''.join(fa.split('\n')[1:]); seqs.append(sq); time.sleep(0.3)
    best=max(seqs,key=len) if seqs else ''
    out.write(f'>{g}|{loc}\n{best}\n')
    log.write(f"{g}\t{loc}\t{s['name']}\t{s.get('otheraliases','')[:60]}\t{s['description'][:60]}\t{len(best)}\n")
    time.sleep(0.3)
