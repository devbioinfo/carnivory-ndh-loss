"""Census of complete plastomes in GenBank for carnivorous genera and their non-carnivorous relatives."""
import urllib.request, urllib.parse, json, time, sys, csv
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(u,p):
    p.update({'tool':'plastome_census'}); 
    for t in range(4):
        try:
            return urllib.request.urlopen(E+u+'?'+urllib.parse.urlencode(p),timeout=60).read()
        except Exception as e: time.sleep(2+t*3)
    raise SystemExit('fail '+u)
def search(taxon, lo=100000):
    term=f'{taxon}[Organism] AND (chloroplast[filter] OR plastid[filter]) AND ("complete genome"[Title] OR "complete sequence"[Title]) AND {lo}:250000[SLEN]'
    r=json.loads(get('esearch.fcgi',{'db':'nuccore','term':term,'retmax':2000,'retmode':'json'}))
    ids=r['esearchresult']['idlist']; out=[]
    for i in range(0,len(ids),200):
        s=json.loads(get('esummary.fcgi',{'db':'nuccore','id':','.join(ids[i:i+200]),'retmode':'json'}))
        for k in s['result']['uids']:
            d=s['result'][k]; out.append((d['accessionversion'],d['organism'],d['slen'],d['title'][:90]))
        time.sleep(0.4)
    time.sleep(0.4); return out
groups=[l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
w=csv.writer(open(sys.argv[2],'w',newline='')); w.writerow(['order','group','role','query','accession','organism','length','title'])
for g in groups:
    order,group,role,taxon=g[:4]; lo=int(g[4]) if len(g)>4 else 100000
    res=search(taxon, lo)
    print(f'{order:14s} {role:13s} {taxon:22s} {len(res):4d} records, {len(set(r[1] for r in res)):3d} species', flush=True)
    for r in res: w.writerow([order,group,role,taxon,*r])
