import csv, os, time, urllib.request, urllib.parse
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
accs=[r['accession'] for r in csv.DictReader(open('taxa_selected.csv'))]
todo=[a for a in accs if not os.path.exists(f'genbank/{a}.gb')]
for i in range(0,len(todo),20):
    batch=todo[i:i+20]
    p=urllib.parse.urlencode({'db':'nuccore','id':','.join(batch),'rettype':'gbwithparts','retmode':'text','tool':'plastome_pan'})
    for t in range(5):
        try: txt=urllib.request.urlopen(E+'?'+p,timeout=180).read().decode(); break
        except Exception as e: print('retry',e); time.sleep(5+5*t)
    recs=[r+'//\n' for r in txt.split('//\n') if r.strip()]
    for r in recs:
        ver=r.split('VERSION')[1].split()[0]
        open(f'genbank/{ver}.gb','w').write(r)
    print(i+len(batch),'/',len(todo),flush=True); time.sleep(0.5)
missing=[a for a in accs if not os.path.exists(f'genbank/{a}.gb')]
print('missing',missing)
