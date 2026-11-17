cd /home/claude/work/pan/nuclear
while read sp acc; do
  [ -s gchk/${sp}_calls.tsv ] && continue
  if [ ! -s genomes/$sp.fna ]; then curl -s -o genomes/$sp.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" && unzip -p genomes/$sp.zip '*.fna' > genomes/$sp.fna; rm -f genomes/$sp.zip; fi
  python3 genome_check.py $sp genomes/$sp.fna && rm -f genomes/$sp.fna
done <<L
Drosophyllum_lusitanicum GCA_032274825.1
Darlingtonia_californica GCA_032270445.1
Pinguicula_primuliflora GCA_025802335.1
Fagopyrum_tataricum GCA_002319775.1
Triantha_glutinosa GCA_051380605.1
Tofieldia_thibetica GCA_056747465.1
Cephalotus_follicularis GCA_001972305.1
Dionaea_muscipula GCA_054771495.1
L
