#!/usr/bin/bash

echo "Castor bean"
cat fasta/UP000008311_3988.fasta | grep -v '^>' | wc -c 

echo "Size of Mouse DB"
cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta | grep -v '^>' | wc -c 

echo "Stripped Castor. We truncated the database until it has the same size as mouse."
cat fasta/UP000008311_3988.fasta | head -n -16952 | grep -v '^>' | wc -c

# Replacing accession numbers, so we quickly can assess the origin of a sequence
cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta | \
    gawk 'BEGIN{i=0}/^>/{i=i+1;print ">Mouse_" i; next}{print}' \
    > mouse.fasta


cat fasta/UP000008311_3988.fasta | head -n -16952 | \
    gawk 'BEGIN{i=0}/^>/{i=i+1;print ">Bean_" i; next}{print}' \
    > bean.fasta

# Reverse the sequences
# First a concatination of the two decoys
cat mouse.fasta bean.fasta | \
    gawk 'BEGIN{str=""}(/^>/ && str==""){print;next}/^>/{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x; x=""; str=""; i=i+1; print $0; next}{str=str $0} END{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x}' \
    > decdec.fasta

# Second a concatination of mous with mouse decoys
cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta > mousemousedec.fasta
cat mouse.fasta | \
    gawk 'BEGIN{str=""}(/^>/ && str==""){print;next}/^>/{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x; x=""; str=""; i=i+1; print $0; next}{str=str $0} END{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x}' \
    >> mousemousedec.fasta

# Third a concatination of mouse with bean decoys
cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta > mousebeandec.fasta
cat bean.fasta | \
    gawk 'BEGIN{str=""}(/^>/ && str==""){print;next}/^>/{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x; x=""; str=""; i=i+1; print $0; next}{str=str $0} END{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x}' \
    >> mousebeandec.fasta

echo "concatenated database sizes"
wc -c *.fasta


./crux-4.1.Linux.x86_64/bin/crux tide-index --overwrite T --decoy-format none  decdec.fasta decdec.idx 
./crux-4.1.Linux.x86_64/bin/crux tide-index --overwrite T --decoy-format none  mousemousedec.fasta mousemousedec.idx 
./crux-4.1.Linux.x86_64/bin/crux tide-index --overwrite T --decoy-format none  mousebeandec.fasta mousebeandec.idx 
./crux-4.1.Linux.x86_64/bin/crux tide-search --output-dir decdec-out --precursor-window 10 --precursor-window-type ppm --overwrite T --top-match 1 ./data/150713IA_Inge_Maria_20.mzML decdec.idx
./crux-4.1.Linux.x86_64/bin/crux tide-search --output-dir mousemousedec-out --precursor-window 10 --precursor-window-type ppm --overwrite T --top-match 1 ./data/150713IA_Inge_Maria_20.mzML mousemousedec.idx
./crux-4.1.Linux.x86_64/bin/crux tide-search --output-dir mousebeandec-out --precursor-window 10 --precursor-window-type ppm --overwrite T --top-match 1 ./data/150713IA_Inge_Maria_20.mzML mousebeandec.idx

echo "Top scoring PSMs among decoy-decoy competitors"
cut -f 15 decdec-out/tide-search.target.txt | gawk -F '_' '{print $1}' | sort | uniq -c 

echo "100 highest scoring Top.scoring PSMs among decoy-decoy competitors"
cut -f 9,15 decdec-out/tide-search.target.txt | tail -n +2 | sort -nr | \
    head -n 100 | cut -f 2 | gawk -F '_' '{print $1}' | sort | uniq -c 

echo "Number of decoys among first 10000 PSMs for TDC. First Mouse then Bean."
cut -f 9,15 mousemousedec-out/tide-search.target.txt | tail -n +2 | sort -nr | head -n 10000 | grep -c 'Mouse_'
cut -f 9,15 mousebeandec-out/tide-search.target.txt | tail -n +2 | sort -nr | head -n 10000 | grep -c 'Bean_'
