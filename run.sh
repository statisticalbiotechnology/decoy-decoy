#!/usr/bin/bash

echo "Castor bean"
cat fasta/UP000008311_3988.fasta | grep -v '^>' | wc -c 

echo "Mouse"
cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta | grep -v '^>' | wc -c 

echo "Stripped Castor"
cat fasta/UP000008311_3988.fasta | head -n -16952 | grep -v '^>' | wc -c

cat fasta/uniprot-proteome_UP000000589+reviewed_yes.fasta | \
    gawk 'BEGIN{i=0}/^>/{i=i+1;print ">Mouse_" i; next}{print}' \
    > mousebean.fasta


cat fasta/UP000008311_3988.fasta | head -n -16952 | \
    gawk 'BEGIN{i=0}/^>/{i=i+1;print ">Bean_" i; next}{print}' \
    >> mousebean.fasta


cat mousebean.fasta | \
    gawk 'BEGIN{str=""}(/^>/ && str==""){print;next}/^>/{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x; x=""; str=""; i=i+1; print $0; next}{str=str $0} END{for(j=length(str);j!=0;j--){x=x substr(str,j,1)}; print x}' \
    > decdec.fasta

echo "concatenated database sizes"
cat mousebean.fasta \
    | grep -v '^>' | wc -c

cat decdec.fasta \
    | grep -v '^>' | wc -c


./crux-4.1.Linux.x86_64/bin/crux tide-index --decoy-format none  decdec.fasta decdec.idx 
./crux-4.1.Linux.x86_64/bin/crux tide-search --precursor-window 10 --precursor-window-type ppm --overwrite T --top-match 1 ./data/150713IA_Inge_Maria_20.mzML decdec.idx

cut -f 15 crux-output/tide-search.target.txt | gawk -F '_' '{print $1}' | sort | uniq -c 
