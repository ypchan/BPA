mkdir 02_mafft
ls 01_data/ | grep fasta | while read a;do echo "mafft --localpair --thread 4 --adjustdirection 01_data/${a} > 02_mafft/${a%.fasta}.mafft.fna";done > mafft.sh
bash mafft.sh
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna

mkdir 03_trimal
ls 01_data/ | grep fasta | sed 's/.fasta//' |while read a;do trimal -in 02_mafft/${a}.mafft.fna -gt 0.5 -out 03_trimal/${a%.mafft.fna}.mafft.trimal.fna;done

mkdir 04_concatenate
input_items=$(ls 03_trimal/*.mafft.trimal.fna)
ftk_concatenate_msa.py -i ${input_items} -p 04_concatenate/concatenate

mkdir 05_modelfinder
iqtree_modelfinder.py -i ${input_items} -o 05_modelfinder

mkdir 06_iqtree
outgroup_label='Dinemasporium_strigosum_MAFF_244355'
nohup iqtree2 -s 04_concatenate/concatenate.fna --seqtype DNA -o ${outgroup_label} --prefix 06_iqtree/iqtree_ml -T AUTO -p 05_modelfinder/best_scheme.txt --ufboot 1000 --alrt 1000 &
