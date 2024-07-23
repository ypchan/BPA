# step0: prepare the taxa table for your group that can be at the genus, family, order... level by comprehensive literature reviews and Indexfungorum.

# step1: download barcode sequences using the script read.GenBank.R, which require that the taxa table must be prepared follow it's rules.

# step2: multiple sequence alignments using mafft
mkdir 02_mafft
ls 01_data/ | grep fasta | while read a;do echo "mafft --localpair --thread 4 --adjustdirection 01_data/${a} > 02_mafft/${a%.fasta}.mafft.fna";done > mafft.sh
bash mafft.sh
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna

# step3: trim the msa using trimal
mkdir 03_trimal
ls 01_data/ | grep fasta | sed 's/.fasta//' |while read a;do trimal -in 02_mafft/${a}.mafft.fna -gt 0.5 -out 03_trimal/${a%.mafft.fna}.mafft.trimal.fna;done

# step4: find the evolutionary models and concatenate the msa together by the identical identifiers
mkdir 04_modelfinder
outgroup_label=Pyrenochaetopsis_uberifbrmis_CBS_142461
iqtree_modelfinder.py -i 03_trimal/Cucurbitaria_ITS.mafft.trimal.fna \
    03_trimal/Cucurbitaria_LSU.mafft.trimal.fna 03_trimal/Cucurbitaria_TEF1.mafft.trimal.fna 03_trimal/Cucurbitaria_RPB2.mafft.trimal.fna \
    -o 04_modelfinder \
    --mrbayes_nexus \
    --outgroup ${outgroup_label}

# step5: maximum likelihood phylogenetic analysis using iqtree
mkdir 05_iqtree
nohup iqtree2 -s 04_modelfinder/concatenate.fna --seqtype DNA -o ${outgroup_label} --prefix 05_iqtree/iqtree_ml -T AUTO -p 05_modelfinder/best_scheme.txt --ufboot 1000 --alrt 1000 &

# step6: Bayesian analysis using Mrbayes
nohup mpirun -n 4 mb < run_mrbayes.sh # if mpirun does not work, please : nohup bash mb < run_mrbayes &
