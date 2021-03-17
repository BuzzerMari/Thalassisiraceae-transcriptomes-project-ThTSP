# Estimation of abundance of Thalassiosirales transcriptomes using Salmon

# Requirements
# Salmon 1.3.0

# First we need to index the Trinity de novo assembly output 

# Salmon
mkdir quants
for file in *-Trinity.fasta
do
index=${file/-Trinity.fasta/_index}
left={file/-Trinity.fasta/.left.fq}
right={file/-Trinity.fasta/.right.fq}
salmon index -t ${file} -i ${index}
salmon quant -i ${index} -l A \
         -1  ${left} \
         -2  ${right} \
         -p 8 --validateMappings -o quants/${file%%}_quant
done


# next step
ls quants/*quant.sf > salmon.sf.file.list
trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method salmon \
    --gene_trans_map  ../trinotate/transcripts.gene2tr.map \
    --out_prefix salmon \
    --name_sample_by_basedir \
    --quant_files salmon.sf.file.list
