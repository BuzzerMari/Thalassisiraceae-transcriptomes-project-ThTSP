# Assembly evaluation of Thalassiosirales transcriptomes using BUSCO and Transrate

# Requirements
# BUSCO 5.0.0, Transrate 1.0.3

# Assembly files generated in the previous step must be in working directory
# Assemblies were moved to a new directory according to used strategy
# for Trinity assemblies go to ~/Trinity_assemblies
# for rnaSPAdes assemblies go to ~/Spades_assemblies
 


# In the NT-server first : conda activate busco

cd ~/Trinity_assemblies

#busco
for file in *-Trinity.fasta
do
out=${file/-Trinity.fasta/}
busco -i ${file}  --auto-lineage-euk -o BUSCO-${out} -m tran
done

# for transrate que need the filtered reads generated (In this case are in the same working directory)

# transrate
for file in *-Trinity.fasta
do
left={file/-Trinity.fasta/.left.fq.paired.fq}
right={file/-Trinity.fasta/.right.fq.paired.fq}
out=${file/-Trinity.fasta/}
transrate --assembly ${file} --left ${left} --right ${right} --output Transrate-${out} --threads 8 
done

cd ../Spades_assemblies


for file in *-Spades.fasta
do
out=${file/.fasta/}
busco -i ${file}  --auto-lineage-euk -o BUSCO-${out} -m tran
done

# for transrate que need the filtered reads generated (In this case are in the same working directory)

# transrate
for file in *-Spades.fasta
do
left={file/-Spades.fasta/.left.fq.paired.fq}
right={file/-Spades.fasta/.right.fq.paired.fq}
out=${file/.fasta/}
transrate --assembly ${file} --left ${left} --right ${right} --output Transrate-${out} --threads 8 
done
