# This is a basic workflow to help you get started with Actions
# pre-processing of Thalassiosirales transcriptomes
# Khmer and trimmomatic should be installed
# fastq files must be in working directory & I have placed the adapter folder in the same directory otherwise change the path below

#Parameters; mismatches=2; palindrome clip threshold= 30; simple clip threshold=10

# trimmomatic

for file in ./*R1*; do
    file2=${file/R1/R2}

    trimmomatic PE ${file} ${file2} ${file%%}_pe ${file%%}_se ${file2%%}_pe ${file2%%}_se \
    ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10
done

mkdir trimmed

# interleave of paired-end output files
for file in ./*R1*_pe; do
    file2=${file/R1/R2}
    outpe=${file/R1.fastq.gz_pe/}

interleave-reads.py ${file} ${file2} | gzip -9c \
    > trimmed/${outpe}.pe.fq.gz

done

#  merge of single-ended files
for file in ../*R1*_se; do
    file2=${file/R1/R2}
    outse=${file/_R1.fastq.gz_se/}
    cat ${file} ${file2}| gzip -9c > trimmed/${outse}.se.fq.gz
done


cd trimmed
chmod u-w *

mkdir Filter

#  quality filter from interleaved sequences
for file in *.fq.gz
do
    
     filterfile=${file%%.fq.gz}.qc.fq.gz
     gunzip -c ${file} | fastq_quality_filter -Q33 -q 25 -p 50 | gzip -9c \
         > Filter/${filterfile}
done

cd Filter
#  extract paired end 

mkdir filtPairs
cd filtPairs
for file in ../*.pe.qc.fq.gz
do
   extract-paired-reads.py ${file}
done

#  simplify names orphan files
for file in *.se
do
  file2=${file%%.pe.qc.fq.gz.se}.se.qc.fq.gz 
  gunzip -c ${file2} > combine
  cat ${file} >> combine
  gzip -c combine > ${file2} # combine single reads together
  rm ${file} combine
done

#  simplify names paired-end files
for file in *.pe
do
   file2=${file%%.pe.qc.fq.gz.pe}.pe.qc.fq
   mv $file $file2
   gzip $file2
done

# here*** correct For optimization now we run the digital normalization from khmer

mkdir digital_normalization
cd digital_normalization

#  First for paired reads

normalize-by-median.py --paired -ksize 20 --cutoff 25 -n_tables 4 \
  --min-tablesize 4e9 --savetable normC25k20.ct \
  ../*.pe.qc.fq.gz

normalize-by-median.py --cutoff 25 --loadtable normC25k20.ct \
  --savetable normC25k20.ct ../*.se.qc.fq.gz


mkdir abundfilt
cd abundfilt
filter-abund.py --variable-coverage ../diginorm/normC20k20.ct \
  --threads ${THREADS:-1} ../diginorm/*.keep

