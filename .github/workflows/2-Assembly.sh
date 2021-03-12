
# Assembly of Thalassiosirales transcriptomes using the Trinity and rnaSPAdes

# Requirements
# Khmer V2.1.1, Trinity V2.12.0,rnaSPAdes V3.13.0

# Files with extension pe.qc.keep.abundfilt.fq.gz generated in the previous step must be in working directory
#

# first we need to split pe.qc.keep.abundfilt.fq.gz files 

for file in *.pe.qc.keep.abundfilt.fq.gz
do
   split-paired-reads.py ${file}
done

for file in *.1
do
file2=${file/.1/.2}
fileout=${file/.pe.qc.keep.abundfilt.fq.gz.1/}
cat ${file} > ${fileout}.left.fq
cat ${file2} > ${fileout}.right.fq
done

for file in *.se.qc.keep.abundfilt.fq.gz
do
fileout=${file/.se.qc.keep.abundfilt.fq.gz/.left.fq}
gunzip -c ${file} >> ${fileout}
done

# Assembly with Trinity 

for fileleft in *.left.fq
do
fileright=${fileleft/.left.fq/.right.fq}
out=${fileleft/.left.fq/}
Trinity --left ${fileleft} \
  --right ${fileright} --seqType fq --max_memory 12G \
  --CPU 10 --no_bowtie --output ${out%%}-Trinity
done

# Assembly with rnaSPAdes

for fileleft in *.left.fq
do
fileright=${fileleft/.left.fq/.right.fq}
out=${fileleft/.left.fq/}
rnaspades.py -1 ${fileleft} -2 ${fileright} -o ${out%%}-SPADES
done
