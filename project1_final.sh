###run code with:
###nohup ./project1_final.sh > project1_test_1.log &
#!/bin/bash
clear
echo "Declare variables for file path"
REF_DIR=/projects/micb405/resources/project_1
PYTHON_DIR=/projects/micb405/resources

mkdir 0_reference
mkdir 1_bam 
mkdir 2_vcf
mkdir 2_vcf_clean
mkdir 3_aln
mkdir 4_tree

echo "Copy reference genome to 0_reference"
cp $REF_DIR/ref_genome.fasta 0_reference
 
echo "Indexing the reference, only required for the first time"
bwa index 0_reference/ref_genome.fasta

echo "Align reads from fastq file to the reference"
echo "LOOP THROUGH ALL THE FILE ENDS WITH "_1.fastq", generate sam file with read 1 and read 2, converted to bam"

for file in ${REF_DIR}/*_1.fastq.gz
do
#STORE FILE NAME WITHOUT SUFFIX, REMOVED THE NAME OF FILE FROM LAST CHARACTER TO THE LAST _, _ INCLUDED
filename=${file%_*} 
#REMOVED THE NAME OF FILE FROM FIRST CHARACTER TO THE LAST /, /INCLUDED 
bwa mem -t 2 0_reference/ref_genome.fasta ${filename}_1.fastq.gz ${filename}_2.fastq.gz | samtools view -@ 2 -b - 1> 1_bam/${filename##*/}.bam
done

echo “Sort BAM & Remove duplicates”
for file in 1_bam/*.bam
do
samtools sort -@ 2 $file | samtools rmdup - ${file%%.*}.sorted.rmdup.bam
done

echo “Index removed duplicates and sorted BAM” 
#INDEX THE SORTED, DUPLICATES REMOVED BAM FILE
for file in 1_bam/*.sorted.rmdup.bam
do
samtools index $file 
done

#BCFTOOLS MPILEUP
echo “Pileup to sum mapping of BAM on reference” 
echo “-q 	skip alignments with mapQ smaller than 30” 
echo “-f 	faidx indexed reference sequence file” 
echo “-O u 	output uncompressed VCF/BCF” 
echo “-I 	indel only” 
#BCFTOOLS CALL
echo “Convert pileup to human-readable format,”
echo “-O v variants only”
echo “-mv ”
for file in 1_bam/*.sorted.rmdup.bam
do
filename=${file##*/} 
bcftools mpileup --threads 10 -q 30 -O u -f 0_reference/ref_genome.fasta $file -I | bcftools call -O v -mv - > 2_vcf/${filename%%.*}.vcf 
done


echo "Variant convert to fasta by the python tool"
echo "vcf_x, exclude snp similar in all variants, since we are comparing to reference0 as a patient, we should use vcf, not vcf_x"
python $PYTHON_DIR/vcf_to_fasta_het.py -x 2_vcf/ project1_vcf

mv 2_vcf/project1_vcf* 3_aln
 
echo "Multiple sequence align by MAFFT, don't do----"
# mafft --auto project1_vcf.fasta > project1_vcf.mafft.mfa

echo "FastTree phylogeny"
FastTree 3_aln/project1_vcf.fasta > 4_tree/project1.nwk


echo "Cleaned up vcf, python - vcf to fasta, fasta to tree"
for file in 2_vcf/*.vcf
do 
filename=${file##*/} 
bcftools filter -e "QUAL < 200" $file > 2_vcf_clean/${filename%%.*}.filtered.vcf
done
python $PYTHON_DIR/vcf_to_fasta_het.py -x 2_vcf_clean/ project1_clean_vcf 
mv 2_vcf_clean/project1_clean_vcf* 3_aln
FastTree 3_aln/project1_clean_vcf.fasta > 4_tree/project1_clean.nwk

echo "DONE"
