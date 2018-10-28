###run code with:
###nohup ./project1_final.sh > project1_test_1.log &
#!/bin/bash
clear

echo `date`
echo "MICB 405 2018W1 Project 1 - Group 4"
echo "Run BASIC  Method & MQUAL200 Method:"

echo "Declare variables for file path & make directories."
REF_DIR=/projects/micb405/resources/project_1
PYTHON_DIR=/projects/micb405/resources

mkdir 0_reference
mkdir 1_bam 
mkdir 2_vcf
mkdir 2_vcf_clean
mkdir 3_aln
mkdir 4_tree

echo "Copy reference genome to dir 0_reference."
cp $REF_DIR/ref_genome.fasta 0_reference
 
echo "BWA INDEX: Indexing the reference, only required for the first time."
bwa index 0_reference/ref_genome.fasta

echo "BWA MEM, SAMTOOLS VIEW, SAMTOOLS SORT, SAMTOOLS RMDUP:"
echo "Align reads from fastq file to the reference, generate sam file with read 1 and read 2, converted to bam, sort BAM, and remove duplicates in BAM."

for file in ${REF_DIR}/*_1.fastq.gz
do
	#REMOVED THE NAME OF FILE FROM LAST CHARACTER TO THE LAST _, _ INCLUDED
	filename=${file%_*} 
	echo "Generating sorted bam with duplicates removed for: ${filename##*/} ######################################" 
	#REMOVED THE NAME OF FILE FROM FIRST CHARACTER TO THE LAST /, /INCLUDED 
	bwa mem -t 2 0_reference/ref_genome.fasta ${filename}_1.fastq.gz ${filename}_2.fastq.gz | samtools view -@ 2 -b - | samtools sort -@ 2 - | samtools rmdup - 1_bam/${filename##*/}.sorted.rmdup.bam
done


echo “SAMTOOLS INDEX: Index removed duplicates and sorted BAM.” 
for file in 1_bam/*.sorted.rmdup.bam
do
	echo "Indexing ${file##*/}."
	samtools index $file 
done


echo “BCFTOOLS MPILEUP: Pileup to sum mapping of BAM on reference” 
echo “-q 	skip alignments with mapQ smaller than 30” 
echo “-f 	faidx indexed reference sequence file” 
echo “-O u 	output uncompressed VCF/BCF” 
echo “BCFTOOLS CALL: Convert pileup to human-readable format,”
echo “-O v	output uncompressed VCF ”
echo “-mv	 alternative model for multiallelic and rare-variant calling, output variant sites only ”

for file in 1_bam/*.sorted.rmdup.bam
do
	filename=${file##*/} 
	echo "Generating human readable pileup to sum mapping for ${filename%%.*} ######################################"
	bcftools mpileup --threads 10 -q 30 -O u -f 0_reference/ref_genome.fasta $file | bcftools call -O v -mv - > 2_vcf/${filename%%.*}.vcf 
done


echo "PYTHON VCF_TO_FASTA_HET.PY:"
echo "Variants convert to fasta by the python tool, also a table of variants."
echo "vcf_x	exclude snp similar in all variants, since we are comparing to reference0 as a patient, we should use vcf, not vcf_x"
python $PYTHON_DIR/vcf_to_fasta_het.py 2_vcf/ project1_vcf

echo "MV: Move variant fasta files to 2_vcf"
mv 2_vcf/project1_vcf* 3_aln

echo "FASTTREE: Generate FastTree phylogeny."
FastTree 3_aln/project1_vcf.fasta > 4_tree/project1.nwk

echo "BCFTOOLS FILTER, PYTHON VCF_TO_FASTA_HET.PY, FASTTREE:"
echo "Discard variants with mapping quality less than 200, python - vcf to fasta, fasta to tree"
for file in 2_vcf/*.vcf
do 
	filename=${file##*/} 
	echo "Discarding MQUAL < 200 for ${filename%%.*} ######################################"
	bcftools filter -e "QUAL < 200" $file > 2_vcf_clean/${filename%%.*}.filtered.vcf
done

python $PYTHON_DIR/vcf_to_fasta_het.py -x 2_vcf_clean/ project1_clean_vcf 
mv 2_vcf_clean/project1_clean_vcf* 3_aln
FastTree 3_aln/project1_clean_vcf.fasta > 4_tree/project1_clean.nwk

echo "Project 1 BASIC and MQUAL200 Completed."
echo `date`