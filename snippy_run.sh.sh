
###run with:
###nohup ./project1_snippy.sh > project1_snippy.out &
#!/bin/bash
clear
echo `date`
echo "MICB 405 2018W1 Project 1 - Group 4"
echo "Run SNIPPY Method:"


echo "Declare variables for file path & make directories."
REF_DIR=/projects/micb405/resources/project_1
mkdir 0_reference

echo "Copy reference genome to dir 0_reference."
cp $REF_DIR/ref_genome.fasta 0_reference
bwa index 0_reference/ref_genome.fasta

echo "Generating input.tab file as input to SNIPPY."
>input.tab
for file in $REF_DIR/*_1.fastq.gz
do
filename=${file%_*}
echo -e "${filename##*/}\t${filename}_1.fastq.gz\t${filename}_2.fastq.gz">> input.tab
done

echo "SNIPPY-MULTI:"
echo "Generating SNIPPY commands to snippy.sh, change permission to 744 to run"
snippy-multi input.tab --cpus 16 --ref 0_reference/ref_genome.fasta > snippy.sh
chmod 744 snippy.sh

echo "Correcting last line of snippy, reference of snippy-core should be the input reference"
sed -i 's/Bat\/ref.fa/0_reference\/ref_genome.fasta/' snippy.sh

echo "RUNNING SNIPPY.SH"
./snippy.sh

echo "SNIPPY-CLEAN_FULL_ALN: Cleanup the random characters in alignment file of snippy."
snippy-clean_full_aln core.aln > clean.core.aln
snippy-clean_full_aln core.full.aln > clean.full.aln

### run for SNP tree (no gubbins)
#run_gubbins.py -p gubbins clean.full.aln
#snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

echo "FASTTREE: Generating Phylogenetic tree using core (variants) only." 
FastTree -gtr -nt clean.core.aln > clean.core.tree

echo "FASTTREE: Generating Phylogenetic tree using full sequence." 
FastTree -gtr -nt clean.full.aln > clean.full.tree

echo "Project 1 SNIPPY Completed."
echo `date`

# --mincov - the minimum number of reads covering a site to be considered (default=10)
# --minfrac - the minimum proportion of those reads which must differ from the reference
# --minqual - the minimum VCF variant call "quality" (default=100)
# --mapqual is the minimum mapping quality to accept in variant calling. BWA MEM using 60 to mean a read is "uniquely mapped".
# --basequal is minimum quality a nucleotide needs to be used in variant calling. We use 13 which corresponds to error probability of ~5%. It is a traditional SAMtools value.
# If you call SNPs for multiple isolates from the same reference, you can produce an alignment of "core SNPs" which can be used to build a high-resolution phylogeny (ignoring possible recombination). A "core site" is a genomic position that is present in all the samples. A core site can have the same nucleotide in every sample ("monomorphic") or some samples can be different ("polymorphic" or "variant"). If we ignore the complications of "ins", "del" variant types, and just use variant sites, these are the "core SNP genome



