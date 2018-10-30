#!/bin/bash
cd ~ 
mkdir ~
for file in /projects/micb405/resources/project_1/*.fastq.gz
do
fastqc -t 2 -o ~/project1_fastqc/ file
echo "Done."
