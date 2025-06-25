#!/bin/bash

set -o pipefail

# read in the STAR module
ml STAR/2.7.10b-GCC-11.3.0

# Define Exon Parent Gene and Transcript Tags

if [ "$ANNOTATION_FORMAT" == "GFF3" ]; then
	echo "Using GFF3 annotation to generate a genome index"
	TRANSCRIPT_TAG="Parent"
	if [[ ! -z "$GENE_PARENT" ]]; then
		GENE_TAG="${GENE_PARENT}"
	else
		GENE_TAG="gene_id"
	fi
elif [ "$ANNOTATION_FORMAT" == "GTF" ]; then
	echo "Using GTF annotation to generate a genome index"
	TRANSCRIPT_TAG="transcript_id"
	GENE_TAG="gene_id"
else
	echo "Please specify whether annotation file is in GTF or GFF3 format"
fi

if [[ "$QUEUE" == "PBS" ]]; then
    echo "PBS is our workload manager/job scheduler."
    ${STAR_FILE} \
--runThreadN $NTHREAD \
--runMode genomeGenerate \
--genomeDir $GEN_DIR \
--genomeFastaFiles $GEN_FASTA \
--sjdbGTFtagExonParentTranscript $TRANSCRIPT_TAG \
--sjdbGTFtagExonParentGene $GENE_TAG \
--sjdbGTFfile $GEN_ANN \
--sjdbOverhang $SPLICE_JUN

elif [[ "${QUEUE}" == "Slurm" ]]; then
    echo "Slurm is our workload manager/job scheduler."	
    STAR \
--runThreadN $NTHREAD \
--runMode genomeGenerate \
--genomeDir $GEN_DIR \
--genomeFastaFiles $GEN_FASTA \
--sjdbGTFtagExonParentTranscript $TRANSCRIPT_TAG \
--sjdbGTFtagExonParentGene $GENE_TAG \
--sjdbGTFfile $GEN_ANN \
--sjdbOverhang $SPLICE_JUN

fi
