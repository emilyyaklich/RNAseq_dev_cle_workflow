#!/bin/bash

set -o pipefail

#   Where is 'dev_rnaseq' located?
DEV_RNASEQ=$(pwd -P)

ROUTINE="$1" # What routine are we running?
CONFIG="$2" # Where is our config file?

   #If the specified config exists
if [[ -f "${CONFIG}" ]]
then
    source "${CONFIG}" # Source it, providing parameters and software
else # If it doesn't
    echo "Please specify a valid config file." >&2 # Print error message
    exit 1 # Exit with non-zero exit status
fi

#   Where do we output the standard error and standard output files?
ERROR="${DEV_RNASEQ}"/ErrorFiles/"${PROJECT}"
mkdir -p "${ERROR}"


#   Run dev_RNAseq
case "${ROUTINE}" in
    1 | Quality_Assessment)
        echo "$(basename $0): Assessing quality..." >&2
        if [[ "$QUEUE" == "PBS" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Quality_Assessment.sh" | qsub -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Assessment
            elif [[ "${QUEUE}" == "Slurm" ]]; then
                echo "Slurm is our workload manager/job scheduler."
                sbatch --job-name=${PROJECT}_Quality_Assessment ${QA_SBATCH} --output=${ERROR}/QA_slurm-%j.out --export=QA_INPUTDIR=${QA_INPUTDIR},SUFFIX=${SUFFIX},QA_OUTPUTDIR=${QA_OUTPUTDIR},QA_TEMP=${QA_TEMP} ${DEV_RNASEQ}/Quality_Assessment.sh
            else
                echo "QUEUE variable in config must be set to PBS or Slurm. Please set to one of the two depending on the workload manager your cluster uses. Exiting..."
                exit 1
        fi
        ;;
    2 | Adapter_Trimming)
        echo "$(basename $0): Trimming Adapters..." >&2
        if [[ -d "$AT_INPUT" ]]; then #if input is a directory
            echo "$AT_INPUT is a directory"
            declare -a files #an array of files
            for f1 in `find $AT_INPUT -name "*$FORWARD_NAMING"`; do
                if [[ -f "$f1" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "Please specify a path to valid files in the config file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$AT_INPUT" ]]; then #if input is a file
            echo "$AT_INPUT is a file"
            Maxarray=$(< $AT_INPUT wc -l)
        else
            echo "Please specify a valid directory or list in the config"
        fi
        echo "Max array index is ${Maxarray}">&2
        if [[ "$QUEUE" == "PBS" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Trimm.sh" | qsub -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Adapter_Trimming -t 1-"${Maxarray}"
        elif [[ "${QUEUE}" == "Slurm" ]]; then
            echo "Slurm is our workload manager/job scheduler."
            Slurm_Maxarray=$(($Maxarray-1))
            sbatch --job-name=${PROJECT}_Adapter_Trimming ${AT_SBATCH} --output=${ERROR}/AT_slurm_%A_%a.out --array=0-${Slurm_Maxarray} \
            --export=QUEUE=${QUEUE},JOB_LOG=${JOB_LOG},AT_INPUT=${AT_INPUT},AT_OUTPUTDIR=${AT_OUTPUTDIR},ADAPTERFILE=${ADAPTERFILE},FORWARD_NAMING=${FORWARD_NAMING},REVERSE_NAMING=${REVERSE_NAMING},SEEDMISMATCH=${SEEDMISMATCH},PALINDROMECLIP=${PALINDROMECLIP},SIMPLECLIP=${SIMPLECLIP},MINADAPTERLEN=${MINADAPTERLEN},KEEPREADS=${KEEPREADS},LEADCUT=${LEADCUT},TRAILCUT=${TRAILCUT},MINLENGTH=${MINLENGTH},PE=${PE} ${DEV_RNASEQ}/Trimm.sh
        else
                    echo "QUEUE variable in config must be set to PBS or Slurm. Please set to one of the two depending on the workload manager your cluster uses. Exiting..."
                    exit 1
        fi
        ;;
    3 | Genome_Index)
        echo "$(basename $0): Generating a genome index from ${ANNOTATION_FORMAT} annotation..." >&2
        if [[ "$QUEUE" == "PBS" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Genome_Index.sh" | qsub -l "${GI_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Genome_Index
        elif [[ "${QUEUE}" == "Slurm" ]]; then
            echo "Slurm is our workload manager/job scheduler."
            sbatch --job-name=${PROJECT}_Genome_Index ${GI_SBATCH} --output=${ERROR}/GI_slurm-%j.out --export=QUEUE=${QUEUE},ANNOTATION_FORMAT=${ANNOTATION_FORMAT},GENE_PARENT=${GENE_PARENT},STAR_FILE=${STAR_FILE},NTHREAD=${NTHREAD},GEN_DIR=${GEN_DIR},GEN_FASTA=${GEN_FASTA},TRANSCRIPT_TAG=${TRANSCRIPT_TAG},GENE_TAG=${GENE_TAG},GEN_ANN=${GEN_ANN},SPLICE_JUN=${SPLICE_JUN} ${DEV_RNASEQ}/Genome_Index.sh
            fi
        ;;
    4 | Collect_Junctions)
        if [[ ! -d "$GEN_DIR" ]]; then #if genome directory is not specified
            echo "Please specify a valid filepath to the genome directory, exiting..."
            exit 1
        fi
        export RM_PASS="first"
        declare -a files #an array of files
        if [[ -d "$RM_INPUT" ]]; then #if input is a directory
            echo "$RM_INPUT is a directory"
            for f in `find $RM_INPUT -name "*$FORWARD"`; do
                if [[ -f "$f" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "$f is not a file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$RM_INPUT" ]]; then #if input is a file
            echo "$RM_INPUT is a file"
            Maxarray=$(< $RM_INPUT wc -l)
        else
            echo "Please specify a valid directory or list in the config"
            exit 1
        fi
        if [ "$PE" == "True" ]; then
            echo "$(basename $0): Running 1st-pass Mapping to Identify Novel Splice Junctions using PE Reads..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Running 1st-pass Mapping to Identify Novel Splice Junctions using SE Reads..." >&2
        else
            echo "Please specify in the config file whether data is PE (True/False), exiting..."
            exit 1
        fi
        echo "Max array index is ${Maxarray}">&2
        if [[ "$QUEUE" == "PBS" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Read_Mapping.sh" | qsub -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Collect_Junctions -V -t 1-"${Maxarray}"
        elif [[ "${QUEUE}" == "Slurm" ]]; then
            echo "Slurm is our workload manager/job scheduler."
            Slurm_Maxarray=$(($Maxarray-1))
	    sbatch --job-name=${PROJECT}_Collect_Junctions ${RM_SBATCH} --output=${ERROR}/CJ_slurm-%j.out --array=0-${Slurm_Maxarray} --export=QUEUE=${QUEUE},RM_INPUT=${RM_INPUT},CJ_JOB_LOG=${CJ_JOB_LOG},RM_NTHREAD=${RM_NTHREAD},FLOWCELL_NAME=${FLOWCELL_NAME},FORWARD=${FORWARD},CJ_OUTPUTDIR=${CJ_OUTPUTDIR},PE=${PE},REVERSE=${REVERSE},GENOMIC_COORDINATE_BAMSORTED=${GENOMIC_COORDINATE_BAMSORTED},RM_PASS=${RM_PASS},MAX_MIS=${MAX_MIS},MAX_N=${MAX_N},UNMAP_F=${UNMAP_F},MINSCORE_READL=${MINSCORE_READL},MINMATCH_READL=${MINMATCH_READL},SEEDSEARCH=${SEEDSEARCH},PLATFORM=${PLATFORM},GEN_DIR=${GEN_DIR} ${DEV_RNASEQ}/Read_Mapping.sh
        else
                echo "QUEUE variable in config must be set to PBS or Slurm. Please set to one of the two depending on the workload manager your cluster uses. Exiting..."
                exit 1
        fi



;;
    5 | Filter_Junctions)
        echo "$(basename $0): Filtering and concatenating junctions for mapping..." >&2
        if [[ "$QUEUE" == "PBS" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Filter_Junctions.sh" | qsub -l "${FJ_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Filter_Junctions
        elif [[ "${QUEUE}" == "Slurm" ]]; then
                echo "Slurm is our workload manager/job scheduler."
                sbatch --job-name=${PROJECT}_Filter_Junctions ${FJ_SBATCH} --output=${ERROR}/FJ_slurm-%j.out --export=JUNCTIONDIR=${JUNCTIONDIR},SJ_LISTNAME=${SJ_LISTNAME},SCAFFOLD_STRING=${SCAFFOLD_STRING},REMOVE_NC_JUNK=${REMOVE_NC_JUNK},UNIQUE_NUM=${UNIQUE_NUM} ${DEV_RNASEQ}/Filter_Junctions.sh
            else
                echo "QUEUE variable in config must be set to PBS or Slurm. Please set to one of the two depending on the workload manager your cluster uses. Exiting..."
                exit 1
fi

;;
    6 | Read_Mapping)
        if [[ ! -d "$GEN_DIR" ]]; then #if genome directory is not specified
            echo "Please specify a valid filepath to the genome directory, exiting..."
            exit 1
        fi
        export RM_PASS="second"
        declare -a files #an array of files
        if [[ -d "$RM_INPUT" ]]; then #if input is a directory
            echo "$RM_INPUT is a directory"
            for f in `find $RM_INPUT -name "*$FORWARD"`; do
                if [[ -f "$f" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "$f is not a file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$RM_INPUT" ]]; then #if input is a file
            echo "$RM_INPUT is a file"
            Maxarray=$(< $RM_INPUT wc -l)
        else
            echo "Please specify a valid directory or list of input files in the config"
            exit 1
        fi
        if [ "$PE" == "True" ]; then
            echo "$(basename $0): Mapping PE Reads..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Mapping SE Reads..." >&2
        else
            echo "Please specify in the config file whether data is PE (True/False), exiting..."
            exit 1
        fi
        if [[ -f "$FILTERED_JUNC_LIST" ]]; then #if junction list variable is set as a file
            file_content=$(head -1 ${FILTERED_JUNC_LIST})
            if [[ -f "${file_content}" ]]; then #if first line is a file
                declare -a junctions ### make an array of filtered junction files
                while read line; do
                    junctions=("${junctions[@]}" "$line")
                done < $FILTERED_JUNC_LIST
                export JUNCTIONS="${junctions[@]}"
                export NUM_JUNCTIONS="${#junctions[@]}"
            else #if file input is not a list of files, but list of junctions
                export JUNCTIONS="${FILTERED_JUNC_LIST}"
                export NUM_JUNCTIONS="one"
            fi
            echo "In second-pass mode using ${NUM_JUNCTIONS} junction files"
        elif [[ -z "$FILTERED_JUNC_LIST" ]]; then #if no input here
            echo "Read Mapping without incorporating un-annotated junctions"
        else
            echo "A junction list is specified in config but is not a valid file, exiting..."
            exit 1
        fi
        echo "Max array index is ${Maxarray}">&2
        if [[ "$QUEUE" == "PBS" ]]; then
		echo "PBS is our workload manager/job scheduler."

            echo "source ${CONFIG} && source ${DEV_RNASEQ}/Read_Mapping.sh" | qsub -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Read_Mapping -V -t 1-"${Maxarray}"
       
        elif [[ "${QUEUE}" == "Slurm" ]]; then
            echo "Slurm is our workload manager/job scheduler."
            Slurm_Maxarray=$(($Maxarray-1))
            sbatch --job-name=${PROJECT}_Read_Mapping ${RM_SBATCH} --output=${ERROR}/RM_slurm-%j.out --array=0-${Slurm_Maxarray} --export=QUEUE=${QUEUE},RM_INPUT=${RM_INPUT},RM_JOB_LOG=${RM_JOB_LOG},RM_NTHREAD=${RM_NTHREAD},FLOWCELL_NAME=${FLOWCELL_NAME},FORWARD=${FORWARD},RM_OUTPUTDIR=${RM_OUTPUTDIR},QUANT=${QUANT},PE=${PE},REVERSE=${REVERSE},GENOMIC_COORDINATE_BAMSORTED=${GENOMIC_COORDINATE_BAMSORTED},RM_PASS=${RM_PASS},MAX_MIS=${MAX_MIS},MAX_N=${MAX_N},UNMAP_F=${UNMAP_F},MINSCORE_READL=${MINSCORE_READL},MINMATCH_READL=${MINMATCH_READL},SEEDSEARCH=${SEEDSEARCH},PLATFORM=${PLATFORM},GEN_DIR=${GEN_DIR} ${DEV_RNASEQ}/Read_Mapping.sh
fi



 ;;
 * )
esac 
