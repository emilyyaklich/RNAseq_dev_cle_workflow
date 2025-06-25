# cle_RNAseq
A pipeline for analyzing expression data. 
Built from Emily Dittmar's sunflower expression pipeline: https://github.com/EDitt/Sunflower_RNAseq


# Part I: Mapping reads and quantifying transcripts

## Programs Used:  
FASTQC: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf  
MultiQC: https://multiqc.info/  
STAR: http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf  

To run `dev_RNAseq`, use the following command, assuming you are in the `dev_RNAseq` directory:  
`./dev_RNAseq.sh <handler> Config`  
Where `<handler>` is one of the handlers listed below, and `Config` is the full file path to the configuration file

## Pre-processing

Start with raw sequence data. Simply copy this data into your working scratch directory.

## Step 1: Quality_Assessment
Run Quality_Assessment on your raw FastQ files. 

To run Quality_Assessment, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Quality_Assessment can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq`)
`./dev_RNAseq.sh Quality_Assessment Config`
where `Config` is the full file path to the configuration file

A directory containing your files to be analyzed must be specified in the config file. It is ok if this is a directory containing sub-directories for each sample (which is the format for raw data as it comes from Basespace).

After quality has been assessed for each sample, the FastQC results will be summarized using MultiQC. These summary statistics will be located in the output directory specified in the config file.


## Step 2: Genome_Index  

This handler will generate a genome index using FASTA and GFF3 or GTF formatted annotations. This step only needs to be performed once for each genome/annotation combination.

If using a GFF3 file for genome indexing rather than the default GTF file, the option `--sjdbGTFtagExonParentTranscript Parent` is added to the script 

To run Genome_Index, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Genome_Index can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq_RNAseq`)  
`./dev_RNAseq.sh Genome_Index Config`  
where `Config` is the full file path to the configuration file.

You will use the contents of the output (directory specified in the Config file) for the next step

## Step 3: Read_Mapping (and transcript quantification)

The Read_Mapping handler uses STAR to map reads to the genome indexed in step 3.

This handler can accept as input EITHER a directory or a text-file list of forward samples (it will find the reverse samples based on the naming suffix specified in the config file)

**We are using -GeneCounts flag of STAR (specified in Config file) which will also quanitfy our transcripts and the output will be in a `.tab` file (one file for each sample) with the rest of the STAR output. These are the files that we will use to analyze our expression data, so this is the most important output from this step.**

### Option for 2-pass mapping
STAR can perform a 2-pass mapping strategy to increase mapping sensitivity around novel splice junctions. This works by running a 1st mapping pass for all samples with the "usual" parameters to identify un-annotated splice junctions. Then mapping is re-run, using the junctions detected in the first pass (in addition to the annotated junctions specified in the genome annotation file). This 2-pass mapping strategy is recommended by GATK and ENCODE best-practices for better alignments around novel splice junctions.

While STAR can perform 2-pass mapping on a per-sample basis, in a study with multiple samples, it is recommended to collect 1st pass junctions from all samples for the second pass. Therefore, the recommended 2-pass procedure is described below: 

#### Step 3a: Collect_Junctions
This step identifies novel junctions during a first read-mapping pass and outputs them as "SJ.out.tab" files for each sample. In the handler used here, only junctions supported by at least 1 uniquely mapped read will be output. This step shares variables for Read_Mapping in the config file, so make sure these are filled out.

Once the variables have been defined, Collect_Junctions can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq`)  
`./dev_RNAseq.sh Collect_Junctions Config`   
where `Config` is the full file path to the configuration file.

#### Step 3b: Filter_Junctions
This step will concatenate the junction files discovered in Step 4a across samples, and then filter them based on user-defined parameters in the config file. This step is not required, as STAR will automatically concatenate the "SJ.out.tab" files before mapping. If you want to skip this step, you can instead just pass a list of all the "SJ.out.tab" files for all samples from Step 4a directly into the `FILTERED_JUNC_LIST` variable for Read_Mapping. (In fact, this is the 2-pass mapping procedure described in the STAR manual). However, filtering junctions is recommended when you have large numbers of samples, and there are several reasons we have implemented this intermediate filtering step: 
1.) Spurious junctions may increase the number of reads that map to multiple places in the genome   
2.) Filtering and concatenating junction files before mapping speeds up the mapping step in 4c - both because of the smaller number of junctions and because STAR then doesn't need to perform the concatenation across large numbers of samples for each sample separately.  
3.) This list or lists can be more easily saved in case one wants to redo the mapping  

Once the variables in the configuration file have been defined, Filter_Junctions can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq`)  
`./dev_RNAseq.sh Filter_Junctions Config`   
where `Config` is the full file path to the configuration file.

#### Step 3c: Read_Mapping  
This step will take all of the novel junctions discovered in the first pass and use them to re-map reads for each sample (in addition to the already-annotated junctions from your annotation file). 

The `FILTERED_JUNC_LIST` variable can be one of two things: 1.) The filepath to the filtered junctions file output from the "Filter_Junctions" handler. Alternatively, if you have multiple outputs from "Filter_Junctions" (if you process samples in batches like we do); this variable can be 2.) A .txt file with a list of the full filepaths to all filtered junction files to be included. STAR can handle multiple junction files as input (and will concatenate before mapping).
As mentioned previously, this list could also just be a list of filepaths to the un-filtered, un-concatenated "SJ.out.tab" files from all samples (if skipping the "Filter_Junctions" step).

To run Read_Mapping, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Read_Mapping can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq`)  
`./dev_RNAseq.sh Read_Mapping Config`   
where `Config` is the full file path to the configuration file.

### Option for 1-pass mapping  
If you want to map your reads without the addition of novel junctions discovered in a first mapping step, you can skip the "Collect_Junctions" and "Filter_Junctions" steps and leave the `FILTERED_JUNC_LIST` variable blank. All other variables need to be specified for Read_Mapping in the configuration file. Once the variables have been defined, Read_Mapping can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `dev_RNAseq`)  
`./dev_RNAseq.sh Read_Mapping Config`   
where `Config` is the full file path to the configuration file.

### Final Notes:  
#### Read groups  
If you have sequence data from the same sample across multiple lanes/runs, the best practice is to map these separately (in order to test for batch effects).

#### Alignment File Output  
Two alignment files are output from STAR with this handler - one alignment file in genomic coordinates and one translated into transcript coordinates (e.g.`*Aligned.toTranscriptome.out.bam`). The default format of the former is an unsorted SAM file. If you plan on using the genomic coordinate alignments for SNP calling, you have the option of getting these output as coordinate-sorted BAM files (similar to the `samtools sort` command) by putting a "yes" for the `GENOMIC_COORDINATE_BAMSORTED` variable in the config. Note that this will add significant computational time and memory.

## Step 4: Prepare for DESEq2

The output of STAR GeneCounts must be properly formatted to read into DESEq2. The file format needs to be changed as the GeneCounts output has 4 columns: 1-GeneID, 2-unstranded counts, 3-forward strand counts, and 4-reverse strand counts. DESeq can only read in files with 2 columns (GeneID and one count column). The ```prepare_DESeq_input.sh``` extracts the first and fourth columns from the files (the data we are using is reverse-stranded).

# Part II: Expression Analysis 


This will mostly be in R using [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for differential expression analysis.

Special note on DGE analysis: There are many, many programs that will run DGE analyses with tons of options for methods within each program (on example being edgeR)...it is crucial to understand each test and how it is handling your data in order to make a decision for which program to use. That being said, the following methods are what I found were best for my data, but there are many ways to do similar analyses. 

## 0. Pre-process data into matrix
Load the GeneCount data into a DESeq dataset and collapse technical replicates (if applicaple)
`load_GC_data.R`
## 1. Visualize data in an MDS plot
`plot_MDS.R`
This will give us an idea of the relationship between our samples. In the next step we will be correcting for varaition, but it is important to visualize at this step.
## 2. Correct for unwanted variation using [ComBat-Seq](https://github.com/zhangyuqing/ComBat-seq) and run DEseq
`run_combatseq.R` and `run_DGE_deseq_CLE.R`
We will visualize the corrected samples in an MDS plot as well (compare with above). The model we are using to run the differential epxression analysis is ~0+dev_stage as we are curious about differential expression wrt dev_stage. 
## 3. Visualize DE in an Upset Plot
`analyze_DGE_deseq_CLE.R`

## 4. Add ortholog info to DESeq output
`add_orthogroup.py`


