# Targeted_Sequencing_Pipeline

# Function
A modified version of the GATK Best Practices pipeline (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) adapted to analyze Haloplex Sequencing reads. With minor modifications the pipeline can be used to analyse any NGS Targeted sequencing. The steps where the pipeline diverges from the most mainstream applications, have been marked in the comments.
The pipeline takes fastq files of multiple samples, performs joint genotyping on them and outputs two separate vcfs: one containing the SNPs of all samples and one containing the indels of all samples. The pipeline already includes hard filtering, which should always be tailored to the dataset at hand. Large datasets, such as WGS or WES data on many individuals can be filtered using VQSR and do not require hard filtering. More information on hard filtering can be found at:
https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations
https://drive.google.com/drive/folders/0BwTg3aXzGxEDczFPaGt0eTVyMjA

# Dependencies:
FastQC, Trimmomatic, Samtools, GATK, Picard,
The script is written to perform genotype calling using the GATK v4.0.1.2, however, as certain walkers of the toolkit were not functioning properly at the time of the creation of the pipeline, variant annotation and variant filtering steps are performed using an older version of GATK.

# Required files:
1. Fastq files of sequencing data (gzipped format is accepted)
2. A FASTA file listing the adapters that were used in the experiments: haloplexadapters.fa
3. A list of the targeted intervals: Target3.interval.list
4. A list file of the paths to where the .vcfs of the joint genotype calling are expected: intervalmap.list
5. A table of the samples and the path to where their .vcfs files are expected: samplemap_all.tsv
Optionally:
6. VCFs of known indels and SNPs in the examined organism

# Divergence from the Best Practices
1. The Haloplex Target Enrichment kit uses restrictive endonucleases for fragmentation, therefore the first 5 nucleotides of each read were trimmed, to avoid reference bias.
2. The Haloplex Target Enrichment kit used in our experiments did not use barcoding for PCR and therefore PCR duplicates are indistinguishable from other fragments that were cut by the same restriction endonuclease. Thus, PCR duplicates are not removed by the pipeline. The pipeline marks PCR duplicates, but downstream programs do not disregard reads flagged as duplicates. For most other applications (except for amplicon sequencing) the removal of PCR duplicates is advisory.
3. The use of hard filtering is only recommended for small datasets (such as WES of < 30 samples or targeted sequencing).
