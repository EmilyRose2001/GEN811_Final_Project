# GEN811_Final_Project: Bacterial Genome Assembly

### By: Emily Persson
## Background
The data used in this project are forward and reverse sequences from two different unidentified bacterial organisms. In total, two forward and two reverse sequences were used. The project aims to identify the two organisms by assembling their genomes from these forward and reverse sequences. 
## Methods
The sequences from the two unidentified bacterial organisms were provided by Kaleb Ducharme, a TA at the University of New Hampshire. They were sequenced using Sanger sequencing and Illumina. Using the university's supercomputer RON and multiple tools, the genomes were assembled.
### Tools
#### FastQC
The purpose of using FastQC was to take a look at the read quality of the sequences. FastQC outputs an HTML and a zip file for each given sequence. There are multiple read quality components in a FastQC HTML file; basic statistics, per base sequence quality, per tile sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, and adapter content. This HTML file was downloaded from RON to examine the FastQC findings.
#### Trimmomatic
Trimmomatic is a tool that removes low-quality reads and adapter content from the sequence. A wrapper script for trimmomatic, created by Joseph Sevigny, was used in this project. The script takes two arguments: the forward read and the reverse read. The output is a directory called 'trimmed-reads' containing trimmed forward, reverse, unpaired forward, and unpaired reverse gzipped reads. FastQC is then performed on the trimmed reads to compare the quality results to the raw reads.
#### SPAdes
SPAdes is a genome assembly tool that is great with most bacterial genomes. SPAdes takes four input files: the trimmed forward, reverse, unpaired forward, and unpaired reverse reads. The output is a directory containing several files. The only files that are important for this assembly are the 'contigs.fasta' and 'spades.log' files. All files except for the ones mentioned were removed from the output spade directory.
#### QUAST
Quast is used to assess the assembly created by SPAdes. It assesses the assembly by taking a look at the quality of the assembly. QUAST only takes one argument, and it is the 'contigs.fasta' file produced by SPAdes. The output of this tool is another directory containing several files. Most of the output files are the same, just in a different format. A list of formats includes txt, gff, csv...
#### BUSCO
BUSCO is used to determine how complete the assembled genome is. BUSCO takes the same input as QUAST, which is the 'contigs.fasta' output file from SPAdes. The output is also similar to QUAST, where there is an output directory containing several files. Two of the main output files are the 'short_summary_busco_output.txt' and the 'full_table_busco_output.tsv'. 
#### PROKKA
PROKKA is a tool used to annotate genomes. The input is also the same as QUAST and BUSCO, which is the 'contigs.fasta' file. The output is also similar, where there is an output directory containing various annotation files. These files include a GFF file for gene annotations, a FFN FASTA file that was used for the annotation, and an FAA FASTA file that was used for the annotation. Using the output GFF file, a count for each gene annotation was also created.
#### 'extract_sequence'
This tool extracts sequences from a PROKKA FFN file based on a specified sequence header. The input is a string representing the sequence header to search for and a PROKKA FFN file to search in. This tool was used to extract the '16S ribosomal RNA' sequences from the output PROKKA FFN file, which were used to identify the mystery organisms.
#### BLAST
BLAST is an alignment tool that aligns a sequence to a reference sequence. BLAST was used multiple times in this project: BLAST the extracted 16S sequences against a custom database created from the 'contigs.fasta' file, BLAST the whole assembly ('contigs.fasta') against the nucleotide database, and BLAST the filtered assembly against the UniVec database. The purpose of BLASTing the 16S sequences is to compare the taxonomic identification to the output of the whole assembly BLAST. To create the 'contigs.fasta' database, the input is the 'contigs.fasta' file, and what database type being created. For this project, a nucleotide database was created, so the input was 'nucl'. The output is multiple 'contigs_db' files. To run blastn on the 16S sequence FASTA file, the input is the FASTA file and the created database. The output is a TSV file. The purpose of BLASTing the whole assembly is to identify which taxonomy each sequence belongs to. The purpose of the UniVec BLAST is to check for contamination. The input is the filtered assembly FASTA, and the output is a file 'genome_vs_univec.6'.
#### blast-ncbi-nt.sh
A tool that will be used further down the pipeline is Blobtools. This tool requires a specific input BLAST file format. To make this easier, a BLAST script created by Joseph Sveigny was used to BLAST the whole assembly against the nucleotide database. This script takes one argument, which is the 'contigs.fasta' file, and outputs a file named 'contigs.fasta.vs.nt.cul5.1e5.megablast.out'. 
#### BWA mem
BWA mem is a read mapping tool that aligns the short sequences of the forward and reverse reads to a reference genome. The input includes the reference genome FASTA file, the forward read FASTQ, and the reverse read FASTQ. The reference genome for this project is the 'contigs.fasta' file. The output is a SAM file. This file contains information on how and where the sequences map to the reference genome.
#### SAMtools
SAMtools is used to convert a SAM file to a BAM file, where a BAM file is the binary version of a SAM file. The input for the conversion is the SAM file created by BWQ mem, and the output is a BAM version of that file.
#### Bedtools
This tool is used to perform genomic arithmetic and analysis. It is used in this project to calculate the per-base coverage of the assembly. The input for this is the BAM file created by SAMtools, and the output is a file named 'coverage.out'.
#### 'gen_input_table.py'
This Python script is used to calculate per-contig coverage. The input for this script is the 'contigs.fasta' file and the 'coverage.out' file created by Bedtools. The output is a file named 'coverage_table.tsv', which is a table with two columns. These columns are the contig header and the average coverage.
#### BlobTools
Blobtools is a genome assembly visualization tool that provides insights into what information should be filtered out of the assembly. The input consists of the 'contigs.fasta' file, the BAM file, and the output file 'contigs.fasta.vs.nt.cul5.1e5.megablast.out' from BLAST. The output is a JSON file, which is used to create a taxonomy table and a plot. These outputs are then used to filter the genome assembly.
#### filter_contigs_by_list.py
This Python script, created by Joseph Seivgny, is designed to filter the assembly based on a list of contigs that have been chosen to remain in the assembly. The input is the 'contigs.fasta' file, the TXT file of the list of contigs, and a desired output name for a FASTA file. The output is the filtered genome assembly.
## Results
![02_fastqc](sample_02_fastqc.png)
#### Figure 1
This figure shows three of the FastqC output plots for the 02_S65_L001 sample. Both the FastQC results for the raw sequence and trimmed sequence are included. The FastQC for the raw sequence takes a FASTQ file of the raw sequence as its input, and the FastQC for the trimmed sequence takes a FASTQ file of the trimmed sequence as its input. For the adapter content plot, the adapter content present in the raw sequence is eliminated in the trimmed plot. The GC content for the sample is shown to be higher than the normal sequence GC content distribution. This is made slightly better by trimming the sequence, however, it is still high. For the sequence length distribution plot, the raw sequence shows a better distribution than the trimmed.
![03_fastqc](sample_03_fastqc.png)
#### Figure 2
![02_sample_blob](blob_out_02_S65_L001.blobDB.json.bestsum.genus.p8.span.100.blobplot.bam0.png)
This figure shows three of the FastQC output plots for the 03_S58_L001 sample. Like Figure 1, the plots shown are adapter content, sequence length distribution, and GC content for both the raw sequence and the trimmed sequence. The findings of these three plots are also very similar to sample 02_S65_L001. The only difference is that 
#### Figure 3
A blobtools plot depicting the greatest taxonomic hits for the unidentified sample 02_S65_L001. As seen by the centralized blue dot, the taxonomy with the highest hit is Clostridioides. The smaller blue dots further down the plot represent contigs that we want to filter out. This plot was created using Blobtools with the output JSON mentioned in the methods section for Blobtools.
![03_sample_blob](blob_out_03_S58_L001.blobDB.json.bestsum.genus.p8.span.100.blobplot.bam0.png)
#### Figure 4
Another blobtools plot showing the greatest taxonomic hits for the other unidentified sample, 03_S58_L001. The taxonomy of this sample is Escherichia. The smaller colorful dots further down the plot represent other taxonomic hits that we want to remove by filtering the genome.
## References
