#! /bin/bash

source activate genomics

#User inputs a directory containing the raw fastq files they want assembled
raw_fastq_dir=$1

#Create a directory for the output fastqc files
mkdir fastqc_raw_reads

#Read assessment

#Loop through the files in the inputed raw fastq directory
for file in $raw_fastq_dir/*.fastq*
do
    #Count the number of reads
    num_reads=$(zgrep "^@" $file | wc -l)

    #Calculate sequencing throughput
    bp_data=$(($num_reads * 250 * 2))

    #Calculate the average coverage
    avg_coverage=$(($bp_data / 7000000))

    #If the calculated average coverage is greater than or equal to 70, fastqc is performed on the file to determine the quality of the reads
    if [[ $avg_coverage -ge 70 ]]
    then
        fastqc $file -o fastqc_raw_reads
    fi
done

#Loop through the forward reads
for forward in $raw_fastq_dir/*_R1_*.fastq*
do
    #Extract the file name of the forward read
    forward_name=$(basename $forward)

    #Splice the first 11 characters from the file name. This is the sample identifier
    forward_identifier=${forward_name:0:11}

    #The forward read is saved to a variable
    forward_read=$forward

    #Loop through the reverse reads
    for reverse in $raw_fastq_dir/*_R2_*.fastq*
    do
        #Extract the file name of the reverse read
        reverse_name=$(basename $reverse)

        #Splice the first 11 characters from the file name. This is the sample identifier
        reverse_identifier=${reverse_name:0:11}

        #The reverse read is saved to a variable
        reverse_read=$reverse

        if [[ $forward_identifier == $reverse_identifier ]]
        then
            #If the sample identifiers for both the forward and reverse reads are the same, then the reads are trimmed
            trim_scriptV2.sh $forward_read $reverse_read
        fi
    done
done

trim_dir="./trimmed-reads"

#Create a directory for the output of the fastqc files
mkdir fastqc_trimmed_reads

#fastqc is performed on the trimmed files
for trim_file in $trim_dir/*.gz
do
    fastqc $trim_file -o fastqc_trimmed_reads
done

#Genome Assembly

#Loop through the forward read files in the trimmed-reads directory
for forward_trim in $trim_dir/*_R1_*.gz
do
    forward_name=$(basename $forward_trim)

    #Splice off the beginning of the forward_name string to identify 'unpaired' files
    unpaired_check=${forward_name:0:8}

    if [[ "unpaired" != $unpaired_check ]]
    then
        forward_identifier=${forward_name:0:11}

        #If the forward read is not 'unpaired', then the file is saved to the forward_read variable
        forward_read=$forward_trim

        #Loop through the reverse read files in the trimmed-reads directory
        for reverse_trim in $trim_dir/*_R2_*.gz
        do
            reverse_name=$(basename $reverse_trim)
            unpaired2_check=${reverse_name:0:8}

            if [[ "unpaired" != $unpaired2_check ]]
            then
                reverse_identifier=${reverse_name:0:11}

                #If the reverse read is not 'unpaired', then the file is saved to the reverse-read variable
                reverse_read=$reverse_trim

                #Loop through the the unpaired forward reads in the trimmed-reads directory
                for forward_unpaired in $trim_dir/unpaired*_R1_*.gz
                do
                    unpaired_for_name=$(basename $forward_unpaired)
                    unpaired_for_identifier=${unpaired_for_name:9:11}

                    #The read is saved to a variable
                    unpaired_for_read=$forward_unpaired

                    #Loop through the unpaired reverse reads in the trimmed-reads directory
                    for reverse_unpaired in $trim_dir/unpaired*_R2_*.gz
                    do
                        unpaired_rev_name=$(basename $reverse_unpaired)
                        unpaired_rev_identifier=${unpaired_rev_name:9:11}

                        #The read is saved to a variable
                        unpaired_rev_read=$reverse_unpaired

                        if [[ ($forward_identifier == $reverse_identifier) && ($unpaired_for_identifier == $unpaired_rev_identifier) && ($forward_identifier == $unpaired_for_identifier) ]]
                        then
                            #Spades is performed on the forward, reverse, unpaired forward, and unpaired reverse reads if their sample identifiers are the same
                            spades.py --isolate -1 $forward_read -2 $reverse_read -s $unpaired_for_read -s $unpaired_rev_read -o spade_assembly_default_$forward_identifier -t 24
                        fi
                    done
                done
            fi
        done
    fi
done


#Genome Assembly Directory Cleanup

#The only files in the spade_assembly_default directories that are needed are the contigs.fasta and spades.log
#Everything else can be removed

#Loop through the spade_assembly_default directories
for dir in ./spade_assembly_default_*
do
    dir_name=$(basename $dir)

    #Create a directory named after the spade_assembly_default directory and save it as the destination_dir variable
    mkdir move_$dir_name
    destination_dir=./move_$dir_name

    #Loop through the files in the directory
    for file in $dir/*
    do
        if [[ ($file == $dir/contigs.fasta) || ($file == $dir/spades.log) ]]
        then
            #If the file is either the contigs.fasta file or the spades.log file, it is moved to the destination directory
            mv $file $destination_dir
        fi
    done
done

#Loop through the spade_assembly_default directories
for rem_dir in ./spade_assembly_default_*
do
    rem_dir_name=$(basename $rem_dir)

    #Remove all files in the directory
    rm -r $rem_dir/*

    #Move the contigs.fasta and spades.log files back into the spade_assembly_default directory
    mv ./move_$rem_dir_name/contigs.fasta ./move_$rem_dir_name/spades.log $rem_dir

    #Remove the directory created for the cleanup
    rm -r ./move_$rem_dir_name
done

#Genome Assessment and Annotation

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the sample's contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Genome assembly assessment
    #Run QUAST on the contigs.fasta file. Output will be a directory named after the sample with QUAST files
    quast.py $contig_file -o quast_results_$sample_name

    #Genome assembly summary
    #Run BUSCO on the contigs.fasta file. Output will be a directory named after the sample with BUSCO files
    busco -i $contig_file -m genome -o busco_results_$sample_name -l bacteria

    #Gene annotation
    #Run PROKKA on the contigs.fasta file. Output will be a directory named after the sample with PROKKA files
    prokka $contig_file --outdir prokka_output_$sample_name --cpus 24 --mincontiglen 200

    #Using the .gff file outputted by PROKKA, an output file with counts for all gene annotations is created
    grep -o "product=.*" prokka_output_$sample_name/PROKKA_*.gff | sed "s/product=//g" | sort | uniq -c | sort -nr > protein_abundances_$sample_name.txt
done

#Organism Identification

for prokka_dir in ./prokka_output_*
do
    prokka_dir_name=$(basename $prokka_dir)
    sample_identity=${prokka_dir_name:14:11}

    #Determine if a 16S sequence exists using the sample's .ffn output file from PROKKA
    grep_output=$(grep -o "16S" $prokka_dir/*.ffn | wc -l)

    if [[ $grep_output -ge 1 ]]
    then
        #If 16S sequences exists (greater than or equal to 1), they are extracted
        extract_sequences "16S ribosomal RNA" $prokka_dir/PROKKA_*.ffn > 16S_sequence_$sample_identity.fasta
    fi
done

#BLAST

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the sample's contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Create a directory named after the sample for all BLAST files belonging to that sample
    mkdir blast_$sample_name

    #16S sequence is BLASTed first for a comparison

    #A database is created from the contigs.fasta file and saved into that sample's blast directory
    makeblastdb -in $contig_file -dbtype nucl -out ./blast_$sample_name/contigs_db_$sample_name
done

#Extract the sample's 16S_sequence.fasta created in the organism identification step
for seq_fasta in ./16S_sequence_*.fasta
do
    seq_fasta_name=$(basename $seq_fasta)
    sample_name=${seq_fasta_name:13:11}

    #BLAST the 16S_sequence.fasta file against the database created in the previous step
    blastn -query $seq_fasta -db ./blast_$sample_name/contigs_db_$sample_name -out ./blast_$sample_name/16S_vs_contigs_6_$sample_name.tsv -outfmt 6
done

#BLAST assembly
for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the sample's contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Change the current working directory to the blast directory named after the sample
    cd ./blast_$sample_name

    #Run BLAST on the contigs.fasta file
    blast-ncbi-nt.sh ../$contig_file

    #Change the current working directory back to where it was
    cd ../
done

#Read Mapping

#Create a directory for the unzipped trimmed forward and reverse reads
mkdir unzipped_trimmed_reads

#Loop through the zipped files in the trimmed-reads directory
for trim_fastq in ./trimmed-reads/*.gz
do
    name=$(basename $trim_fastq)

    #Splice off the beginning of the string to identify 'unpaired' files
    unpaired_check=${name:0:8}

    if [[ "unpaired" != $unpaired_check ]]
    then
        #Splice the sample identifier from the file name
        identifier=${name:0:18}

        #If the trimmed file is not an unpaired file, then the file is unzipped and saved into the unzipped_trimmed_reads directory
        zcat $trim_fastq > ./unzipped_trimmed_reads/$identifier.fastq
    fi
done

#Create a directory for all files involved in read mapping
mkdir read_mapping

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the sample's contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Extract the sample's unzipped forward trimmed read
    for forward_trim in ./unzipped_trimmed_reads/*_R1_*.fastq
    do
        forward_name=$(basename $forward_trim)
        forward_identifier=${forward_name:0:11}
        forward_read=$forward_trim

        #Extract the sample's unzipped reverse trimmed read
        for reverse_trim in ./unzipped_trimmed_reads/*_R2_*.fastq
        do
            reverse_name=$(basename $reverse_trim)
            reverse_identifier=${reverse_name:0:11}
            reverse_read=$reverse_trim

            if [[ ($sample_name == $forward_identifier) && ($forward_identifier == $reverse_identifier) ]]
            then
                #Index the contigs.fasta file. This is the reference genome
                bwa index $contig_file

                #Map the forward and reverse reads to the reference genome. The output is a .sam file
                bwa mem -t 24 $contig_file $forward_trim $reverse_trim > ./read_mapping/raw_mapped_$sample_name.sam
            fi
        done
    done
done

#Activate the samtools environment for access to all samtools commands
source activate samtools

#Extract the sample's .sam file
for sam_file in ./read_mapping/raw_mapped_*.sam
do
    sam_name=$(basename $sam_file)
    sam_identifier=${sam_name:11:11}

    #Reads that didn't match to the assembly are removed and a .bam file is outputted
    samtools view -@ 24 -Sb $sam_file | samtools sort -@ 24 -o ./read_mapping/sorted_mapped_$sam_identifier.bam
done

#Extract the sample's .bam file
for bam_file in ./read_mapping/sorted_mapped_*.bam
do
    bam_name=$(basename $bam_file)
    bam_identifier=${bam_name:14:11}

    #Determine how many reads mapped to the assembly
    samtools flagstat $bam_file

    #Index the .bam file
    samtools index $bam_file
done

#Activate the genomics environment again
source activate genomics

#Extract the sample's .bam file
for bam_file in ./read_mapping/sorted_mapped_*.bam
do
    bam_name=$(basename $bam_file)
    bam_identifier=${bam_name:14:11}

    #Use bedtools on the .bam file to calculate per base coverage
    bedtools genomecov -ibam $bam_file > ./read_mapping/coverage_$bam_identifier.out
done

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the sample's contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Extract the sample's coverage.out file
    for coverage_file in ./read_mapping/coverage_*.out
    do
        coverage_file_name=$(basename $coverage_file)
        coverage_identifier=${coverage_file_name:9:11}

        if [[ $sample_name == $coverage_identifier ]]
        then
            #Use the contigs.fasta and coverage.out files to calculate per contig coverage
            /tmp/genome_back/gen_input_table.py --isbedfiles $contig_file $coverage_file > ./read_mapping/coverage_table_$coverage_identifier.csv
        fi
    done
done
