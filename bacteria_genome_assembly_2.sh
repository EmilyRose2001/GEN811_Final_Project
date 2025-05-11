#! /bin/bash

source activate genomics

#Non-Target Contig Removal

#Create a directory for all output files of blobtools
mkdir blobtools_contig_removal

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Extract the .bam files
    for bam_file in ./read_mapping/sorted_mapped_*.bam
    do
        bam_name=$(basename $bam_file)
        bam_identifier=${bam_name:14:11}

        for blast_dir in ./blast_*
        do
            blast_dir_name=$(basename $blast_dir)
            blast_identifier=${blast_dir_name:6:11}

            #Extract the blast output megablast.out file
            blast_file=$blast_dir/contigs.fasta.vs.nt.cul5.maxt10.1e5.megablast.out

            if [[ ($sample_name == $bam_identifier) && ($bam_identifier == $blast_identifier) ]]
            then
                #Create a blobtools lookup table with the extracted contigs.fasta, .bam, and megablast.out files
                blobtools create -i $contig_file -b $bam_file -t $blast_file -o ./blobtools_contig_removal/blob_out_$sample_name

                #Create a blobtools output table useing the output lookup table
                blobtools view -i ./blobtools_contig_removal/blob_out_$sample_name.blobDB.json -r all -o ./blobtools_contig_removal/blob_taxonomy_$sample_name

                #Create a blobtools plot using the output lookup table
                blobtools plot -i ./blobtools_contig_removal/blob_out_$sample_name.blobDB.json -r genus

                #Output plot files are moved into the blobtools directory. The other files were created in the blobtools directory and don't need to be moved
                mv blob_out_* ./blobtools_contig_removal
            fi
        done
    done
done

#Filter Genome Assembly

#Create a directory to store all filtered files
mkdir ./filtered_assembly

for tax_file in ./blobtools_contig_removal/blob_taxonomy*.txt
do
    #Copy the taxonomy file outputted from blobtools to the created filtered_assembly directory for each sample
    cp $tax_file ./filtered_assembly
done

#Extract the copied taxonomy blob files
for tax_file in ./filtered_assembly/blob_taxonomy*.txt
do
    tax_file_name=$(basename $tax_file)
    file_identifier=${tax_file_name:14:11}

    #Filter by 500 bp length
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' > ./filtered_assembly/length_filter_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | wc -l >> ./filtered_assembly/length_filter_$file_identifier.txt

    #Check what contigs are lost
    grep -v '#' $tax_file | awk -F'\t' '$2 < 500' > ./filtered_assembly/removed_length_filter_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 < 500' | wc -l >> ./filtered_assembly/removed_length_filter_$file_identifier.txt

    #Filter by coverage starting at 5
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5' > ./filtered_assembly/coverage_filter_5_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5' | wc -l >> ./filtered_assembly/coverage_filter_5_$file_identifier.txt

    #Check what contigs are lost
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 5' > ./filtered_assembly/removed_coverage_filter_5_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 5' | wc -l >> ./filtered_assembly/removed_coverage_filter_5_$file_identifier.txt

    #Filter by 15 coverage
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5' > ./filtered_assembly/coverage_filter_15_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5' | wc -l >> ./filtered_assembly/coverage_filter_15_$file_identifier.txt

    #Check what contigs are lost
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 5' > ./filtered_assembly/removed_coverage_filter_15_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 5' | wc -l >> ./filtered_assembly/removed_coverage_filter_15_$file_identifier.txt

    #Filter by 20 coverage
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' > ./filtered_assembly/coverage_filter_20_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | wc -l >> ./filtered_assembly/coverage_filter_20_$file_identifier.txt

    #Check what contigs are lost
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 20' > ./filtered_assembly/removed_coverage_filter_20_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 20' | wc -l >> ./filtered_assembly/removed_coverage_filter_20_$file_identifier.txt

    #Filter by 25 coverage
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 25' > ./filtered_assembly/coverage_filter_25_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 25' | wc -l >> ./filtered_assembly/coverage_filter_25_$file_identifier.txt

    #Check what contigs are lost
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 25' > ./filtered_assembly/removed_coverage_filter_25_$file_identifier.txt
    grep -v '#' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 25' | wc -l >> ./filtered_assembly/removed_coverage_filter_25_$file_identifier.txt

    #Create list of contigs to keep
    grep -v '##' $tax_file | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_len500_cov20_$file_identifier.txt
done

for spades_dir in ./spade_assembly_default_*
do
    spades_dir_name=$(basename $spades_dir)
    sample_name=${spades_dir_name:23:11}

    #Extract the contigs.fasta file
    contig_file=$spades_dir/contigs.fasta

    #Extract the file that contains the list of contigs that will be kept
    for list_contig in ./list_of_contigs_to_keep_len500_cov20_*.txt
    do
        list_contig_name=$(basename $list_contig)
        list_name=${list_contig_name:37:11}

        if [[ $sample_name == $list_name ]]
        then
            #Create a filtered genome based on the contigs that are being kept
            filter_contigs_by_list.py $contig_file $list_contig filtered_$sample_name.fasta
        fi
    done
done

#Average Coverage Using Blobtools

#Extract the taxonomy table outputted by blobtools
for tax_file in ./filtered_assembly/blob_taxonomy*.txt
do
    tax_file_name=$(basename $tax_file)
    file_identifier=${tax_file_name:14:11}

    #Extract the file that contains the list of contigs that will be kept
    for list_contig in ./list_of_contigs_to_keep_len500_cov20_*.txt
    do
        list_contig_name=$(basename $list_contig)
        list_name=${list_contig_name:37:11}

        if [[ $file_identifier == $list_name ]]
        then
            #Calculate the average coverage using blobtools and save the output into a file named after the sample
            grep -f $list_contig $tax_file | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}' > average_coverage_filtered_$file_identifier.txt
        fi
    done
done

#Check for contamination by BLASTing against UniVec

#Download the UniVec database
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"

for genome_fasta in ./filtered_*.fasta
do
    genome_fasta_name=$(basename $genome_fasta)
    genome_name=${genome_fasta_name:9:11}

    #BLAST all filtered genome fasta files against UniVec
    blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query $genome_fasta -subject UniVec -outfmt 6 -out genome_vs_univec_$genome_name.6
done
