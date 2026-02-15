############################################################
#  Spike-in ChIP protocol
############################################################

# Count the number of reads that map to each genome and print the ratio
rule spikein_stats:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.dedup.bam"
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.spikein.stats"
    shell:
        """
        primary_reads=$(samtools idxstats {input} | grep -E '^{index_primary}_' | \
                            grep -v -E '(_random|_alt|_fix|_decoy|chrM|chrY|chrUn)' | \
                                awk '{{sum += $3}} END {{print sum}}')
        spikein_reads=$(samtools idxstats {input} | grep -E '^{index_spikein}_' | \
                            grep -v -E '(_random|_alt|_fix|_decoy|chrM|chrY|chrUn)' | \
                                awk '{{sum += $3}} END {{print sum}}')
        echo -e "{index_primary}_reads=$primary_reads" >> {output} 
        echo "{index_spikein}_reads=$spikein_reads" >> {output} 
        echo "ratio of {index_primary} to {index_spikein} reads is" >> {output} 
        echo "scale=2; $primary_reads/$spikein_reads" | bc >> {output} 
        """

# Use grep to create a bam file with only one species reads
rule separate_reads:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.dedup.bam"
    params:
        primary_bam_tmp1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}1.bam",
        primary_bam_tmp2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}2.bam",
        primary_header_tmp = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}.newHeader.txt",
        spikein_bam_tmp1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_spikein}1.bam",
        spikein_bam_tmp2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_spikein}2.bam",
        spikein_header_tmp = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_spikein}.newHeader.txt",
        
    threads: process
    output:
        output1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}.sort.bam",
        output2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_spikein}.sort.bam"
    shell:
        """
        samtools view -h {input} | grep '.*{index_primary}.*' | samtools view -b -> {params.primary_bam_tmp1}
        samtools view -H {params.primary_bam_tmp1} | sed 's/{index_primary}_chr/chr/' > {params.primary_header_tmp}
        samtools reheader {params.primary_header_tmp} {params.primary_bam_tmp1} > {params.primary_bam_tmp2}
        samtools sort {params.primary_bam_tmp2} -o {output.output1}
        samtools index --threads {process} {output.output1}
        rm {params.primary_bam_tmp1} {params.primary_bam_tmp2} {params.primary_header_tmp}

        samtools view -h {input} | grep '.*{index_spikein}.*' | samtools view -b -> {params.spikein_bam_tmp1}
        samtools view -H {params.spikein_bam_tmp1} | sed 's/{index_spikein}_chr/chr/' > {params.spikein_header_tmp}
        samtools reheader {params.spikein_header_tmp} {params.spikein_bam_tmp1} > {params.spikein_bam_tmp2}
        samtools sort {params.spikein_bam_tmp2} -o {output.output2}
        samtools index --threads {process} {output.output2}
        rm {params.spikein_bam_tmp1} {params.spikein_bam_tmp2} {params.spikein_header_tmp}
        """