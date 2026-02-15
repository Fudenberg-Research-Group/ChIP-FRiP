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
    threads: process
    output:
        output1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}.sort.bam",
        output2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_spikein}.sort.bam"
    shell:
        """
        samtools view -h {input} | \
            awk '$1 ~ /^@/ || $3 ~ /^{index_primary}_/' | \
            samtools view -b - | \
            samtools reheader -c 'sed "s/{index_primary}_chr/chr/"' - | \
            samtools sort -@ {threads} -o {output.output1} -
        samtools index -@ {threads} {output.output1}

        samtools view -h {input} | \
            awk '$1 ~ /^@/ || $3 ~ /^{index_spikein}_/' | \
            samtools view -b - | \
            samtools reheader -c 'sed "s/{index_spikein}_chr/chr/"' - | \
            samtools sort -@ {threads} -o {output.output2} -
        samtools index -@ {threads} {output.output2}
        """