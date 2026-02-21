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
        primary_reads=$(samtools idxstats {input} | grep -E '^{primary_assembly}_' | \
                            grep -v -E '(_random|_alt|_fix|_decoy|chrM|chrY|chrUn)' | \
                                awk '{{sum += $3}} END {{print sum}}')
        spikein_reads=$(samtools idxstats {input} | grep -E '^{spikein_assembly}_' | \
                            grep -v -E '(_random|_alt|_fix|_decoy|chrM|chrY|chrUn)' | \
                                awk '{{sum += $3}} END {{print sum}}')
        echo -e "{primary_assembly}_reads=$primary_reads" >> {output} 
        echo "{spikein_assembly}_reads=$spikein_reads" >> {output} 
        echo "ratio of {primary_assembly} to {spikein_assembly} reads is" >> {output} 
        echo "scale=2; $primary_reads/$spikein_reads" | bc >> {output} 
        """

# Use grep to create a bam file with only one species reads
rule separate_reads:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.dedup.bam"
    threads: num_processes
    output:
        output1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{primary_assembly}.dedup.bam",
        output2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{spikein_assembly}.dedup.bam"
    shell:
        """
        samtools view -h {input} | \
            awk '$1 ~ /^@/ || $3 ~ /^{primary_assembly}_/' | \
            samtools view -b - | \
            samtools reheader -c 'sed "s/{primary_assembly}_chr/chr/"' - | \
            samtools sort -@ {threads} -o {output.output1} -
        samtools index -@ {threads} {output.output1}

        samtools view -h {input} | \
            awk '$1 ~ /^@/ || $3 ~ /^{spikein_assembly}_/' | \
            samtools view -b - | \
            samtools reheader -c 'sed "s/{spikein_assembly}_chr/chr/"' - | \
            samtools sort -@ {threads} -o {output.output2} -
        samtools index -@ {threads} {output.output2}
        """