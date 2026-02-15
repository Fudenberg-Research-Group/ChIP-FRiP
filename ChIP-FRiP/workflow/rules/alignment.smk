##############################
# Alignment
##############################

# automatic trim and alignment
rule bowtie2_map:
    input:
        files = get_fastq_paths
    params:
        trimmed_fastq_path1 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.trimmed.fastq",
        trimmed_fastq_path2 = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}_2.trimmed.fastq"
    wildcard_constraints:
        sample_name = "[a-zA-Z0-9_-]+"
    threads: process
    output:
        temp(f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.sam")
    shell:
        """
        files=({input.files})
        num_files=${{#files[@]}}

        if [ "$num_files" -eq 1 ]; then
            fastp --thread {process} -i {input.files} -o {params.trimmed_fastq_path1}
            bowtie2 -p {process} -x {index} -U {params.trimmed_fastq_path1} -S {output}
            rm {params.trimmed_fastq_path1}
        elif [ "$num_files" -eq 2 ]; then
            fastp --thread {process} -i {input.files[0]} -I {input.files[1]} -o {params.trimmed_fastq_path1} -O {params.trimmed_fastq_path2}
            bowtie2 -p {process} -x {index} -1 {params.trimmed_fastq_path1} -2 {params.trimmed_fastq_path2}  -S {output}
            rm {params.trimmed_fastq_path1} {params.trimmed_fastq_path2}
        fi
        """

# skip alignments with MAPQ smaller than INT {quality} and type of reads include in argument -F 1804, which includes read unmapped, mate unmapped, not primary alignment, Read fails platform/vendor quality checks, Read is PCR or optical duplicate
rule samtools_filter:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.sam"
    wildcard_constraints:
        sample_name="[^.]+"
    threads: process
    output:
        temp(f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}{{end_type}}.q{quality}.bam")
    shell:
        "samtools view -F 1804 --threads {process} -h -q {quality} {input} > {output}"

# This is for paired-end fastq file only. Sort by reads name and then correctly remove secondary and unmapped reads
rule samtools_fixmate:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.pairedend.q{quality}.bam"
    params:
        sorted_name = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.sorted.bam"
    wildcard_constraints:
        sample_name="[^.]+"
    threads: process
    output:
        temp(f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.bam")
    shell:
        """
        samtools sort -n -@{process} {input} -o {params.sorted_name}
        samtools fixmate -m -r {params.sorted_name} {output}
        rm {params.sorted_name}
        """

# back to sorting by genomic coordinates for markdup
rule samtools_sort:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.bam"
    wildcard_constraints:
        sample_name="[^.]+"
    threads: process
    output:
        temp(f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.sort.bam")
    shell:
        """
        samtools sort -@{process} {input} -o {output}
        samtools index --threads {process} {output}
        """

# remove duplicate reads
rule samtools_markdup:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.sort.bam"
    params:
        stats_file_path = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.markdup.stats"
    wildcard_constraints:
        sample_name="[^.]+"
    threads: process
    output:
        bam = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.dedup.bam"
    shell:
        """
        samtools markdup -f {params.stats_file_path} -r -d {distance} {input} {output.bam}
        samtools index --threads {process} {output.bam}
        """