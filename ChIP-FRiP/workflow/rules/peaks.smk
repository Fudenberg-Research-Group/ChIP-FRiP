##############################
# Calling peaks
##############################

rule call_peaks:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{{type}}.bam"
    params:
        output_prefix = f"{{sample_name}}.q{quality}.{{type}}", 
        output_dir = f"{{pathway_to_folder}}/{{sample_name}}/",
        fastq_files = lambda wildcards: get_fastq_paths(wildcards) # q: can this be done as get_fastq_paths(wildcards, chip_fastqs_input)?
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{{type}}_peaks.narrowPeak"
    shell:
        """
        files=({params.fastq_files})
        num_files=${{#files[@]}}

        if [ "$num_files" -eq 1 ]; then
            macs2 callpeak {broad_peaks_option} -t {input} -n {params.output_prefix} --outdir {params.output_dir} --gsize {effective_genome_size} -q 0.05
        elif [ "$num_files" -eq 2 ]; then
            macs2 callpeak {broad_peaks_option} -f BAMPE -t {input} -n {params.output_prefix} --outdir {params.output_dir} --gsize {effective_genome_size} -q 0.05
        fi
        """

rule call_peaks_with_input:
    input:
        bam = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{{type}}.bam",
        chip_input_file = lambda wildcards: get_chipseq_input_path(wildcards),
    params:
        output_prefix = f"{{sample_name}}.q{quality}.{{type}}.withinput", 
        output_dir = f"{{pathway_to_folder}}/{{sample_name}}/",
        fastq_files = lambda wildcards: get_fastq_paths(wildcards)
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{{type}}.withinput_peaks.narrowPeak"
    shell:
        """
        files=({params.fastq_files})
        num_files=${{#files[@]}}

        if [ "$num_files" -eq 1 ]; then
            macs2 callpeak {broad_peaks_option} -t {input.bam} -c {input.chip_input_file} -n {params.output_prefix} --outdir {params.output_dir} --gsize {effective_genome_size} -q 0.05
        elif [ "$num_files" -eq 2 ]; then
            macs2 callpeak {broad_peaks_option} -f BAMPE -t {input.bam} -c {input.chip_input_file} -n {params.output_prefix} --outdir {params.output_dir} --gsize {effective_genome_size} -q 0.05
        fi
        """