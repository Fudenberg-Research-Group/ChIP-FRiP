############################################################
#  Create bigwig files
############################################################

# create an unscaled bigwig file
rule bw:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.dedup.bam"
    wildcard_constraints:
        sample_name="[^.]+"
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.bw"
    params:
        blacklist=lambda wildcards: f'--blackListFileName {config["bamcoverage_params"]["blacklist"]}' if config["bamcoverage_params"]["blacklist"] else "",
        extendReads=lambda wildcards: f'--extendReads {config["bamcoverage_params"]["extendReads"]}' if config["bamcoverage_params"]["extendReads"] > 0 else ""
    shell:
        "bamCoverage {params.blacklist} {params.extendReads} -b {input} -o {output} -of bigwig --binSize {binsize} --effectiveGenomeSize {effective_genome_size}"

# create an unscaled bigwig file for spike-in
rule spikein_bw:
    input:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}.sort.bam"
    wildcard_constraints:
        sample_name="[^.]+"
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.{index_primary}.bw"
    params:
        blacklist=lambda wildcards: f'--blackListFileName {config["bamcoverage_params"]["blacklist"]}' if config["bamcoverage_params"]["blacklist"] else "",
        extendReads=lambda wildcards: f'--extendReads {config["bamcoverage_params"]["extendReads"]}' if config["bamcoverage_params"]["extendReads"] > 0 else ""
    shell:
        "bamCoverage {params.blacklist} {params.extendReads} -b {input} -o {output} -of bigwig --binSize {binsize} --effectiveGenomeSize {effective_genome_size}"

# Spike-in Normalization
rule rescaling:
    input:
        chip_stats = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.spikein.stats",
        chip_input_stats = lambda wildcards: get_chipseq_input_spikestats_path(wildcards),
        chip_bam = f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.q{quality}.{index_primary}.sort.bam"
    output:
        f"{{pathway_to_folder}}/{{sample_name}}/{{sample_name}}.{index_primary}.rescale.bw"
    params:
        blacklist=lambda wildcards: f'--blackListFileName {config["bamcoverage_params"]["blacklist"]}' if config["bamcoverage_params"]["blacklist"] else "",
        extendReads=lambda wildcards: f'--extendReads {config["bamcoverage_params"]["extendReads"]}' if config["bamcoverage_params"]["extendReads"] > 0 else ""
    shell:
        """
        primary_reads=$(awk -F'=' '/{index_primary}_reads/{{print $2}}' {input.chip_stats}) 
        spikein_reads=$(awk -F'=' '/{index_spikein}_reads/{{print $2}}' {input.chip_stats}) 
        chip_input_primary_reads=$(awk -F'=' '/{index_primary}_reads/{{print $2}}' {input.chip_input_stats})  
        chip_input_spikein_reads=$(awk -F'=' '/{index_spikein}_reads/{{print $2}}' {input.chip_input_stats})
        factor=`echo "scale=20; $chip_input_spikein_reads / $chip_input_primary_reads / $spikein_reads * 15000000" | bc`
        echo "Scaling_factor=$factor" >> {input.chip_stats}
        bamCoverage {params.blacklist} {params.extendReads} -b {input.chip_bam} -o {output} -of bigwig --binSize {binsize} --scaleFactor $factor --effectiveGenomeSize {effective_genome_size}
        """