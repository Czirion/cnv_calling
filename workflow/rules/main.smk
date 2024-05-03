rule links:
    input:
        REFS / "{lineage}.fasta"
    output:
        REFS_DIR / "{lineage}" / "{lineage}.fasta"
    log:
        "logs/references/links/{lineage}.log"
    shell:
        "ln -s -r -f {input} {output} &> {log}"

# Run RepeatModeler for each reference genome
rule repeat_modeler:
    input:
        rules.links.output
    output:
        known = REFS_DIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        unknown = REFS_DIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    params:
        repdir = "repeats"
    threads:
        config["repeats"]["threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/references/repeats/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {input} {params.repdir} &> {log}"

# Run RepeatMasker for each reference genome. Obtain a BED file with the location of the repetitive sequences
rule repeat_masker:
    input:
        database = config["repeats"]["repeats_database"],
        fasta = rules.links.output,
        known = rules.repeat_modeler.output.known,
        unknown = rules.repeat_modeler.output.unknown
    output:
        REFS_DIR / "{lineage}" / "repeats" / "{lineage}_repeats.bed"
    threads:
        config["repeats"]["threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/references/repeats/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {input.database} {input.fasta} {input.known} {input.unknown} {output} &> {log}"


# Get average read depth per window of each sample 
rule mosdepth:
    input:
        bam = BAMS / "{sample}.bam",
        bai = BAMS / "{sample}.bam.bai"
    output:
        bed = SAMPLES_DIR / "mosdepth" / "{sample}" / "depth.regions.bed.gz"
    params:
        window_size = config["mosdepth"]["window_size"],
        min_mapq = config["mosdepth"]["min_mapq"],
        extra = config["mosdepth"]["extra"],
        outdir = SAMPLES_DIR / "mosdepth"
    threads:
       config["mosdepth"]["threads"]   
    conda:
        "../envs/depth.yaml"
    log:
        "logs/samples/mosdepth/mosdepth_good_{sample}.log"
    shell:
        "mosdepth -n "
        "--by {params.window_size} "
        "--mapq {params.min_mapq} "
        "-t {threads} "
        "{params.extra} "
        "{params.outdir}/{wildcards.sample}/depth "
        "{input.bam} "
        "&> {log} "

def cnv_calling_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": SAMPLES_DIR / "mosdepth" / s["sample"] / "depth.regions.bed.gz" ,
        "repeats": REFS_DIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule cnv_calling:
    input:
        unpack(cnv_calling_input)
    output:
        chromosome = SAMPLES_DIR / "cnv" / "{sample}" / "chromosome_depth.tsv",
        window = SAMPLES_DIR / "cnv" / "{sample}" / "window_depth.tsv",
        cnv = SAMPLES_DIR / "cnv" / "{sample}" / "copy_number_variants.tsv"
    params:
        window = config["mosdepth"]["window_size"],
        smooth = config["cnv_calling"]["smoothing_size"],
        depth_threshold = config["cnv_calling"]["depth_threshold"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samples/mosdepth/depth_{sample}.log"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh "
        "-di {input.depth} "
        "-ri {input.repeats} "
        "-co {output.chromosome} "
        "-wo {output.window} "
        "-vo {output.cnv} "
        "-np {wildcards.sample} "
        "-wp {params.window} "
        "-sp {params.smooth} "
        "-dp {params.depth_threshold} "
        "&> {log}"

rule dataset_cnv:
    input:
        expand(SAMPLES_DIR / "cnv" / "{sample}" / "copy_number_variants.tsv", sample = SAMPLES)
    output:
        DATASET_DIR / "copy_number_variants_dataset.tsv"
    run:
        cnv = pd.concat([pd.read_csv(f, sep="\t") for f in input])
        cnv.to_csv(output[0], sep="\t", index=False)
