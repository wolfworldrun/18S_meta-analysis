# Snakemake file - batch download and process multiple SRA studies
## Modified for batch processing multiple projects
## Original adapted from Sarah Hu https://forum.qiime2.org/t/qiime2-snakemake-workflow-tutorial-18s-16s-tag-sequencing/11334

configfile: "config.yml"

#-----SET VARIABLES-----#

# Get list of projects to process
PROJECTS = config["projects"]

# Shared paths
RAW_DATA = config["raw_data"]
OUTPUT = config["output"]
VISUALIZATION = config["visualization"]
REFFASTA = config["taxonomy_fasta"]
CLASSIFIER = config["trained_ref_database"]

# Helper function to get project-specific or global settings
def get_setting(wildcards, key, default=None):
    """
    Get a setting for a specific project.
    First checks project_settings[project][key], then falls back to config[key]
    """
    if "project_settings" in config:
        if wildcards.project in config["project_settings"]:
            if key in config["project_settings"][wildcards.project]:
                return config["project_settings"][wildcards.project][key]
    return config.get(key, default)


rule all:
    input:
        # Download and import - for all projects
        rawreads = expand(RAW_DATA + "{project}/reads", project=PROJECTS),
        q2_import = expand(OUTPUT + "{project}/{project}-PE-demux.qza", project=PROJECTS),
        
        # Processing - for all projects
        q2_primerRM = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer.qza", project=PROJECTS),
        q2_merged = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza", project=PROJECTS),
        q2_qualfilter = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza", project=PROJECTS),
        
        # Core outputs - for all projects
        table = expand(OUTPUT + "{project}/{project}-table.qza", project=PROJECTS),
        rep = expand(OUTPUT + "{project}/{project}-rep-seqs.qza", project=PROJECTS),
        stats = expand(OUTPUT + "{project}/{project}-deblur-stats.qza", project=PROJECTS),
        sklearn = expand(OUTPUT + "{project}/{project}-tax_sklearn.qza", project=PROJECTS),
        
        # Exported files - for all projects
        table_biom = expand(OUTPUT + "{project}/export/table/feature-table.biom", project=PROJECTS),
        table_tsv = expand(OUTPUT + "{project}/export/table/feature-table.tsv", project=PROJECTS),
        table_tax = expand(OUTPUT + "{project}/export/taxonomy.tsv", project=PROJECTS),
        rep_seqs_fasta = expand(OUTPUT + "{project}/export/rep-seqs.fasta", project=PROJECTS),
        
        # Visualizations - for all projects
        raw = expand(VISUALIZATION + "{project}/{project}-PE-demux.qzv", project=PROJECTS),
        primer = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer.qzv", project=PROJECTS),
        merged = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged.qzv", project=PROJECTS),
        filtered = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged-filtered.qzv", project=PROJECTS),
        deblur_stats = expand(VISUALIZATION + "{project}/{project}-deblur-stats.qzv", project=PROJECTS)


rule download_runinfo:
    output:
        runinfo = RAW_DATA + "{project}/runinfo.csv"
    conda:
        "envs/parallelfastqdump.yaml"
    params:
        sra_id = lambda wildcards: get_setting(wildcards, "SRAid", wildcards.project)
    shell:
        """
        mkdir -p $(dirname {output.runinfo})
        esearch -db sra -query {params.sra_id} | efetch -format runinfo > {output.runinfo}
        """
        

rule get_runs:
    input:
        runinfo = RAW_DATA + "{project}/runinfo.csv"
    output:
        SRRnumbers = RAW_DATA + "{project}/SRR.numbers"
    shell:
        """
        cat {input.runinfo} | cut -d ',' -f 1 > {output.SRRnumbers}
        sed -i '1d' {output.SRRnumbers}
        """
        

rule fasterq_dump:
    input: 
        SRRnumbers = RAW_DATA + "{project}/SRR.numbers"
    output:
        rawreads = directory(RAW_DATA + "{project}/reads")
    conda:
        "envs/parallelfastqdump.yaml"
    threads: 10
    shell:
        """
        mkdir -p {output.rawreads}
        for srr in $(cat {input.SRRnumbers})
        do
            parallel-fastq-dump --tmpdir . --outdir {output.rawreads} --threads {threads} --gzip --split-files --sra-id $srr
        done
        """
        

rule rename_files_import:
    input:
        rawreads = RAW_DATA + "{project}/reads"
    output:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza"
    log:
        OUTPUT + "{project}/logs/{project}_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell: 
        """
        cd {input.rawreads}
        rename 's/_/_00_L001_/g' * 2>/dev/null || true
        rename 's/.fastq.gz/_001.fastq.gz/g' * 2>/dev/null || true
        rename 's/_1/_R1/g' * 2>/dev/null || true
        rename 's/_2/_R2/g' * 2>/dev/null || true
        
        mkdir -p $(dirname {output.q2_import})
        mkdir -p $(dirname {log})
        
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input.rawreads} \
            --output-path {output.q2_import} \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
        """


rule rm_primers:
    input:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza"
    output:
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza"
    log:
        OUTPUT + "{project}/logs/{project}_primer_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        primerF = lambda wildcards: get_setting(wildcards, "primerF"),
        primerR = lambda wildcards: get_setting(wildcards, "primerR"),
        min_length = lambda wildcards: get_setting(wildcards, "p_min_length", 10),
        primer_err = lambda wildcards: get_setting(wildcards, "primer_err", 0.1)
    shell:
        """
        mkdir -p $(dirname {log})
        
        qiime cutadapt trim-paired \
           --i-demultiplexed-sequences {input.q2_import} \
           --p-front-f {params.primerF} \
           --p-front-r {params.primerR} \
           --p-minimum-length {params.min_length} \
           --p-error-rate {params.primer_err} \
           --p-discard-untrimmed \
           --p-match-adapter-wildcards \
           --p-match-read-wildcards \
           --verbose \
           --o-trimmed-sequences {output.q2_primerRM}
        """


rule merge:
    input:
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza"
    output:
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza"
    log:
        OUTPUT + "{project}/logs/{project}_merge_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        max_diffs = lambda wildcards: get_setting(wildcards, "max_diffs_merge", 40),
        minovlen = lambda wildcards: get_setting(wildcards, "minovlen_merge", 40)
    threads: 4
    shell:
        """
        mkdir -p $(dirname {log})
        
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs {input.q2_primerRM} \
            --o-merged-sequences {output.q2_merged} \
            --p-maxdiffs {params.max_diffs} \
            --p-allowmergestagger \
            --p-minovlen {params.minovlen} \
            --p-threads {threads} \
            --verbose
        """


rule quality_filter:
    input:
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza"
    output: 
        q2_qualfilter = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza",
        q2_qualfilter_stats = OUTPUT + "{project}/{project}-PE-demux-filter-stats.qza",
        q2_qualfilter_stats_viz = VISUALIZATION + "{project}/{project}-PE-demux-filter-stats.qzv"
    log:
        OUTPUT + "{project}/logs/{project}_filter_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        min_quality = lambda wildcards: get_setting(wildcards, "filter_minquality", 20)
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.q2_qualfilter_stats_viz})
        
        qiime quality-filter q-score \
            --i-demux {input.q2_merged} \
            --p-min-quality {params.min_quality} \
            --o-filtered-sequences {output.q2_qualfilter} \
            --o-filter-stats {output.q2_qualfilter_stats}
        
        qiime metadata tabulate \
            --m-input-file {output.q2_qualfilter_stats} \
            --o-visualization {output.q2_qualfilter_stats_viz}
        """


rule deblur: 
    input:
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza",
        ref_database_fasta = REFFASTA
    output:
        table = OUTPUT + "{project}/{project}-table.qza",
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        stats = OUTPUT + "{project}/{project}-deblur-stats.qza"
    log:
        OUTPUT + "{project}/logs/{project}-deblur.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        trim_length = lambda wildcards: get_setting(wildcards, "p_trim_length", 325)
    threads: 16
    shell:
        """
        mkdir -p $(dirname {log})
        
        qiime deblur denoise-other \
            --i-demultiplexed-seqs {input.q2_merged} \
            --i-reference-seqs {input.ref_database_fasta} \
            --p-trim-length {params.trim_length} \
            --p-no-hashed-feature-ids \
            --o-representative-sequences {output.rep} \
            --o-table {output.table} \
            --p-sample-stats \
            --p-jobs-to-start {threads} \
            --o-stats {output.stats} \
            --p-min-reads 1
        """


rule get_stats:
    input:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza",
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza",
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza",
        q2_qualfilter = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza",
        stats = OUTPUT + "{project}/{project}-deblur-stats.qza"
    output:
        raw = VISUALIZATION + "{project}/{project}-PE-demux.qzv",
        primer = VISUALIZATION + "{project}/{project}-PE-demux-noprimer.qzv",
        merged = VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged.qzv",
        filtered = VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged-filtered.qzv",
        deblur_stats = VISUALIZATION + "{project}/{project}-deblur-stats.qzv"
    log:
        OUTPUT + "{project}/logs/{project}_getviz_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.raw})
        
        qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
        qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
        qiime demux summarize --i-data {input.q2_merged} --o-visualization {output.merged}
        qiime demux summarize --i-data {input.q2_qualfilter} --o-visualization {output.filtered}
        qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output.deblur_stats}
        """


rule assign_tax:
    input:
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        classifier = CLASSIFIER 
    output:
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    log:
        OUTPUT + "{project}/logs/{project}-sklearn_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.rep} \
            --o-classification {output.sklearn}
        """


rule gen_table:
    input:
        table = OUTPUT + "{project}/{project}-table.qza"
    output:
        table_biom = OUTPUT + "{project}/export/table/feature-table.biom"
    log:
        OUTPUT + "{project}/logs/{project}_exportBIOM_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        outdir = OUTPUT + "{project}/export/table"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        
        qiime tools export --input-path {input.table} --output-path {params.outdir}
        """


rule convert:
    input:
        table_biom = OUTPUT + "{project}/export/table/feature-table.biom"
    output:
        table_tsv = OUTPUT + "{project}/export/table/feature-table.tsv"
    log:
        OUTPUT + "{project}/logs/{project}_exportTSV_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        
        biom convert -i {input.table_biom} -o {output.table_tsv} --to-tsv
        """


rule gen_tax:
    input:
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    output:
        table_tax = OUTPUT + "{project}/export/taxonomy.tsv"
    log:
        OUTPUT + "{project}/logs/{project}_exportTAXTSV_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        outdir = OUTPUT + "{project}/export/"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        
        qiime tools export --input-path {input.sklearn} --output-path {params.outdir}
        """
    

rule convert_seqs:
    input:
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza"
    output: 
        rep_seqs_fasta = OUTPUT + "{project}/export/rep-seqs.fasta"
    log:
        OUTPUT + "{project}/logs/{project}_exportFASTA_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        outdir = OUTPUT + "{project}/export/"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        
        qiime tools export --input-path {input.rep} --output-path {params.outdir}
        """
