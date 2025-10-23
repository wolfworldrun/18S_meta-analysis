# Snakemake file - download raw fastq files from a study on ncbi SRA 
## Adapted from Sarah Hu https://forum.qiime2.org/t/qiime2-snakemake-workflow-tutorial-18s-16s-tag-sequencing/11334

configfile: "config.yaml"

#-----SET VARIABLES-----#

SRAID = config["SRAid"] # accession ID for study in SRA
PROJ = config["proj_name"] #name of project directory as named in raw_data_files folder
INPUTDIR = config["raw_data"] #path to raw fastq files
OUTPUT = config["output"]
VISUALIZATION = config["visualization"]
REFFASTA = config["taxonomy_fasta"]
CLASSIFIER = config["trained_ref_database"]


rule all:
  input:
    # download
    rawreads = directory(INPUTDIR+"reads"),
    # fastqc
    # qiime2 
    q2_import = OUTPUT + PROJ + "-PE-demux.qza",
    q2_primerRM = OUTPUT + PROJ + "-PE-demux-noprimer.qza",
    #q2_primerRM_readthrough = OUTPUT + PROJ + "-PE-demux-noprimer_readthrough.qza,
    q2_merged = OUTPUT + PROJ + "-PE-demux-noprimer-merged.qza",
    q2_qualfilter = OUTPUT + PROJ + "-PE-demux-noprimer-merged-filtered.qza",
    table = OUTPUT + PROJ + "-table.qza",
    rep = OUTPUT + PROJ + "-rep-seqs.qza",
    stats = OUTPUT + PROJ + "-deblur-stats.qza",
    sklearn = OUTPUT + PROJ + "-tax_sklearn.qza",
    table_biom = OUTPUT + "export/table/feature-table.biom",
    table_tsv = OUTPUT + "export/table/feature-table.biom",
    table_tax = OUTPUT + "export/taxonomy.tsv",
    rep_seqs_fasta = OUTPUT + "export/rep-seqs.fasta",
    # visualizations
    q2_qualfilter_stats = OUTPUT + PROJ + "-PE-demux-filter-stats.qza",
    raw = VISUALIZATION + PROJ + "-PE-demux.qzv",
    primer = VISUALIZATION + PROJ + "-PE-demux-noprimer.qzv",
    merged = VISUALIZATION + PROJ + "-PE-demux-noprimer-merged.qzv",
    filtered = VISUALIZATION + PROJ + "-PE-demux-noprimer-merged-filtered.qzv"



rule download_runinfo:
    output:
        runinfo = INPUTDIR + "runinfo.csv"
    conda:
        "envs/parallelfastqdump.yaml"
    shell:
        """
        esearch -db sra -query {config[SRAid]} | efetch -format runinfo > {output.runinfo}
        """
        
rule get_runs:
    input:
        runinfo = INPUTDIR + "runinfo.csv"
    output:
        SRRnumbers = INPUTDIR + "SRR.numbers"
    shell:
        """
        cat {input.runinfo} | cut -d ',' -f 1 > {output.SRRnumbers}
        sed -i '1d' {output.SRRnumbers}
        """
        
rule fasterq_dump:
    input: 
        SRRnumbers = INPUTDIR + "SRR.numbers"
    output:
        rawreads = directory(INPUTDIR+"reads")
    conda:
        "envs/parallelfastqdump.yaml"
    shell:
        """
        mkdir -p {output.rawreads}
        for srr in $(cat {input.SRRnumbers})
            do
            parallel-fastq-dump --tmpdir . --outdir {output.rawreads} --threads 10 --gzip --split-files --sra-id srr
            done
        """
        
rule rename_files_import:
    input:
        rawreads = directory(INPUTDIR+"reads")
    output:
        q2_import = OUTPUT + PROJ + "-PE-demux.qza"
    log:
        OUTPUT + "qiime2/logs/" + PROJ + "_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell: 
        """
        cd {input.rawreads}
        rename 's/_/_00_L001_/g' * /
        rename 's/.fastq.gz/_001.fastq.gz/g' * /
        rename 's/_1/_R1/g' * /
        rename 's/_2/_R2/g' *
        
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input.rawreads} \
            --output-path {output.q2_import} \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
        """

# uncomment if you know there is primer read-through in your reads (read length > amplicon length)
# if read length < amplicon length then there will be no read-through
#rule rm_primers_readthrough:
  #input:
    #q2_import = OUTPUT + PROJ + "-PE-demux.qza"
  #output:
    #q2_primerRM_readthrough = OUTPUT + PROJ + "-PE-demux-noprimer_readthrough.qza"
  #log:
    #OUTPUT + "qiime2/logs/" + PROJ + "_primer_q2.log"
  #conda:
    #"envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  #shell:
    #"""
    #qiime cutadapt trim-paired \
       #--i-demultiplexed-sequences {input.q2_import} \
       #--p-adapter-f {config[adapterF]} \
       #--p-adapter-r {config[adapterR]} \
       #--p-match-adapter-wildcards \
       #--p-match-read-wildcards \
       #--verbose \
       #--o-trimmed-sequences {output.q2_primerRM}
    #"""
       
rule rm_primers:
  input:
    q2_import = OUTPUT + PROJ + "-PE-demux.qza"
    # q2_primerRM_readthrough = OUTPUT + PROJ + "-PE-demux-noprimer_readthrough.qza" # use this one if doing trim step above
  output:
    q2_primerRM = OUTPUT + PROJ + "-PE-demux-noprimer.qza"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "_primer_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  shell:
    """
    qiime cutadapt trim-paired \
       --i-demultiplexed-sequences {input.q2_import} \
       --p-front-f {config[primerF]} \
       --p-front-r {config[primerR]} \
       --p-minimum-length {config[p_min_length]} \
       --p-error-rate {config[primer_err]} \
       --p-discard-untrimmed \
       --p-match-adapter-wildcards \
       --p-match-read-wildcards \
       --verbose \
       --o-trimmed-sequences {output.q2_primerRM}
    """


rule merge:
    input:
        q2_primerRM = OUTPUT + PROJ + "-PE-demux-noprimer.qza"
    output:
        q2_merged = OUTPUT + PROJ + "-PE-demux-noprimer-merged.qza"
    log:
        OUTPUT + "qiime2/logs/" + PROJ + "_merge_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell: 
        "qiime vsearch merge-pairs \
            --i-demultiplexed-seqs {input.q2_primerRM} \
            --o-merged-sequences {output.q2_merged} \
            --p-maxdiffs {config[max_diffs_merge]} \
            --p-allowmergestagger \
            --p-minovlen {config[minovlen_merge]} \
            --p-threads 4 \
            --verbose"

rule quality_filter:
    input:
        q2_merged = OUTPUT + PROJ + "-PE-demux-noprimer-merged.qza"
    output: 
        q2_qualfilter = OUTPUT + PROJ + "-PE-demux-noprimer-merged-filtered.qza",
        q2_qualfilter_stats = OUTPUT + PROJ + "-PE-demux-filter-stats.qza",
        q2_qualfilter_stats_viz = VISUALIZATION + PROJ + "-PE-demux-filter-stats.qzv"
    log:
        OUTPUT + "qiime2/logs/" + PROJ + "filter_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    shell:
        """
        qiime quality-filter q-score \
            --i-demux {input.q2_merged}] \
            --p-min-quality {config[filter_minquality]} \
            --o-filtered-sequences {output.q2_qualfilter} \
            --o-filter-stats {output.q2_qualfilter_stats}
        
        qiime quality-filter visualize-stats \
            --i-filter-stats {output.q2_qualfilter_stats} \
            --o-visualization {output.q2_qualfilter_stats_viz}
        """

rule deblur: 
  input:
    q2_merged = OUTPUT + PROJ + "-PE-demux-noprimer-merged.qza",
    ref_database_fasta = REFFASTA
  output:
    table = OUTPUT + PROJ + "-table.qza",
    rep = OUTPUT + PROJ + "-rep-seqs.qza",
    stats = OUTPUT + PROJ + "-deblur-stats.qza"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "-deblur.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  shell:
    """
    qiime deblur denoise-other \
        --i-demultiplexed-seqs {input.q2_merged} \
        --i-reference-seqs {input.ref_database_fasta} \
        --p-trim-length {config[p_trim_length]} \
        --p-no-hashed-feature-ids \
        --o-representative-sequences {output.rep} \
        --o-table {output.table} \
        --p-sample-stats \
        --p-jobs-to-start 16 \
        --o-stats {output.stats} \
        --p-min-reads=1
    """

rule get_stats:
  input:
    q2_import = OUTPUT + PROJ + "-PE-demux.qza",
    q2_primerRM = OUTPUT + PROJ + "-PE-demux-noprimer.qza",
    q2_merged = OUTPUT + PROJ + "-PE-demux-noprimer-merged.qza",
    q2_qualfilter = OUTPUT + PROJ + "-PE-demux-noprimer-merged-filtered.qza",
    stats = OUTPUT + PROJ + "-deblur-stats.qza"
  output:
    raw = VISUALIZATION + PROJ + "-PE-demux.qzv",
    primer = VISUALIZATION + PROJ + "-PE-demux-noprimer.qzv",
    merged = VISUALIZATION + PROJ + "-PE-demux-noprimer-merged.qzv",
    filtered = VISUALIZATION + PROJ + "-PE-demux-noprimer-merged-filtered.qzv",
    deblur_stats = VISUALIZATION + PROJ + "-deblur-stats.qzv"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "_getviz_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  shell:
    """
    qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
    qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    qiime demux summarize --i-data {input.q2_merged} --o-visualization {output.merged}
    qiime demux summarize --i-data {input.q2_qualfilter} --o-visualization {output.filtered}
    qiime demux summarize --i-data {input.stats} --o-visualization {output.stats}
    """

rule assign_tax:
  input:
    rep = OUTPUT + PROJ + "-rep-seqs.qza",
    classifier = CLASSIFIER 
  output:
    sklearn = OUTPUT + PROJ + "-tax_sklearn.qza"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "-sklearn_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  shell:
    "qiime feature-classifier classify-sklearn \
	--i-classifier {input.classifier} \
	--i-reads {input.rep} \
	--o-classification {output.sklearn}"

rule gen_table:
  input:
    table = OUTPUT + PROJ + "-table.qza"
  output:
    table_biom = OUTPUT + "export/table/feature-table.biom"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "_exportBIOM_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  params:
    directory(OUTPUT + "export/table")
  shell:
    "qiime tools export --input-path {input.table_biom} --output-path {params}"

rule convert:
  input:
    table_biom = OUTPUT + "export/table/feature-table.biom"
  output:
    table_tsv = OUTPUT + "export/table/feature-table.biom"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "_exportTSV_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  shell:
    "biom convert -i {input} -o {output} --to-tsv"

rule gen_tax:
  input:
    sklearn = OUTPUT + PROJ + "-tax_sklearn.qza"
  output:
    table_tax = OUTPUT + "export/taxonomy.tsv"
  log:
    OUTPUT + "qiime2/logs/" + PROJ + "_exportTAXTSV_q2.log"
  conda:
    "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
  params:
    directory(OUTPUT + "export/")
  shell:
    "qiime tools export --input-path {input.sklearn} --output-path {params}"
    
rule convert_seqs:
    input:
        rep = OUTPUT + PROJ + "-rep-seqs.qza"
    output: 
        rep_seqs_fasta = OUTPUT + "export/rep-seqs.fasta"
    log:
        OUTPUT + "qiime2/logs/" + PROJ + "_exportFASTA_q2.log"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        directory(OUTPUT + "export/")
    conda:
        "envs/qiime2-2019.4.yaml"
    shell: "qiime tools export --input-path {input.rep} --output-path {params}"

