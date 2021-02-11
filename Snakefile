rule correct_batch:
    input:
        "batch_correction/data/{dataset}/{dataset}_expr.txt.gz",
        "batch_correction/data/{dataset}/{dataset}_metadata.txt.gz",
    params:
        output_prefix="{dataset}",
        output_dir="batch_correction/{software}/output/"
    output:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    shell:
        "Rscript batch_correction/{wildcards.software}/run_{wildcards.software}.R {input} '{params.output_prefix}' '{params.output_dir}'| tee '{params.output_dir}/{wildcards.dataset}_log'"
       
#nPC_ft=lambda wildcards: 50 if wildcards.software == "SIMBA" else 20
nPC_ft = 20

rule evaluate_ARI:
    input:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_ARI_PC{nPCs}.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/",
    shell:
        "Rscript batch_correction/metrics/run_ARISampled.R {wildcards.software} {input} '{params.output_dir}' {wildcards.nPCs}"

rule evaluate_ASW:
    input:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_ASW_PC{nPCs}.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/",
    shell:
        "python batch_correction/metrics/run_ASW.py {wildcards.software} {input} {params.output_dir} {wildcards.nPCs}"

SOFTWARES=["Raw_nopp", "Raw_PCA", "Seurat3", "Harmony", "LIGER", "SIMBA", "SIMBA-new", "SIMBA-anchor"]
rule evaluate_LISI:
    input:
        expand("batch_correction/{software}/output/{{dataset}}_{software}_pca.txt", software=SOFTWARES)
    output:
        "batch_correction/metrics/{dataset}/LISI_40_PC{nPCs}.txt" 
    params:
        output_dir="batch_correction/metrics/{dataset}/",
    shell: 
        "Rscript batch_correction/metrics/run_LISI.R {input} {params.output_dir} {wildcards.nPCs}"

METRICS=["ARI", "ASW"]
rule collect_metrics_for_method:
    input:
        expand("batch_correction/metrics/{{dataset}}/{{software}}_{metric}_PC{{nPCs}}.txt", metric=METRICS),
        "batch_correction/metrics/{dataset}/LISI_40_PC{nPCs}.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_metrics_collected_PC{nPCs}.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/"
    shell:
        "Rscript batch_correction/metrics/collect_metrics_for_method.R {wildcards.software} {input} {params.output_dir} {wildcards.nPCs}"

rule collect_metrics:
    input: 
        expand("batch_correction/metrics/{{dataset}}/{software}_metrics_collected_PC{{nPCs}}.txt", software=SOFTWARES)
    output:
        "batch_correction/metrics/{dataset}/metrics_collected_PC{nPCs}.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/"
    shell:
        "Rscript batch_correction/metrics/collect_metrics_all.R {input} {params.output_dir} {wildcards.nPCs}"
