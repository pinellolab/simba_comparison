rule batch_correction:
    input:
        "batch_correction/data/{dataset}/{dataset}_expr.txt.gz",
        "batch_correction/data/{dataset}/{dataset}_metadata.txt.gz",
    params:
        output_prefix="{dataset}",
        output_dir="batch_correction/{software}/output/"
    output:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    shell:
        "Rscript batch_correction/{wildcards.software}/run_{wildcards.software}.R {input} '{params.output_prefix}' '{params.output_dir}'"
       
rule evaluate_ARI:
    input:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_ARI.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/"
    shell:
        "Rscript batch_correction/metrics/run_ARISampled.R {wildcards.software} {input} '{params.output_dir}'"

rule evaluate_ASW:
    input:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_ASW.txt"
    params:
        output_dir="batch_correction/metrics/{dataset}/"
    shell:
        "python batch_correction/metrics/run_ASW.py {wildcards.software} {input} {params.output_dir}"

rule evaluate_LISI:
    input:
        "batch_correction/{software}/output/{dataset}_{software}_pca.txt"
    output:
        "batch_correction/metrics/{dataset}/{software}_LISI_40.txt" 
    params:
        output_dir="batch_correction/metrics/{dataset}/"
    shell: 
        "Rscript batch_correction/metrics/run_LISI.R {wildcards.software} {input} {params.output_dir}"
