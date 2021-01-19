rule batch_correction:
    input:
        "batch_correction/data/{dataset}/{dataset}_expr.txt.gz",
        "batch_correction/data/{dataset}/{dataset}_metadata.txt.gz",
        {dataset},
        "batch_correction/{software}/output"
    output:
        "batch_correction/{software}/output/"
    shell:
        "Rscript batch_correction/{software}/run_{software}.R {input}"
        
