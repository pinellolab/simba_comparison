import pandas as pd

metric_files = map(float, input.strip('[]').split(','))

for metric_file in metric_files:
    if metric_file.endswith("ARI.txt"):
        df = pd.read_csv(metric_file, sep = "\t")
        df.loc([df.use_case.endswith("_median"), ]
