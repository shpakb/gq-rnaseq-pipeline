import pandas as pd

df = pd.read_csv("/home/boris/gq-rnaseq-pipeline/files/srr_to_gsm.tsv", sep="\t")
print(df[df.gsm == "GSM261958"]["srr"].tolist())
