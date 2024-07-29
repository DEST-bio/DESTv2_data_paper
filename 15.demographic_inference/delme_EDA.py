import pandas as pd

df = pd.read_csv("output/opt_est_confidence_intervals.tsv", sep="\t")
df = df[df['region'] == "Europe"]
df = df.iloc[:3, [0, 1, 8, 9, 10, 11]]

print(df)