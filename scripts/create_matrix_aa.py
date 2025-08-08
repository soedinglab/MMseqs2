import os

import pandas as pd

base_remote_path = "ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/"
base_local_path = "../matrices/"
names = [
    "BLOSUM45",
    "BLOSUM50",
    "BLOSUM62",
    "BLOSUM80",
    "BLOSUM90",
    "PAM100",
    "PAM120",
    "PAM160",
    "PAM200",
    "PAM250"
]

os.makedirs(base_local_path, exist_ok = True)

for name in names:
    path = base_remote_path + name
    df = pd.read_csv(path, delim_whitespace = True, comment = "#")
    df = df.drop(index = "*")
    df = df.drop(columns = "*")

    min_val = -128

    # note: '[' = 'A' + 26
    for i in range(27):
        c = chr(ord("A") + i)
        if not c in df.index:
            df.loc[c, :] = min_val

    for i in range(32):
        c = chr(ord("A") + i)
        if not c in df.columns:
            df.loc[:, c] = min_val

    # alphabetically sort the amino acids
    df = df.sort_index(axis = 0)
    df = df.sort_index(axis = 1)

    for col in df.columns:
        df[col] = df[col].astype(int)

    res = "[" + df.to_csv(index = False, header = False) + "]"
    res = res.replace("\n", ",\n")
    res_path = base_local_path + name

    print(name)
    print(df)
    print()

    with open(res_path, "w") as f:
        f.write(res)
