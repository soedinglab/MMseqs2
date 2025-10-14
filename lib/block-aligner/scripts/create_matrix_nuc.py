import os
import sys

import numpy as np

matrices = [
    {
        "name": "NW1",
        "scores": {
            ("A", "A"): 1,
            ("A", "C"): -1,
            ("A", "G"): -1,
            ("A", "N"): -1,
            ("A", "T"): -1,
            ("C", "C"): 1,
            ("C", "G"): -1,
            ("C", "N"): -1,
            ("C", "T"): -1,
            ("G", "G"): 1,
            ("G", "N"): -1,
            ("G", "T"): -1,
            ("N", "N"): 1,
            ("N", "T"): -1,
            ("T", "T"): 1
        }
    },
]

base_local_path = "../matrices/"

os.makedirs(base_local_path, exist_ok = True)

for matrix in matrices:
    res = np.full((8, 16), -128)

    for key, val in matrix["scores"].items():
        res[ord(key[0]) & 0b111][ord(key[1]) & 0b1111] = val

        if key[0] != key[1]:
            res[ord(key[1]) & 0b111][ord(key[0]) & 0b1111] = val

    res = str(res.flatten().tolist()).replace(" ", "")

    res_path = base_local_path + matrix["name"]

    print(matrix["name"])
    print(res)
    print()

    with open(res_path, "w") as f:
        f.write(res)
