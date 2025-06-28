import random

lookup_path = "../data/scop/scop_lookup.fix.tsv"
pssm_path = "../data/scop/scop_mmseqs_pssm.pssm"
seq_path = "../data/scop/scop.fasta"
res_path = "../data/scop/pairs.pssm"

seq_to_scop = {}
families = {}

def process_scop_id(scop_id):
    #return scop_id[:scop_id.rindex(".")]
    return scop_id

with open(lookup_path) as f:
    for line in f:
        seq_id, scop_id = line.strip().split()
        scop_id = process_scop_id(scop_id)
        seq_to_scop[seq_id] = scop_id
        families[scop_id] = ([], [])

seq_lines = {}
seq_seqs = {}

with open(pssm_path) as f:
    curr_seq = None
    for line in f:
        if line.startswith("#"):
            curr_seq = line[1:].strip()
            seq_lines[curr_seq] = []
        else:
            seq_lines[curr_seq].append(line.strip())

with open(seq_path) as f:
    curr_seq = None
    for line in f:
        if line.startswith(">"):
            curr_seq = line[1:].strip()
            seq_seqs[curr_seq] = ""
        else:
            seq_seqs[curr_seq] += line.strip()

for seq_id, pssm in seq_lines.items():
    scop_id = seq_to_scop[seq_id]
    families[scop_id][0].append(pssm)

for seq_id, seq in seq_seqs.items():
    if not seq_id in seq_to_scop:
        continue
    scop_id = seq_to_scop[seq_id]
    families[scop_id][1].append(seq)

seq_pssm_pairs = []

def consensus_seq(lines):
    return "".join([s.split()[1] for s in lines[1:]])

for _, (pssm_family, seq_family) in families.items():
    random.shuffle(pssm_family)
    random.shuffle(seq_family)
    for i in range(min(len(pssm_family), len(seq_family))):
        pssm = pssm_family[i]
        seq = seq_family[i]
        seq_pssm_pairs.append((seq, consensus_seq(pssm), pssm))

print("Number of seq-pssm pairs:", len(seq_pssm_pairs))

with open(res_path, "w") as f:
    for seq, cns, pssm in seq_pssm_pairs:
        f.write("#" + seq + "\n")
        f.write("#" + cns + "\n")
        f.write("\n".join(pssm) + "\n")
