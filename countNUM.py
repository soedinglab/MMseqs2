#Counting the unique target ID from result2.m8 & reclass_result.txt
from collections import defaultdict
import pdb

PATH1 = '/home/yakim/benchmark/result_mmseq.m8'    # mmseq
PATH2 = '/home/yakim/benchmark/reclass_bit.m8'     # after reclassification - bit score
PATH3 = '/home/yakim/benchmark/reclass_seqid.m8'  # after reclassification - seq id
PATH4 = '/home/yakim/benchmark/reclass0922.m8'  # after reclassification - C++ (seqid)


# MMseq2 result count
qandt = []
hit_num = 0
redund1 = []
with open(PATH1, 'r') as f:
    for line in f:
        a = line.split()
        qandt.append([a[0],a[1]])

for q, t in qandt:
    query = q.split('_')[-1]  #query
    qnum = q.split('_')[0]
    if qnum not in redund1:
      redund1.append(qnum)
      if query == t:
          hit_num += 1
print('Hit # of mmseqs result:', hit_num, '&& Total #:', len(redund1))


# Mmseq2 result count (count also for every top 10)
'''
hit_num1 = 0
grouped = defaultdict(list)

with open(PATH1, 'r') as f:
  for line in f:
      a = line.split()
      query_full = a[0]      # query705_from_UniRef90_A0AA46S9G1
      target = a[1]          # A0AA46S9G1
      seq_identity = round(float(a[2]), 3)
      query_num = query_full.split('_')[0]      # query705
      query_last_id = query_full.split('_')[-1] # A0AA46S9G1

        # query 번호별로 그룹핑
      grouped[query_num].append({
            'query_last_id': query_last_id,
            'target': target,
            'seq_identity': seq_identity
        })        # row[0]                       row[1]                          row[2]
        #{query1: (A0AA46S9G1, A0AA46S9R1, 0.9), (A0AA46S9G1, A0AA46S9G1, 0.9), (A0AA46S9G1, A0AA34S9T1, 0.85)}

for query_num, rows in grouped.items():
    if not rows:
        continue
    top_seq_identity = rows[0]['seq_identity']    # 0.9
    top_query_last_id = rows[0]['query_last_id']  # A0AA46S9G1

    count = sum(
        (row['query_last_id'] == row['target']) and
        (row['seq_identity'] == top_seq_identity)
        for row in rows
    )
    #print(f'{query_num}: {count}')   # 각 그룹별로 개수 출력
    hit_num1 += count

print('Hit # of mmseqs result(top 10):', hit_num1, '&& Total #:', len(redund))
'''

# After reclassification result (bit score)
qandt2 = []
hit_num2 = 0
redund2 = []
with open(PATH2, 'r') as f:
    for line in f:
        a = line.split()
        qandt2.append([a[0],a[1]])

for q, t in qandt2:
    query = q.split('_')[-1]  #query
    qnum = q.split('_')[0]
    if qnum not in redund2:
      redund2.append(qnum)
      if query == t:
          hit_num2 += 1
print('Hit # of reclassification result(bit score):', hit_num2, '&& Total #:', len(redund2))


# After reclassification result (seq id)
qandt3 = []
hit_num3 = 0
redund3 = []
with open(PATH3, 'r') as f:
    for line in f:
        a = line.split()
        qandt3.append([a[0],a[1]])

for q, t in qandt3:
    query = q.split('_')[-1]  #query
    qnum = q.split('_')[0]
    if qnum not in redund3:
      redund3.append(qnum)
      if query == t:
          hit_num3 += 1
print('Hit # of reclassification result(seqid):', hit_num3, '&& Total #:', len(redund3))

# after reclass result count (C++)
qandt4 = []
hit_num4 = 0
redund4 = []
with open(PATH4, 'r') as f:
    for line in f:
        a = line.split()
        qandt4.append([a[0],a[1]])

for q, t in qandt4:
    query = q.split('_')[-1]  #query
    qnum = q.split('_')[0]
    if qnum not in redund4:
      redund4.append(qnum)
      if query == t:
          hit_num4 += 1
print('Hit # of reclassification result(C++):', hit_num4, '&& Total #:', len(redund4))

