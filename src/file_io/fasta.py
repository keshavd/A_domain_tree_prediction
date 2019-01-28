from Bio import SeqIO
import numpy as np
from sklearn.model_selection import train_test_split

fasta = "/mnt/storage/grid/home/keshav/projects/A_domain_tree_prediction/platinum_sequences/platinum_a_domains.faa"

records = {record.id:record for record in list(SeqIO.parse(fasta, 'fasta'))}

entries = dict()
for record in records.keys():
    correct_1 = record.split('~')[0]
    correct_2 = record.split('~')[1]
    for correct in [correct_1, correct_2]:
        try:
            entries[correct].append(record)
        except KeyError:
            entries[correct] = []
            entries[correct].append(record)

embedded ={k:i for i,k in enumerate(sorted(entries.keys()))}

bad_keys = ['unused',  'not_used', 'Unused']
for k in bad_keys:
    del entries[k]

small_keys = [k for k in entries.keys() if len(entries[k]) < 4]
for k in small_keys:
    del entries[k]

y = []
X = []
for aa, record_list in entries.items():
    embed = embedded[aa]
    for record in record_list:
        y.append(embed)
        X.append(record)

y = np.array(y)
X = np.array(X)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=4206969)

# Get Basal Accuracy Measure from Training Files

train_fasta = "/mnt/storage/grid/home/keshav/projects/cdhit_domain_scorer/model/train.fasta"

train_records = [records[x] for x in X_train]

SeqIO.write(train_records, train_fasta, "fasta")