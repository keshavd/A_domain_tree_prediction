import numpy as np
from subprocess import Popen, PIPE
import shlex
from tqdm import tqdm

# SETTINGS
cdhit = "/mnt/storage/grid/home/keshav/bin/cd-hit"
num_cores = 2
cdhit_train = "/mnt/storage/grid/home/keshav/projects/cdhit_domain_scorer/model/train.fasta"
decimal_places = 3

# GENERATE PARAMETERS
ws_5 = np.arange(0.7, 1.0, 0.01)
ws_4 = np.arange(0.6, 0.7, 0.01)
ws_3 = np.arange(0.5, 0.6, 0.01)
ws_2 = np.arange(0.4, 0.5, 0.01)
parameter_combos = {2: ws_2,
                    3: ws_3,
                    4: ws_4,
                    5: ws_5}

arg_combos = []

for k, v in parameter_combos.items():
    arg_1 = k
    for arg_2 in v:
        arg_combos.append((arg_1, arg_2))

# PERFORM CDHIT CLUSTERING

# Train
output_train = []
for word_count, threshold in tqdm(arg_combos, desc="CDHIT Analysis"):
    cdhit_out = "%s_%s_%s.cdhit.out" % (cdhit_train, str(word_count), str(threshold))
    cmd = "%s -i %s -o %s -c %s -d 0 -n %d -T %d" % (cdhit, cdhit_train,
            cdhit_out, str(round(threshold, decimal_places)), word_count,
            num_cores)
    args=shlex.split(cmd)
    p = Popen(args, stderr=PIPE, stdout=PIPE)
    out,err =p.communicate()
    if len(err) > 0:
        print("###################################")
        print(out.decode('utf-8'))
        print("###################################")
        print(err)
    else:
        output_train.append(cdhit_out)
