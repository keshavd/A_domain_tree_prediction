import re
import sys
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def rename(x):
    rn = '{}~{}~{}~{}~{}'.format(x['validated_monomer'],x['orig_prediction'],x['module_num'],x['orf_idx'],x['seq_name'])
    return re.sub('(?<!fasta\|unknown) ','_',rn)

recs,seenseqs = [],[]
d = json.load(open(sys.argv[1])) #final_domains.json
for x in d:
    if x['seq'] not in seenseqs:
        recs.append(SeqRecord(Seq(x['seq']),id=rename(x)))
        seenseqs.append(x['seq'])
SeqIO.write(recs,'platinum_a_domains.faa','fasta')
