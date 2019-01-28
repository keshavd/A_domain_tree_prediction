#!/bin/bash

refaln=$1 #reference alingment
qseqs=$2 #fasta of querys to fit to tree
raxtree=$3 #raxml tree
raxinfo=$4 #RAxML_info file generated with the tree

bname=$(basename $refaln)
hmm="$bname".hmm

hmmbuild $refaln $hmm
hmmalign -o query_and_reference.sto --mapali $refaln $hmm $qseqs
pplacer -t $raxtree -s $raxinfo query_and_reference.sto
