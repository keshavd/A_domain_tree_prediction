import os
import sys
import time
import pandas as pd
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_label(label):
    conversion = {
                      'MeAsp':'asp',
                      'MeGlu':'glu',
                      'MePhe':'phe',
                      'MePro':'pro',
                      'me-pro':'pro',
                      'Pro|MePro':'pro',
                      'OHAsn':'asn',
                      'OHAsp':'asp',
                      'OHLeu':'leu',
                      'OHOrn':'orn',
                      'OHPhe':'phe',
                      'OHTyr':'tyr',
                      'OHTyr|Tyr':'tyr',
                      'OHVal':'val',
                      '&beta;-Ala':'bala',
                      '5-Hydroxy-L-leucine':'leu',
                      '3-me-glu|glu':'glu',
                      '4-5-dehydro-arg|me-tyr':'arg|tyr',
                      '6cl-trp|trp':'trp',
                      'Pro|MePro':'pro',
                      'ala|d-ala':'ala',
                      'alpha-me-ser':'aser',
                      'aoh-gly|gly':'gly',
                      'asn|d-asn':'asn',
                      'asp|bme-asp':'asp',
                      'boh-d-leu|d-leu':'leu',
                      'boh-d-leu|d-leu|leu':'leu',
                      'd-cyclo-horn|horn':'horn',
                      'd-dab|dab':'dab',
                      'd-gln':'gln',
                      'd-pro':'pro',
                      'd-ser|ser':'ser',
                      'd-val|val':'val',
                      'glu|3-me-glu':'glu',
                      'goh-val|val':'val',
                      'me-pro|pro':'pro',
                      'MeAsp': 'asp',
                      'MePhe': 'phe',
                      'MePro': 'pro',
                      'OHAsn': 'asn',
                      'OHAsp': 'asp',
                      'OHLeu': 'leu',
                      'OHOrn': 'orn',
                      'OHTyr': 'tyr',
                      'OHVal': 'val',
                      'lys-beta':'blys',
                      'ala-beta':'bala',
                     }
    nullabels = ['not used','not_used','notused','Unused','unused','no_call','no_confident_result','no_force_needed','unk']
    if label in nullabels:
        return 'null'
    elif label in list(conversion.keys()):
        return conversion[label].lower()
    else:
        return label.lower()

def csvtofasta(csvfile):
    df= pd.read_csv(sys.argv[1])
    recs = [(SeqRecord(Seq(row['seq']),description='',id='{}~{}'.format(clean_label(row['label']),row['seq_id']))) for idx,row in df.iterrows()]
    name = os.path.basename(sys.argv[1]).split('.')[0]
    outname = '{}.faa'.format(name)
    SeqIO.write(recs,outname,'fasta')
    return outname

def clean_fasta_header(header):
#expect the format label~label~idstuff
    bh = header.split('~')
    if len(bh) > 2:
        l1 = clean_label(bh[0])
        l2 = clean_label(bh[1])
        return '{}~{}~{}.'.format(l1,l2,'~'.join(bh[2:]))
    elif len(bh) == 2:
        l1 = clean_label(bh[0])
        return '{}~{}.'.format(l1,bh[1])
    else:
        return header

def sequence_cleaner(fasta_file):
    sequences={}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(seq_record.seq).upper()
        if sequence not in sequences:
            sequences[sequence] = clean_fasta_header(seq_record.id)
        else:
            sequences[sequence] += "_" + clean_fasta_header(seq_record.id)
    recs = [SeqRecord(Seq(key),description='',id=value) for key,value in sequences.items()]
    name = os.path.basename(fasta_file).split('.')[0]
    outname = '{}.deduped.faa'.format(name)
    SeqIO.write(recs,outname,'fasta')
    return outname

def callmafft(fasta_file,tophylip=False):
    fname = '{}.mafft.aln.faa'.format(os.path.basename(fasta_file).split('.')[0])
    msa_cmd = 'mafft --auto {} >| {}'.format(fasta_file,fname)
    os.system(msa_cmd)
    if tophylip:
        while not os.path.exists(fname):
            time.sleep(1)
        if os.path.isfile(fname):
            pname = '{}.mafft.aln.phylip'.format(os.path.basename(fasta_file).split('.')[0])
            AlignIO.convert(fname,'fasta',pname,'phylip-relaxed')
        print('finished')
        return pname
    else:
        return fname

if __name__ == '__main__':
    inputfile = sys.argv[1]
    fileext = os.path.basename(inputfile).split('.')[1]
    if fileext ==  'csv':
        fasta = csvtofasta(inputfile)
        ddfasta = sequence_cleaner(fasta)
        msa_phylip = callmafft(ddfasta)
        print('phylip msa file is at {}'.format(msa_phylip))
    elif fileext in ['faa','fasta','fna']:
        ddfasta = sequence_cleaner(inputfile)
        msa_phylip = callmafft(ddfasta,False)
        print('msa file is at {}'.format(msa_phylip))


