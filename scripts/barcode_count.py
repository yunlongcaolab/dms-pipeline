# import dms_variants.codonvarianttable
import argparse
import sys, os
import pandas as pd
import numpy as np
from Bio import SeqIO
import datetime
import yaml
import gzip

parser = argparse.ArgumentParser(description="Count barcodes for fastq file.")
parser.add_argument('--sample', type=str, help="Sample name.")
parser.add_argument('--batch', type=str, help="Batch name.")
parser.add_argument('--library', type=str, help="Library name.")
parser.add_argument('--input_fastq', '-i', nargs='+', type=str, help="Input fastq file(s).")
parser.add_argument('--output_dir', '-o', type=str)
parser.add_argument('--bclen', '-b', type=int)
parser.add_argument('--table', '-t', type=str)
parser.add_argument('--allowed_lowq', type=int)
parser.add_argument('--lowq', type=int)
parser.add_argument('--max_distance', type=int, help="Max allowed nearest Hamming distance for barcode mapping using a BK-Tree.")
parser.add_argument('--min_dist_diff', type=int, help="Min required difference in distance between the query barcode and its nearset and second nearset barcode.")

# def build_variant_table(table_path, wt_seq_path):
#     wt_seqrecord = SeqIO.read(wt_seq_path, 'fasta')
#     geneseq = str(wt_seqrecord.seq)
#     primary_target = wt_seqrecord.name
#     variants = dms_variants.codonvarianttable.CodonVariantTable(
#                geneseq = geneseq,
#                barcode_variant_file = table_path,
#                substitutions_are_codon = True,
#                substitutions_col = 'codon_substitutions',
#                primary_target = primary_target)
#     return variants, primary_target

def count_variants(variants, lib, fastq_data, bclen, max_dist, min_dist_diff, allowed_lowq, lowq=20):
    valid_barcodes = list(pd.unique(variants.query('library == @lib')['barcode']))
    
    sys.stdout.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+f' | N barcodes: {len(valid_barcodes)}\n')
    sys.stdout.flush()
    
    if max_dist >= 1:
        from bktree import BKTree
        bktree = BKTree(bclen)
        for bc in valid_barcodes:
            bktree.insert(bc)
        sys.stdout.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+' | BK-Tree constructed.\n')
        sys.stdout.flush()
    
    fates = {
        'processed': 0,
        'valid barcode': 0, # including recovered
        'invalid barcode (unmatched)': 0, 
        'invalid barcode (ambiguous)': 0, 
        'low quality barcode': 0,
        'unparseable barcode': 0, # length of read < bclen
        'containing N': 0,
        'recovered': 0, # approximate match
    }
    counts = {}
    for bc in valid_barcodes:
        counts[bc] = 0
    
    for fastq_file in fastq_data:
        if fastq_file[-3:] == '.gz':
            f = gzip.open(fastq_file, 'rt')
        else:
            f = open(fastq_file, 'rt')
        for entry in SeqIO.parse(f, 'fastq'):
            fates['processed'] += 1
            if (fates['processed'] % 500000 == 0):
                sys.stdout.write(
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+f' | %d processed, %d valid, %d recovered.\n'%(fates['processed'], fates['valid barcode'], fates['recovered']))
                sys.stdout.flush()
        
            if bclen > len(entry):
                fates['unparseable barcode'] += 1
                continue
            entry = entry[0:bclen]
            
            if 'N' in entry.seq:
                fates['containing N'] += 1
                continue
                
            if (sum(np.array(entry.letter_annotations['phred_quality']) < lowq) > allowed_lowq):
                fates['low quality barcode'] += 1
                continue
            
            if entry.seq in counts:
                counts[entry.seq] += 1
                fates['valid barcode'] += 1
                continue
            
            if max_dist < 1:
                fates['invalid barcode (unmatched)'] += 1
            else:
                _bc = str(entry.seq)
                _res = bktree.find_nearest(_bc, max_dist)
                if _res[1] == "NotFound":
                    fates['invalid barcode (unmatched)'] += 1
                    continue
                if min_dist_diff > 0:    
                    _res = bktree.find2(_bc, max_dist+min_dist_diff)
                    if (_res[1][0]-_res[0][0] < min_dist_diff):
                        fates['invalid barcode (ambiguous)'] += 1
                        continue
                
                fates['valid barcode'] += 1
                fates['recovered'] += 1
                counts[_res[0][1]] += 1
        f.close()

    counts = pd.DataFrame(list(counts.items()), columns=['barcode', 'count']).sort_values(['count', 'barcode'], ascending=[False, True]).reset_index(drop=True).query("count > 0")

    fates['valid ratio'] = fates['valid barcode'] / fates['processed']
    fates['detected barcodes'] = len(counts)
    
    return counts, fates


def main():
    args = parser.parse_args()
    
    # variants, primary_target = build_variant_table(args.table, args.wildtype)
    variants = pd.read_csv(args.table)
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    sys.stdout.write(f'variants table: {args.table}\n')

    sys.stdout.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+' | Start parsing barcodes...\n')
    sys.stdout.write(f'    fastq: {args.input_fastq}\n')
    sys.stdout.write(f'    max nearest distance: {args.max_distance}\n')
    sys.stdout.write(f'    min distance difference: {args.min_dist_diff}\n')
    sys.stdout.write(f'    barcode length: {args.bclen}\n')
    sys.stdout.write(f'    allowed lowQ bases: {args.allowed_lowq}\n')
    sys.stdout.flush()
    
    counts, fates = count_variants(variants, args.library, fastq_data = args.input_fastq, 
                                   bclen=args.bclen, 
                                   max_dist=args.max_distance, min_dist_diff=args.min_dist_diff, 
                                   allowed_lowq=args.allowed_lowq, lowq=args.lowq)
    sys.stdout.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+' | Done.\n')
    
    counts.to_csv(os.path.join(output_dir, 'counts.csv'), index=None)

    info = {
        'sample': args.sample,
        'batch': args.batch,
        'library': args.library,
        'fastq': args.input_fastq,
        'max_distance': args.max_distance,
        'min_dist_diff': args.min_dist_diff,
        'barcode_length': args.bclen,
        'allowed_lowq': args.allowed_lowq,
        'lowq': args.lowq,
        'stat': fates,
    }
    yaml.dump(info, open(os.path.join(output_dir, 'barcode_count_info.yaml'), 'w'))

if __name__ == '__main__':
    main()
