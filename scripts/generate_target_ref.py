from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from pathlib import Path
import re, yaml

target_dir = Path(snakemake.input[0]).parent
target = Path(snakemake.input[0]).parent.name
plasmid = open(snakemake.input[0], 'r').read().strip().upper()

specs = snakemake.config['parse_specs']['default']
if target in snakemake.config['parse_specs']:
    specs.update(snakemake.config['parse_specs'][target])

if (target_dir / 'template.txt').exists():
    print(f"Using custom template for {target}")
    template = re.compile(open(target_dir / 'template.txt', 'r').read().strip().upper())
else:
    print(f"Using default template for {target}")
    template = re.compile(open(snakemake.config['default_template'], 'r').read().strip().upper())

search = template.search(plasmid) # 5 groups, corresponding to "termini5", "gene", "spacer", "barcode", "termini3"

if not search:
    raise ValueError(f"Plasmid sequence does not match template for {target}")

termini5, gene, spacer, barcode, termini3 = search.groups()

ref_seq = Seq(termini5 + gene + spacer + 'N'*len(barcode) + termini3)

with open(target_dir / f'{target}.fasta', 'w') as f:
    f.write(f'>{target}\n{ref_seq}\n')

gb = SeqRecord(ref_seq, id=target, name=target, description=f'{target} yeast display library reference.')

# Add molecule type
gb.annotations['molecule_type'] = 'DNA'

st = search.start(1)
gb.features = [
    SeqFeature(FeatureLocation(0, search.end(1)-st, strand=1), type='termini5', qualifiers={'locus_tag': 'termini5', 'label': 'termini5'}),
    SeqFeature(FeatureLocation(search.end(1)-st, search.end(2)-st, strand=1), type='gene', qualifiers={'locus_tag': 'gene', 'label': 'gene'}),
    SeqFeature(FeatureLocation(search.end(2)-st, search.end(3)-st, strand=1), type='spacer', qualifiers={'locus_tag': 'spacer', 'label': 'spacer'}),
    SeqFeature(FeatureLocation(search.end(3)-st, search.end(4)-st, strand=1), type='barcode', qualifiers={'locus_tag': 'barcode', 'label': 'barcode'}),
    SeqFeature(FeatureLocation(search.end(4)-st, len(ref_seq), strand=1), type='termini3', qualifiers={'locus_tag': 'termini3', 'label': 'termini3'})
]

with open(target_dir / f'pacbio_amplicons.gb', 'w') as f:
    SeqIO.write(gb, f, 'genbank')

with open(target_dir / f'parse_specs.yaml', 'w') as f:
    yaml.dump({target: specs}, f)