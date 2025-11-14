
import itertools
import json
from pathlib import Path

# Standard DNA nucleotides
nucleotides = ['A', 'C', 'G', 'T']  

nonuplets = [''.join(t) for t in itertools.product(nucleotides, repeat=9)]
nonuplets_map = {triplet: {"count": 0, "indices": []} for triplet in nonuplets}

out_path = Path(__file__).with_name('nucleotide_nonuplets_9.json')
with out_path.open('w', encoding='utf-8') as f:
    json.dump(nonuplets_map, f, indent=2, sort_keys=True)

print(f"Wrote {len(nonuplets_map)} triplets to {out_path}")