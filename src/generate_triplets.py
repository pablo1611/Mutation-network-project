import itertools
import json
from pathlib import Path

# 21 amino acids (20 standard + selenocysteine 'U')
amino_acids = [
    'A','R','N','D','C','E','Q','G','H','I',
    'L','K','M','F','P','S','T','W','Y','V','U'
]

triplets = [''.join(t) for t in itertools.product(amino_acids, repeat=3)]
triplet_map = {triplet: {"count": 0, "indices": []} for triplet in triplets}

out_path = Path(__file__).with_name('amino_acid_triplets.json')
with out_path.open('w', encoding='utf-8') as f:
    json.dump(triplet_map, f, indent=2, sort_keys=True)

print(f"Wrote {len(triplet_map)} triplets to {out_path}")

