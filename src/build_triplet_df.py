import pandas as pd
from src.load_clones import load_clones


def build_triplet_df(clones):
    rows = []
    for clone_id, clone in clones.items():
        clone.extract_nines()
        clone.translate_nines()
        for idx, aa_triplet in clone.nine_aa_triplets:
            if aa_triplet is not None:
                rows.append({
                    'triplet': aa_triplet,
                    'clone_id': clone_id,
                    'index': idx
                })
    return pd.DataFrame(rows)


