import pandas as pd
from src.load_clones import load_clones


def build_triplet_df(clones, region: str = None):
    """Build occurrence-level DataFrame of amino-acid triplets for all clones.

    If region is provided ('r1' or 'r2'), only include triplet occurrences whose
    nucleotide start index falls into that region. Regions are defined by 0-based
    nucleotide start index: r1 => 0..220 inclusive, r2 => 221+.

    The DataFrame columns are: ['triplet', 'clone_id', 'index'] where index is the
    1-based start position (keeps existing behavior).
    """
    rows = []
    for clone_id, clone in clones.items():
        clone.extract_nines()
        clone.translate_nines()
        for idx, aa_triplet in clone.nine_aa_triplets:
            if aa_triplet is None:
                continue
            # idx is 1-based start index from Clone.extract_nines()
            if region is not None:
                start0 = idx - 1
                if region == 'r1' and not (0 <= start0 <= 220):
                    continue
                if region == 'r2' and not (start0 > 220):
                    continue
            rows.append({
                'triplet': aa_triplet,
                'clone_id': clone_id,
                'index': idx
            })
    return pd.DataFrame(rows)


