# src/codon_mapper.py
from pathlib import Path
import json

class CodonMapper:
    """
    Maps every 9-nucleotide window (nonuplet) in a DNA sequence
    to its corresponding amino-acid triplet translation.
    """

    def __init__(self):
        # Load amino-acid triplets file (structure only, optional)
        triplet_path = Path(__file__).with_name("amino_acid_triplets.json")
        if triplet_path.exists():
            with triplet_path.open("r", encoding="utf-8") as f:
                self.amino_acid_triplets = json.load(f)
        else:
            self.amino_acid_triplets = {}

        # Standard genetic code
        self.codon_table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        }

    def translate_nine_mer(self, nine_mer: str) -> str:
        """Translate a 9-base DNA fragment into its amino-acid triplet."""
        nine_mer = nine_mer.upper()
        if len(nine_mer) != 9:
            return "XXX"

        codons = [nine_mer[i:i+3] for i in range(0, 9, 3)]
        aa_triplet = "".join(self.codon_table.get(c, "X") for c in codons)
        return aa_triplet

    def map_sequence(self, dna_sequence: str):
        """
        Slide a 9-base window **one base at a time** across the DNA sequence.
        Yields (start_index, nine_mer, translated_triplet).
        """
        dna_sequence = dna_sequence.upper()
        for i in range(0, len(dna_sequence) - 9 + 1):
            nine_mer = dna_sequence[i:i+9]
            yield i, nine_mer, self.translate_nine_mer(nine_mer)
