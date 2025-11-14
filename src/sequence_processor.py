# src/sequence_processor.py
import json, os
from src.codon_mapper import CodonMapper
from src.generate_triplets import generate_triplets
from src.generate_nonuplets import generate_nonuplets

class SequenceProcessor:
    """
    Reads a nucleotide dataset, extracts 9-mers, translates to amino acid triplets,
    updates counters and index arrays, and populates both JSONs.
    """

    def __init__(self, dataset_path: str):
        self.dataset_path = dataset_path
        self.mapper = CodonMapper()
        self.triplets_data = generate_triplets()
        self.nonuplets_data = generate_nonuplets()

    def _read_dataset(self) -> str:
        """Read the sequence dataset (assumes plain text or FASTA-like format)."""
        with open(self.dataset_path, "r") as f:
            lines = [line.strip() for line in f if not line.startswith(">")]
        return "".join(lines).upper()

    def process(self):
        """Main entry point: reads dataset and populates both JSONs."""
        sequence = self._read_dataset()
        n = len(sequence)

        for i in range(0, n - 9 + 1, 3):
            nine_mer = sequence[i:i+9]
            if len(nine_mer) < 9:
                continue

            amino_triplet = self.mapper.translate(nine_mer)

            # Update nonuplet stats
            if nine_mer not in self.nonuplets_data:
                self.nonuplets_data[nine_mer] = {"count": 0, "indices": []}
            self.nonuplets_data[nine_mer]["count"] += 1
            self.nonuplets_data[nine_mer]["indices"].append(i)

            # Update triplet stats
            if amino_triplet not in self.triplets_data:
                self.triplets_data[amino_triplet] = {"count": 0, "indices": []}
            self.triplets_data[amino_triplet]["count"] += 1
            self.triplets_data[amino_triplet]["indices"].append(i)

        return self.triplets_data, self.nonuplets_data

    def save_results(self, output_dir="data"):
        """Save populated JSONs."""
        os.makedirs(output_dir, exist_ok=True)
        trip_path = os.path.join(output_dir, "populated_triplets.json")
        nonup_path = os.path.join(output_dir, "populated_nonuplets.json")

        with open(trip_path, "w") as f:
            json.dump(self.triplets_data, f, indent=4)
        with open(nonup_path, "w") as f:
            json.dump(self.nonuplets_data, f, indent=4)

        return trip_path, nonup_path
