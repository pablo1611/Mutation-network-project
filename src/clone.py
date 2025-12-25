
from src.codon_mapper import CodonMapper

class Clone:
    def __init__(self, seq_id=None, ai=None, sample_=None, subject_=None, clone_id=None, function=None, copy_nu=None, cdr3_aa=None, sequence=None, germline=None, **extra_fields):
        # First 10 stable attributes (always present in this order)
        self.seq_id = seq_id
        self.ai = ai
        self.sample_ = sample_
        self.subject_ = subject_
        self.clone_id = clone_id
        self.function = function
        self.copy_nu = copy_nu
        self.cdr3_aa = cdr3_aa
        self.sequence = sequence
        self.germline = germline

        # Dynamic attributes from dataset (anything after the first 10)
        self.extra_fields = {}
        for k, v in (extra_fields or {}).items():
            setattr(self, k, v)
            self.extra_fields[k] = v

        self.nines = []  # List to store extracted nonuplets (nines)
        self.nine_aa_triplets = []  # List to store (index, translated amino acid triplet or None)


    def extract_nines(self):
        """
        Extracts all possible nonuplets (sliding window of size 9, step 1) from the sequence and stores (index, nine) tuples in self.nines.
        """
        if self.sequence:
            self.nines = [(i+1, self.sequence[i:i+9]) for i in range(len(self.sequence) - 8)]
        else:
            self.nines = []

    def translate_nines(self):
        """
        Translates each nonuplet in self.nines to an amino acid triplet using CodonMapper.
        If a nonuplet contains '-' or 'N', stores None at that index.
        Stores results as a list of (index, aa_triplet or None).
        """
        codon_mapper = CodonMapper()
        self.nine_aa_triplets = []
        for idx, nine in self.nines:
            if '-' in nine or 'N' in nine:
                self.nine_aa_triplets.append((idx, None))
            else:
                aa_triplet = codon_mapper.translate_nine_mer(nine)
                self.nine_aa_triplets.append((idx, aa_triplet))

    def __repr__(self):
        base = (
            f"Clone(seq_id={self.seq_id}, ai={self.ai}, sample_={self.sample_}, subject_={self.subject_}, "
            f"clone_id={self.clone_id}, function={self.function}, copy_nu={self.copy_nu}, cdr3_aa={self.cdr3_aa}, "
            f"sequence={self.sequence}, germline={self.germline}"
        )
        if self.extra_fields:
            extras = ", ".join(f"{k}={v!r}" for k, v in self.extra_fields.items())
            return base + ", " + extras + ")"
        return base + ")"
