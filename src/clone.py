
from src.codon_mapper import CodonMapper

class Clone:
    def __init__(self, seq_id=None, ai=None, sample_=None, subject_=None, clone_id=None, function=None, copy_nu=None, cdr3_aa=None, sequence=None, germline=None, ab_target=None, time_po=None):
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
        self.ab_target = ab_target
        self.time_po = time_po
        self.nines = []  # List to store extracted nonuplets (nines)
        self.nine_aa_triplets = []  # List to store translated amino acid triplets


    def extract_nines(self):
        """
        Extracts all possible nonuplets (sliding window of size 9, step 1) from the sequence and stores them in self.nines.
        """
        if self.sequence:
            self.nines = [self.sequence[i:i+9] for i in range(len(self.sequence) - 8)]
        else:
            self.nines = []

    def translate_nines(self):
        """
        Translates each nonuplet in self.nines to an amino acid triplet using CodonMapper.
        If a nonuplet contains '-' or 'N', stores None at that index.
        """
        codon_mapper = CodonMapper()
        self.nine_aa_triplets = []
        for nine in self.nines:
            if '-' in nine or 'N' in nine:
                self.nine_aa_triplets.append(None)
            else:
                aa_triplet = codon_mapper.translate_nine_mer(nine)
                self.nine_aa_triplets.append(aa_triplet)

    def __repr__(self):
        return (
            f"Clone(seq_id={self.seq_id}, ai={self.ai}, sample_={self.sample_}, subject_={self.subject_}, "
            f"clone_id={self.clone_id}, function={self.function}, copy_nu={self.copy_nu}, cdr3_aa={self.cdr3_aa}, "
            f"sequence={self.sequence}, germline={self.germline}, ab_target={self.ab_target}, time_po={self.time_po})"
        )
