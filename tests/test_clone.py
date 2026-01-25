"""Unit tests for Clone class"""
import unittest
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.clone import Clone


class TestClone(unittest.TestCase):
    
    def test_clone_initialization(self):
        """Test basic clone initialization with required fields"""
        clone = Clone(
            seq_id="seq001",
            ai=1,
            sample_="sample1",
            subject_="subject1",
            clone_id=100,
            function="productive",
            copy_nu=5,
            cdr3_aa="CARDRGYW",
            sequence="ATGCGATCGATCG",
            germline="IGHV1-1"
        )
        
        self.assertEqual(clone.seq_id, "seq001")
        self.assertEqual(clone.clone_id, 100)
        self.assertEqual(clone.sequence, "ATGCGATCGATCG")
        self.assertEqual(clone.cdr3_aa, "CARDRGYW")
    
    def test_clone_extra_fields(self):
        """Test that extra fields are properly stored"""
        clone = Clone(
            seq_id="seq001",
            clone_id=100,
            sequence="ATGCGATCGATCG",
            experiment="exp1",
            timepoint="baseline"
        )
        
        self.assertEqual(clone.experiment, "exp1")
        self.assertEqual(clone.timepoint, "baseline")
        self.assertIn("experiment", clone.extra_fields)
        self.assertIn("timepoint", clone.extra_fields)
    
    def test_extract_nines_basic(self):
        """Test nonuplet extraction from sequence"""
        clone = Clone(clone_id=1, sequence="ATGCGATCG")
        clone.extract_nines()
        
        # Sequence of length 9 should produce 1 nonuplet
        self.assertEqual(len(clone.nines), 1)
        self.assertEqual(clone.nines[0], (1, "ATGCGATCG"))
    
    def test_extract_nines_sliding_window(self):
        """Test sliding window extraction"""
        clone = Clone(clone_id=1, sequence="ATGCGATCGATCG")
        clone.extract_nines()
        
        # Length 13, should produce 5 nonuplets
        self.assertEqual(len(clone.nines), 5)
        self.assertEqual(clone.nines[0][1], "ATGCGATCG")
        self.assertEqual(clone.nines[1][1], "TGCGATCGA")
        self.assertEqual(clone.nines[4][1], "GATCGATCG")
    
    def test_extract_nines_empty_sequence(self):
        """Test extraction with empty or None sequence"""
        clone = Clone(clone_id=1, sequence=None)
        clone.extract_nines()
        self.assertEqual(len(clone.nines), 0)
        
        clone2 = Clone(clone_id=2, sequence="")
        clone2.extract_nines()
        self.assertEqual(len(clone2.nines), 0)
    
    def test_extract_nines_too_short(self):
        """Test extraction with sequence shorter than 9"""
        clone = Clone(clone_id=1, sequence="ATGC")
        clone.extract_nines()
        self.assertEqual(len(clone.nines), 0)
    
    def test_translate_nines_valid(self):
        """Test translation of valid nonuplets"""
        # ATG = M, CGT = R, ATC = I -> MRI
        clone = Clone(clone_id=1, sequence="ATGCGTATC")
        clone.extract_nines()
        clone.translate_nines()
        
        self.assertEqual(len(clone.nine_aa_triplets), 1)
        idx, aa_triplet = clone.nine_aa_triplets[0]
        self.assertEqual(idx, 1)
        self.assertEqual(aa_triplet, "MRI")
    
    def test_translate_nines_with_gaps(self):
        """Test translation handles gaps (N or -)"""
        clone = Clone(clone_id=1, sequence="ATGCGN-TC")
        clone.extract_nines()
        clone.translate_nines()
        
        # Should produce 1 nonuplet but translation should be None due to N and -
        self.assertEqual(len(clone.nine_aa_triplets), 1)
        idx, aa_triplet = clone.nine_aa_triplets[0]
        self.assertIsNone(aa_triplet)
    
    def test_translate_nines_multiple(self):
        """Test translation of multiple nonuplets"""
        # Create sequence with multiple valid nonuplets
        clone = Clone(clone_id=1, sequence="ATGCGTATCGATCGA")
        clone.extract_nines()
        clone.translate_nines()
        
        # Should have multiple translations
        self.assertGreater(len(clone.nine_aa_triplets), 1)
        
        # Check first one is valid
        idx, aa_triplet = clone.nine_aa_triplets[0]
        self.assertIsNotNone(aa_triplet)
        self.assertEqual(len(aa_triplet), 3)


if __name__ == '__main__':
    unittest.main()
