"""Unit tests for CodonMapper class"""
import unittest
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.codon_mapper import CodonMapper


class TestCodonMapper(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.mapper = CodonMapper()
    
    def test_translate_single_codon(self):
        """Test translation of single codons via codon_table"""
        self.assertEqual(self.mapper.codon_table["ATG"], "M")
        self.assertEqual(self.mapper.codon_table["TGG"], "W")
        self.assertEqual(self.mapper.codon_table["TAA"], "*")  # Stop codon
        self.assertEqual(self.mapper.codon_table["TAG"], "*")  # Stop codon
        self.assertEqual(self.mapper.codon_table["TGA"], "*")  # Stop codon
    
    def test_translate_codon_lowercase(self):
        """Test that lowercase input works via nine_mer"""
        # Test via nine_mer method which handles uppercase
        result = self.mapper.translate_nine_mer("atgcgtatc")
        self.assertEqual(result, "MRI")
    
    def test_translate_codon_invalid(self):
        """Test invalid codon handling returns X"""
        # CodonMapper uses X for invalid codons, not None
        # Test via nine_mer with invalid codon
        result = self.mapper.translate_nine_mer("ATNGGGTTT")
        self.assertIn("X", result)  # Contains X due to ATN
    
    def test_translate_nine_mer_valid(self):
        """Test translation of valid 9-mers"""
        # ATG CGT ATC = M R I -> "MRI"
        result = self.mapper.translate_nine_mer("ATGCGTATC")
        self.assertEqual(result, "MRI")
        
        # TTT CCC AAA = F P K -> "FPK"
        result = self.mapper.translate_nine_mer("TTTCCCAAA")
        self.assertEqual(result, "FPK")
    
    def test_translate_nine_mer_with_stop(self):
        """Test 9-mer containing stop codon"""
        # ATG TAA CCC = M * P -> "M*P"
        result = self.mapper.translate_nine_mer("ATGTAACCC")
        self.assertEqual(result, "M*P")
    
    def test_translate_nine_mer_invalid(self):
        """Test invalid 9-mer handling"""
        # Wrong length - returns XXX
        result = self.mapper.translate_nine_mer("ATGCGT")
        self.assertEqual(result, "XXX")
        
        # Contains N - translates to X
        result = self.mapper.translate_nine_mer("ATGNCGATC")
        self.assertIn("X", result)
        
        # Contains - - translates to X
        result = self.mapper.translate_nine_mer("ATG-GTATC")
        self.assertIn("X", result)
    
    def test_all_standard_codons(self):
        """Test that all 64 standard codons are handled"""
        # Just verify the codon table has correct number of entries
        self.assertEqual(len(self.mapper.codon_table), 64)
        
        # Test a few from each amino acid
        test_cases = [
            ("GGG", "G"), ("GGA", "G"), ("GGC", "G"), ("GGT", "G"),  # Glycine
            ("TTT", "F"), ("TTC", "F"),  # Phenylalanine
            ("ATG", "M"),  # Methionine (start)
            ("TGG", "W"),  # Tryptophan
        ]
        
        for codon, expected_aa in test_cases:
            self.assertEqual(
                self.mapper.codon_table[codon], 
                expected_aa,
                f"Codon {codon} should translate to {expected_aa}"
            )


if __name__ == '__main__':
    unittest.main()
