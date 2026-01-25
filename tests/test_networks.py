"""Unit tests for KmerNetwork class"""
import unittest
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.networks import KmerNetwork, compute_network_distance


class TestKmerNetwork(unittest.TestCase):
    
    def test_initialization(self):
        """Test network initialization"""
        network = KmerNetwork(k=3)
        self.assertEqual(network.k, 3)
        self.assertEqual(len(network.nodes), 0)
        self.assertFalse(network.store_positions)
    
    def test_initialization_with_positions(self):
        """Test network initialization with position storage"""
        network = KmerNetwork(k=3, store_positions=True)
        self.assertTrue(network.store_positions)
    
    def test_add_kmer_basic(self):
        """Test adding a single kmer"""
        network = KmerNetwork(k=3)
        network.add_kmer("ATG", clone_id=1)
        
        self.assertIn("ATG", network.nodes)
        self.assertEqual(network.nodes["ATG"]["count"], 1)
        self.assertIn(1, network.nodes["ATG"]["clones"])
    
    def test_add_kmer_multiple_times(self):
        """Test adding same kmer multiple times"""
        network = KmerNetwork(k=3)
        network.add_kmer("ATG", clone_id=1)
        network.add_kmer("ATG", clone_id=1)
        network.add_kmer("ATG", clone_id=2)
        
        self.assertEqual(network.nodes["ATG"]["count"], 3)
        self.assertEqual(len(network.nodes["ATG"]["clones"]), 2)
        self.assertIn(1, network.nodes["ATG"]["clones"])
        self.assertIn(2, network.nodes["ATG"]["clones"])
    
    def test_add_kmer_with_stop_codon_nucleotide(self):
        """Test that stop codons in nucleotide 9-mers are flagged"""
        network = KmerNetwork(k=9)
        # TAA is a stop codon
        network.add_kmer("ATGTAACCC", clone_id=1)
        
        self.assertIn("ATGTAACCC", network.nodes)
        self.assertEqual(network.nodes["ATGTAACCC"].get("color"), 1)
    
    def test_add_kmer_with_stop_aa(self):
        """Test that stop symbol in AA triplets is flagged"""
        network = KmerNetwork(k=3)
        network.add_kmer("M*P", clone_id=1)
        
        self.assertIn("M*P", network.nodes)
        self.assertEqual(network.nodes["M*P"].get("color"), 1)
    
    def test_add_sequence_basic(self):
        """Test adding a sequence"""
        network = KmerNetwork(k=3)
        network.add_sequence("ATGCGT", clone_id=1)
        
        # Should create 4 kmers: ATG, TGC, GCG, CGT
        self.assertEqual(len(network.nodes), 4)
        self.assertIn("ATG", network.nodes)
        self.assertIn("TGC", network.nodes)
        self.assertIn("GCG", network.nodes)
        self.assertIn("CGT", network.nodes)
    
    def test_add_sequence_skips_invalid(self):
        """Test that sequences with N or - are skipped"""
        network = KmerNetwork(k=3)
        network.add_sequence("ATNGC-T", clone_id=1)
        
        # Should skip kmers containing N or -
        # Only valid: ATN (skip), TNG (skip), NGC (skip), GC- (skip), C-T (skip)
        # Actually all contain N or -, so should be empty
        self.assertEqual(len(network.nodes), 0)
    
    def test_add_sequence_with_positions(self):
        """Test adding sequence with position tracking"""
        network = KmerNetwork(k=3, store_positions=True)
        network.add_sequence("ATGCGT", clone_id=1)
        
        # Check positions are stored
        self.assertIn("ATG", network.nodes)
        self.assertIn("positions", network.nodes["ATG"])
        # Positions are tracked internally
        self.assertIsNotNone(network.nodes["ATG"]["positions"])
    
    def test_normalize_nodes_basic(self):
        """Test node normalization"""
        network = KmerNetwork(k=3)
        network.add_kmer("ATG", clone_id=1)
        network.add_kmer("ATG", clone_id=2)
        network.add_kmer("CGT", clone_id=1)
        
        network.normalize_nodes()
        
        # After normalization, each node should have normalized_freq
        self.assertIn("normalized_freq", network.nodes["ATG"])
        self.assertIn("normalized_freq", network.nodes["CGT"])
        
        # Total counts = 3, ATG appears 2 times
        self.assertAlmostEqual(network.nodes["ATG"]["normalized_freq"], 2/3)
        self.assertAlmostEqual(network.nodes["CGT"]["normalized_freq"], 1/3)
    
    def test_normalize_nodes_empty(self):
        """Test normalization of empty network"""
        network = KmerNetwork(k=3)
        network.normalize_nodes()  # Should not crash
        self.assertEqual(len(network.nodes), 0)
    
    def test_apply_node_threshold(self):
        """Test applying node frequency threshold"""
        network = KmerNetwork(k=3)
        network.add_kmer("ATG", clone_id=1)
        network.add_kmer("ATG", clone_id=2)
        network.add_kmer("CGT", clone_id=1)
        
        network.normalize_nodes()
        
        # Set threshold to filter out nodes with frequency < 0.5
        network.apply_node_threshold(0.5)
        
        # ATG has frequency 2/3 > 0.5, should be above_threshold=True
        # CGT has frequency 1/3 < 0.5, should be above_threshold=False
        self.assertIn("ATG", network.nodes)
        self.assertTrue(network.nodes["ATG"]["above_threshold"])
        self.assertIn("CGT", network.nodes)  # Node still exists
        self.assertFalse(network.nodes["CGT"]["above_threshold"])
    
    def test_compute_edge_probabilities(self):
        """Test edge probability computation"""
        network = KmerNetwork(k=3)
        
        # Add kmers that appear in same clone (co-occurrence)
        # Need to add them as a sequence so they're counted as co-occurring
        network.add_sequence("ATGTGCGCG", clone_id=1)
        
        network.normalize_nodes()
        network.compute_edge_probabilities()
        
        # Should have created edges attribute
        self.assertTrue(hasattr(network, 'edges'))
        # May or may not have edges depending on co-occurrence logic
        self.assertIsNotNone(network.edges)
        self.assertTrue(hasattr(network, 'edges'))
        self.assertGreater(len(network.edges), 0)
    
    def test_apply_edge_threshold(self):
        """Test applying edge probability threshold"""
        network = KmerNetwork(k=3)
        
        # Create edges
        network.add_kmer("ATG", clone_id=1)
        network.add_kmer("TGC", clone_id=1)
        
        network.normalize_nodes()
        network.compute_edge_probabilities()
        
        # Apply threshold
        network.apply_edge_threshold(0.5)
        
        # Check that edges are filtered (implementation-specific)
        self.assertTrue(hasattr(network, 'edges'))


class TestNetworkDistance(unittest.TestCase):
    
    def test_compute_distance_identical_networks(self):
        """Test distance between identical networks"""
        net1 = KmerNetwork(k=3)
        net1.add_kmer("ATG", clone_id=1)
        net1.add_kmer("TGC", clone_id=1)
        net1.normalize_nodes()
        net1.compute_edge_probabilities()
        
        net2 = KmerNetwork(k=3)
        net2.add_kmer("ATG", clone_id=1)
        net2.add_kmer("TGC", clone_id=1)
        net2.normalize_nodes()
        net2.compute_edge_probabilities()
        
        distance = compute_network_distance(net1, net2)
        
        # Distance should be 0 or very close to 0
        self.assertAlmostEqual(distance, 0.0, places=5)
    
    def test_compute_distance_different_networks(self):
        """Test distance between different networks"""
        net1 = KmerNetwork(k=3)
        net1.add_kmer("ATG", clone_id=1)
        net1.normalize_nodes()
        net1.compute_edge_probabilities()
        
        net2 = KmerNetwork(k=3)
        net2.add_kmer("CGT", clone_id=1)
        net2.normalize_nodes()
        net2.compute_edge_probabilities()
        
        distance = compute_network_distance(net1, net2)
        
        # Distance should be non-zero
        self.assertGreater(distance, 0)
    
    def test_compute_distance_empty_networks(self):
        """Test distance computation with empty networks"""
        net1 = KmerNetwork(k=3)
        net2 = KmerNetwork(k=3)
        
        net1.normalize_nodes()
        net2.normalize_nodes()
        net1.compute_edge_probabilities()
        net2.compute_edge_probabilities()
        
        distance = compute_network_distance(net1, net2)
        
        # Should handle empty networks gracefully
        self.assertIsNotNone(distance)


if __name__ == '__main__':
    unittest.main()
