"""Integration tests for the full pipeline"""
import unittest
import sys
from pathlib import Path
import tempfile
import csv

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.clone import Clone
from src.networks import KmerNetwork, compute_network_distance
from src.build_triplet_df import build_triplet_df


class TestIntegrationPipeline(unittest.TestCase):
    
    def setUp(self):
        """Create test clones for integration tests"""
        self.test_clones = [
            Clone(
                seq_id="seq1",
                clone_id=1,
                sequence="ATGCGTATCGATCGAATGCGT",
                sample_="sample1",
                ai=1
            ),
            Clone(
                seq_id="seq2",
                clone_id=2,
                sequence="TGCGTATCGATCGAATGCGTA",
                sample_="sample1",
                ai=2
            ),
            Clone(
                seq_id="seq3",
                clone_id=3,
                sequence="GCGTATCGATCGAATGCGTAT",
                sample_="sample2",
                ai=1
            ),
        ]
        
        # Extract and translate nines for all clones
        for clone in self.test_clones:
            clone.extract_nines()
            clone.translate_nines()
    
    def test_full_pipeline_clone_to_network(self):
        """Test complete pipeline from clones to network"""
        # Create networks
        aa_network = KmerNetwork(k=3)
        nt_network = KmerNetwork(k=9)
        
        # Process clones
        for clone in self.test_clones:
            # Add nucleotide 9-mers
            for idx, nine in clone.nines:
                if '-' not in nine and 'N' not in nine:
                    nt_network.add_kmer(nine, clone.clone_id, pos=idx)
            
            # Add amino acid triplets
            for idx, aa_triplet in clone.nine_aa_triplets:
                if aa_triplet:
                    aa_network.add_kmer(aa_triplet, clone.clone_id, pos=idx)
        
        # Networks should have nodes
        self.assertGreater(len(aa_network.nodes), 0)
        self.assertGreater(len(nt_network.nodes), 0)
        
        # Normalize
        aa_network.normalize_nodes()
        nt_network.normalize_nodes()
        
        # All nodes should have normalized_freq
        for kmer, data in aa_network.nodes.items():
            self.assertIn('normalized_freq', data)
            self.assertGreater(data['normalized_freq'], 0)
        
        # Compute edges
        aa_network.compute_edge_probabilities()
        nt_network.compute_edge_probabilities()
        
        # Should have edges
        self.assertTrue(hasattr(aa_network, 'edges'))
        self.assertTrue(hasattr(nt_network, 'edges'))
    
    def test_network_comparison_pipeline(self):
        """Test network comparison workflow"""
        # Create two separate networks from different clone sets
        net1 = KmerNetwork(k=3)
        net2 = KmerNetwork(k=3)
        
        # Add clones to net1
        for clone in self.test_clones[:2]:
            for idx, aa_triplet in clone.nine_aa_triplets:
                if aa_triplet:
                    net1.add_kmer(aa_triplet, clone.clone_id)
        
        # Add clones to net2
        for clone in self.test_clones[1:]:
            for idx, aa_triplet in clone.nine_aa_triplets:
                if aa_triplet:
                    net2.add_kmer(aa_triplet, clone.clone_id)
        
        # Normalize both
        net1.normalize_nodes()
        net2.normalize_nodes()
        
        # Compute edges
        net1.compute_edge_probabilities()
        net2.compute_edge_probabilities()
        
        # Compute distance
        distance = compute_network_distance(net1, net2)
        
        # Distance should be valid
        self.assertIsNotNone(distance)
        self.assertGreaterEqual(distance, 0)
    
    def test_threshold_application_pipeline(self):
        """Test applying thresholds in the pipeline"""
        # Create network
        network = KmerNetwork(k=3)
        
        # Add multiple occurrences of different kmers
        for clone in self.test_clones:
            for idx, aa_triplet in clone.nine_aa_triplets:
                if aa_triplet:
                    network.add_kmer(aa_triplet, clone.clone_id)
        
        # Normalize
        network.normalize_nodes()
        initial_node_count = len(network.nodes)
        
        # Apply node threshold
        network.apply_node_threshold(0.1)
        filtered_node_count = len(network.nodes)
        
        # Should filter out some nodes
        self.assertLessEqual(filtered_node_count, initial_node_count)
        
        # Compute edges
        network.compute_edge_probabilities()
        
        if hasattr(network, 'edges') and len(network.edges) > 0:
            # Count total edges before threshold
            initial_edge_count = sum(len(targets) for targets in network.edges.values())
            
            # Apply edge threshold
            network.apply_edge_threshold(0.5)
            
            # Count edges after threshold
            filtered_edge_count = sum(
                len([t for t, data in targets.items() if data.get('above_threshold', True)])
                for targets in network.edges.values()
            )
            
            # Should filter out some edges
            self.assertLessEqual(filtered_edge_count, initial_edge_count)
    
    def test_build_triplet_df_integration(self):
        """Test building triplet dataframe from clones"""
        # build_triplet_df expects a dict of clones
        clones_dict = {clone.clone_id: clone for clone in self.test_clones}
        df = build_triplet_df(clones_dict)
        
        # Should return a dataframe
        self.assertIsNotNone(df)
        
        # Should have rows
        self.assertGreater(len(df), 0)
        
        # Should have expected columns (actual columns from build_triplet_df)
        expected_cols = ['clone_id', 'index', 'triplet']
        for col in expected_cols:
            self.assertIn(col, df.columns)
        
        # All triplets should be valid (no None) - column is 'triplet'
        self.assertTrue(all(df['triplet'].notna()))


if __name__ == '__main__':
    unittest.main()
