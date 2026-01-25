"""Unit tests for UI utility functions"""
import unittest
import sys
from pathlib import Path
import tempfile
import os

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from ui.main_window import save_edges_to_csv, get_output_dir
from src.networks import KmerNetwork


class TestUIFunctions(unittest.TestCase):
    
    def test_get_output_dir_creates_directory(self):
        """Test that get_output_dir creates directory"""
        output_dir = get_output_dir()
        
        # Should return a valid path
        self.assertIsNotNone(output_dir)
        self.assertIsInstance(output_dir, str)
        
        # Directory should exist
        self.assertTrue(os.path.exists(output_dir))
        self.assertTrue(os.path.isdir(output_dir))
    
    def test_save_edges_to_csv_basic(self):
        """Test saving edges to CSV"""
        # Create a simple network with edges
        network = KmerNetwork(k=3)
        # Add as sequence to create co-occurrences
        network.add_sequence("ATGTGCGCG", clone_id=1)
        network.add_sequence("ATGTGCGCG", clone_id=2)
        
        network.normalize_nodes()
        network.compute_edge_probabilities()
        
        # Create a temp file for output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            output_path = f.name
        
        try:
            # Save edges
            result = save_edges_to_csv(
                network,
                output_path,
                dataset_name="test_dataset",
                analysis_type="test"
            )
            
            # If network has edges, should return path
            if result is not None:
                self.assertEqual(result, output_path)
                # File should exist
                self.assertTrue(os.path.exists(output_path))
                
                # File should have content
                self.assertGreater(os.path.getsize(output_path), 0)
                
                # Should be readable as CSV
                import pandas as pd
                df = pd.read_csv(output_path)
                
                # Should have expected columns
                expected_columns = ['source_node', 'target_node', 'weight', 
                                  'above_threshold', 'dataset', 'analysis_type', 'timestamp']
                for col in expected_columns:
                    self.assertIn(col, df.columns)
                
                # Should have rows
                # Should have rows
            self.assertGreater(len(df), 0)
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_save_edges_to_csv_empty_network(self):
        """Test saving edges from empty network"""
        network = KmerNetwork(k=3)
        network.normalize_nodes()
        network.compute_edge_probabilities()
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            output_path = f.name
        
        try:
            # Should handle empty network gracefully
            result = save_edges_to_csv(network, output_path)
            
            # Should return None for empty network
            self.assertIsNone(result)
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_save_edges_csv_content(self):
        """Test CSV content is correct"""
        network = KmerNetwork(k=3)
        # Add as sequence to create co-occurrences
        network.add_sequence("ATGTGCGCGATGTGC", clone_id=1)
        network.add_sequence("ATGTGCGCGATGTGC", clone_id=2)
        
        network.normalize_nodes()
        network.compute_edge_probabilities()
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            output_path = f.name
        
        try:
            result = save_edges_to_csv(
                network,
                output_path,
                dataset_name="TestDataset",
                analysis_type="single_dataset"
            )
            
            if result is not None:
                import pandas as pd
                df = pd.read_csv(output_path)
                
                # Check dataset name is correct
                self.assertTrue(all(df['dataset'] == "TestDataset"))
                
                # Check analysis type
                self.assertTrue(all(df['analysis_type'] == "single_dataset"))
                
                # Check weights are numeric
                self.assertTrue(pd.api.types.is_numeric_dtype(df['weight']))
                
                # Weights should be probabilities (0-1)
                self.assertTrue(all(df['weight'] >= 0))
                self.assertTrue(all(df['weight'] >= 0))
            self.assertTrue(all(df['weight'] <= 1))
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)


if __name__ == '__main__':
    unittest.main()
