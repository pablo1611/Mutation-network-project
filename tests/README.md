# Test Suite

This directory contains unit tests and integration tests for the Antibody Sequence Loader project.

## Test Files

- **test_clone.py** - Tests for the Clone class
  - Initialization with required and extra fields
  - Nonuplet extraction (sliding window)
  - Translation of nonuplets to amino acid triplets
  - Handling of invalid sequences (gaps, N's)

- **test_codon_mapper.py** - Tests for the CodonMapper class
  - Single codon translation
  - 9-mer translation
  - Stop codon handling
  - Invalid input handling

- **test_networks.py** - Tests for KmerNetwork and network distance
  - Network initialization
  - Adding kmers and sequences
  - Node normalization
  - Edge probability computation
  - Threshold application
  - Network distance calculation

- **test_ui_functions.py** - Tests for UI utility functions
  - Output directory creation
  - CSV edge export functionality
  - File content validation

- **test_integration.py** - End-to-end integration tests
  - Full pipeline from clones to networks
  - Network comparison workflow
  - Threshold application in complete pipeline
  - Triplet dataframe generation

## Running Tests

### Run all tests:
```bash
python tests/run_all_tests.py
```

### Run specific test file:
```bash
python -m unittest tests/test_clone.py
```

### Run with pytest (if installed):
```bash
pytest tests/
pytest tests/ -v  # verbose output
pytest tests/test_clone.py  # specific file
```

### Run from project root:
```bash
cd /path/to/Mutation-network-project
python -m unittest discover -s tests -p "test_*.py"
```

## Test Coverage

The test suite covers:
- ✅ Core data structures (Clone, KmerNetwork)
- ✅ Sequence processing (extraction, translation)
- ✅ Network operations (normalization, thresholding)
- ✅ Edge computation and filtering
- ✅ Network comparison metrics
- ✅ UI utility functions
- ✅ File I/O operations
- ✅ End-to-end workflows

## Adding New Tests

To add new tests:
1. Create a new file named `test_*.py` in this directory
2. Import unittest and the modules you want to test
3. Create test classes that inherit from `unittest.TestCase`
4. Write test methods starting with `test_`
5. Run the test suite to verify

Example:
```python
import unittest
from src.your_module import YourClass

class TestYourClass(unittest.TestCase):
    def test_something(self):
        obj = YourClass()
        self.assertEqual(obj.method(), expected_value)
```
