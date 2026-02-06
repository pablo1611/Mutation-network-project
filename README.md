# Mutation Network Project

A Python-based application for analyzing antibody sequence mutation networks. This tool loads antibody clone data from CSV files, builds triplet-based networks, computes distance metrics between networks, and provides interactive visualizations. It features a user-friendly GUI for dataset management, analysis, and comparison.

## Features

- **Dataset Loading**: Load and preview antibody clone data from CSV files with metadata extraction (subjects, samples, time points, etc.).
- **Network Analysis**: Build k-mer networks (triplets and nonuplets) from sequences and compute dissimilarity distances.
- **Visualization**: Interactive Plotly-based AA3 network plots for exploring mutation patterns.
- **GUI Interface**: CustomTkinter-based app for easy dataset upload, analysis, and comparison (supports up to 2 datasets).
- **Comparison Tools**: Compare networks between datasets or within a dataset (e.g., R1 vs. R2 regions).
- **Export Functionality**: Save triplet occurrence data as CSV files and open them automatically.
- **Automatic Edge Export**: Every analysis automatically saves network edges with their calculated probabilities (weights) to CSV files for further analysis.
- **Cross-Platform**: Works on macOS and Windows; package as native executables using PyInstaller.


## Installation

### Prerequisites
- Python 3.12 or higher
- Virtual environment tool (venv)

### Setup
1. Clone or download the project:
   ```
   git clone <repository-url>
   cd mutation-network-project
   ```

2. Create and activate a virtual environment:
   ```
   python -m venv .venv
   source .venv/bin/activate  # On macOS/Linux
   # or on Windows: .venv\Scripts\activate
   ```

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

### Running the Application
```bash
python ui/main_window.py
```

The GUI allows you to:
1. **Upload Dataset(s)**: Select CSV files containing antibody clone data
2. **View Metadata**: See extracted information about subjects, samples, and time points
3. **Execute Triplets**: Generate and export triplet occurrence data
4. **Analyze Network**: Visualize amino acid triplet networks with interactive plots
5. **Compare**: Analyze differences between datasets or regions (R1 vs R2)

### Output Files

All output files are saved to the `output/` directory:

- **Triplet Data**: `triplets_dataset_1.csv`, `triplets_dataset_2.csv` - K-mer occurrence data
- **Network Plots**: `aa3_network_dataset_1.html`, `aa3_network_dataset_2.html` - Interactive visualizations
- **Edge Datasets** (automatic): CSV files containing all network edges with calculated probabilities:
  - `edges_dataset_1.csv` - Edges from single dataset analysis
  - `edges_dataset_2.csv` - Edges from second dataset analysis
  - `edges_dataset1_comparison.csv`, `edges_dataset2_comparison.csv` - Edges from dataset comparison
  - `edges_region1_comparison.csv`, `edges_region2_comparison.csv` - Edges from region comparison

### Edge Dataset Structure

Each edge CSV file contains:
- **source_node** & **target_node**: Connected amino acid triplets
- **weight** / **probability**: Calculated edge probability (0.0 to 1.0)
- **above_threshold**: Whether edge passes the threshold filter
- **dataset**: Dataset or region identifier
- **analysis_type**: Type of analysis performed
- **timestamp**: When the analysis was run

Example usage with pandas:
```python
import pandas as pd
edges = pd.read_csv('output/edges_dataset_1.csv')
top_edges = edges.nlargest(10, 'weight')  # Get top 10 strongest connections
```

For detailed documentation on the edge dataset feature, see [EDGES_DATASET_DOCUMENTATION.md](EDGES_DATASET_DOCUMENTATION.md).

## Building the Application

To create standalone executables for distribution:

### On macOS
```
pyinstaller --windowed --onefile ui/main_window.py
```
- Generates `dist/main_window.app` (native macOS app).

### On Windows
On a Windows machine:
```
pyinstaller --windowed --onefile ui/main_window.py
```
- Generates `dist/main_window.exe` (Windows executable).

**Note**: Builds must be done on the target OS. The `--windowed` flag hides the console for GUI apps.

## Data Format

CSV files should contain columns like:
- `seq_id`, `sample_id`, `subject_id`, `clone_id`, `functional`, `copy_number`, `cdr3_aa`, `sequence`, `germline`, `ab_target`, `time_point`

The app filters to the highest `copy_number` per `clone_id` for analysis.

## Note on Analysis Duration

Important: analyses may take a while depending on dataset size and selected options. Large datasets or extended comparisons can run from several minutes up to hours. Do not close the application, terminal, or interrupt the process while an analysis is running; wait for it to finish. Output files will be written to the `output/` directory when the run completes.
