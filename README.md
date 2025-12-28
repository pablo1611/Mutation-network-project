# Mutation Network Project

A Python-based application for analyzing antibody sequence mutation networks. This tool loads antibody clone data from CSV files, builds triplet-based networks, computes distance metrics between networks, and provides interactive visualizations. It features a user-friendly GUI for dataset management, analysis, and comparison.

## Features

- **Dataset Loading**: Load and preview antibody clone data from CSV files with metadata extraction (subjects, samples, time points, etc.).
- **Network Analysis**: Build k-mer networks (triplets and nonuplets) from sequences and compute dissimilarity distances.
- **Visualization**: Interactive Plotly-based AA3 network plots for exploring mutation patterns.
- **GUI Interface**: CustomTkinter-based app for easy dataset upload, analysis, and comparison (supports up to 2 datasets).
- **Comparison Tools**: Compare networks between datasets or within a dataset (e.g., R1 vs. R2 regions).
- **Export Functionality**: Save triplet occurrence data as CSV files and open them automatically.
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
