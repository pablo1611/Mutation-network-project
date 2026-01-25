# Project Diagrams

This directory contains UML diagrams documenting the architecture and behavior of the Antibody Sequence Loader application.

## Diagrams

### 1. Class Diagram (`class_diagram.puml`)
Shows the main classes and their relationships:
- **Clone**: Represents antibody clone with sequence data
- **CodonMapper**: Translates DNA sequences to amino acids
- **KmerNetwork**: Core data structure for k-mer networks
- **NetworkCollection**: Manages multiple network types
- **AntibodySequenceLoaderApp**: Main UI application
- Utility functions for data processing

**Key Relationships:**
- Clone uses CodonMapper for translation
- NetworkCollection contains KmerNetwork instances
- App orchestrates all components

### 2. Use Case Diagram (`usecase_diagram.puml`)
Illustrates user interactions with the system:

**Main Use Cases:**
- **Data Loading**: Browse and load CSV datasets
- **Network Analysis**: Process sequences with configurable thresholds
- **Visualization**: Generate interactive network plots
- **Data Export**: Save triplets and edge data
- **Comparison**: Compare datasets or regions
- **Multi-Dataset Management**: Handle multiple datasets

**Actors:**
- Researcher: Primary user for analysis
- Data Analyst: Focus on export and comparison

### 3. Sequence Diagram (`sequence_diagram.puml`)
Details the flow of dataset analysis:

**Process Flow:**
1. **Data Loading**: CSV → Clones with metadata
2. **Extraction**: Extract nonuplets and translate to triplets
3. **Network Creation**: Build k-mer networks (AA and NT)
4. **Processing**: Normalize, apply thresholds, compute edges
5. **Export**: Save edges to CSV, generate HTML visualization
6. **Comparison**: Calculate distance between networks (optional)

## Viewing the Diagrams

### Online Viewers
Upload `.puml` files to:
- [PlantUML Online Editor](http://www.plantuml.com/plantuml/uml/)
- [PlantText](https://www.planttext.com/)

### VS Code
Install the PlantUML extension:
```
code --install-extension jebbs.plantuml
```
Then open any `.puml` file and press `Alt+D` to preview.

### Command Line
If you have PlantUML installed:
```bash
# Generate PNG
plantuml diagrams/class_diagram.puml

# Generate SVG
plantuml -tsvg diagrams/*.puml

# Generate all formats
plantuml -tpng -tsvg diagrams/*.puml
```

### Installation
**macOS:**
```bash
brew install plantuml
```

**Linux:**
```bash
sudo apt-get install plantuml
```

**Windows:**
Download from [PlantUML website](https://plantuml.com/download)

## Diagram Formats

All diagrams use PlantUML syntax (`.puml` files):
- Text-based and version control friendly
- Easy to modify and maintain
- Can generate PNG, SVG, PDF, and more
- Cross-platform compatibility

## Updating Diagrams

When modifying the code architecture:
1. Update the corresponding `.puml` file
2. Regenerate the images if needed
3. Commit both the `.puml` source and generated images

## Architecture Notes

### Key Design Decisions
- **Separation of Concerns**: UI, data processing, and visualization are separate
- **Flexible Thresholds**: Node and edge thresholds configurable per dataset
- **Multiple K-mer Types**: Supports both amino acid (k=3) and nucleotide (k=9) networks
- **Region Comparison**: Can analyze CDR1 vs CDR2 regions within same dataset
- **Export Everything**: All intermediate data (triplets, edges) exportable

### Data Flow
```
CSV → Clones → Nonuplets → Triplets → K-mer Networks → Thresholded Networks → Visualization/Export
```

### Network Distance Metric
Uses Jaccard distance on node sets combined with frequency comparison for quantitative network similarity.
