# Edge Dataset Auto-Save Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                    USER ACTIONS                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ├──────┐
                              │      │
          ┌───────────────────┤      ├───────────────────┐
          │                   │      │                   │
          ▼                   ▼      ▼                   ▼
    ┌─────────┐        ┌──────────────────┐      ┌──────────────┐
    │ Analyze │        │ Compare Between  │      │   Compare    │
    │ Network │        │    Datasets      │      │    Areas     │
    └─────────┘        └──────────────────┘      └──────────────┘
          │                   │                          │
          │                   │                          │
┌─────────────────────────────────────────────────────────────────┐
│               AUTOMATIC PROCESSING                               │
├─────────────────────────────────────────────────────────────────┤
│  1. Load network data                                            │
│  2. Normalize node frequencies                                   │
│  3. Compute edge probabilities                                   │
│  4. Apply thresholds                                             │
│  5. Call save_edges_to_csv()  ◄── NEW FEATURE                   │
│  6. Generate visualization                                       │
└─────────────────────────────────────────────────────────────────┘
          │                   │                          │
          ▼                   ▼                          ▼
    ┌─────────┐        ┌──────────────────┐      ┌──────────────┐
    │ Dataset │        │  Dataset 1 & 2   │      │  Region 1&2  │
    │ Edges   │        │     Edges        │      │    Edges     │
    └─────────┘        └──────────────────┘      └──────────────┘
          │                   │                          │
          ▼                   ▼                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                    OUTPUT FILES                                  │
├─────────────────────────────────────────────────────────────────┤
│  edges_dataset_1.csv                                             │
│  edges_dataset_2.csv                                             │
│  edges_dataset1_comparison.csv                                   │
│  edges_dataset2_comparison.csv                                   │
│  edges_region1_comparison.csv                                    │
│  edges_region2_comparison.csv                                    │
└─────────────────────────────────────────────────────────────────┘
```

## CSV File Structure

```
┌────────────┬────────────┬────────┬────────────┬────────────────┬───────────┬───────────────┬───────────────┐
│ source_node│ target_node│ weight │ probability│ above_threshold│  dataset  │ analysis_type │  timestamp    │
├────────────┼────────────┼────────┼────────────┼────────────────┼───────────┼───────────────┼───────────────┤
│    AAA     │    AAC     │ 0.523  │   0.523    │     True       │ Dataset_1 │ single_...    │ 2026-01-24... │
│    AAC     │    AAA     │ 0.478  │   0.478    │     True       │ Dataset_1 │ single_...    │ 2026-01-24... │
│    AAA     │    ACA     │ 0.234  │   0.234    │     False      │ Dataset_1 │ single_...    │ 2026-01-24... │
│    ...     │    ...     │  ...   │    ...     │      ...       │    ...    │     ...       │     ...       │
└────────────┴────────────┴────────┴────────────┴────────────────┴───────────┴───────────────┴───────────────┘
                                    ▲
                                    │
                            This is the WEIGHT
                        (probability of edge connection)
```

## Edge Probability Calculation

For each node `i` and its neighbor `j`:

```
         frequency(j)
P(j→i) = ────────────────────────────
         Σ frequency(all neighbors)


Example:
  Node AAA has neighbors: AAC (freq=50), AAT (freq=30), ACA (freq=20)
  
  P(AAC→AAA) = 50 / (50+30+20) = 50/100 = 0.50
  P(AAT→AAA) = 30 / (50+30+20) = 30/100 = 0.30
  P(ACA→AAA) = 20 / (50+30+20) = 20/100 = 0.20
```

## Data Flow

```
Clone Data (CSV)
      ↓
Load Clones
      ↓
Build Networks
      ↓
Calculate K-mers
      ↓
Compute Node Frequencies
      ↓
Compute Edge Probabilities  ◄── Each edge gets a weight
      ↓
Apply Thresholds
      ↓
Save to Pandas DataFrame    ◄── NEW: Automatic save
      ↓
Export as CSV
      ↓
output/edges_*.csv
```
