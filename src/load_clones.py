import pandas as pd
from src.clone import Clone



def load_clones(csv_path):
    """
    Loads clone data from a CSV file, filters to keep only the row with the maximum copy_number for each clone_id,
    and returns a dictionary mapping clone_id to Clone objects for fast lookup.

    Parameters:
        csv_path (str): Path to the CSV file containing clone data. The file must have columns:
            seq_id, sample_id, subject_id, clone_id, functional, copy_number, cdr3_aa, sequence, germline, ab_target, time_point

    Returns:
        dict: Dictionary where keys are clone_id values and values are Clone objects representing each unique clone_id
              with the highest copy_number in the dataset.

    Example:
        clones_dict = load_clones('data/clones.csv')
        clone = clones_dict['646357']  # Access Clone object by clone_id
    """
    df = pd.read_csv(csv_path)
    # Filter to keep only the row with the max copy_number for each clone_id
    filtered_df = df.loc[df.groupby('clone_id')['copy_number'].idxmax()]
    clones_dict = {}
    for _, row in filtered_df.iterrows():
        clone = Clone(
            seq_id=row.get('seq_id'),
            ai=row.get('ai'),
            sample_=row.get('sample_id'),
            subject_=row.get('subject_id'),
            clone_id=row.get('clone_id'),
            function=row.get('functional'),
            copy_nu=row.get('copy_number'),
            cdr3_aa=row.get('cdr3_aa'),
            sequence=row.get('sequence'),
            germline=row.get('germline'),
            ab_target=row.get('ab_target'),
            time_po=row.get('time_point')
        )
        clones_dict[row.get('clone_id')] = clone
    return clones_dict

