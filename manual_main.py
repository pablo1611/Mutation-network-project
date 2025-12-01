from src.load_clones import load_clones
from src.build_triplet_df import build_triplet_df
import time
import os
    

if __name__ == "__main__":
    

    # load .env from project root (where this script is run)
    #load_dotenv()

    #csv_path = os.getenv('CLONES_CSV')
    #if not csv_path:
     #   csv_path = "Users/leepotashnik/Desktop/study/semester 8+9/Final project/cleaned_seqs.csv"
    csv_path = r"C:\Studies\Semester 2 24-25\project\cleaned_seqs.csv"
    print("Starting clone loading...")
    start_time = time.time()
    result = load_clones(csv_path, True)
    # load_clones returns either clones_dict or (clones_dict, networks) when build_networks=True
    if isinstance(result, tuple):
        clones, networks = result
    else:
        clones = result
        networks = None

    print(f"Loaded {len(clones)} clones:")
    print(f"took {time.time() - start_time:.2f} seconds.")
    #for clone in clones:
    #    print(clone)
    # Build triplet DataFrame for all clones
    start_time = time.time()
    triplet_df = build_triplet_df(clones)
    print(f"Building triplet DataFrame took {time.time() - start_time:.2f} seconds.")
    print(triplet_df.head())

    # Verifier: sum of non-None triplets for all clones should equal DataFrame length
    total_triplets = sum(
        sum(1 for idx, aa_triplet in clone.nine_aa_triplets if aa_triplet is not None)
        for clone in clones.values()
    )
    print(f"Verifier: sum of non-None triplets across all clones = {total_triplets}")
    print(f"Triplet DataFrame length = {len(triplet_df)}")
    print(f"Match: {total_triplets == len(triplet_df)}")


    clone_id =711418
    clone = clones.get(clone_id)
    if clone:
        clone.extract_nines()
        clone.translate_nines()
        print(f"Nines for clone_id {clone_id}:")
        print(clone.nines)
        print("and matching triplet: ")
        print(clone.nine_aa_triplets)
        """
        # If networks were built, print aa3 network nodes and counts
        if networks:
            aa_net = networks.get_aa_network()
            if aa_net:
                print("\nAA network nodes (key, count):")
                try:
                    total = aa_net.node_count()
                except Exception:
                    total = None
                print(f"Total aa nodes: {total}")
                # iterate deterministically
                for kmer in sorted(aa_net.nodes.keys()):
                    node = aa_net.nodes[kmer]
                    print(kmer, node.get('count'))
        """
    else:
        print(f"Clone with clone_id {clone_id} not found.")
    print(f"Finished in {time.time() - start_time:.2f} seconds.")
