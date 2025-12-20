from src.load_clones import load_clones
from src.build_triplet_df import build_triplet_df
from src.networks import compute_network_distance
import time
import os
from dotenv import load_dotenv
    

if __name__ == "__main__":
    

    # load .env from project root
    load_dotenv()

    csv_path = os.getenv('CLONES_CSV')
    if not csv_path:
        csv_path = "/Users/leepotashnik/Desktop/study/semester 8+9/Final project/cleaned_seqs.csv"
    print("Starting clone loading...")
    start_time1 = time.time()
    result = load_clones(csv_path, True)
    # load_clones returns either clones_dict or (clones_dict, networks) when build_networks=True
    if isinstance(result, tuple):
        clones, networks = result
    else:
        clones = result
        networks = None

    print(f"Loaded {len(clones)} clones:")
    print(f"took {time.time() - start_time1:.2f} seconds.")
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
    print(f"Finished cloning in {time.time() - start_time1:.2f} seconds.")

    # --- EXPORT: write occurrence-level triplet DFs and aggregated network CSVs ---
# Uncomment and edit `out_dir` before running.

    out_dir = os.path.expanduser("/Users/leepotashnik/Desktop/study/semester 8+9/Final project")
    os.makedirs(out_dir, exist_ok=True)

    # 1) occurrence-level AA triplet dataframes (exact shape as build_triplet_df)
    triplets_r1 = build_triplet_df(clones, region='r1')   # triplet, clone_id, index
    triplets_r2 = build_triplet_df(clones, region='r2')
    triplets_r1.to_csv(os.path.join(out_dir, "triplets_r1_occurrences.csv"), index=False)
    triplets_r2.to_csv(os.path.join(out_dir, "triplets_r2_occurrences.csv"), index=False)

    # Compute distance between R1 and R2 AA triplet networks
    if networks:
        aa_net_r1 = networks.aa_network_region['r1']
        aa_net_r2 = networks.aa_network_region['r2']
        distance = compute_network_distance(aa_net_r1, aa_net_r2)
        print(f"Distance between R1 and R2 AA triplet networks: {distance:.4f}")

    # 2) aggregated node-level CSVs for networks (kmer, count, clones)
    # Uses the compact aggregated CSV writer added to NetworksManager
    # Writes:
    #   base.aa3.r1.aggregated.csv, base.aa3.r2.aggregated.csv,
    #   base.k9.r1.aggregated.csv, base.k9.r2.aggregated.csv
    if networks:
        # choose base filename (no extension)
        base = os.path.join(out_dir, "networks_base")
        networks.dump_region_aggregated_csvs(base)

    print("Exports written to:", out_dir)
    print("total time:", time.time() - start_time1, "seconds.")
