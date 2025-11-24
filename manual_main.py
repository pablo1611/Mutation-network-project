from src.load_clones import load_clones

if __name__ == "__main__":
    import time
    import os
    from dotenv import load_dotenv

    # load .env from project root (where this script is run)
    load_dotenv()

    csv_path = os.getenv('CLONES_CSV')
    if not csv_path:
        # FALLBACK: raise error if not set
        csv_path = "Users/leepotashnik/Desktop/study/semester 8+9/Final project/cleaned_seqs.csv"

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
    clone_id =646357
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
