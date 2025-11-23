from src.load_clones import load_clones

if __name__ == "__main__":
    clones = load_clones(r"C:\Studies\Semester 2 24-25\project\cleaned_seqs.csv")
    print(f"Loaded {len(clones)} clones:")
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
    else:
        print(f"Clone with clone_id {clone_id} not found.")
