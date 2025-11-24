from collections import defaultdict
import json
import heapq
from typing import Iterable, Tuple, Optional

from src.codon_mapper import CodonMapper

# Keep a module-level CodonMapper so it's available when needed but not
# attached to the NetworksManager instance returned to callers.
_CODON_MAPPER: Optional[CodonMapper] = None


class KmerNetwork:
    """
    Track counts and occurrences for a given k (k-mer network).

    Internal structure:
      nodes: dict mapping kmer -> {
          'count': int,            # total occurrences across all clones
          'clones': set(int),      # set of clone_id values where kmer appeared
          'positions': { clone_id: [pos, ...], ... } # positions per clone
      }

    Note: sets are used internally for fast deduplication; use `to_serializable` to convert to JSON-friendly types.
    """

    def __init__(self, k, store_positions=False):
        """
        store_positions: if True, store positions per clone (memory heavy). Default False.
        """
        self.k = k
        self.store_positions = store_positions
        self.nodes = {}

    def add_kmer(self, kmer, clone_id, pos=None):
        node = self.nodes.get(kmer)
        if node is None:
            # only create positions container when requested
            node = {'count': 0, 'clones': set()}
            if self.store_positions:
                node['positions'] = defaultdict(list)
            self.nodes[kmer] = node
        node['count'] += 1
        node['clones'].add(clone_id)
        if self.store_positions and pos is not None:
            node['positions'][clone_id].append(pos)

    def add_sequence(self, sequence, clone_id):
        """Add all k-mers from sequence (sliding window) to the network.

        If positions are enabled, records positions, otherwise only updates counts and clone sets.
        """
        if not sequence or len(sequence) < self.k:
            return
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            # skip ambiguous/reserved characters
            if '-' in kmer or 'N' in kmer:
                continue
            # pass pos only when storing positions to keep node signature simple
            if self.store_positions:
                self.add_kmer(kmer, clone_id, i)
            else:
                self.add_kmer(kmer, clone_id)

    def get_node(self, kmer):
        """Return the raw node structure for a kmer, or None if missing."""
        return self.nodes.get(kmer)

    def node_count(self) -> int:
        return len(self.nodes)

    def total_occurrences(self) -> int:
        return sum(node['count'] for node in self.nodes.values())

    def top_kmers(self, n: int = 20) -> Iterable[Tuple[str, int, int]]:
        """Return top-n kmers by count: (kmer, count, clone_count)"""
        if not self.nodes:
            return []
        it = ((node['count'], kmer) for kmer, node in self.nodes.items())
        largest = heapq.nlargest(n, it)
        result = []
        for count, kmer in largest:
            node = self.nodes[kmer]
            clone_count = len(node.get('clones', ()))
            result.append((kmer, count, clone_count))
        return result

    def sample_kmers(self, n: int = 20):
        items = list(self.nodes.items())[:n]
        return [(kmer, node['count']) for kmer, node in items]

    def to_serializable(self):
        """Return JSON-serializable representation (converts sets and defaultdicts)."""
        out = {}
        for kmer, node in self.nodes.items():
            entry = {
                'count': node['count'],
                'clones': sorted(node['clones'])
            }
            if self.store_positions:
                entry['positions'] = {str(cid): pos_list[:] for cid, pos_list in node['positions'].items()}
            out[kmer] = entry
        return out

    def dump_json(self, filepath):
        data = self.to_serializable()
        with open(filepath, 'w') as fh:
            json.dump({'k': self.k, 'nodes': data}, fh)


class NetworksManager:
    """Manage networks.

    Design change: nucleotide 3-mers are NOT created. Instead, we maintain:
      - nucleotide k-mer networks for requested k values (default only k=9)
      - an optional amino-acid triplet network (aa3) derived from sliding 9-mers
    """

    def __init__(self, ks=(9,), store_positions=False):
        """Create networks.

        ks: iterable of nucleotide k values to build (3 will be ignored).
        """
        # ignore nucleotide k==3 (we no longer track nucleotide triplets)
        ks_filtered = [k for k in ks if k != 3]
        self.nucleotide_network = {k: KmerNetwork(k, store_positions=store_positions) for k in ks_filtered}
        self.store_positions = store_positions

        # amino-acid triplet network (k=3 in amino-acid alphabet)
        self.aa_network: Optional[KmerNetwork] = KmerNetwork(3, store_positions=store_positions)

    def add_sequence(self, sequence, clone_id):
        # nucleotide networks
        for k, net in self.nucleotide_network.items():
            net.add_sequence(sequence, clone_id)

        # amino-acid triplet network built from sliding 9-nt windows
        if self.aa_network:
            global _CODON_MAPPER
            if _CODON_MAPPER is None:
                _CODON_MAPPER = CodonMapper()
            for i, nine, aa_triplet in _CODON_MAPPER.map_sequence(sequence):
                # skip ambiguous translations
                if 'X' in aa_triplet:
                    continue
                # record aa triplet; pos is nucleotide index i
                self.aa_network.add_kmer(aa_triplet, clone_id, pos=i)

    def get_network(self, k):
        return self.nucleotide_network.get(k)

    def get_aa_network(self):
        return self.aa_network

    def to_serializable(self):
        out = {k: net.to_serializable() for k, net in self.nucleotide_network.items()}
        if self.aa_network:
            out['aa3'] = self.aa_network.to_serializable()
        return out

    def dump_json(self, basepath):
        for k, net in self.nucleotide_network.items():
            path = f"{basepath}.k{k}.json"
            net.dump_json(path)
        if self.aa_network:
            self.aa_network.dump_json(f"{basepath}.aa3.json")

    def summary(self) -> dict:
        out = {}
        for k, net in self.nucleotide_network.items():
            out[k] = {'nodes': net.node_count(), 'occurrences': net.total_occurrences()}
        if self.aa_network:
            out['aa3'] = {'nodes': self.aa_network.node_count(), 'occurrences': self.aa_network.total_occurrences()}
        return out
