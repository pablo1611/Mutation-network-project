from collections import defaultdict
import json
import heapq
from typing import Iterable, Tuple, Optional
import csv
import pandas as pd

from src.codon_mapper import CodonMapper

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
        # mark nucleotide 9-mers that contain stop codons with a color flag
        if self.k == 9:
            # ensure uppercase
            kmer_up = kmer.upper()
            try:
                codons = [kmer_up[j:j+3] for j in (0, 3, 6)]
            except Exception:
                codons = []
            stops = {"TAA", "TAG", "TGA"}
            if any(c in stops for c in codons):
                # set color flag to 1 (only set once)
                if 'color' not in node:
                    node['color'] = 1
        # mark amino-acid triplets that contain stop symbol '*' with color flag
        if self.k == 3:
            # for aa triplets, '*' denotes a stop residue
            if '*' in kmer:
                if 'color' not in node:
                    node['color'] = 1

    def add_sequence(self, sequence, clone_id):
        """
        Add all k-mers from sequence (sliding window) to the network.
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

    def to_rows(self):
        """Return list of rows representing nodes suitable for CSV writing.

        Each row: { 'kmer', 'count', 'clone_count', 'clones', 'positions', 'color' }
        """
        rows = []
        for kmer, node in self.nodes.items():
            clones = sorted(node.get('clones', []))
            clone_count = len(clones)
            positions = None
            if self.store_positions:
                positions = {str(cid): pos_list[:] for cid, pos_list in node.get('positions', {}).items()}
            row = {
                'kmer': kmer,
                'count': node.get('count', 0),
                'clone_count': clone_count,
                'clones': ';'.join(str(c) for c in clones),
                'positions': json.dumps(positions) if positions is not None else '',
            }
            if 'color' in node:
                row['color'] = node['color']
            else:
                row['color'] = ''
            rows.append(row)
        return rows

    def to_aggregated_rows(self):
        """Return compact aggregated rows: kmer, count, clones (semi-colon list)."""
        rows = []
        for kmer, node in self.nodes.items():
            clones = sorted(node.get('clones', []))
            row = {
                'kmer': kmer,
                'count': node.get('count', 0),
                'clones': ';'.join(str(c) for c in clones)
            }
            rows.append(row)
        return rows

    def normalize_nodes(self):
        """Normalize node frequencies by dividing each node's count by the total sum of all counts.
        
        Stores the normalized value as 'normalized_freq' field on each node.
        This represents the proportion of the total clone population containing this k-mer.
        """
        total_count = sum(node.get('count', 0) for node in self.nodes.values())
        if total_count == 0:
            return
        
        for kmer, node in self.nodes.items():
            count = node.get('count', 0)
            node['normalized_freq'] = count / total_count
            # Initialize threshold flag as True (included by default)
            node['above_threshold'] = True

    def apply_node_threshold(self, threshold: float):
        """Apply frequency threshold to nodes.
        
        Nodes with normalized_freq below threshold will have 'above_threshold' = False.
        This is logical filtering - nodes remain in the structure but are marked as filtered.
        
        Args:
            threshold: minimum normalized frequency (0.0 to 1.0)
        """
        for kmer, node in self.nodes.items():
            normalized_freq = node.get('normalized_freq', 0.0)
            node['above_threshold'] = normalized_freq >= threshold

    def compute_edge_probabilities(self):
        """Compute edge probabilities for the network.
        
        For each node i, computes the probability that it originated from each neighbor j as:
        P(j -> i) = frequency(j) / sum of frequencies of all neighbors of i
        
        Stores edge probabilities in 'edges' dict where:
        edges[node_i][node_j] = probability that node_i came from node_j
        
        This creates a directed, weighted graph.
        """
        # Initialize edges structure if not present
        if not hasattr(self, 'edges') or self.edges is None:
            self.edges = {}
        
        # For k=3 (triplets), neighbors differ by exactly one position
        for node_i, node_data_i in self.nodes.items():
            if len(node_i) != self.k:
                continue
                
            # Find all neighbors (differ by exactly one character)
            neighbors = []
            for node_j, node_data_j in self.nodes.items():
                if node_i == node_j:
                    continue
                # Count differences
                diffs = sum(1 for a, b in zip(node_i, node_j) if a != b)
                if diffs == 1:
                    neighbors.append((node_j, node_data_j.get('count', 0)))
            
            # Compute probabilities
            if neighbors:
                total_neighbor_freq = sum(freq for _, freq in neighbors)
                if total_neighbor_freq > 0:
                    self.edges[node_i] = {}
                    for neighbor_id, neighbor_freq in neighbors:
                        probability = neighbor_freq / total_neighbor_freq
                        self.edges[node_i][neighbor_id] = {
                            'probability': probability,
                            'above_threshold': True  # Default to True, will be updated by threshold
                        }

    def apply_edge_threshold(self, threshold: float):
        """Apply probability threshold to edges.
        
        Edges with probability below threshold will have 'above_threshold' = False.
        This is logical filtering - edges remain but are marked as filtered.
        
        Args:
            threshold: minimum edge probability (0.0 to 1.0)
        """
        if not hasattr(self, 'edges') or self.edges is None:
            return
            
        for node_i, edges_dict in self.edges.items():
            for node_j, edge_data in edges_dict.items():
                probability = edge_data.get('probability', 0.0)
                edge_data['above_threshold'] = probability >= threshold

    def get_frequency_distribution(self):
        """Get distribution of node frequencies for histogram.
        
        Returns:
            dict: {'frequencies': list of counts, 'normalized_frequencies': list of normalized values}
        """
        frequencies = [node.get('count', 0) for node in self.nodes.values()]
        normalized = [node.get('normalized_freq', 0.0) for node in self.nodes.values()]
        return {
            'frequencies': frequencies,
            'normalized_frequencies': normalized
        }

    def get_edge_probability_distribution(self):
        """Get distribution of edge probabilities for histogram.
        
        Returns:
            list: list of all edge probabilities
        """
        if not hasattr(self, 'edges') or self.edges is None:
            return []
        
        probabilities = []
        for edges_dict in self.edges.values():
            for edge_data in edges_dict.values():
                probabilities.append(edge_data.get('probability', 0.0))
        return probabilities


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
        # create per-region nucleotide networks: r1 = indices 0..220, r2 = 221..end
        self.nucleotide_network_region = {
            k: {'r1': KmerNetwork(k, store_positions=store_positions), 'r2': KmerNetwork(k, store_positions=store_positions)}
            for k in ks_filtered
        }
        # per-region amino-acid triplet networks (aa3)
        self.aa_network_region = {'r1': KmerNetwork(3, store_positions=store_positions), 'r2': KmerNetwork(3, store_positions=store_positions)}

    def add_sequence(self, sequence, clone_id):
        # nucleotide networks
        # manually slide windows so we can split into regions
        for k, net in self.nucleotide_network.items():
            L = len(sequence)
            if L < k:
                continue
            for i in range(0, L - k + 1):
                kmer = sequence[i:i+k]
                if '-' in kmer or 'N' in kmer:
                    continue
                # add to global nucleotide network
                net.add_kmer(kmer, clone_id, pos=i if net.store_positions else None)
                # add to appropriate region network
                region = 'r1' if i <= 220 else 'r2'
                region_net = self.nucleotide_network_region[k][region]
                region_net.add_kmer(kmer, clone_id, pos=i if region_net.store_positions else None)

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
                # also add to aa region split
                region = 'r1' if i <= 220 else 'r2'
                self.aa_network_region[region].add_kmer(aa_triplet, clone_id, pos=i if self.aa_network_region[region].store_positions else None)

    def get_network(self, k):
        return self.nucleotide_network.get(k)

    def get_network_region(self, k, region: str):
        """Return the nucleotide network for a specific region ('r1' or 'r2')."""
        rn = self.nucleotide_network_region.get(k)
        if not rn:
            return None
        return rn.get(region)

    def get_aa_network(self):
        return self.aa_network

    def get_aa_network_region(self, region: str):
        """Return the amino-acid triplet network for a specific region ('r1' or 'r2')."""
        return self.aa_network_region.get(region)

    def to_serializable(self):
        out = {k: net.to_serializable() for k, net in self.nucleotide_network.items()}
        # include per-region nucleotide networks
        for k, regions in self.nucleotide_network_region.items():
            out[f"k{k}_r1"] = regions['r1'].to_serializable()
            out[f"k{k}_r2"] = regions['r2'].to_serializable()
        # include per-region aa networks
        out['aa3_r1'] = self.aa_network_region['r1'].to_serializable()
        out['aa3_r2'] = self.aa_network_region['r2'].to_serializable()
        if self.aa_network:
            out['aa3'] = self.aa_network.to_serializable()
        return out

    def dump_json(self, basepath):
        for k, net in self.nucleotide_network.items():
            path = f"{basepath}.k{k}.json"
            net.dump_json(path)
            # dump region splits
            r1_path = f"{basepath}.k{k}.r1.json"
            r2_path = f"{basepath}.k{k}.r2.json"
            self.nucleotide_network_region[k]['r1'].dump_json(r1_path)
            self.nucleotide_network_region[k]['r2'].dump_json(r2_path)
        # dump aa region splits
        self.aa_network_region['r1'].dump_json(f"{basepath}.aa3.r1.json")
        self.aa_network_region['r2'].dump_json(f"{basepath}.aa3.r2.json")
        if self.aa_network:
            self.aa_network.dump_json(f"{basepath}.aa3.json")

    def dump_region_csvs(self, basepath):
        """Dump per-region CSVs for nucleotide and amino-acid networks.

        Writes files:
          {basepath}.aa3.r1.csv, {basepath}.aa3.r2.csv,
          {basepath}.k{K}.r1.csv, {basepath}.k{K}.r2.csv for each nucleotide k
        """
        # write aa region CSVs
        for region in ('r1', 'r2'):
            aa_net = self.aa_network_region.get(region)
            if aa_net:
                path = f"{basepath}.aa3.{region}.csv"
                rows = aa_net.to_rows()
                if rows:
                    with open(path, 'w', newline='') as fh:
                        writer = csv.DictWriter(fh, fieldnames=['kmer', 'count', 'clone_count', 'clones', 'positions', 'color'])
                        writer.writeheader()
                        for r in rows:
                            writer.writerow(r)

        # write nucleotide region CSVs
        for k, regions in self.nucleotide_network_region.items():
            for region in ('r1', 'r2'):
                net = regions.get(region)
                if not net:
                    continue
                path = f"{basepath}.k{k}.{region}.csv"
                rows = net.to_rows()
                if rows:
                    with open(path, 'w', newline='') as fh:
                        writer = csv.DictWriter(fh, fieldnames=['kmer', 'count', 'clone_count', 'clones', 'positions', 'color'])
                        writer.writeheader()
                        for r in rows:
                            writer.writerow(r)

    def dump_region_aggregated_csvs(self, basepath):
        """Dump compact aggregated CSVs for each region (kmer, count, clones).

        Writes files:
          {basepath}.aa3.r1.aggregated.csv, {basepath}.aa3.r2.aggregated.csv,
          {basepath}.k{K}.r1.aggregated.csv, {basepath}.k{K}.r2.aggregated.csv
        """
        # aa aggregated
        for region in ('r1', 'r2'):
            aa_net = self.aa_network_region.get(region)
            if aa_net:
                path = f"{basepath}.aa3.{region}.aggregated.csv"
                rows = aa_net.to_aggregated_rows()
                if rows:
                    df = pd.DataFrame(rows)
                    df.to_csv(path, index=False)

        # nucleotide aggregated
        for k, regions in self.nucleotide_network_region.items():
            for region in ('r1', 'r2'):
                net = regions.get(region)
                if not net:
                    continue
                path = f"{basepath}.k{k}.{region}.aggregated.csv"
                rows = net.to_aggregated_rows()
                if rows:
                    df = pd.DataFrame(rows)
                    df.to_csv(path, index=False)

    def summary(self) -> dict:
        out = {}
        for k, net in self.nucleotide_network.items():
            out[k] = {'nodes': net.node_count(), 'occurrences': net.total_occurrences()}
            # include region summaries
            regions = self.nucleotide_network_region.get(k)
            if regions:
                out[f"{k}_r1"] = {'nodes': regions['r1'].node_count(), 'occurrences': regions['r1'].total_occurrences()}
                out[f"{k}_r2"] = {'nodes': regions['r2'].node_count(), 'occurrences': regions['r2'].total_occurrences()}
        # aa region summaries
        out['aa3_r1'] = {'nodes': self.aa_network_region['r1'].node_count(), 'occurrences': self.aa_network_region['r1'].total_occurrences()}
        out['aa3_r2'] = {'nodes': self.aa_network_region['r2'].node_count(), 'occurrences': self.aa_network_region['r2'].total_occurrences()}
        if self.aa_network:
            out['aa3'] = {'nodes': self.aa_network.node_count(), 'occurrences': self.aa_network.total_occurrences()}
        return out


def compute_network_distance(network1: KmerNetwork, network2: KmerNetwork) -> float:
    """
    Compute the distance (dissimilarity) between two KmerNetwork instances based on shared triplet clones.

    The distance R is the average over all unique kmers of:
    R_kmer = 1 - (ab / (A + B))
    where A is number of clones for kmer in network1,
    B in network2, ab is intersection.

    Returns the average R across all kmers present in at least one network.
    """
    all_kmers = set(network1.nodes.keys()) | set(network2.nodes.keys())
    if not all_kmers:
        return 0.0  # no kmers, distance 0

    total_r = 0.0
    count = 0
    for kmer in all_kmers:
        clones1 = network1.nodes.get(kmer, {}).get('clones', set())
        clones2 = network2.nodes.get(kmer, {}).get('clones', set())
        if clones1 == clones2:
            continue  # identical clone sets, skip
        A = len(clones1)
        B = len(clones2)
        ab = len(clones1 & clones2)
        if A + B > 0:
            r = 1 - (ab / (A + B))
        else:
            r = 0.0  # shouldn't happen
        total_r += r
        count += 1
    return total_r / count if count > 0 else 0.0


def compute_degree_distribution(network: KmerNetwork) -> dict:
    """Compute degree distribution for the network.
    
    Returns a dictionary with degree statistics:
    - in_degree: dict mapping node -> number of incoming edges
    - out_degree: dict mapping node -> number of outgoing edges
    - degree_dist: distribution of degrees
    """
    if not hasattr(network, 'edges') or network.edges is None:
        network.compute_edge_probabilities()
    
    in_degree = {}
    out_degree = {}
    
    # Initialize all nodes with 0 degree
    for node in network.nodes.keys():
        in_degree[node] = 0
        out_degree[node] = 0
    
    # Count edges (only those above threshold)
    for node_i, edges_dict in network.edges.items():
        for node_j, edge_data in edges_dict.items():
            if edge_data.get('above_threshold', True):
                out_degree[node_i] += 1
                in_degree[node_j] += 1
    
    # Compute distribution
    all_degrees = list(in_degree.values()) + list(out_degree.values())
    degree_counts = {}
    for deg in all_degrees:
        degree_counts[deg] = degree_counts.get(deg, 0) + 1
    
    return {
        'in_degree': in_degree,
        'out_degree': out_degree,
        'degree_distribution': degree_counts,
        'avg_in_degree': sum(in_degree.values()) / len(in_degree) if in_degree else 0,
        'avg_out_degree': sum(out_degree.values()) / len(out_degree) if out_degree else 0
    }


def compute_degree_centrality(network: KmerNetwork) -> dict:
    """Compute degree centrality for all nodes.
    
    Degree centrality is the fraction of nodes that a node is connected to.
    For directed graphs, we compute both in-degree and out-degree centrality.
    
    Returns:
        dict: {
            'in_centrality': dict mapping node -> in-degree centrality,
            'out_centrality': dict mapping node -> out-degree centrality,
            'top_central_nodes': list of (node, centrality) tuples sorted by total centrality
        }
    """
    degree_stats = compute_degree_distribution(network)
    n = len(network.nodes)
    
    if n <= 1:
        return {
            'in_centrality': {},
            'out_centrality': {},
            'top_central_nodes': []
        }
    
    # Normalize by (n-1) to get centrality scores
    in_centrality = {node: deg / (n - 1) for node, deg in degree_stats['in_degree'].items()}
    out_centrality = {node: deg / (n - 1) for node, deg in degree_stats['out_degree'].items()}
    
    # Compute total centrality (in + out)
    total_centrality = {node: in_centrality[node] + out_centrality[node] for node in network.nodes.keys()}
    
    # Get top central nodes
    top_nodes = sorted(total_centrality.items(), key=lambda x: x[1], reverse=True)[:20]
    
    return {
        'in_centrality': in_centrality,
        'out_centrality': out_centrality,
        'top_central_nodes': top_nodes
    }


def compare_networks(network1: KmerNetwork, network2: KmerNetwork) -> dict:
    """Compare two networks using multiple metrics.
    
    Returns a comprehensive comparison including:
    - Basic statistics (node count, edge count)
    - Degree distributions
    - Centrality measures
    - Structural distance
    
    Args:
        network1: First KmerNetwork
        network2: Second KmerNetwork
    
    Returns:
        dict: Comprehensive comparison metrics
    """
    # Ensure both networks have edge computations
    if not hasattr(network1, 'edges') or network1.edges is None:
        network1.compute_edge_probabilities()
    if not hasattr(network2, 'edges') or network2.edges is None:
        network2.compute_edge_probabilities()
    
    # Basic statistics
    stats1 = {
        'node_count': len(network1.nodes),
        'edge_count': sum(len(edges) for edges in network1.edges.values()),
        'total_frequency': sum(node.get('count', 0) for node in network1.nodes.values())
    }
    stats2 = {
        'node_count': len(network2.nodes),
        'edge_count': sum(len(edges) for edges in network2.edges.values()),
        'total_frequency': sum(node.get('count', 0) for node in network2.nodes.values())
    }
    
    # Degree distributions
    degree_dist1 = compute_degree_distribution(network1)
    degree_dist2 = compute_degree_distribution(network2)
    
    # Centrality measures
    centrality1 = compute_degree_centrality(network1)
    centrality2 = compute_degree_centrality(network2)
    
    # Structural distance
    distance = compute_network_distance(network1, network2)
    
    return {
        'network1': {
            'stats': stats1,
            'degree': degree_dist1,
            'centrality': centrality1
        },
        'network2': {
            'stats': stats2,
            'degree': degree_dist2,
            'centrality': centrality2
        },
        'distance': distance
    }

