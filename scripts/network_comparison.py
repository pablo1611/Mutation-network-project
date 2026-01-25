"""Network comparison and analysis visualization module.

This module provides functions to compare two networks and visualize
their structural differences using various metrics.
"""
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import os


def compute_degree_distribution(aa_network):
    """Compute degree distribution for a network.
    
    Args:
        aa_network: KmerNetwork instance with edges
        
    Returns:
        dict: {'in_degree': list, 'out_degree': list, 'total_degree': list}
    """
    if not hasattr(aa_network, 'edges') or not aa_network.edges:
        return {'in_degree': [], 'out_degree': [], 'total_degree': []}
    
    # Count in-degrees and out-degrees
    in_degree = {}
    out_degree = {}
    
    # Initialize all nodes with 0 degree
    for node in aa_network.nodes:
        in_degree[node] = 0
        out_degree[node] = 0
    
    # Count edges
    for source, targets in aa_network.edges.items():
        out_degree[source] = len([t for t, data in targets.items() if data.get('above_threshold', True)])
        for target, edge_data in targets.items():
            if edge_data.get('above_threshold', True):
                in_degree[target] = in_degree.get(target, 0) + 1
    
    return {
        'in_degree': list(in_degree.values()),
        'out_degree': list(out_degree.values()),
        'total_degree': [in_degree[n] + out_degree[n] for n in aa_network.nodes]
    }


def compute_degree_centrality(aa_network):
    """Compute degree centrality for all nodes.
    
    Degree centrality = degree / (n - 1) where n is number of nodes
    
    Args:
        aa_network: KmerNetwork instance with edges
        
    Returns:
        dict: node_id -> centrality score
    """
    n = len(aa_network.nodes)
    if n <= 1:
        return {node: 0.0 for node in aa_network.nodes}
    
    degree_dist = compute_degree_distribution(aa_network)
    total_degrees = degree_dist['total_degree']
    
    centrality = {}
    for i, node in enumerate(aa_network.nodes):
        centrality[node] = total_degrees[i] / (n - 1)
    
    return centrality


def compare_networks(net1, net2):
    """Compare two networks structurally.
    
    Returns a dictionary with comparison metrics:
    - Common nodes, unique nodes
    - Degree distributions
    - Network density
    - Average centrality
    
    Args:
        net1: First KmerNetwork
        net2: Second KmerNetwork
        
    Returns:
        dict: Comparison metrics
    """
    nodes1 = set(net1.nodes.keys())
    nodes2 = set(net2.nodes.keys())
    
    common_nodes = nodes1 & nodes2
    unique_to_1 = nodes1 - nodes2
    unique_to_2 = nodes2 - nodes1
    
    # Compute degree distributions
    deg_dist1 = compute_degree_distribution(net1)
    deg_dist2 = compute_degree_distribution(net2)
    
    # Compute centralities
    centrality1 = compute_degree_centrality(net1)
    centrality2 = compute_degree_centrality(net2)
    
    return {
        'common_nodes': len(common_nodes),
        'unique_to_1': len(unique_to_1),
        'unique_to_2': len(unique_to_2),
        'total_nodes_1': len(nodes1),
        'total_nodes_2': len(nodes2),
        'degree_dist_1': deg_dist1,
        'degree_dist_2': deg_dist2,
        'avg_centrality_1': np.mean(list(centrality1.values())) if centrality1 else 0,
        'avg_centrality_2': np.mean(list(centrality2.values())) if centrality2 else 0,
        'avg_in_degree_1': np.mean(deg_dist1['in_degree']) if deg_dist1['in_degree'] else 0,
        'avg_in_degree_2': np.mean(deg_dist2['in_degree']) if deg_dist2['in_degree'] else 0,
        'avg_out_degree_1': np.mean(deg_dist1['out_degree']) if deg_dist1['out_degree'] else 0,
        'avg_out_degree_2': np.mean(deg_dist2['out_degree']) if deg_dist2['out_degree'] else 0
    }


def plot_network_comparison(net1, net2, name1="Network 1", name2="Network 2", out_html=None, return_html: bool = False):
    """Create comprehensive comparison visualization.
    
    Args:
        net1: First KmerNetwork
        net2: Second KmerNetwork
        name1: Name for first network
        name2: Name for second network
        out_html: Optional path to save HTML
        return_html: If True and out_html is not provided, return the full HTML string
    """
    comparison = compare_networks(net1, net2)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Node Overlap',
            'Degree Distribution Comparison',
            'In-Degree vs Out-Degree',
            'Centrality Comparison'
        ),
        specs=[
            [{'type': 'pie'}, {'type': 'bar'}],
            [{'type': 'scatter'}, {'type': 'bar'}]
        ],
        vertical_spacing=0.12,
        horizontal_spacing=0.12
    )
    
    # 1. Node overlap pie chart
    labels = [f'Common', f'Unique to {name1}', f'Unique to {name2}']
    values = [
        comparison['common_nodes'],
        comparison['unique_to_1'],
        comparison['unique_to_2']
    ]
    colors = ['#4CAF50', '#2196F3', '#FF9800']
    
    fig.add_trace(
        go.Pie(labels=labels, values=values, marker=dict(colors=colors)),
        row=1, col=1
    )
    
    # 2. Degree distribution comparison
    deg1 = comparison['degree_dist_1']['total_degree']
    deg2 = comparison['degree_dist_2']['total_degree']
    
    if deg1 and deg2:
        max_degree = max(max(deg1), max(deg2))
        bins = list(range(0, max_degree + 2))
        
        hist1, _ = np.histogram(deg1, bins=bins)
        hist2, _ = np.histogram(deg2, bins=bins)
        
        fig.add_trace(
            go.Bar(x=bins[:-1], y=hist1, name=name1, marker_color='#2196F3'),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(x=bins[:-1], y=hist2, name=name2, marker_color='#FF9800'),
            row=1, col=2
        )
    
    # 3. In-degree vs Out-degree scatter
    in_deg1 = comparison['degree_dist_1']['in_degree']
    out_deg1 = comparison['degree_dist_1']['out_degree']
    in_deg2 = comparison['degree_dist_2']['in_degree']
    out_deg2 = comparison['degree_dist_2']['out_degree']
    
    if in_deg1 and out_deg1:
        fig.add_trace(
            go.Scatter(x=in_deg1, y=out_deg1, mode='markers', name=name1,
                      marker=dict(color='#2196F3', size=8, opacity=0.6)),
            row=2, col=1
        )
    if in_deg2 and out_deg2:
        fig.add_trace(
            go.Scatter(x=in_deg2, y=out_deg2, mode='markers', name=name2,
                      marker=dict(color='#FF9800', size=8, opacity=0.6)),
            row=2, col=1
        )
    
    # 4. Centrality comparison
    metrics = ['Avg Centrality', 'Avg In-Degree', 'Avg Out-Degree']
    values1 = [
        comparison['avg_centrality_1'],
        comparison['avg_in_degree_1'],
        comparison['avg_out_degree_1']
    ]
    values2 = [
        comparison['avg_centrality_2'],
        comparison['avg_in_degree_2'],
        comparison['avg_out_degree_2']
    ]
    
    fig.add_trace(
        go.Bar(x=metrics, y=values1, name=name1, marker_color='#2196F3'),
        row=2, col=2
    )
    fig.add_trace(
        go.Bar(x=metrics, y=values2, name=name2, marker_color='#FF9800'),
        row=2, col=2
    )
    
    # Update axes labels
    fig.update_xaxes(title_text="Degree", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    
    fig.update_xaxes(title_text="In-Degree", row=2, col=1)
    fig.update_yaxes(title_text="Out-Degree", row=2, col=1)
    
    fig.update_xaxes(title_text="Metric", row=2, col=2)
    fig.update_yaxes(title_text="Value", row=2, col=2)
    
    # Update layout
    fig.update_layout(
        title_text=f"Network Comparison: {name1} vs {name2}",
        showlegend=True,
        height=900,
        hovermode='closest'
    )
    
    if out_html:
        fig.write_html(out_html)
        print(f"Comparison visualization saved to: {out_html}")
        return comparison

    if return_html:
        html = fig.to_html(full_html=True, include_plotlyjs='cdn')
        return comparison, html

    fig.show()
    return comparison


def generate_comparison_report(
    net1,
    net2,
    name1="Network 1",
    name2="Network 2",
    output_dir=None,
    save_to_disk: bool = True,
):
    """Generate a comprehensive comparison report with visualizations and statistics.
    
    Args:
        net1: First KmerNetwork
        net2: Second KmerNetwork
        name1: Name for first network
        name2: Name for second network
        output_dir: Directory to save outputs
        save_to_disk: When False, returns an in-memory HTML report + text summary (no files written)
        
    Returns:
        If save_to_disk=True:
            str: Path to output directory
        If save_to_disk=False:
            dict: {'comparison': dict, 'html': str, 'text': str}
    """
    def _build_text_report(cmp):
        lines = []
        lines.append("Network Comparison Report")
        lines.append("=" * 60)
        lines.append("")
        lines.append(f"{name1} vs {name2}")
        lines.append("")

        lines.append("Node Statistics:")
        lines.append(f"  {name1}: {cmp['total_nodes_1']} nodes")
        lines.append(f"  {name2}: {cmp['total_nodes_2']} nodes")
        lines.append(f"  Common nodes: {cmp['common_nodes']}")
        lines.append(f"  Unique to {name1}: {cmp['unique_to_1']}")
        lines.append(f"  Unique to {name2}: {cmp['unique_to_2']}")
        lines.append("")

        lines.append("Degree Statistics:")
        lines.append(f"  {name1}:")
        lines.append(f"    Avg In-Degree: {cmp['avg_in_degree_1']:.4f}")
        lines.append(f"    Avg Out-Degree: {cmp['avg_out_degree_1']:.4f}")
        lines.append(f"  {name2}:")
        lines.append(f"    Avg In-Degree: {cmp['avg_in_degree_2']:.4f}")
        lines.append(f"    Avg Out-Degree: {cmp['avg_out_degree_2']:.4f}")
        lines.append("")

        lines.append("Centrality:")
        lines.append(f"  {name1}: {cmp['avg_centrality_1']:.6f}")
        lines.append(f"  {name2}: {cmp['avg_centrality_2']:.6f}")
        lines.append("")
        return "\n".join(lines)

    if not save_to_disk:
        comparison, html = plot_network_comparison(net1, net2, name1, name2, out_html=None, return_html=True)
        text = _build_text_report(comparison)
        return {"comparison": comparison, "html": html, "text": text}

    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "output", "comparisons")

    os.makedirs(output_dir, exist_ok=True)

    # Generate visualization
    viz_path = os.path.join(output_dir, f"comparison_{name1}_vs_{name2}.html")
    comparison = plot_network_comparison(net1, net2, name1, name2, out_html=viz_path)

    # Generate text report
    report_path = os.path.join(output_dir, f"comparison_{name1}_vs_{name2}_report.txt")
    with open(report_path, 'w') as f:
        f.write(_build_text_report(comparison))

    print(f"\nComparison report generated in: {output_dir}")
    return output_dir
