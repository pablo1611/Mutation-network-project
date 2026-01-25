"""Plot AA-triplet network (Plotly prototype)

Usage:
  python scripts/plot_aa3_plotly.py --csv PATH_TO_CLONES_CSV --target ABC

This script builds `NetworksManager` from the clones CSV (using
`src.load_clones.load_clones(..., build_networks=True)`), extracts the
amino-acid triplet network and renders a Plotly 3D view showing all
nodes as a faint point cloud and a highlighted "star" of neighbors for
the chosen target node. Neighbor arcs are colored by which position
(1/2/3) differs.

"""
import argparse
import math
from typing import List, Tuple

import numpy as np
import plotly.graph_objects as go
import json

from src.load_clones import load_clones


def compute_pos_change(a: str, b: str) -> int:
    """Return 0-based index of differing position if exactly one differs, else -1."""
    if not (isinstance(a, str) and isinstance(b, str) and len(a) == len(b) == 3):
        return -1
    diffs = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return diffs[0] if len(diffs) == 1 else -1


def quad_bezier(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, steps: int = 24) -> Tuple[List[float], List[float], List[float]]:
    t = np.linspace(0, 1, steps)
    pts = (1 - t)[:, None] ** 2 * p0 + 2 * (1 - t)[:, None] * t[:, None] * p1 + t[:, None] ** 2 * p2
    return pts[:, 0].tolist(), pts[:, 1].tolist(), pts[:, 2].tolist()


def build_coords_map(nodes: List[str]) -> Tuple[dict, List[str]]:
    # derive sorted alphabet from observed letters (should be 21 letters -> 0..20)
    letters = sorted({c for k in nodes for c in k})
    mapping = {ch: i for i, ch in enumerate(letters)}
    coords = {k: (mapping[k[0]], mapping[k[1]], mapping[k[2]]) for k in nodes}
    return coords, letters


def plot_aa_network(nodes: List[str], coords_map: dict, letters: List[str], target: str = None, out_html: str = None, counts_map: dict = None, aa_network=None, node_threshold: float = 0.0, edge_threshold: float = 0.0):
    """Plot the AA network with optional filtering by thresholds.
    
    Args:
        nodes: List of node IDs
        coords_map: Mapping of node IDs to 3D coordinates
        letters: Sorted alphabet of amino acids
        target: Optional initial target node
        out_html: Optional path to save HTML
        counts_map: Mapping of node IDs to metadata (count, color)
        aa_network: KmerNetwork instance (optional, for edge probability filtering)
        node_threshold: Minimum normalized frequency for nodes (0.0 to 1.0)
        edge_threshold: Minimum probability for edges (0.0 to 1.0)
    """
    # Filter nodes based on threshold
    if aa_network:
        filtered_nodes = [
            n for n in nodes
            if aa_network.nodes.get(n, {}).get('above_threshold', True)
        ]
    else:
        filtered_nodes = nodes
    
    all_x = [coords_map[n][0] for n in filtered_nodes]
    all_y = [coords_map[n][1] for n in filtered_nodes]
    all_z = [coords_map[n][2] for n in filtered_nodes]

    fig = go.Figure()

    # full node cloud (faint) with hover showing triplet and count
    hover_texts = []
    sizes = []
    node_colors = []
    for n in filtered_nodes:
        meta = counts_map.get(n, {}) if counts_map else {}
        cnt = meta.get('count', 0) if isinstance(meta, dict) else meta
        is_stop = bool(meta.get('color')) if isinstance(meta, dict) else False
        
        # Add normalized frequency to hover text if available
        if aa_network:
            node_data = aa_network.nodes.get(n, {})
            norm_freq = node_data.get('normalized_freq', 0.0)
            hover_texts.append(f"{n}<br>count: {cnt}<br>norm_freq: {norm_freq:.6f}")
        else:
            hover_texts.append(f"{n}<br>count: {cnt}")
        
        # scale marker size slightly by log(count)
        try:
            sizes.append(max(2, math.log1p(cnt) * 3))
        except Exception:
            sizes.append(2)
        # color stop-containing triplets red, others purple
        node_colors.append('red' if is_stop else 'rgba(111,45,168,0.8)')

    fig.add_trace(go.Scatter3d(
        x=all_x, y=all_y, z=all_z,
        mode='markers',
        marker=dict(size=sizes, color=node_colors, opacity=0.9),
        name='All nodes', showlegend=False,
        hoverinfo='text',
        text=hover_texts,
        # store raw node id so JS can reliably identify the clicked node (one per point)
        customdata=[[nid] for nid in filtered_nodes]
    ))

    # Do not force a server-side default target: let the user pick one on the page.
    # If a target is provided, validate it; otherwise operate without an initial highlight.
    if target is not None and target not in coords_map:
        raise ValueError(f"Target node '{target}' not found in network nodes")

    colors = ['red', 'green', 'blue']
    # If no initial target was provided, don't draw a server-side highlighted target.
    # We'll add an empty placeholder trace for the target marker so the client can populate it.
    if target is not None:
        # draw initial target marker and its static neighbor arcs (kept for backwards compatibility)
        tx, ty, tz = coords_map[target]
        target_meta = counts_map.get(target, {}) if counts_map else {}
        target_count = target_meta.get('count', 0) if isinstance(target_meta, dict) else target_meta
        
        # Add normalized frequency to target display if available
        if aa_network:
            target_node_data = aa_network.nodes.get(target, {})
            target_norm_freq = target_node_data.get('normalized_freq', 0.0)
            target_text = f"<b>{target}</b><br>count: {target_count}<br>norm_freq: {target_norm_freq:.6f}"
        else:
            target_text = f"<b>{target}</b><br>count: {target_count}"
        
        fig.add_trace(go.Scatter3d(
            x=[tx], y=[ty], z=[tz],
            mode='markers+text',
            marker=dict(size=10 + (math.log1p(target_count) * 2), color='black', symbol='diamond', line=dict(color='white', width=1)),
            text=[target_text], textposition='top center',
            name='Target', showlegend=False,
            hoverinfo='text'
        ))

        nodes_set = set(filtered_nodes)
        max_dim = max(len(letters) - 1, 1)
        offset_mag = max(3, max_dim * 0.2)

        # For each neighbor differing by one position, draw a curved arc with midpoint offset along axis
        # Use edge probabilities if available and respect edge threshold
        for dim in range(3):
            for alt in letters:
                if alt == target[dim]:
                    continue
                neigh = list(target)
                neigh[dim] = alt
                neigh = ''.join(neigh)
                if neigh not in nodes_set:
                    continue

                # Check if edge passes threshold (if aa_network is provided)
                if aa_network and hasattr(aa_network, 'edges') and aa_network.edges:
                    # Check both directions
                    # P(target | neigh) means "target came from neigh", stored in edges[target][neigh]
                    outgoing_data = aa_network.edges.get(target, {}).get(neigh, {})  # target came from neigh
                    # P(neigh | target) means "neigh came from target", stored in edges[neigh][target]  
                    incoming_data = aa_network.edges.get(neigh, {}).get(target, {})  # neigh came from target
                    
                    # Skip if neither direction passes threshold
                    if not (outgoing_data.get('above_threshold', False) or incoming_data.get('above_threshold', False)):
                        continue
                    
                    edge_label = f"Edge: {target} ↔ {neigh}<br>Position {dim+1} differs"
                else:
                    edge_label = f"Edge: {target} ↔ {neigh}<br>Position {dim+1} differs"

                sx, sy, sz = coords_map[target]
                tx2, ty2, tz2 = coords_map[neigh]
                p0 = np.array([sx, sy, sz], dtype=float)
                p2 = np.array([tx2, ty2, tz2], dtype=float)
                mid = (p0 + p2) / 2.0
                # offset control point along axis corresponding to dim
                axis = np.zeros(3)
                axis[dim] = 1.0
                p1 = mid + axis * offset_mag

                xs, ys, zs = quad_bezier(p0, p1, p2, steps=30)
                fig.add_trace(go.Scatter3d(
                    x=xs, y=ys, z=zs,
                    mode='lines',
                    line=dict(color=colors[dim], width=3),
                    hoverinfo='text',
                    text=[edge_label] * len(xs),
                    showlegend=False
                ))

    # Calculate total normalized frequency of filtered nodes
    total_norm_freq = 0.0
    total_edges = 0
    if aa_network:
        for n in filtered_nodes:
            node_data = aa_network.nodes.get(n, {})
            total_norm_freq += node_data.get('normalized_freq', 0.0)
        
        # Count edges passing threshold
        if hasattr(aa_network, 'edges') and aa_network.edges:
            for node_i, edges_dict in aa_network.edges.items():
                if node_i not in filtered_nodes:
                    continue
                for node_j, edge_data in edges_dict.items():
                    if node_j not in filtered_nodes:
                        continue
                    if edge_data.get('above_threshold', True):
                        total_edges += 1

    fig.update_layout(
        title=dict(
            text=f"AA-triplet network<br><sub>Displaying {len(filtered_nodes)} nodes, {total_edges} edges | Total norm_freq: {total_norm_freq:.4f} ({total_norm_freq*100:.2f}%)</sub>",
            x=0.5,
            xanchor='center'
        ),
        scene=dict(aspectmode='cube'),
        # give more top margin so overlaid controls don't cover the title
        margin=dict(l=0, r=0, b=0, t=80),
        # include simple Zoom Out / Zoom In buttons positioned to the right
        updatemenus=[dict(type='buttons', direction='left', buttons=[
            dict(label='Zoom Out', method='relayout', args=[{"scene.xaxis.range": [ -1, len(letters)],
                                                               "scene.yaxis.range": [ -1, len(letters)],
                                                               "scene.zaxis.range": [ -1, len(letters)],
                                                               "scene.camera.eye": {"x": 1.5, "y": 1.5, "z": 1.5}}]),
            dict(label='Zoom In', method='relayout', args=[{"scene.xaxis.range": [0, 1],
                                                           "scene.yaxis.range": [0, 1],
                                                           "scene.zaxis.range": [0, 1],
                                                           "scene.camera.eye": {"x": 0.6, "y": 0.6, "z": 0.6}}])
        ], pad={"r": 10, "t": 10}, showactive=True, x=0.98, xanchor='right', y=1.05, yanchor='top')]
    )
    # Add three empty placeholder traces for hovered-neighbor arcs (pos1,pos2,pos3).
    hover_trace_idxs = []
    for i, col in enumerate(colors):
        tr = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color=col, width=4), hoverinfo='text', text=[], name=f'neighbors-pos{i+1}')
        fig.add_trace(tr)
        hover_trace_idxs.append(len(fig.data) - 1)

    # add a placeholder trace for the dynamic target marker (client will populate it)
    target_trace = go.Scatter3d(x=[], y=[], z=[], mode='markers+text', marker=dict(size=10, color='black', symbol='diamond', line=dict(color='white', width=1)), text=[], textposition='top center', name='Target', showlegend=False, hoverinfo='text')
    fig.add_trace(target_trace)
    target_trace_idx = len(fig.data) - 1

    # If an output path is provided, write an interactive HTML file and print location.
    if out_html:
        # prepare embedded data for JS: coords_map, neighbors, letters
        payload = {
            'coords': coords_map,
            'neighbors': {},
            'letters': letters,
            'colors': colors,
            # only expose an initial target if explicitly provided by the caller
            'initial_target': target if target is not None else '',
            'hover_trace_idxs': hover_trace_idxs,
            'target_trace_idx': target_trace_idx,
            'counts': {k: (counts_map.get(k, {}).get('count', 0) if isinstance(counts_map.get(k, {}), dict) else counts_map.get(k, 0)) for k in filtered_nodes},
            'normalized_freqs': {},
            'edge_probabilities': {},
            # Only include nodes above threshold (keep JSON-serializable)
            'filtered_nodes': list(filtered_nodes)
        }
        
        # Add normalized frequencies if available
        if aa_network:
            for n in filtered_nodes:
                node_data = aa_network.nodes.get(n, {})
                payload['normalized_freqs'][n] = node_data.get('normalized_freq', 0.0)
            
            # Add edge probabilities if available
            if hasattr(aa_network, 'edges') and aa_network.edges:
                for node_i, edges_dict in aa_network.edges.items():
                    if node_i not in filtered_nodes:
                        continue
                    payload['edge_probabilities'][node_i] = {}
                    for node_j, edge_data in edges_dict.items():
                        if node_j not in filtered_nodes:
                            continue
                        if edge_data.get('above_threshold', True):
                            payload['edge_probabilities'][node_i][node_j] = edge_data.get('probability', 0.0)
        
        # build neighbors mapping - include if EITHER direction passes threshold
        nodes_set = set(filtered_nodes)
        for n in filtered_nodes:
            neighs = []
            for d in range(3):
                for alt in letters:
                    if alt == n[d]:
                        continue
                    m = list(n)
                    m[d] = alt
                    m = ''.join(m)
                    if m in nodes_set:
                        # Check if edge passes threshold in EITHER direction
                        include_edge = False
                        if aa_network and hasattr(aa_network, 'edges') and aa_network.edges:
                            # P(n | m): n came from m, stored in edges[n][m]
                            n_from_m_passes = aa_network.edges.get(n, {}).get(m, {}).get('above_threshold', False)
                            # P(m | n): m came from n, stored in edges[m][n]
                            m_from_n_passes = aa_network.edges.get(m, {}).get(n, {}).get('above_threshold', False)
                            # Include if either direction passes
                            include_edge = n_from_m_passes or m_from_n_passes
                        else:
                            include_edge = True
                        
                        if include_edge:
                            neighs.append({'id': m, 'pos': d})
            payload['neighbors'][n] = neighs

        html_str = fig.to_html(include_plotlyjs='cdn', full_html=True)
        # append JS to handle hover/click events and draw neighbor arcs dynamically
        extra_js = """
<script>
const payload = """ + json.dumps(payload) + """;
function quadBezierPoints(p0, p1, p2, steps) {
    const pts = [];
    for (let i=0;i<steps;i++){
        const t = i/(steps-1);
        const x = (1-t)*(1-t)*p0[0] + 2*(1-t)*t*p1[0] + t*t*p2[0];
        const y = (1-t)*(1-t)*p0[1] + 2*(1-t)*t*p1[1] + t*t*p2[1];
        const z = (1-t)*(1-t)*p0[2] + 2*(1-t)*t*p1[2] + t*t*p2[2];
        pts.push([x,y,z]);
    }
    return pts;
}

function buildLineSegments(nodeId) {
    const neighs = payload.neighbors[nodeId] || [];
    const coords = payload.coords;
    const letters = payload.letters;
    const max_dim = Math.max(letters.length-1, 1);
    const offset_mag = Math.max(3, max_dim * 0.2);
    const steps = 30;
    const perDim = [[],[],[]];
    
    // Build edges for each neighbor
    for (const nb of neighs){
        const dim = nb.pos;
        const p0 = coords[nodeId];
        const p2 = coords[nb.id];
        const mid = [(p0[0]+p2[0])/2, (p0[1]+p2[1])/2, (p0[2]+p2[2])/2];
        const axis = [0,0,0]; axis[dim]=1;
        const p1 = [mid[0]+axis[0]*offset_mag, mid[1]+axis[1]*offset_mag, mid[2]+axis[2]*offset_mag];
        const pts = quadBezierPoints(p0,p1,p2,steps);
        const xs = [], ys = [], zs = [];
        for (let i=0;i<pts.length;i++){ xs.push(pts[i][0]); ys.push(pts[i][1]); zs.push(pts[i][2]); }
        xs.push(null); ys.push(null); zs.push(null);
        
        // Keep edge hover text minimal (no numeric probability/percent details)
        let hoverText = `Edge: ${nodeId} ↔ ${nb.id}<br>Position ${dim+1} differs`;
        
        perDim[dim].push({x: xs, y: ys, z: zs, text: hoverText});
    }
    const out = [ {x:[], y:[], z:[], text:[]}, {x:[], y:[], z:[], text:[]}, {x:[], y:[], z:[], text:[]} ];
    for (let d=0; d<3; d++){
        for (const seg of perDim[d]){ 
            out[d].x = out[d].x.concat(seg.x); 
            out[d].y = out[d].y.concat(seg.y); 
            out[d].z = out[d].z.concat(seg.z);
            // Add hover text for each point in the segment
            for (let i=0; i<seg.x.length; i++) {
                out[d].text.push(seg.text);
            }
        }
    }
    return out;
}

document.addEventListener('DOMContentLoaded', function() {
    const gd = document.querySelectorAll('.plotly-graph-div')[0];

    function extractNodeId(point) {
        // Try customdata first (most reliable for clicks)
        if (point && point.customdata) {
            try {
                const cd = point.customdata;
                if (Array.isArray(cd) && cd.length > 0) {
                    return String(cd[0]).trim().toUpperCase();
                }
                return String(cd).trim().toUpperCase();
            } catch (e) {
                console.warn('Error extracting customdata:', e);
            }
        }
        // Fallback: parse from hover text
        if (point && point.text) {
            try {
                const parts = String(point.text).split('<br>');
                if (parts.length > 0) {
                    return parts[0].replace(/<[^>]*>/g, '').trim().toUpperCase();
                }
            } catch (e) {
                console.warn('Error parsing text:', e);
            }
        }
        return null;
    }

    // Debug status panel - shows click events and locked node
    const debugPanel = document.createElement('div');
    debugPanel.id = 'debug-panel';
    debugPanel.style.position = 'fixed';
    debugPanel.style.top = '12px';
    debugPanel.style.left = '12px';
    debugPanel.style.width = '280px';
    debugPanel.style.zIndex = '2147483647';
    debugPanel.style.background = 'rgba(50,50,50,0.9)';
    debugPanel.style.color = 'white';
    debugPanel.style.borderRadius = '6px';
    debugPanel.style.padding = '8px 10px';
    debugPanel.style.fontFamily = 'monospace';
    debugPanel.style.fontSize = '11px';
    debugPanel.innerHTML = '<div style="font-weight:600; margin-bottom:4px;">DEBUG STATUS</div>' +
                          '<div id="debug-status">Waiting for interaction...</div>';
    try { document.body.appendChild(debugPanel); } catch(e) {}

    function updateDebugStatus(msg) {
        const el = document.getElementById('debug-status');
        if (el) {
            const timestamp = new Date().toLocaleTimeString();
            el.innerHTML = `[${timestamp}] ${msg}`;
        }
    }

    // Side panel: show outgoing edges + probabilities for clicked node
    const edgePanel = document.createElement('div');
    edgePanel.id = 'edge-panel';
    edgePanel.style.position = 'fixed';
    edgePanel.style.top = '64px';
    edgePanel.style.right = '12px';
    edgePanel.style.width = '340px';
    edgePanel.style.maxHeight = '80vh';
    edgePanel.style.overflow = 'auto';
    edgePanel.style.zIndex = '2147483647';
    edgePanel.style.background = 'rgba(255,255,255,0.95)';
    edgePanel.style.borderRadius = '6px';
    edgePanel.style.boxShadow = '0 2px 8px rgba(0,0,0,0.12)';
    edgePanel.style.padding = '8px 10px';
    edgePanel.style.fontFamily = 'system-ui, -apple-system, Segoe UI, Roboto, sans-serif';
    edgePanel.style.fontSize = '12px';
    edgePanel.innerHTML = '<div style="font-weight:600; margin-bottom:6px;">Outgoing edges</div>' +
                          '<div id="edge-panel-body" style="white-space:pre-wrap;">Click a node to list outgoing edges.</div>';
    try { document.body.appendChild(edgePanel); } catch(e) {}

    function renderOutgoingEdges(nodeId) {
        const body = document.getElementById('edge-panel-body');
        if (!body) return;
        const outgoing = (payload.edge_probabilities && payload.edge_probabilities[nodeId]) ? payload.edge_probabilities[nodeId] : null;
        if (!outgoing || Object.keys(outgoing).length === 0) {
            body.textContent = nodeId ? (nodeId + "\n(no outgoing edges above threshold)") : 'Click a node to list outgoing edges.';
            return;
        }
        const rows = Object.entries(outgoing)
            .sort((a, b) => (b[1] || 0) - (a[1] || 0))
            .map(([dst, prob]) => `${nodeId} → ${dst}: ${Number(prob).toFixed(6)}`);
        body.textContent = nodeId + "\n" + rows.join("\n");
    }
    // add input + button to allow dynamic zoom to typed node id
    const controls = document.createElement('div');
    controls.style.padding = '6px';
    controls.style.display = 'flex';
    controls.style.gap = '6px';
    controls.style.alignItems = 'center';
    const input = document.createElement('input');
    input.placeholder = 'Enter node id (e.g. ABC)';
    input.style.padding = '4px';
    input.style.width = '140px';
    // attach a datalist so users can pick an existing node from the site
    const dataList = document.createElement('datalist');
    dataList.id = 'node-list';
    for (const id of Object.keys(payload.coords)){
        const opt = document.createElement('option'); opt.value = id; dataList.appendChild(opt);
    }
    // link input to datalist and set initial value if provided
    input.setAttribute('list', 'node-list');
    if (payload.initial_target) input.value = payload.initial_target;
    const btn = document.createElement('button');
    btn.textContent = 'Zoom & Highlight';
    // hide the small auxiliary button to avoid confusion; use the explicit Zoom In button below
    btn.style.display = 'none';
    btn.style.padding = '4px 8px';
    controls.appendChild(input);
    controls.appendChild(dataList);
    controls.appendChild(btn);
    if (gd && gd.parentNode) {
        try { document.body.appendChild(controls); } catch(e) { gd.parentNode.insertBefore(controls, gd); }
    }
    let lockedNode = null;
    let hoverNode = null;
    // preserve original title to restore when clearing selection
    const originalTitle = (gd && gd.layout && gd.layout.title && (gd.layout.title.text || gd.layout.title)) || document.title || '';
    // style controls to avoid overlapping plot title and add a separate title span
    try {
        // Use fixed positioning so controls remain visible regardless of plot container
        controls.style.position = 'fixed';
        controls.style.top = '12px';
        controls.style.right = '12px';
        controls.style.zIndex = '2147483647';
        controls.style.background = 'rgba(255,255,255,0.95)';
        controls.style.borderRadius = '6px';
        controls.style.boxShadow = '0 2px 8px rgba(0,0,0,0.12)';
        controls.style.padding = '6px';
    } catch(e) {}
    const titleSpan = document.createElement('div');
    titleSpan.style.fontSize = '14px';
    titleSpan.style.fontWeight = '600';
    titleSpan.style.marginRight = '12px';
    titleSpan.textContent = originalTitle || 'AA-triplet network';
    // show initial target in the small title span if provided
    try { if (payload.initial_target) titleSpan.textContent = `AA-triplet network: target ${payload.initial_target}`; } catch(e) {}
    try { controls.insertBefore(titleSpan, controls.firstChild); } catch(e) {}
    // keep built-in Plotly updatemenus (we position them to the right in layout)
    // add explicit Zoom Out / Zoom In buttons to controls (placed to the right)
    const zoomOutBtn = document.createElement('button');
    zoomOutBtn.textContent = 'Zoom Out';
    zoomOutBtn.style.padding = '4px 8px';
    zoomOutBtn.style.marginLeft = '6px';
    const zoomInBtn = document.createElement('button');
    zoomInBtn.textContent = 'Zoom In (input)';
    zoomInBtn.style.padding = '4px 8px';
    zoomInBtn.style.marginLeft = '6px';
    // hide auxiliary highlight button
    btn.style.display = 'none';
    // add a visible Search button that will zoom to the typed triplet
    const searchBtn = document.createElement('button');
    searchBtn.textContent = 'Search';
    searchBtn.style.padding = '4px 8px';
    searchBtn.style.marginLeft = '6px';
    controls.appendChild(searchBtn);
    controls.appendChild(zoomOutBtn);
    controls.appendChild(zoomInBtn);
    let zoomInTarget = null;
    let builtZoomInBtn = null;
    let builtZoomOutBtn = null;
    function zoomToNode(nodeId) {
        if(!payload.coords[nodeId]) { alert('Node "' + nodeId + '" not found'); return; }
        const c = payload.coords[nodeId];
        const range = 3;
        const xmin = c[0]-range, xmax = c[0]+range;
        const ymin = c[1]-range, ymax = c[1]+range;
        const zmin = c[2]-range, zmax = c[2]+range;
        Plotly.relayout(gd, {"scene.xaxis.range": [xmin, xmax], "scene.yaxis.range": [ymin, ymax], "scene.zaxis.range": [zmin, zmax], "scene.camera.eye": {"x":0.6,"y":0.6,"z":0.6}});
        // highlight neighbours and lock
        if(!payload.neighbors[nodeId]) return;
        lockedNode = nodeId;
        const segs = buildLineSegments(nodeId);
        for(let d=0; d<3; d++){ const idx = payload.hover_trace_idxs[d]; Plotly.restyle(gd, {x: [segs[d].x], y: [segs[d].y], z: [segs[d].z], text: [segs[d].text]}, [idx]); }
        // update plot title to reflect current target
        try { Plotly.relayout(gd, {"title.text": `AA-triplet network: target ${nodeId}`}); } catch(e) { console.warn('Could not update title', e); }
        // update dynamic title span (controls label)
        try { if (titleSpan) titleSpan.textContent = `AA-triplet network: target ${nodeId}`; } catch(e) {}
        // update Zoom In button label and remember target
        try { zoomInTarget = nodeId; zoomInBtn.textContent = `Zoom In (${nodeId})`; } catch(e) {}
        try { if (builtZoomInBtn) builtZoomInBtn.textContent = `Zoom In (${nodeId})`; } catch(e) {}
        // update the placeholder target trace so the selected node is highlighted
        try {
            const cnt = payload.counts && payload.counts[nodeId] ? payload.counts[nodeId] : 0;
            const normFreq = payload.normalized_freqs && payload.normalized_freqs[nodeId] ? payload.normalized_freqs[nodeId] : 0;
            const sz = 10 + (Math.log(1 + cnt) * 2);
            let displayText = `<b>${nodeId}</b><br>count: ${cnt}`;
            if (normFreq > 0) displayText += `<br>norm_freq: ${normFreq.toFixed(6)}`;
            Plotly.restyle(gd, {x: [[c[0]]], y: [[c[1]]], z: [[c[2]]], text: [[displayText]], 'marker.size': [[sz]]}, [payload.target_trace_idx]);
        } catch(e) { console.warn('Could not update target trace', e); }
    }
    btn.addEventListener('click', function() { zoomToNode(input.value.trim().toUpperCase()); });
    // allow pressing Enter in the input to trigger the same zoom/highlight
    input.addEventListener('keydown', function(e) { if (e.key === 'Enter') { e.preventDefault();
            // perform search on Enter
            const val = input.value.trim().toUpperCase();
            if(!val){ alert('Please enter a triplet to search'); return; }
            if(!payload.coords[val]){ alert('Triplet "' + val + '" not found in network'); return; }
            zoomToNode(val);
        } });

    // Search button click: validate and zoom to typed triplet
    searchBtn.addEventListener('click', function(){
        const val = (input.value||'').trim().toUpperCase();
        if(!val){ alert('Please enter a triplet to search'); return; }
        if(!payload.coords[val]){ alert('Triplet "' + val + '" not found in network'); return; }
        zoomToNode(val);
    });
    // Zoom Out button: restore full ranges
    zoomOutBtn.addEventListener('click', function(){
        const full = Math.max((payload.letters?payload.letters.length:0)-1, 1);
        try { Plotly.relayout(gd, {"scene.xaxis.range": [-1, full], "scene.yaxis.range": [-1, full], "scene.zaxis.range": [-1, full], "scene.camera.eye": {"x":1.5, "y":1.5, "z":1.5}}); } catch(e) {}
        zoomInTarget = null; zoomInBtn.textContent = 'Zoom In';
        try { if (titleSpan) titleSpan.textContent = originalTitle; } catch(e) {}
    });
    // Zoom In button: prioritize typed input then fallback to last target or initial
    zoomInBtn.addEventListener('click', function(){
        consif (!nodeId) {
                console.log('Hover: no node ID extracted');
                return;
            }
            console.log('Hover on node:', nodeId, 'Locked node:', lockedNode);
            hoverNode = nodeId;
            // If a different node is locked, don't show hover edges
            if(lockedNode && lockedNode !== nodeId) {
                updateDebugStatus(`Hover: ${nodeId} (ignored, ${lockedNode} locked)`);
                return;
            }
            // If this node has no neighbors, skip
            if(!payload.neighbors[nodeId]) {
                console.log('Node has no neighbors:', nodeId);
                return;
            }
            // Draw edges for this node'); return; }
        zoomToNode(target);
    });
    // remove Plotly's built-in updatemenus so our on-page controls are the authoritative Zoom buttons
    try {
        Plotly.relayout(gd, {updatemenus: []});
    } catch(e) {
        console.warn('Could not remove built-in updatemenus', e);
    }
    gd.on('plotly_hover', function(evt){
        try{
            const point = evt.points[0];
            const nodeId = extractNodeId(point);
            hoverNode = nodeId;
        console.log('Click event received:', evt);
        updateDebugStatus('CLICK event received!');
        if(evt && evt.points && evt.points.length>0){
            const point = evt.points[0];
            console.log('Click point data:', point);
            const nodeId = extractNodeId(point);
            console.log('Extracted node ID from click:', nodeId);
            updateDebugStatus(`Click: ${nodeId || 'NO ID'}`);
            
            if(!nodeId) {
                console.warn('Click: Could not extract node ID');
                updateDebugStatus('Click: Failed to extract node ID');
                return;
            }
            
            if(!payload.neighbors[nodeId]) {
                console.log('Click: Node has no neighbors:', nodeId);
                updateDebugStatus(`Click: ${nodeId} has no neighbors`);
                return;
            }
            
            // Toggle: if clicking the same locked node, unlock it
            if(lockedNode === nodeId){
                console.log('Unlocking node:', nodeId);
                updateDebugStatus(`UNLOCKED: ${nodeId}`);
                lockedNode = null;
                // If still hovering over a node, show its edges; otherwise clear
                if(hoverNode && payload.neighbors[hoverNode]){
                    const segs = buildLineSegments(hoverNode);
                    for(let d=0; d<3; d++){
                        const idx = payload.hover_trace_idxs[d];
                        Plotly.restyle(gd, {x: [segs[d].x], y: [segs[d].y], z: [segs[d].z], text: [segs[d].text]}, [idx]);
                    }
                } else {
                    for(let d=0; d<3; d++){
                        const idx = payload.hover_trace_idxs[d];
                        Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]]}, [idx]);
                    }
                }
                // restore original title
                try { Plotly.relayout(gd, {"title.text": originalTitle}); } catch(e) {}
                try { if (titleSpan) titleSpan.textContent = originalTitle; } catch(e) {}
                // clear dynamic target marker
                try { Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]], text: [[]]}, [payload.target_trace_idx]); } catch(e) {}
                // clear side panel
                try { renderOutgoingEdges(null); } catch(e) {}
            } else {
                // Lock this node: draw and persist its edges
                updateDebugStatus(`🔒 LOCKED: ${nodeId}`);
                console.log('Locking node:', nodeId);
                lockedNode = nodeId;
                const segs = buildLineSegments(nodeId);
                for(let d=0; d<3; d++){
                    const idx = payload.hover_trace_idxs[d];
                    Plotly.restyle(gd, {x: [segs[d].x], y: [segs[d].y], z: [segs[d].z], text: [segs[d].text]}, [idx]);
                }
                // update title to show locked node
                try { Plotly.relayout(gd, {"title.text": `AA-triplet network: locked on ${nodeId}`}); } catch(e) {}
                try { if (titleSpan) titleSpan.textContent = `AA-triplet network: locked on ${nodeId}`; } catch(e) {}
                // highlight the locked node with a marker
                try {
                    const c = payload.coords[nodeId];
                    const cnt = payload.counts && payload.counts[nodeId] ? payload.counts[nodeId] : 0;
                    const sz = 10 + (Math.log(1 + cnt) * 2);
                    Plotly.restyle(gd, {x: [[c[0]]], y: [[c[1]]], z: [[c[2]]], text: [[`<b>${nodeId}</b> (LOCKED)<br>count: ${cnt}`]], 'marker.size': [[sz]]}, [payload.target_trace_idx]);
                } catch(e) { console.warn('Could not update target trace on click', e); }
                // update side panel with outgoing edges
                try { renderOutgoingEdges(nodeId); } catch(e) {}
            }
        } else {
            // Click on empty space: clear lock
            updateDebugStatus('Click on empty space - cleared');
            console.log('Click on empty space - clearing lock');
            lockedNode = null;
            for(let d=0; d<3; d++){
                const idx = payload.hover_trace_idxs[d];
                Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]]}, [idx]);
            }
            try { Plotly.relayout(gd, {"title.text": originalTitle}); } catch(e) {}
            try { if (titleSpan) titleSpan.textContent = originalTitle; } catch(e) {}
        console.log('Unhover event - locked node:', lockedNode);
        hoverNode = null;
        // If a node is locked, keep showing its edges (don't clear on unhover)
        if(lockedNode) {
            console.log('Node is locked, keeping edges visible');
            return;
        }
        // Otherwise, clear the edge traces
        console.log('Clearing edges on unhover');
            for(let d=0; d<3; d++){
                const idx = payload.hover_trace_idxs[d];
                Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]]}, [idx]);
            }
            try { Plotly.relayout(gd, {"title.text": originalTitle}); } catch(e) {}
            try { if (titleSpan) titleSpan.textContent = originalTitle; } catch(e) {}
            // clear dynamic target trace
            try { Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]], text: [[]]}, [payload.target_trace_idx]); } catch(e) {}

            // clear side panel when clicking empty space
            try { renderOutgoingEdges(null); } catch(e) {}
        }
    });

    gd.on('plotly_unhover', function(evt){
        hoverNode = null;
        if(lockedNode) return;
        for(let d=0; d<3; d++){
            const idx = payload.hover_trace_idxs[d];
            Plotly.restyle(gd, {x: [[]], y: [[]], z: [[]]} , [idx]);
        }
    });
});
</script>
"""
        # Ensure the injected JS ends up inside the HTML document (before </body>)
        if '</body>' in html_str:
            html_str = html_str.replace('</body>', extra_js + '\n</body>')
        else:
            html_str = html_str + extra_js

        with open(out_html, 'w', encoding='utf-8') as fh:
            fh.write(html_str)
        print(f"Wrote interactive plot to: {out_html}")
    else:
        fig.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', required=True, help='Path to clones CSV (same shape used by load_clones)')
    parser.add_argument('--target', required=False, help='AA triplet to focus (e.g. "ARN")')
    parser.add_argument('--out', required=False, help='Write interactive HTML to this path instead of calling fig.show()')
    args = parser.parse_args()

    print('Loading clones and building networks (this may take a while)')
    result = load_clones(args.csv, build_networks=True)
    if isinstance(result, tuple):
        clones, networks = result
    else:
        raise RuntimeError('load_clones did not return networks; run with build_networks=True')

    aa_net = networks.get_aa_network()
    if aa_net is None:
        raise RuntimeError('No aa_network present in NetworksManager')

    nodes = sorted(aa_net.nodes.keys())
    coords_map, letters = build_coords_map(nodes)
    # build counts + color map from network node metadata
    counts_map = {k: {'count': aa_net.nodes[k].get('count', 0), 'color': aa_net.nodes[k].get('color', 0)} for k in nodes}

    # do not force a default target — let the page input control selection
    target = args.target if args.target else None
    plot_aa_network(nodes, coords_map, letters, target=target, out_html=args.out, counts_map=counts_map)


if __name__ == '__main__':
    main()


def plot_from_aa_network(aa_net, target: str = None, out_html: str = None, node_threshold: float = 0.0, edge_threshold: float = 0.0):
    """Convenience wrapper to plot directly from a NetworksManager.aa_network

    Usage from Python:
        from scripts.plot_aa3_plotly import plot_from_aa_network
        plot_from_aa_network(networks.get_aa_network(), target='ABC', node_threshold=0.001, edge_threshold=0.05)
    
    Args:
        aa_net: KmerNetwork instance
        target: Optional initial target node
        out_html: Optional path to save HTML
        node_threshold: Minimum normalized frequency for nodes (0.0 to 1.0)
        edge_threshold: Minimum probability for edges (0.0 to 1.0)
    """
    if aa_net is None:
        raise RuntimeError('aa_net is None')
    
    # Apply normalization and thresholds before plotting
    aa_net.normalize_nodes()
    aa_net.apply_node_threshold(node_threshold)
    aa_net.compute_edge_probabilities()
    aa_net.apply_edge_threshold(edge_threshold)
    
    nodes = sorted(aa_net.nodes.keys())
    coords_map, letters = build_coords_map(nodes)
    tgt = target
    counts_map = {k: {'count': aa_net.nodes[k].get('count', 0), 'color': aa_net.nodes[k].get('color', 0)} for k in nodes}
    plot_aa_network(nodes, coords_map, letters, target=tgt, out_html=out_html, counts_map=counts_map, 
                   aa_network=aa_net, node_threshold=node_threshold, edge_threshold=edge_threshold)
