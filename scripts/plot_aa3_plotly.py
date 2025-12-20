"""Plot AA-triplet network (Plotly prototype)

Usage:
  python scripts/plot_aa3_plotly.py --csv PATH_TO_CLONES_CSV --target ABC

This script builds `NetworksManager` from the clones CSV (using
`src.load_clones.load_clones(..., build_networks=True)`), extracts the
amino-acid triplet network and renders a Plotly 3D view showing all
nodes as a faint point cloud and a highlighted "star" of neighbors for
the chosen target node. Neighbor arcs are colored by which position
(1/2/3) differs.

Note: rendering all pairwise edges for 9,261 nodes is very heavy; this
prototype intentionally renders only the node cloud and a focused star
of neighbors for interactivity and quick validation.
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


def plot_aa_network(nodes: List[str], coords_map: dict, letters: List[str], target: str = None, out_html: str = None, counts_map: dict = None):
    all_x = [coords_map[n][0] for n in nodes]
    all_y = [coords_map[n][1] for n in nodes]
    all_z = [coords_map[n][2] for n in nodes]

    fig = go.Figure()

    # full node cloud (faint) with hover showing triplet and count
    hover_texts = []
    sizes = []
    node_colors = []
    for n in nodes:
        meta = counts_map.get(n, {}) if counts_map else {}
        cnt = meta.get('count', 0) if isinstance(meta, dict) else meta
        is_stop = bool(meta.get('color')) if isinstance(meta, dict) else False
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
        text=hover_texts
    ))

    if not target:
        target = nodes[0]

    if target not in coords_map:
        raise ValueError(f"Target node '{target}' not found in network nodes")

    colors = ['red', 'green', 'blue']
    # draw target marker
    tx, ty, tz = coords_map[target]
    target_meta = counts_map.get(target, {}) if counts_map else {}
    target_count = target_meta.get('count', 0) if isinstance(target_meta, dict) else target_meta
    fig.add_trace(go.Scatter3d(
        x=[tx], y=[ty], z=[tz],
        mode='markers+text',
        marker=dict(size=10 + (math.log1p(target_count) * 2), color='black', symbol='diamond', line=dict(color='white', width=1)),
        text=[f"<b>{target}</b><br>count: {target_count}"], textposition='top center',
        name='Target', showlegend=False,
        hoverinfo='text'
    ))

    nodes_set = set(nodes)
    max_dim = max(len(letters) - 1, 1)
    offset_mag = max(3, max_dim * 0.2)

    # For each neighbor differing by one position, draw a curved arc with midpoint offset along axis
    for dim in range(3):
        for alt in letters:
            if alt == target[dim]:
                continue
            neigh = list(target)
            neigh[dim] = alt
            neigh = ''.join(neigh)
            if neigh not in nodes_set:
                continue

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
                text=[f"{target} → {neigh} (pos {dim+1})"] * len(xs),
                showlegend=False
            ))

    fig.update_layout(
        title=f"AA-triplet network: target {target}",
        scene=dict(aspectmode='cube'),
        margin=dict(l=0, r=0, b=0, t=50),
        updatemenus=[dict(type='buttons', direction='left', buttons=[
            dict(label='Zoom Out', method='relayout', args=[{"scene.xaxis.range": [ -1, len(letters)],
                                                               "scene.yaxis.range": [ -1, len(letters)],
                                                               "scene.zaxis.range": [ -1, len(letters)],
                                                               "scene.camera.eye": {"x": 1.5, "y": 1.5, "z": 1.5}}]),
            dict(label=f'Zoom In ({target})', method='relayout', args=[{"scene.xaxis.range": [tx-3, tx+3],
                                                                         "scene.yaxis.range": [ty-3, ty+3],
                                                                         "scene.zaxis.range": [tz-3, tz+3],
                                                                         "scene.camera.eye": {"x": 0.6, "y": 0.6, "z": 0.6}}])
        ], pad={"r": 10, "t": 10}, showactive=True, x=0.1, xanchor='left', y=1.05, yanchor='top')]
    )
    # Add three empty placeholder traces for hovered-neighbor arcs (pos1,pos2,pos3).
    hover_trace_idxs = []
    for i, col in enumerate(colors):
        tr = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color=col, width=4), hoverinfo='none', name=f'neighbors-pos{i+1}')
        fig.add_trace(tr)
        hover_trace_idxs.append(len(fig.data) - 1)

    # If an output path is provided, write an interactive HTML file and print location.
    if out_html:
        # prepare embedded data for JS: coords_map, neighbors, letters
        payload = {
            'coords': coords_map,
            'neighbors': {},
            'letters': letters,
            'colors': colors,
            'hover_trace_idxs': hover_trace_idxs
        }
        # build neighbors mapping now
        nodes_set = set(nodes)
        for n in nodes:
            neighs = []
            for d in range(3):
                for alt in letters:
                    if alt == n[d]:
                        continue
                    m = list(n)
                    m[d] = alt
                    m = ''.join(m)
                    if m in nodes_set:
                        neighs.append({'id': m, 'pos': d})
            payload['neighbors'][n] = neighs

        html_str = fig.to_html(include_plotlyjs='cdn', full_html=True)
        # append JS to handle hover/click events and draw neighbor arcs dynamically
        extra_js = f"""
<script>
const payload = {json.dumps(payload)};
function quadBezierPoints(p0, p1, p2, steps) {{
    const pts = [];
    for (let i=0;i<steps;i++){{
        const t = i/(steps-1);
        const x = (1-t)*(1-t)*p0[0] + 2*(1-t)*t*p1[0] + t*t*p2[0];
        const y = (1-t)*(1-t)*p0[1] + 2*(1-t)*t*p1[1] + t*t*p2[1];
        const z = (1-t)*(1-t)*p0[2] + 2*(1-t)*t*p1[2] + t*t*p2[2];
        pts.push([x,y,z]);
    }}
    return pts;
}}

function buildLineSegments(nodeId) {{
    const neighs = payload.neighbors[nodeId] || [];
    const coords = payload.coords;
    const letters = payload.letters;
    const max_dim = Math.max(letters.length-1, 1);
    const offset_mag = Math.max(3, max_dim * 0.2);
    const steps = 30;
    const perDim = [[],[],[]];
    for (const nb of neighs){{
        const dim = nb.pos;
        const p0 = coords[nodeId];
        const p2 = coords[nb.id];
        const mid = [(p0[0]+p2[0])/2, (p0[1]+p2[1])/2, (p0[2]+p2[2])/2];
        const axis = [0,0,0]; axis[dim]=1;
        const p1 = [mid[0]+axis[0]*offset_mag, mid[1]+axis[1]*offset_mag, mid[2]+axis[2]*offset_mag];
        const pts = quadBezierPoints(p0,p1,p2,steps);
        const xs = [], ys = [], zs = [];
        for (let i=0;i<pts.length;i++){{ xs.push(pts[i][0]); ys.push(pts[i][1]); zs.push(pts[i][2]); }}
        xs.push(null); ys.push(null); zs.push(null);
        perDim[dim].push({x: xs, y: ys, z: zs});
    }}
    const out = [ {x:[], y:[], z:[]}, {x:[], y:[], z:[]}, {x:[], y:[], z:[]} ];
    for (let d=0; d<3; d++){{
        for (const seg of perDim[d]){{ out[d].x = out[d].x.concat(seg.x); out[d].y = out[d].y.concat(seg.y); out[d].z = out[d].z.concat(seg.z); }}
    }}
    return out;
}}

document.addEventListener('DOMContentLoaded', function() {{
    const gd = document.querySelectorAll('.plotly-graph-div')[0];
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
    const btn = document.createElement('button');
    btn.textContent = 'Zoom & Highlight';
    btn.style.padding = '4px 8px';
    controls.appendChild(input);
    controls.appendChild(btn);
    if (gd && gd.parentNode) gd.parentNode.insertBefore(controls, gd);
    let lockedNode = null;
    let hoverNode = null;
    function zoomToNode(nodeId) {{
        if(!payload.coords[nodeId]) {{ alert('Node "' + nodeId + '" not found'); return; }}
        const c = payload.coords[nodeId];
        const range = 3;
        const xmin = c[0]-range, xmax = c[0]+range;
        const ymin = c[1]-range, ymax = c[1]+range;
        const zmin = c[2]-range, zmax = c[2]+range;
        Plotly.relayout(gd, {{"scene.xaxis.range": [xmin, xmax], "scene.yaxis.range": [ymin, ymax], "scene.zaxis.range": [zmin, zmax], "scene.camera.eye": {{"x":0.6,"y":0.6,"z":0.6}}}});
        // highlight neighbours and lock
        if(!payload.neighbors[nodeId]) return;
        lockedNode = nodeId;
        const segs = buildLineSegments(nodeId);
        for(let d=0; d<3; d++){{ const idx = payload.hover_trace_idxs[d]; Plotly.restyle(gd, {{x: [segs[d].x], y: [segs[d].y], z: [segs[d].z]}}, [idx]); }}
    }}
    btn.addEventListener('click', function() {{ zoomToNode(input.value.trim()); }});
    gd.on('plotly_hover', function(evt){{
        try{{
            const point = evt.points[0];
            let hoverText = point.text || point.customdata || '';
            const nodeId = String(hoverText).split('<br>')[0];
            hoverNode = nodeId;
            if(lockedNode && lockedNode !== nodeId) return;
            if(!payload.neighbors[nodeId]) return;
            const segs = buildLineSegments(nodeId);
            for(let d=0; d<3; d++){{
                const idx = payload.hover_trace_idxs[d];
                Plotly.restyle(gd, {{x: [segs[d].x], y: [segs[d].y], z: [segs[d].z]}}, [idx]);
            }}
        }} catch(e){{ console.error('hover handler error', e); }}
    }});

    gd.on('plotly_click', function(evt){{
        if(evt && evt.points && evt.points.length>0){{
            const point = evt.points[0];
            const nodeId = String((point.text||'').split('<br>')[0]);
            if(!payload.neighbors[nodeId]) return;
            if(lockedNode === nodeId){{
                lockedNode = null;
                if(hoverNode){{
                    const segs = buildLineSegments(hoverNode);
                    for(let d=0; d<3; d++){{
                        const idx = payload.hover_trace_idxs[d];
                        Plotly.restyle(gd, {{x: [segs[d].x], y: [segs[d].y], z: [segs[d].z]}}, [idx]);
                    }}
                }} else {{
                    for(let d=0; d<3; d++){{
                        const idx = payload.hover_trace_idxs[d];
                        Plotly.restyle(gd, {{x: [[]], y: [[]], z: [[]]}}, [idx]);
                    }}
                }}
            }} else {{
                lockedNode = nodeId;
                const segs = buildLineSegments(nodeId);
                for(let d=0; d<3; d++){{
                    const idx = payload.hover_trace_idxs[d];
                    Plotly.restyle(gd, {{x: [segs[d].x], y: [segs[d].y], z: [segs[d].z]}}, [idx]);
                }}
            }}
        }} else {{
            lockedNode = null;
            for(let d=0; d<3; d++){{
                const idx = payload.hover_trace_idxs[d];
                Plotly.restyle(gd, {{x: [[]], y: [[]], z: [[]]}}, [idx]);
            }}
        }}
    }});

    gd.on('plotly_unhover', function(evt){{
        hoverNode = null;
        if(lockedNode) return;
        for(let d=0; d<3; d++){{
            const idx = payload.hover_trace_idxs[d];
            Plotly.restyle(gd, {{x: [[]], y: [[]], z: [[]]}}, [idx]);
        }}
    }});
}});
</script>
"""
        with open(out_html, 'w', encoding='utf-8') as fh:
            fh.write(html_str)
            fh.write(extra_js)
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

    target = args.target if args.target else nodes[0]
    plot_aa_network(nodes, coords_map, letters, target=target, out_html=args.out, counts_map=counts_map)


if __name__ == '__main__':
    main()


def plot_from_aa_network(aa_net, target: str = None, out_html: str = None):
    """Convenience wrapper to plot directly from a NetworksManager.aa_network

    Usage from Python:
        from scripts.plot_aa3_plotly import plot_from_aa_network
        plot_from_aa_network(networks.get_aa_network(), target='ABC')
    """
    if aa_net is None:
        raise RuntimeError('aa_net is None')
    nodes = sorted(aa_net.nodes.keys())
    coords_map, letters = build_coords_map(nodes)
    tgt = target if target else nodes[0]
    counts_map = {k: {'count': aa_net.nodes[k].get('count', 0), 'color': aa_net.nodes[k].get('color', 0)} for k in nodes}
    plot_aa_network(nodes, coords_map, letters, target=tgt, out_html=out_html, counts_map=counts_map)
