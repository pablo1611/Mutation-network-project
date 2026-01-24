import customtkinter as ctk
from tkinter import filedialog, messagebox
import os
import subprocess
import sys
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import socket
import pandas as pd
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.load_clones import load_clones
from src.build_triplet_df import build_triplet_df
from src.networks import compute_network_distance
from scripts.plot_aa3_plotly import plot_from_aa_network


# Example imports from your logic layer
# from src.generate_triplets import generate_triplets  # Not used
# from src.generate_nonuplets import generate_nonuplets  # Not used

def save_edges_to_csv(network, output_path, dataset_name="", analysis_type=""):
    """Save all edges with their probabilities (weights) to a CSV file.
    
    Args:
        network: KmerNetwork instance with edges computed
        output_path: Path where the CSV file will be saved
        dataset_name: Optional name/identifier for the dataset
        analysis_type: Optional type of analysis (e.g., 'single_dataset', 'compare_datasets', 'compare_regions')
    
    Returns:
        str: Path to the saved CSV file
    """
    if not hasattr(network, 'edges') or not network.edges:
        print(f"Warning: No edges found in network for {dataset_name}")
        return None
    
    # Collect all edges with their data
    edges_data = []
    for source_node, targets in network.edges.items():
        for target_node, edge_info in targets.items():
            probability = edge_info.get('probability', 0.0)
            above_threshold = edge_info.get('above_threshold', True)
            
            edges_data.append({
                'source_node': source_node,
                'target_node': target_node,
                'weight': probability,
                'above_threshold': above_threshold,
                'dataset': dataset_name,
                'analysis_type': analysis_type,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            })
    
    # Create DataFrame
    df = pd.DataFrame(edges_data)
    
    # Sort by probability/weight in descending order
    df = df.sort_values('weight', ascending=False)
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} edges to {output_path}")
    
    return output_path

class AntibodySequenceLoaderApp(ctk.CTk):
    def __init__(self):
        # Start a persistent local HTTP server to serve generated plot HTML files.
        # This avoids starting per-plot servers and makes the behavior identical
        # to running `app.py` from source where Plotly uses a localhost URL.
        try:
            serve_dir = os.getcwd()
            class _Handler(SimpleHTTPRequestHandler):
                def __init__(self, *args, **kwargs):
                    super().__init__(*args, directory=serve_dir, **kwargs)

            # find a free port and start server in background
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as _s:
                _s.bind(('127.0.0.1', 0))
                _port = _s.getsockname()[1]

            try:
                httpd = ThreadingHTTPServer(('127.0.0.1', _port), _Handler)
                def _serve():
                    try:
                        httpd.serve_forever()
                    except Exception:
                        pass
                t = threading.Thread(target=_serve, daemon=True)
                t.start()
                # expose the server port to other modules via env var
                os.environ['PLOT_SERVER_PORT'] = str(_port)
            except Exception:
                # if server fails, don't block app startup
                pass
        except Exception:
            pass
        super().__init__()

        # Cursor for hover on clickable widgets
        self._hand_cursor = "pointinghand" if sys.platform == "darwin" else "hand2"

        # Window config
        self.title("Antibody Sequence Loader")
        self.geometry("800x700")  # Increased size for two datasets
        ctk.set_appearance_mode("light")
        ctk.set_default_color_theme("blue")

        # Internal state for datasets (dynamic, up to 2)
        self.datasets = []  
        self.dataset_frames = []  

        # Make the window scrollable
        self.scrollable_frame = ctk.CTkScrollableFrame(self, width=800, height=600)
        self.scrollable_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Header
        self.label_title = ctk.CTkLabel(self.scrollable_frame, text="Antibody Sequence Loader", font=ctk.CTkFont(size=20, weight="bold"))
        self.label_title.pack(pady=(25, 15))

        # Container for datasets
        self.datasets_container = ctk.CTkFrame(self.scrollable_frame)
        self.datasets_container.pack(pady=10, fill="x", padx=20)

        # Create initial frame for Dataset 1
        self.add_dataset_frame()

        # Buttons frame
        self.buttons_frame = ctk.CTkFrame(self.scrollable_frame)
        self.buttons_frame.pack(pady=10)
        # Add buttons
        self.btn_add_dataset = ctk.CTkButton(self.buttons_frame, text="Add Another Dataset", command=self.add_another_dataset)
        self._apply_click_cursor(self.btn_add_dataset)
        self.btn_add_dataset.grid(row=0, column=0, padx=10, pady=5)
        self.btn_compare_areas = ctk.CTkButton(self.buttons_frame, text="Compare Areas", command=self.compare_areas)
        self.btn_compare_datasets = ctk.CTkButton(self.buttons_frame, text="Compare Between Datasets", command=self.compare_datasets)
        self.btn_clear = ctk.CTkButton(self.buttons_frame, text="Clear All", fg_color="#e5e7eb", text_color="black", command=self.clear_form)
        self._apply_click_cursor(self.btn_compare_areas)
        self._apply_click_cursor(self.btn_compare_datasets)
        self._apply_click_cursor(self.btn_clear)
        self.btn_clear.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
        self.update_buttons()

        # Store clones for use in execute_triplets
        self.clones_r1 = None  
        self.clones_r2 = None
        self.triplet_df_r1 = None
        self.triplet_df_r2 = None
        self.networks_r1 = None
        self.networks_r2 = None
        
        # Threshold settings
        self.node_threshold = 0.0  # Default: no filtering
        self.edge_threshold = 0.0  # Default: no filtering

    def _apply_click_cursor(self, widget):
        """Make hover cursor indicate clickability (hand/finger)."""
        try:
            widget.configure(cursor=self._hand_cursor)
        except Exception:
            # Some Tk builds only support 'hand2'
            try:
                widget.configure(cursor="hand2")
            except Exception:
                pass

    def add_dataset_frame(self):
        index = len(self.datasets)
        title = f"Dataset {index + 1}"
        frame = ctk.CTkFrame(self.datasets_container)
        frame.pack(pady=5, fill="x", padx=10)

        # Title
        title_label = ctk.CTkLabel(frame, text=title, font=ctk.CTkFont(size=16, weight="bold"))
        title_label.pack(pady=(10, 5))

        # Dropzone
        drop_frame = ctk.CTkFrame(frame, width=300, height=60, corner_radius=10)
        drop_frame.pack(pady=5)
        drop_label = ctk.CTkLabel(drop_frame, text="Select your Dataset", text_color="#555")
        drop_label.place(relx=0.5, rely=0.5, anchor="center")

        # Browse button
        btn_browse = ctk.CTkButton(frame, text="Browse", command=lambda idx=index: self.browse_file(idx), width=100)
        self._apply_click_cursor(btn_browse)
        btn_browse.pack(pady=5)

        # Metadata frame
        metadata_frame = ctk.CTkFrame(frame)
        metadata_frame.pack(pady=5, fill="x", padx=10)

        meta_title = ctk.CTkLabel(metadata_frame, text="Metadata Preview", font=ctk.CTkFont(size=14, weight="bold"))
        meta_title.pack(anchor="w", pady=(5, 2))

        label_experiment = ctk.CTkLabel(metadata_frame, text="Experiment: –")
        label_experiment.pack(anchor="w", padx=10)
        label_group = ctk.CTkLabel(metadata_frame, text="Group ID: –")
        label_group.pack(anchor="w", padx=10)
        label_source = ctk.CTkLabel(metadata_frame, text="Sample Source: –")
        label_source.pack(anchor="w", padx=10)
        label_clones = ctk.CTkLabel(metadata_frame, text="Total Clones: –")
        label_clones.pack(anchor="w", padx=10)
        label_time = ctk.CTkLabel(metadata_frame, text="Time Points: –")
        label_time.pack(anchor="w", padx=10)

        # Dataset buttons frame
        dataset_buttons_frame = ctk.CTkFrame(frame)
        # Don't pack initially

        # Threshold controls frame
        threshold_frame = ctk.CTkFrame(frame)
        threshold_frame.pack(pady=5, fill="x", padx=10)
        
        threshold_title = ctk.CTkLabel(threshold_frame, text="Network Filtering Thresholds", font=ctk.CTkFont(size=12, weight="bold"))
        threshold_title.pack(anchor="w", pady=(5, 2))
        
        # Node threshold slider
        node_threshold_label = ctk.CTkLabel(threshold_frame, text="Node Frequency Threshold: 0.000")
        node_threshold_label.pack(anchor="w", padx=10)
        node_threshold_slider = ctk.CTkSlider(threshold_frame, from_=0.0, to=0.01, number_of_steps=100, 
                                             command=lambda v, lbl=node_threshold_label, idx=index: self.update_node_threshold(v, lbl, idx))
        node_threshold_slider.set(0.0)
        node_threshold_slider.pack(fill="x", padx=10, pady=2)
        
        # Edge threshold slider
        edge_threshold_label = ctk.CTkLabel(threshold_frame, text="Edge Probability Threshold: 0.00")
        edge_threshold_label.pack(anchor="w", padx=10)
        edge_threshold_slider = ctk.CTkSlider(threshold_frame, from_=0.0, to=1.0, number_of_steps=100,
                                             command=lambda v, lbl=edge_threshold_label, idx=index: self.update_edge_threshold(v, lbl, idx))
        edge_threshold_slider.set(0.0)
        edge_threshold_slider.pack(fill="x", padx=10, pady=2)

        btn_execute = ctk.CTkButton(dataset_buttons_frame, text="Execute Triplets Occurrences", command=lambda idx=index: self.execute_triplets(idx))
        self._apply_click_cursor(btn_execute)
        # Don't pack initially

        btn_analyze = ctk.CTkButton(dataset_buttons_frame, text="Analyze Dataset", command=lambda idx=index: self.analyze_dataset(idx))
        self._apply_click_cursor(btn_analyze)
        # Don't pack initially

        # Store in lists
        self.datasets.append({
            'path': None,
            'metadata': {},
            'extraction_success': False,
            'drop_label': drop_label,
            'metadata_labels': {
                'experiment': label_experiment,
                'group': label_group,
                'source': label_source,
                'clones': label_clones,
                'time': label_time
            },
            'buttons': {
                'execute': btn_execute,
                'analyze': btn_analyze,
                'frame': dataset_buttons_frame
            },
            'thresholds': {
                'node_threshold': 0.0,
                'edge_threshold': 0.0,
                'node_slider': node_threshold_slider,
                'edge_slider': edge_threshold_slider,
                'node_label': node_threshold_label,
                'edge_label': edge_threshold_label
            }
        })
        self.dataset_frames.append(frame)

    def add_another_dataset(self):
        self.add_dataset_frame()
        self.btn_add_dataset.grid_forget()

    def update_buttons(self):
        # Clear existing buttons except clear and add
        for widget in self.buttons_frame.winfo_children():
            if widget not in [self.btn_clear, self.btn_add_dataset, self.btn_compare_areas, self.btn_compare_datasets]:
                widget.destroy()

        uploaded_count = sum(1 for d in self.datasets if d['path'] is not None)

        if len(self.dataset_frames) < 2:
            self.btn_add_dataset.grid(row=0, column=0, padx=10, pady=5)
            if uploaded_count == 1:
                self.btn_compare_areas.grid(row=0, column=1, padx=10, pady=5)
            else:
                self.btn_compare_areas.grid_forget()
            self.btn_compare_datasets.grid_forget()
        else:
            self.btn_add_dataset.grid_forget()
            self.btn_compare_areas.grid_forget()
            if uploaded_count == 2:
                self.btn_compare_datasets.grid(row=0, column=0, padx=10, pady=5)
            else:
                self.btn_compare_datasets.grid_forget()
        

    def browse_file(self, index):
        if all(d['path'] is not None for d in self.datasets):
            messagebox.showinfo("All datasets uploaded", "All available datasets have been uploaded.")
            return
        file_path = filedialog.askopenfilename(
            title="Select a CSV File",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if file_path:
            self.datasets[index]['path'] = file_path
            self.datasets[index]['drop_label'].configure(text=os.path.basename(file_path))
            # Process the file to extract metadata
            self.parse_dataset_and_generate(index)
            self.update_buttons()

    def parse_dataset_and_generate(self, index):
        """Parse dataset and populate metadata"""
        try:
            import pandas as pd
            df = pd.read_csv(self.datasets[index]['path'])
            # Extract metadata: e.g., unique subjects, samples, etc.
            experiment = df['subject_id'].unique().tolist() if 'subject_id' in df else ['Unknown']
            group_id = df['sample_id'].unique().tolist() if 'sample_id' in df else ['Unknown']
            sample_source = df['ab_target'].unique().tolist() if 'ab_target' in df else ['Unknown']
            time_points = df['time_point'].unique().tolist() if 'time_point' in df else ['Unknown']
            total_clones = len(df)

            self.datasets[index]['metadata'] = {
                "experiment": ", ".join(map(str, experiment[:3])),  # Limit to first 3
                "group_id": ", ".join(map(str, group_id[:3])),
                "sample_source": ", ".join(map(str, sample_source[:3])),
                "total_clones": str(total_clones),
                "time_points": ", ".join(map(str, time_points[:3]))
            }

            # Update UI labels
            labels = self.datasets[index]['metadata_labels']
            labels['experiment'].configure(text=f"Experiment: {self.datasets[index]['metadata']['experiment']}")
            labels['group'].configure(text=f"Group ID: {self.datasets[index]['metadata']['group_id']}")
            labels['source'].configure(text=f"Sample Source: {self.datasets[index]['metadata']['sample_source']}")
            labels['clones'].configure(text=f"Total Clones: {self.datasets[index]['metadata']['total_clones']}")
            labels['time'].configure(text=f"Time Points: {self.datasets[index]['metadata']['time_points']}")

            self.datasets[index]['extraction_success'] = True

            # Show buttons frame and buttons
            self.datasets[index]['buttons']['frame'].pack(pady=5, fill="x", padx=10)
            self.datasets[index]['buttons']['execute'].pack(side="left", padx=5, pady=5)
            self.datasets[index]['buttons']['analyze'].pack(side="left", padx=5, pady=5)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse dataset {index+1}: {e}")
            self.datasets[index]['extraction_success'] = False
            # Hide buttons frame
            self.datasets[index]['buttons']['frame'].pack_forget()

    def analyze_dataset(self, index):
        """Analyze the dataset and show visualiization"""
        print(f"Debug: Starting analyze_dataset for index {index}")
        if not self.datasets[index]['extraction_success']:
            messagebox.showerror("Error", f"Dataset {index+1} not processed successfully.")
            return
        try:
            print(f"Debug: Loading clones from {self.datasets[index]['path']}")
            # Load clones
            result = load_clones(self.datasets[index]['path'], build_networks=True)
            print(f"Debug: Load result type: {type(result)}")
            if isinstance(result, tuple):
                clones, networks = result
                print(f"Debug: Clones loaded: {len(clones) if clones else 0}, Networks: {networks is not None}")
            else:
                clones = result
                networks = None
                print(f"Debug: Clones loaded (no networks): {len(clones) if clones else 0}")
            # Build triplet DataFrame
            print("Debug: Building triplet_df")
            triplet_df = build_triplet_df(clones)
            print(f"Debug: Triplet_df shape: {triplet_df.shape}")
            if index == 0:
                self.clones_r1 = clones
                self.triplet_df_r1 = triplet_df
                self.networks_r1 = networks if networks else None
                print("Debug: Stored in r1")
            else:
                self.clones_r2 = clones
                self.triplet_df_r2 = triplet_df
                self.networks_r2 = networks if networks else None
                print("Debug: Stored in r2")
            if networks:
                # Get threshold values for this dataset
                node_thresh = self.datasets[index]['thresholds']['node_threshold']
                edge_thresh = self.datasets[index]['thresholds']['edge_threshold']
                print(f"Debug: Using thresholds - node: {node_thresh}, edge: {edge_thresh}")
                out_dir = os.path.join(os.getcwd(), "output")
                os.makedirs(out_dir, exist_ok=True)
                out_html = os.path.join(out_dir, f"aa3_network_dataset_{index+1}.html")
                
                # Get the network and apply thresholds
                aa_net = networks.get_aa_network()
                aa_net.normalize_nodes()
                aa_net.compute_edge_probabilities()
                aa_net.apply_node_threshold(node_thresh)
                aa_net.apply_edge_threshold(edge_thresh)
                
                # Save edges to CSV automatically
                edges_csv_path = os.path.join(out_dir, f"edges_dataset_{index+1}.csv")
                save_edges_to_csv(
                    aa_net,
                    edges_csv_path,
                    dataset_name=f"Dataset_{index+1}",
                    analysis_type="single_dataset_analysis"
                )
                
                plot_from_aa_network(
                    aa_net,
                    target=None,
                    out_html=out_html,
                    node_threshold=node_thresh,
                    edge_threshold=edge_thresh,
                )

                # Open the generated HTML via the local server (so embedded JS runs reliably)
                port = os.environ.get('PLOT_SERVER_PORT')
                if port:
                    rel_path = os.path.relpath(out_html, os.getcwd()).replace(os.sep, '/')
                    url = f"http://127.0.0.1:{port}/{rel_path}"
                    subprocess.run(['open', url])
                else:
                    subprocess.run(['open', out_html])
        except Exception as e:
            print(f"Debug: Exception in analyze_dataset: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to analyze dataset {index+1}: {e}")
            return

    def update_node_threshold(self, value, label, index):
        """Update node threshold slider value"""
        self.datasets[index]['thresholds']['node_threshold'] = float(value)
        label.configure(text=f"Node Frequency Threshold: {float(value):.6f}")
    
    def update_edge_threshold(self, value, label, index):
        """Update edge threshold slider value"""
        self.datasets[index]['thresholds']['edge_threshold'] = float(value)
        label.configure(text=f"Edge Probability Threshold: {float(value):.4f}")
    
    def execute_triplets(self, reigon_index):
        if self.clones_r1 is None and self.clones_r2 is None:
            self.analyze_dataset(reigon_index)
        # save the triplet_df
        if reigon_index == 0 and self.triplet_df_r1 is not None:
            out_dir = os.path.join(os.getcwd(), "output")
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"triplets_dataset_1.csv")
            self.triplet_df_r1.to_csv(out_path, index=False)
            subprocess.run(['open', out_path])
        elif reigon_index == 1 and self.triplet_df_r2 is not None:
            out_dir = os.path.join(os.getcwd(), "output")
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"triplets_dataset_2.csv")
            self.triplet_df_r2.to_csv(out_path, index=False)
            subprocess.run(['open', out_path])
        else:
            messagebox.showerror("Error", "Triplets not available. Please analyze the dataset first.")

    def compare_datasets(self):
        """Compare between the two datasets"""
        if not all(d['extraction_success'] for d in self.datasets):
            messagebox.showerror("Error", "Both datasets must be processed successfully.")
            return
        if self.networks_r1 and self.networks_r2:
            try:
                # Compute distance for AA triplet networks in R1 and R2
                aa_net_r1 = self.networks_r1.aa_network
                aa_net_r2 = self.networks_r2.aa_network
                
                # Apply normalization and thresholds
                aa_net_r1.normalize_nodes()
                aa_net_r1.compute_edge_probabilities()
                aa_net_r1.apply_node_threshold(self.datasets[0]['thresholds']['node_threshold'])
                aa_net_r1.apply_edge_threshold(self.datasets[0]['thresholds']['edge_threshold'])
                
                aa_net_r2.normalize_nodes()
                aa_net_r2.compute_edge_probabilities()
                aa_net_r2.apply_node_threshold(self.datasets[1]['thresholds']['node_threshold'])
                aa_net_r2.apply_edge_threshold(self.datasets[1]['thresholds']['edge_threshold'])
                
                # Compute basic distance metric
                distance = compute_network_distance(aa_net_r1, aa_net_r2)
                
                # Save edges for both datasets
                out_dir = os.path.join(os.getcwd(), "output")
                os.makedirs(out_dir, exist_ok=True)
                
                edges_csv_r1 = os.path.join(out_dir, "edges_dataset1_comparison.csv")
                save_edges_to_csv(
                    aa_net_r1,
                    edges_csv_r1,
                    dataset_name="Dataset_1",
                    analysis_type="compare_datasets"
                )
                
                edges_csv_r2 = os.path.join(out_dir, "edges_dataset2_comparison.csv")
                save_edges_to_csv(
                    aa_net_r2,
                    edges_csv_r2,
                    dataset_name="Dataset_2",
                    analysis_type="compare_datasets"
                )
                
                summary = f"Dataset Comparison Results:\n\n"
                summary += f"Dataset 1: {len(aa_net_r1.nodes)} nodes\n"
                summary += f"Dataset 2: {len(aa_net_r2.nodes)} nodes\n\n"
                summary += f"Network Distance: {distance:.6f}\n\n"
                summary += f"\nEdges saved to:\n- {edges_csv_r1}\n- {edges_csv_r2}"
                
                messagebox.showinfo("Comparison Result", summary)

                # Do not open or save comparison HTML
                    
            except Exception as e:
                import traceback
                traceback.print_exc()
                messagebox.showerror("Error", f"Failed to compare datasets: {e}")
        else:
            self.analyze_dataset(0)
            self.analyze_dataset(1)
            self.compare_datasets()  # Retry after analysis

    def compare_areas(self):
        """Compare between R1 and R2 areas in the uploaded dataset"""
        if not self.networks_r1:
            self.analyze_dataset(0)
            self.compare_areas()  # Retry after analysis
            return
        try:
            r1_net = self.networks_r1.aa_network_region['r1']
            r2_net = self.networks_r1.aa_network_region['r2']
            
            # Apply normalization and thresholds
            r1_net.normalize_nodes()
            r1_net.compute_edge_probabilities()
            node_thresh = self.datasets[0]['thresholds']['node_threshold']
            edge_thresh = self.datasets[0]['thresholds']['edge_threshold']
            r1_net.apply_node_threshold(node_thresh)
            r1_net.apply_edge_threshold(edge_thresh)
            
            r2_net.normalize_nodes()
            r2_net.compute_edge_probabilities()
            r2_net.apply_node_threshold(node_thresh)
            r2_net.apply_edge_threshold(edge_thresh)
            
            # Compute distance
            distance = compute_network_distance(r1_net, r2_net)
            
            # Save edges for both regions
            out_dir = os.path.join(os.getcwd(), "output")
            os.makedirs(out_dir, exist_ok=True)
            
            edges_csv_r1 = os.path.join(out_dir, "edges_region1_comparison.csv")
            save_edges_to_csv(
                r1_net,
                edges_csv_r1,
                dataset_name="Region_1",
                analysis_type="compare_regions"
            )
            
            edges_csv_r2 = os.path.join(out_dir, "edges_region2_comparison.csv")
            save_edges_to_csv(
                r2_net,
                edges_csv_r2,
                dataset_name="Region_2",
                analysis_type="compare_regions"
            )
            
            summary = f"Region Comparison Results (R1 vs R2):\n\n"
            summary += f"Region 1: {len(r1_net.nodes)} nodes\n"
            summary += f"Region 2: {len(r2_net.nodes)} nodes\n\n"
            summary += f"Network Distance: {distance:.6f}\n\n"
            summary += f"\nEdges saved to:\n- {edges_csv_r1}\n- {edges_csv_r2}"
            
            messagebox.showinfo("Compare Areas", summary)

            # Do not open or save comparison HTML
                
        except Exception as e:
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to compare areas: {e}")

    def clear_form(self):
        """Reset all datasets and UI."""
        # Remove extra frames
        while len(self.dataset_frames) > 1:
            frame = self.dataset_frames.pop()
            frame.destroy()
        self.datasets = self.datasets[:1]  # Keep only first
        # Reset first dataset
        self.datasets[0]['path'] = None
        self.datasets[0]['metadata'] = {}
        self.datasets[0]['extraction_success'] = False
        self.datasets[0]['drop_label'].configure(text="Select your Dataset")
        labels = self.datasets[0]['metadata_labels']
        labels['experiment'].configure(text="Experiment: –")
        labels['group'].configure(text="Group ID: –")
        labels['source'].configure(text="Sample Source: –")
        labels['clones'].configure(text="Total Clones: –")
        labels['time'].configure(text="Time Points: –")        # Hide buttons frame
        self.datasets[0]['buttons']['frame'].pack_forget()        
        self.btn_add_dataset.grid(row=0, column=0, padx=10, pady=5)  # Re-show add button
        self.update_buttons()
        self.clones_r1 = None  
        self.clones_r2 = None
        self.triplet_df_r1 = None
        self.triplet_df_r2 = None
        self.networks_r1 = None
        self.networks_r2 = None
        

if __name__ == "__main__":
    app = AntibodySequenceLoaderApp()
    app.mainloop()