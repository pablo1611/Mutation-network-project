import customtkinter as ctk
from tkinter import filedialog, messagebox
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.load_clones import load_clones
from src.build_triplet_df import build_triplet_df
from src.networks import compute_network_distance
from scripts.plot_aa3_plotly import plot_from_aa_network


# Example imports from your logic layer
# from src.generate_triplets import generate_triplets  # Not used
# from src.generate_nonuplets import generate_nonuplets  # Not used

class AntibodySequenceLoaderApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Window config
        self.title("Antibody Sequence Loader")
        self.geometry("800x700")  # Increased size for two datasets
        ctk.set_appearance_mode("light")
        ctk.set_default_color_theme("blue")

        # Internal state for datasets (dynamic, up to 2)
        self.datasets = []  # Will hold dicts for each dataset
        self.dataset_frames = []  # List of frames

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
        self.btn_add_dataset.grid(row=0, column=0, padx=10, pady=5)
        self.btn_compare_areas = ctk.CTkButton(self.buttons_frame, text="Compare Areas", command=self.compare_areas)
        self.btn_compare_datasets = ctk.CTkButton(self.buttons_frame, text="Compare Between Datasets", command=self.compare_datasets)
        self.btn_clear = ctk.CTkButton(self.buttons_frame, text="Clear All", fg_color="#e5e7eb", text_color="black", command=self.clear_form)
        self.btn_clear.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
        self.update_buttons()

        # Store clones for use in execute_triplets
        self.clones_r1 = None  
        self.clones_r2 = None
        self.triplet_df_r1 = None
        self.triplet_df_r2 = None
        self.networks_r1 = None
        self.networks_r2 = None

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

        btn_execute = ctk.CTkButton(dataset_buttons_frame, text="Execute Triplets Occurrences", command=lambda idx=index: self.execute_triplets(idx))
        # Don't pack initially

        btn_analyze = ctk.CTkButton(dataset_buttons_frame, text="Analyze Dataset", command=lambda idx=index: self.analyze_dataset(idx))
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
                plot_from_aa_network(networks.get_aa_network())
        except Exception as e:
            print(f"Debug: Exception in analyze_dataset: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to analyze dataset {index+1}: {e}")
            return
        
    
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
            messagebox.showinfo("Success", f"Triplets saved and opened: {out_path}")
        elif reigon_index == 1 and self.triplet_df_r2 is not None:
            out_dir = os.path.join(os.getcwd(), "output")
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"triplets_dataset_2.csv")
            self.triplet_df_r2.to_csv(out_path, index=False)
            subprocess.run(['open', out_path])
            messagebox.showinfo("Success", f"Triplets saved and opened: {out_path}")
        else:
            messagebox.showerror("Error", "Triplets not available. Please analyze the dataset first.")

    def compare_datasets(self):
        """Compare between the two datasets"""
        if not all(d['extraction_success'] for d in self.datasets):
            messagebox.showerror("Error", "Both datasets must be processed successfully.")
            return
        if self.networks_r1 and self.networks_r2:
            # Compute distance for AA triplet networks in R1 and R2, average them
            aa_net_r1 = self.networks_r1.aa_network
            aa_net_r2 = self.networks_r2.aa_network
            distance = compute_network_distance(aa_net_r1, aa_net_r2)
            messagebox.showinfo("Comparison Result", f"Average distance between datasets: {distance:.4f}")
        else:
            self.analyze_dataset(0)
            self.analyze_dataset(1)
            self.compare_datasets()  # Retry after analysis
            messagebox.showerror("Error", "Failed to build networks for comparison.")

    def compare_areas(self):
        """Compare between R1 and R2 areas in the uploaded dataset"""
        if not self.networks_r1:
            self.analyze_dataset(0)
            self.compare_areas()  # Retry after analysis
            return
        try:
            r1_net = self.networks_r1.aa_network_region['r1']
            r2_net = self.networks_r1.aa_network_region['r2']
            distance = compute_network_distance(r1_net, r2_net)
            messagebox.showinfo("Compare Areas", f"Distance between R1 and R2 in Dataset 1: {distance:.4f}")
        except Exception as e:
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