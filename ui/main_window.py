import customtkinter as ctk
from tkinter import filedialog, messagebox
from src.sequence_processor import SequenceProcessor
import json
import os

# Example imports from your logic layer
from src.generate_triplets import generate_triplets
from src.generate_nonuplets import generate_nonuplets

class AntibodySequenceLoaderApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Window config
        self.title("Antibody Sequence Loader")
        self.geometry("560x450")
        ctk.set_appearance_mode("light")
        ctk.set_default_color_theme("blue")

        # Internal state
        self.file_path = None
        self.metadata = {}
        self.extraction_success = False

        # Header
        self.label_title = ctk.CTkLabel(self, text="Antibody Sequence Loader", font=ctk.CTkFont(size=20, weight="bold"))
        self.label_title.pack(pady=(25, 15))

        # Dropzone / Browse button
        self.drop_frame = ctk.CTkFrame(self, width=400, height=80, corner_radius=10)
        self.drop_frame.pack(pady=10)
        self.drop_label = ctk.CTkLabel(self.drop_frame, text="Drag & Drop your file or click Browse", text_color="#555")
        self.drop_label.place(relx=0.5, rely=0.5, anchor="center")

        self.btn_browse = ctk.CTkButton(self, text="Browse", command=self.browse_file, width=120)
        self.btn_browse.pack(pady=(5, 10))

        # Metadata frame
        self.metadata_frame = ctk.CTkFrame(self)
        self.metadata_frame.pack(pady=10, fill="x", padx=20)

        self.label_meta_title = ctk.CTkLabel(
            self.metadata_frame,
            text="Metadata Preview",
            font=ctk.CTkFont(size=14, weight="bold"),
        )
        self.label_meta_title.pack(anchor="w", pady=(10, 2), padx=10)

        self.label_experiment = ctk.CTkLabel(self.metadata_frame, text="Experiment: –")
        self.label_experiment.pack(anchor="w", padx=20)
        self.label_group = ctk.CTkLabel(self.metadata_frame, text="Group ID: –")
        self.label_group.pack(anchor="w", padx=20)
        self.label_source = ctk.CTkLabel(self.metadata_frame, text="Sample Source: –")
        self.label_source.pack(anchor="w", padx=20)

        # Status label
        self.status_label = ctk.CTkLabel(self, text="Waiting for dataset upload...", text_color="#555")
        self.status_label.pack(pady=(15, 10))

        # Buttons
        self.frame_buttons = ctk.CTkFrame(self)
        self.frame_buttons.pack(pady=10)

        self.btn_clear = ctk.CTkButton(self.frame_buttons, text="Clear Form", fg_color="#e5e7eb", text_color="black", command=self.clear_form)
        self.btn_clear.grid(row=0, column=0, padx=10)
        self.btn_next = ctk.CTkButton(self.frame_buttons, text="Next", command=self.process_next)
        self.btn_next.grid(row=0, column=1, padx=10)

    def browse_file(self):
        file_path = filedialog.askopenfilename(
            title="Select a Sequence File",
            filetypes=[("Text Files", "*.txt"), ("JSON Files", "*.json"), ("All Files", "*.*")]
        )
        if file_path:
            self.file_path = file_path
            self.drop_label.configure(text=os.path.basename(file_path))
            self.status_label.configure(text="File loaded successfully", text_color="green")

    def clear_form(self):
        """Reset form and UI."""
        self.file_path = None
        self.metadata = {}
        self.extraction_success = False
        self.drop_label.configure(text="Drag & Drop your dataset or click Browse")
        self.status_label.configure(text="Form cleared", text_color="#555")
        self.label_experiment.configure(text="Experiment: –")
        self.label_group.configure(text="Group ID: –")
        self.label_source.configure(text="Sample Source: –")
        self.btn_next.configure(state="disabled")

    def process_next(self):
        """Triggered when user clicks Next"""
        if not self.extraction_success:
            messagebox.showerror("Error", "Extraction not successful yet.")
            return
        messagebox.showinfo("Next Step", "Proceeding to next step of analysis...")

    def parse_dataset_and_generate(self):
        """Parse dataset and populate metadata + generate structures"""
        try:
            # Temporary placeholder to simulate metadata extraction
            self.metadata = {
                "experiment": "Trial A",
                "group_id": "12345",
                "sample_source": "Peripheral Blood",
            }

            # Update UI labels
            self.label_experiment.configure(text=f"Experiment: {self.metadata['experiment']}")
            self.label_group.configure(text=f"Group ID: {self.metadata['group_id']}")
            self.label_source.configure(text=f"Sample Source: {self.metadata['sample_source']}")

            # Call generator functions
            from src.sequence_processor import SequenceProcessor  # local import to avoid circular dependency
            processor = SequenceProcessor(self.file_path)

            triplets_data, nonuplets_data = processor.process()
            triplet_path, nonup_path = processor.save_results()

            # Save outputs
            # TODO: think how we want to handle output paths
            print(f"✅ Triplets saved at: {triplet_path}")
            print(f"✅ Nonuplets saved at: {nonup_path}")

            # success
            self.status_label.configure(
                text="Triplets extracted successfully", text_color="green"
            )
            self.btn_next.configure(state="normal")  # enable Next
            self.extraction_success = True

        except Exception as e:
            # failure
            self.status_label.configure(
                text=f"Extraction failed: {e}", text_color="red"
            )
            self.btn_next.configure(state="disabled")
            self.extraction_success = False


 