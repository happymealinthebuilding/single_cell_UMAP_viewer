import os
import zipfile
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import scanpy
import pandas as pd
from shiny import App, reactive, render, ui

# Define the UMAP processing function
def umap_process(data):
    adata = scanpy.read(data)
    scanpy.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata.raw = adata
    scanpy.pp.filter_cells(adata, min_genes=200)
    scanpy.pp.filter_cells(adata, max_genes=8000)
    scanpy.pp.filter_genes(adata, min_cells=3)
    scanpy.pp.normalize_total(adata, target_sum=1e4)
    scanpy.pp.log1p(adata)
    scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    scanpy.pp.scale(adata, max_value=10)
    scanpy.pp.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, color="Cell type", show=False)
    return umap_visualization()

def umap_visualization():
    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    plt.savefig(temp_file.name, bbox_inches='tight')
    plt.close()
    return temp_file.name

app_ui = ui.page_fluid(
    ui.input_file("file", None, accept=".zip"),
    ui.output_image("umap_plot"),
    ui.output_text("error_message") # added to display error messages
)

def server(input, output, session):
    @reactive.Calc
    def uploaded_file():
        files = input.file()
        if not files:
            return None

        if isinstance(files, list) and files:
            file = files[0]
        else:
            return {"error": "No file uploaded."}

        if not file['name'].endswith('.zip'):
            return {"error": "Invalid file type. Please upload a .zip file."}

        filepath = file['datapath']
        temp_dir = tempfile.TemporaryDirectory()
        zip_filepath = os.path.join(temp_dir.name, file['name'])
        os.rename(filepath, zip_filepath)

        try:
            with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
                zip_ref.extractall(temp_dir.name)

            print(f"Contents of unzipped directory: {os.listdir(temp_dir.name)}") # for debugging

            h5ad_file = next((f for f in os.listdir(temp_dir.name) if f.endswith('.h5ad')), None)
            if h5ad_file:
                return os.path.join(temp_dir.name, h5ad_file)
            else:
                return {"error": "No .h5ad file found in the ZIP."}
        except Exception as e:
            return {"error": f"Error processing ZIP file: {e}"}

    @output
    @render.image
    def umap_plot():
        filepath = uploaded_file()
        if isinstance(filepath, dict) and "error" in filepath:
            output.error_message.set(filepath["error"]) # display error
            return None
        else:
            output.error_message.set("") # clear previous error
        if filepath:
            return {"src": umap_process(filepath), "alt": "UMAP Plot"}

    @output
    @render.text
    def error_message():
        return "" # initial error message
app = App(app_ui, server)