{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
}
import shiny
from shiny import App, render, ui, reactive
import scanpy as sc
import matplotlib.pyplot as plt
import zipfile
import os
import shutil  # For removing directories

app_ui = ui.page_fluid(
    ui.input_file("file", "Upload a ZIP file containing an H5AD file", accept=".zip"),
    ui.output_plot("umap_plot"),
    ui.output_text("error_message")  # For displaying error messages
)

def server(input, output, session):
    temp_dir = "extracted_data"  # Define a temporary directory

    @reactive.Effect
    @reactive.observe
    def cleanup_temp_dir(): # Clean up on app start and when a new file is uploaded
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)  # Remove directory and all its contents
        os.makedirs(temp_dir, exist_ok=True) # Recreate the directory

    @output
    @render.plot
    def umap_plot():
        uploaded_file = input.file()
        if uploaded_file is None:
            return "Please upload a file"

        error_message = ""  # Initialize error message

        try:
            # Unzip the file
            with zipfile.ZipFile(uploaded_file["path"], 'r') as zip_ref:
                zip_ref.extractall(temp_dir)

            # Find the .h5ad file
            h5ad_file = None
            for file in os.listdir(temp_dir):
                if file.endswith(".h5ad"):
                    h5ad_file = file
                    break

            if h5ad_file is None:
                raise FileNotFoundError("No .h5ad file found in the ZIP.")

            # Load and plot the data
            adata = sc.read_h5ad(os.path.join(temp_dir, h5ad_file)) # Use os.path.join
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            fig, ax = plt.subplots()
            sc.pl.umap(adata, ax=ax)
            return fig

        except zipfile.BadZipFile as e:
            error_message = f"Invalid ZIP file: {e}"
        except FileNotFoundError as e:
            error_message = str(e)  # Convert exception to string for display
        except ValueError as e: # Handle potential scanpy errors
            error_message = f"Error processing data: {e}. Check if the .h5ad file is valid."
        except Exception as e: # Catch other potential errors
            error_message = f"An unexpected error occurred: {e}"
        finally:
            if error_message: # If there was an error, print and display it
                print(error_message) # Print to console for debugging
                return error_message  # Return the error message to display in the UI

    @output
    @render.text
    def error_message(): # This output displays the error message in the UI
        uploaded_file = input.file()
        if uploaded_file is None:
            return "" # Return empty string if no file uploaded yet

        # Call the umap_plot function to trigger error handling
        umap_plot()
        return "" # Return empty string as errors are now handled directly in umap_plot

app = App(app_ui, server)

# To run in Jupyter Notebook:
# app.run()

# To run as a standalone app:
if __name__ == "__main__":
    app.run()
    
