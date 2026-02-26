import os

# NCBI Entrez requires an email address to track API usage.
ENTREZ_EMAIL = "example@gmail.com"

# Define base directory (assuming config.py is in the root of your project)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Define Data Directories
DATA_DIR = os.path.join(BASE_DIR, "data")
RAW_DATA_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DATA_DIR = os.path.join(DATA_DIR, "processed")
STRUCTURES_DIR = os.path.join(DATA_DIR, "structures")

# Define Results Directories
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ALIGNMENTS_DIR = os.path.join(RESULTS_DIR, "alignments")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
METRICS_DIR = os.path.join(RESULTS_DIR, "metrics")

# Automatically create directories if they don't exist
DIRECTORIES = [
    RAW_DATA_DIR, 
    PROCESSED_DATA_DIR, 
    STRUCTURES_DIR, 
    ALIGNMENTS_DIR, 
    FIGURES_DIR, 
    METRICS_DIR
]

for dir_path in DIRECTORIES:
    os.makedirs(dir_path, exist_ok=True)

print("Configuration loaded and directory structure verified.")