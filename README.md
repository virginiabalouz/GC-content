# [`GCanner`](https://gcanner.streamlit.app/)
### An Interactive Tool for Genome-Wide GC Analysis  
**Authors:** Virginia Balouz and Carlos A. Buscaglia  

GCanner is a Python-based interactive tool designed for the genome-wide analysis of GC composition. It enables users to identify and classify genomic compartments in *Trypanosoma cruzi* and other organisms based on their GC content.

## 🌍 Online Version: GCanner Web App
The user-friendly GCanner app is available [here](https://gcanner.streamlit.app/)

## 💡 Quick Tutorial

### 1️⃣ Upload a FASTA File
- Click the **'Browse files'** button to upload a FASTA file for analysis.
- The default page allows uploading files up to **200 MB**.
- To modify this restriction, run the following command (available only for the locally-run option) :

  ```bash
  streamlit run GC-content-streamlit.py --server.maxUploadSize 2000
  ```
    this setting will allow the user to upload files up to 2GB
### 2️⃣ Modify Analysis Parameters
GC scanning and smoothing parameters are initially set to default values to determine core and disruptive compartments in *T. cruzi* but can be adjusted:
  1. On the **left panel**, locate the **'General Settings'** section.
  2. Adjust the parameters according to your needs:
     - Minimum sequence length for analysis.
     - Window size.
     - Step size.
     - Number of points for Lowess smoothing.
     - Threshold cutoff for GC-low and GC-rich regions.
     - Number of GC content plots to display.

### 3️⃣ GOI-based search
To identify whether a **Gene of Interest (GOI)** is located in a core or disruptive compartment in addition to uploading the FASTA file (as described in 1️⃣) follow these steps:

  1. Toggle the button for **GOI-based analysis**.
  2. Upload a GFF file with annotated features. **Note:**Ensure the sequence names in both files match exactly for accurate mapping.
  3. Enter the GOI keyword: In the box **Search by gene of interest (keyword)**, specify the keyword associated with your GOI to search in the GFF file.   

## 📊 Outputs
### 1️⃣ Table of Mapped Regions
This table contains six columns: an index column enumerating each row, the 'Sequence_name' column for the sequence's name, 'Start' and 'End' columns indicating the beginning and end of each region, and the 'Region type' and 'Region Length' columns, which specify whether the regions were classified as core or disruptive and their corresponding lengths.
The table can be downloaded as a CSV file by clicking the 📥 in the top right corner.

### 2️⃣ Compartment Proportion Plot
This barplot summarizes the proportion of bases classified as core (green) or disruptive (pink). The total count of bases for each compartment (in Kb) and their proportion are displayed at the top of each column.

### 3️⃣ GC Content Plots
These plots display the GC content (blue lines), Lowess smoothing (red line), and the predicted regions shaded in green for core regions and pink for disruptive regions along the DNA sequence. Each plot corresponds to one sequence in the FASTA file and is saved as **'plot_sequence_name.svg'**.

### 3️⃣ Keyword Match Compartment Proportion Plot (for GOI-based searches)
A bar plot summarizing the frequency of GOI occurrences within each compartment. Matches at compartment boundaries will be categorized as "Edge."

### 📥 Downloading Results
By clicking the **‘Download Figures as CSV’** button at the bottom part of the screen, all figures seen in the page will be saved separately in a compressed file named ‘plots_svg.zip’ to the default folder.

---

## 🚀 How to Run GCanner Locally
1. Download the script **`GC-content-streamlit.py`** and save it.
2. Open a terminal and navigate to the folder where the script is saved.
3. Run the following command:

   ```bash
   streamlit run GC-content-streamlit.py
   ```

The message *'You can now view your Streamlit app in your browser.'* will be printed on the command line, and the application will open in a new browser tab.

### 📦 Dependencies

- Python 3.7+
- Biopython
- NumPy and math
- Pandas
- Matplotlib and Seaborn
- Statsmodels
- Streamlit
- Empfile, Zipfile, and OS modules for temporary file handling and exporting outputs
---
✨ Enjoy using GCanner for your genome-wide GC analysis! 🚀
📖 Please cite Balouz et al.(2025)
