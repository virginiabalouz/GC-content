import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from io import StringIO
import tempfile
import zipfile
import os

# Importación tardía para diagnosticar si el error persiste
try:
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess
except ImportError as e:
    st.error(f"Error cargando statsmodels: {e}. Revisa el archivo requirements.txt")

# --- FUNCIONES ---
def f_smoothcal(dataframe, column, f_smooth_points=100):
    points = len(dataframe[column])
    if points == 0: return 1
    f_val = (f_smooth_points / points)
    return min(f_val, 1)

def calculate_gc_content(sequence, window, step):
    gc_content = []
    sequence_length = len(sequence)
    sequence = sequence.upper()
    for i in range(0, sequence_length, step):
        start = i
        end = min(i + window, sequence_length)
        window_sequence = sequence[start:end]
        if len(window_sequence) > 0:
            gc_count = window_sequence.count('G') + window_sequence.count('C')
            gc_percent = gc_count / len(window_sequence)
            gc_content.append((start + 1, end, gc_percent))
    return gc_content

def map_regions(df, fasta_len, cutoff):
    if df.empty or 'Smoothed_GC' not in df.columns: return []
    df['Core'] = np.where(df["Smoothed_GC"] < cutoff, 1, 0)
    df['Changes'] = df['Core'].diff().fillna(0).abs()
    hits = df[df['Changes'] == 1].index.tolist()
    regions = []
    
    cutoff_status = "Disruptive" if df["Smoothed_GC"].iloc[0] >= cutoff else "Core"
    last_pos = 1
    for change in hits:
        end = df['Start'].iloc[change]
        regions.append((last_pos, end, cutoff_status))
        last_pos = end + 1
        cutoff_status = "Core" if cutoff_status == "Disruptive" else "Disruptive"
    
    regions.append((last_pos, fasta_len, cutoff_status))
    return regions

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    if not handles: return
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc='upper right')

# --- SIDEBAR ---
with st.sidebar:
    st.subheader("Select analysis")
    analysis_type = st.toggle("GOI-based analysis")
    GOI = ""
    if analysis_type:
        GOI = st.text_input("Search by gene (keyword)", value="MASP").upper()
    
    st.subheader("General settings")
    min_len = st.number_input("Min length to analyze:", min_value=0, value=50000)
    window_size = st.number_input("Window size:", min_value=1, value=500)
    step_size = st.number_input("Step size:", min_value=1, value=300)
    smooth_f = st.slider("Smoothing factor:", 0, 200, 100)
    cutoff_value = st.slider("Core/Disruptive Cutoff:", 0.1, 1.0, 0.51)
    num_plots = st.slider("Max plots to display:", 0, 20, 5)
    
    st.subheader("Graph settings")
    ymin_graph = st.number_input("Y min:", value=0.20, step=0.01)
    ymax_graph = st.number_input("Y max:", value=0.80, step=0.01)

# --- APP ---
st.markdown("<h1 style='text-align: center; color: red;'>GCanner</h1>", unsafe_allow_html=True)
st.markdown("<h5 style='text-align: center; color: grey;'>Interactive genome-wide GC analysis</h5>", unsafe_allow_html=True)
st.subheader("", divider="red")

fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna"])

all_figs = []
all_regions = []
filtered_Contigs_to_analyze = []

if fasta_file:
    # Lógica de GFF integrada
    all_GOI = []
    set_contigs_GOI = set()
    if analysis_type:
        st.info("Please upload the GFF file to search for " + GOI)
        GFF_file = st.file_uploader("Upload GFF", type=["gff"])
        if GFF_file:
            gff_text = GFF_file.getvalue().decode("utf-8")
            for line in gff_text.splitlines():
                if not line.startswith("#") and GOI in line.upper():
                    parts = line.strip().split("\t")
                    if len(parts) > 8:
                        all_GOI.append({
                            "Sequence name": parts[0].strip(),
                            "pbi": int(parts[3]), "pbf": int(parts[4]),
                            "Description": parts[8]
                        })
                        set_contigs_GOI.add(parts[0].strip())

    # Procesar FASTA
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = list(SeqIO.parse(fasta_content, "fasta"))
    plot_count = 0
    
    for record in sequences:
        if len(record.seq) < min_len: continue
        if analysis_type and record.id not in set_contigs_GOI: continue
        
        filtered_Contigs_to_analyze.append(record.id)
        gc_data = calculate_gc_content(str(record.seq), window_size, step_size)
        df_gc = pd.DataFrame(gc_data, columns=["Start", "End", "GC_Content"])
        
        # Aplicar Lowess
        frac = f_smoothcal(df_gc, "GC_Content", smooth_f)
        df_gc['Smoothed_GC'] = lowess(df_gc["GC_Content"], df_gc["Start"], frac=frac)[:, 1]
        
        # Regiones
        current_regions = map_regions(df_gc, len(record.seq), cutoff_value)
        for r in current_regions:
            all_regions.append({"Sequence name": record.id, "Start": r[0], "End": r[1], "Region Type": r[2]})

        # Gráfico
        if plot_count < num_plots:
            fig, ax = plt.subplots(figsize=(10, 4))
            sns.lineplot(data=df_gc, x="Start", y="GC_Content", ax=ax, alpha=0.3, label="GC Content")
            sns.lineplot(data=df_gc, x="Start", y="Smoothed_GC", ax=ax, color="red", label="Lowess")
            ax.axhline(y=cutoff_value, color="white", linestyle="-")
            
            for r in current_regions:
                ax.axvspan(r[0], r[1], facecolor=('springgreen' if r[2]=="Core" else 'salmon'), alpha=0.2)
            
            if analysis_type:
                for g in [x for x in all_GOI if x["Sequence name"] == record.id]:
                    ax.axvspan(g["pbi"], g["pbf"], color="black", alpha=0.8, ymin=0, ymax=0.05, label=GOI)
                legend_without_duplicate_labels(ax)
                
            ax.set_ylim(ymin_graph, ymax_graph)
            ax.set_title(f"Sequence: {record.id}")
            st.pyplot(fig)
            all_figs.append((record.id, fig))
            plot_count += 1

    # --- RESULTADOS FINALES ---
    if all_regions:
        regions_df = pd.DataFrame(all_regions)
        regions_df["Region Length"] = regions_df["End"] - regions_df["Start"] + 1
        
        # Gráfico de resumen
        c1, c2 = st.columns(2)
        with c1:
            st.write("### Genome Proportions")
            stats = regions_df.groupby("Region Type")["Region Length"].sum().reset_index()
            fig_stat, ax_stat = plt.subplots()
            sns.barplot(data=stats, x="Region Type", y="Region Length", palette=["green", "red"], ax=ax_stat)
            st.pyplot(fig_stat)
            all_figs.append(("Summary_Stats", fig_stat))

        # Descarga
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "plots.zip")
            with zipfile.ZipFile(zip_path, "w") as z:
                for name, f in all_figs:
                    p = os.path.join(tmpdir, f"{name}.svg")
                    f.savefig(p, format="svg")
                    z.write(p, arcname=f"{name}.svg")
            with open(zip_path, "rb") as f_zip:
                st.download_button("Download SVG Plots", f_zip, "results.zip")

        st.write("### Data Table")
        st.dataframe(regions_df)

st.caption("Disruptomics project - Balouz V, Buscaglia CA. 2026")
