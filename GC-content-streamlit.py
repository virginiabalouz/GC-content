import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from Bio import SeqIO
from io import StringIO
import tempfile
import zipfile
import os

# --- FUNCIONES CORE ---
def f_smoothcal(dataframe, column, f_smooth_points=100):
    points = len(dataframe[column])
    if points == 0: return 1
    f_val = (f_smooth_points / points)
    return min(f_val, 1)

def calculate_gc_content(sequence, window, step):
    gc_content = []
    sequence_length = len(sequence)
    sequence = sequence.upper() # Maneja case-insensitivity de una vez
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
    if df.empty: return []
    df['Core'] = np.where(df["Smoothed_GC"] < cutoff, 1, 0)
    df['Changes'] = df['Core'].diff().fillna(0).abs()
    hits = df[df['Changes'] == 1].index.tolist()
    regions = []
    
    # Estado inicial
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

# --- CONFIGURACIÓN SIDEBAR ---
with st.sidebar:
    st.subheader("Select analysis")
    analysis_type = st.toggle("GOI-based analysis")
    GOI = ""
    if analysis_type:
        GOI = st.text_input("Search by gene of interest (keyword)", value="MASP").upper()
    
    st.subheader("General settings")
    min_len = st.number_input("Minimum sequence length:", value=50000)
    window_size = st.number_input("Window size:", value=500)
    step_size = st.number_input("Step size:", value=300)
    smooth_f = st.slider("Smoothing factor (Lowess):", 0, 200, 100)
    cutoff_value = st.slider("Core/Disruptive Cutoff:", 0.1, 1.0, 0.51)
    num_plots = st.slider("Number of plots to display:", 0, 20, 1)
    
    st.subheader("Graph settings")
    ymin_graph = st.number_input("Y min:", value=0.20, step=0.01)
    ymax_graph = st.number_input("Y max:", value=0.80, step=0.01)

# --- CUERPO DE LA APP ---
st.markdown("<h1 style='text-align: center; color: red;'>GCanner</h1>", unsafe_allow_html=True)
st.markdown("<h5 style='text-align: center; color: grey;'>Interactive tool for genome-wide GC analysis</h5>", unsafe_allow_html=True)
st.subheader("", divider="red")

fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna"])
st.subheader("", divider="red")

# Inicialización de estructuras de datos
all_figs = []
all_regions = []
all_GOI = []
set_contigs_GOI = set()
filtered_Contigs_to_analyze = []

if fasta_file:
    # 1. Carga de GFF si el análisis está activado
    if analysis_type:
        GFF_file = st.file_uploader("Upload GFF file", type=["gff"])
        if GFF_file:
            gff_content = StringIO(GFF_file.getvalue().decode("utf-8"))
            for line in gff_content:
                if not line.startswith("#") and GOI in line.upper():
                    parts = line.strip().split("\t")
                    if len(parts) > 8:
                        contig = parts[0].strip()
                        all_GOI.append({
                            "Sequence name": contig, "pbi": int(parts[3]), 
                            "pbf": int(parts[4]), "Description": parts[8]
                        })
                        set_contigs_GOI.add(contig)
            if not all_GOI:
                st.warning(f"No genes found with keyword '{GOI}'")

    # 2. Procesamiento de FASTA
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = list(SeqIO.parse(fasta_content, "fasta"))
    plot_count = 0
    
    st.subheader("Results")
    
    for record in sequences:
        if len(record.seq) < min_len: continue
        if analysis_type and record.id not in set_contigs_GOI: continue
        
        filtered_Contigs_to_analyze.append(record.id)
        gc_data = calculate_gc_content(str(record.seq), window_size, step_size)
        df_gc = pd.DataFrame(gc_data, columns=["Start", "End", "GC_Content"])
        
        # Suavizado Lowess (Importado directamente)
        sf = f_smoothcal(df_gc, "GC_Content", smooth_f)
        df_gc['Smoothed_GC'] = lowess(df_gc["GC_Content"], df_gc["Start"], frac=sf)[:, 1]
        
        # Mapeo de regiones
        regions = map_regions(df_gc, len(record.seq), cutoff_value)
        for r in regions:
            all_regions.append({"Sequence name": record.id, "Start": r[0], "End": r[1], "Region Type": r[2]})

        # Renderizado de gráficos (hasta el límite de num_plots)
        if plot_count < num_plots:
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=df_gc, x="Start", y="GC_Content", ax=ax, alpha=0.4, label="Raw GC")
            sns.lineplot(data=df_gc, x="Start", y="Smoothed_GC", color="red", ax=ax, label="Lowess")
            ax.axhline(y=cutoff_value, color="black", linestyle="--", alpha=0.3)
            
            for r in regions:
                color = 'springgreen' if r[2] == "Core" else 'salmon'
                ax.axvspan(r[0], r[1], facecolor=color, alpha=0.2)
            
            if analysis_type:
                for gene in [g for g in all_GOI if g["Sequence name"] == record.id]:
                    ax.axvspan(gene["pbi"], gene["pbf"], color="black", alpha=0.7, ymin=0, ymax=0.05, label=GOI)
                legend_without_duplicate_labels(ax)
            
            ax.set_ylim(ymin_graph, ymax_graph)
            ax.set_title(f"Sequence: {record.id}")
            st.pyplot(fig)
            all_figs.append((record.id, fig))
            plot_count += 1

    # 3. Estadísticas y Resúmenes Finales
    if all_regions:
        regions_df = pd.DataFrame(all_regions)
        regions_df["Region Length"] = regions_df["End"] - regions_df["Start"] + 1
        
        # Gráficos de Proporción
        dis_total = regions_df[regions_df["Region Type"] == "Disruptive"]["Region Length"].sum()
        core_total = regions_df[regions_df["Region Type"] == "Core"]["Region Length"].sum()
        total_bp = regions_df["Region Length"].sum()
        
        summary_data = pd.DataFrame({
            'Compartment': ['Core', 'Disruptive'],
            '%': [core_total/total_bp*100, dis_total/total_bp*100],
            'bases': [core_total, dis_total]
        })

        col1, col2 = st.columns(2)
        with col1:
            fig_sum, ax_sum = plt.subplots(figsize=(4, 5))
            sns.barplot(data=summary_data, x='Compartment', y='%', palette=["green", "red"], alpha=0.3, ax=ax_sum)
            for i, p in enumerate(ax_sum.patches):
                val = summary_data.iloc[i]
                ax_sum.annotate(f'{round(val["bases"]/1000,1)}Kb\n{round(val["%"],1)}%', 
                                (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom')
            ax_sum.set_ylim(0, 115)
            st.pyplot(fig_sum)
            all_figs.append(("Proportions_Summary", fig_sum))

        # Lógica GOI y tabla final
        if analysis_type and all_GOI:
            # Mapeo de genes a regiones para la tabla final
            for gene in all_GOI:
                match = regions_df[(regions_df["Sequence name"] == gene["Sequence name"]) & 
                                   (gene["pbi"] >= regions_df["Start"]) & 
                                   (gene["pbi"] <= regions_df["End"])]
                gene["Region Type"] = match.iloc[0]["Region Type"] if not match.empty else "Edge/Unknown"

            filtered_GOI_df = pd.DataFrame(all_GOI)
            with col2:
                if not filtered_GOI_df.empty:
                    fig_goi, ax_goi = plt.subplots(figsize=(4, 5))
                    counts = filtered_GOI_df["Region Type"].value_counts()
                    counts.plot(kind='bar', color=['green', 'red', 'grey'], alpha=0.3, ax=ax_goi)
                    ax_goi.set_title(f"Distribution of {GOI}")
                    st.pyplot(fig_goi)
                    all_figs.append((f"{GOI}_Stats", fig_goi))

        # Descarga de resultados (ZIP de SVGs)
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "plots.zip")
            with zipfile.ZipFile(zip_path, "w") as z:
                for name, f in all_figs:
                    p = os.path.join(tmpdir, f"{name}.svg")
                    f.savefig(p, format="svg")
                    z.write(p, arcname=f"{name}.svg")
            with open(zip_path, "rb") as fz:
                st.download_button("Download all figures (SVG)", fz, "GCanner_results.zip")

        st.write("### Regions Table")
        st.dataframe(regions_df)
        if analysis_type and all_GOI:
            st.write(f"### {GOI} Mapping Table")
            st.dataframe(pd.DataFrame(all_GOI))

st.caption("Disruptomics project. Balouz V, Buscaglia CA. GCanner, a Genome-Wide GC Composition Tool for the Unbiased Assessment of Trypanosoma cruzi Genomic Compartments. Methods Mol Biol. 2026;2982:31-45. doi: 10.1007/978-1-0716-4848-3_3. PMID: 41182609.")
