import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from Bio import SeqIO
from io import StringIO
import tempfile
import zipfile
import os

# --- FUNCIONES ---
def f_smoothcal(dataframe, column, f_smooth_points=100):
    points = len(dataframe[column])
    if points == 0: return 1
    f_val = (f_smooth_points / points)
    return min(f_val, 1)

def calculate_gc_content(sequence, window, step):
    gc_content = []
    sequence_length = len(sequence)
    # Optimizamos para manejar mayúsculas y minúsculas de una vez
    for i in range(0, sequence_length, step):
        start = i
        end = min(i + window, sequence_length)
        window_sequence = sequence[start:end].upper()
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
    
    # Estado inicial basado en el primer punto suavizado
    cutoff_status = "Disruptive" if df["Smoothed_GC"].iloc[0] >= cutoff else "Core"

    last_pos = 1
    for change in hits:
        end = df['Start'].iloc[change]
        regions.append((last_pos, end, cutoff_status))
        last_pos = end + 1
        cutoff_status = "Core" if cutoff_status == "Disruptive" else "Disruptive"
    
    # Añadir la última región hasta el final del FASTA
    regions.append((last_pos, fasta_len, cutoff_status))
    return regions

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    if unique:
        ax.legend(*zip(*unique), loc='upper right')

# --- CONFIGURACIÓN DE PÁGINA ---
st.set_page_config(page_title="GCanner", layout="wide")

# --- SIDEBAR ---
with st.sidebar:
    st.subheader("Select analysis")
    analysis_type = st.toggle("GOI-based analysis")
    GOI = ""
    if analysis_type:
        GOI = st.text_input("Search by gene of interest (keyword)", value="MASP", max_chars=100).upper()
    
    st.subheader("General settings")
    min_len = st.number_input("Minimum sequence length to analyze:", min_value=0, value=50000)
    window_size = st.number_input("Window size (0-1000):", min_value=1, max_value=1000, value=500)
    step_size = st.number_input("Step size (0-500):", min_value=1, max_value=500, value=300)
    smooth_f = st.slider("Smoothing factor (Lowess):", 0, 200, 100)
    cutoff_value = st.slider("Core/Disruptive Cutoff:", 0.1, 1.0, 0.51)
    num_plots = st.slider("Number of GC content plots to display:", 0, 20, 1)
    
    st.subheader("Graph settings")
    ymin_graph = st.number_input("Y min:", min_value=0.00, max_value=0.50, value=0.20, step=0.01)
    ymax_graph = st.number_input("Y max:", min_value=0.50, max_value=1.00, value=0.80, step=0.01)

# --- CUERPO PRINCIPAL ---
st.markdown("<h1 style='text-align: center; color: red; line-height: 0.1'>GCanner</h1><h5 style='text-align: center; color: grey;line-height:1 '> An interactive tool for genome-wide GC analysis and isochore identification</h5>", unsafe_allow_html=True)
st.subheader("", divider="red")
st.markdown("#### Upload FASTA file")
fasta_file = st.file_uploader("", type=["fasta", "fa", "fna"])
st.subheader("", divider="red")

# Inicialización de variables para evitar NameErrors
all_figs = []
all_regions = []
filtered_Contigs_to_analyze = []
regions_df = pd.DataFrame()
filtered_DF_GOI = pd.DataFrame()
GFF_file = None

# --- LÓGICA DE CARGA GFF (Si aplica) ---
all_GOI = []
set_contigs_GOI = set()
if fasta_file and analysis_type:
    st.markdown("#### Upload GFF file")
    GFF_file = st.file_uploader("", type=["gff"])
    st.subheader("", divider="red")
    if GFF_file:
        GFF_content = StringIO(GFF_file.getvalue().decode("utf-8"))
        for line in GFF_content:
            if not line.startswith("#") and GOI in line.upper():
                parts = line.strip().split("\t")
                if len(parts) > 8:
                    contig = parts[0].strip()
                    all_GOI.append({
                        "Sequence name": contig, 
                        "pbi": int(parts[3]), 
                        "pbf": int(parts[4]), 
                        "Description": parts[8]
                    })
                    set_contigs_GOI.add(contig)
        
        if not all_GOI:
            st.error(f"No genes found with keyword: {GOI}")
            st.stop()

# --- PROCESAMIENTO ---
if fasta_file:
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = list(SeqIO.parse(fasta_content, "fasta"))
    plot_count = 0
    
    status_text = "Results: Core and Disruptive regions" if not analysis_type else "Results: GOI-based search"
    st.subheader(status_text)

    for record in sequences:
        seq_len = len(record.seq)
        if seq_len < min_len:
            continue
        
        # Si es GOI-based, filtrar por contigs que tengan el gen
        if analysis_type and record.id not in set_contigs_GOI:
            continue

        filtered_Contigs_to_analyze.append(record.id)
        
        # Calcular GC
        gc_data = calculate_gc_content(str(record.seq), window_size, step_size)
        df_gc = pd.DataFrame(gc_data, columns=["Start", "End", "GC_Content"])
        
        # Suavizado Lowess
        frac = f_smoothcal(df_gc, "GC_Content", smooth_f)
        df_gc['Smoothed_GC'] = sm.nonparametric.lowess(df_gc["GC_Content"], df_gc["Start"], frac=frac)[:, 1]
        
        # Mapear Regiones
        current_regions = map_regions(df_gc, seq_len, cutoff_value)
        for r_start, r_end, r_type in current_regions:
            all_regions.append({
                "Sequence name": record.id, 
                "Start": r_start, 
                "End": r_end, 
                "Region Type": r_type
            })

        # Gráfico
        if plot_count < num_plots:
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=df_gc, x="Start", y="GC_Content", ax=ax, label="GC-Content", alpha=0.4)
            sns.lineplot(data=df_gc, x="Start", y="Smoothed_GC", color="red", ax=ax, label="Lowess smoothing")
            ax.axhline(y=cutoff_value, color="black", linestyle="--", alpha=0.5)
            
            for r in current_regions:
                color = 'springgreen' if r[2] == "Core" else 'salmon'
                ax.axvspan(r[0], r[1], facecolor=color, alpha=0.2)
            
            # Dibujar marcas de GOI si aplica
            if analysis_type:
                contig_genes = [g for g in all_GOI if g["Sequence name"] == record.id]
                for gene in contig_genes:
                    ax.axvspan(gene["pbi"], gene["pbf"], color="black", alpha=0.6, ymin=0, ymax=0.05, label=GOI)
                legend_without_duplicate_labels(ax)

            ax.set_ylim(ymin_graph, ymax_graph)
            ax.set_title(f"Sequence: {record.id} | Window: {window_size} | Cutoff: {cutoff_value}")
            st.pyplot(fig)
            all_figs.append((record.id, fig))
            plot_count += 1

    # --- RESUMEN Y EXPORTACIÓN ---
    if not all_regions:
        st.warning(f"No sequences longer than {min_len} bp found.")
    else:
        regions_df = pd.DataFrame(all_regions)
        regions_df["Region Length"] = regions_df["End"] - regions_df["Start"] + 1
        
        # Estadísticas Globales
        dis_len = regions_df[regions_df["Region Type"] == "Disruptive"]["Region Length"].sum()
        core_len = regions_df[regions_df["Region Type"] == "Core"]["Region Length"].sum()
        total_len = regions_df["Region Length"].sum()
        
        summary_df = pd.DataFrame({
            'Compartment': ['Core', 'Disruptive'],
            '%': [core_len/total_len*100, dis_len/total_len*100],
            'bases': [core_len, dis_len]
        })

        col6, col7 = st.columns(2)
        with col6:
            fig_sum, ax_sum = plt.subplots(figsize=(4, 4))
            sns.barplot(data=summary_df, x='Compartment', y='%', palette=["green", "red"], alpha=0.3, ax=ax_sum)
            for i, p in enumerate(ax_sum.patches):
                txt = f'{round(summary_df.iloc[i]["bases"]/1000,1)}Kb\n{round(summary_df.iloc[i]["%"],1)}%'
                ax_sum.annotate(txt, (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom')
            ax_sum.set_ylim(0, 110)
            st.pyplot(fig_sum)
            all_figs.append(("Summary_Proportions", fig_sum))

        # Lógica adicional para GOI
        if analysis_type and all_GOI:
            # Mapear cada gen a su tipo de región
            for gene in all_GOI:
                # Buscar en qué región cae el inicio del gen
                match = regions_df[(regions_df["Sequence name"] == gene["Sequence name"]) & 
                                   (gene["pbi"] >= regions_df["Start"]) & 
                                   (gene["pbi"] <= regions_df["End"])]
                if not match.empty:
                    gene["Region Type"] = match.iloc[0]["Region Type"]
                    gene["Contig-len"] = regions_df[regions_df["Sequence name"] == gene["Sequence name"]]["End"].max()
                else:
                    gene["Region Type"] = "Unknown"
            
            df_GOI_final = pd.DataFrame(all_GOI)
            filtered_DF_GOI = df_GOI_final[df_GOI_final.get("Contig-len", 0) > min_len]
            
            with col7:
                if not filtered_DF_GOI.empty:
                    counts = filtered_DF_GOI["Region Type"].value_counts(normalize=True) * 100
                    fig_goi, ax_goi = plt.subplots(figsize=(4, 4))
                    counts.plot(kind='bar', ax=ax_goi, color=['green', 'red', 'grey'], alpha=0.3)
                    ax_goi.set_title(f"Location of {GOI}")
                    ax_goi.set_ylabel("% Genes")
                    st.pyplot(fig_goi)
                    all_figs.append((f"{GOI}_Distribution", fig_goi))

        # Botón de Descarga
        if all_figs:
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = os.path.join(tmpdir, "plots_svg.zip")
                with zipfile.ZipFile(zip_path, "w") as zipf:
                    for name, f in all_figs:
                        f_path = os.path.join(tmpdir, f"{name}.svg")
                        f.savefig(f_path, format="svg")
                        zipf.write(f_path, arcname=f"{name}.svg")
                
                with open(zip_path, "rb") as f_zip:
                    st.download_button("Download figures as SVG", data=f_zip, file_name="plots_svg.zip", mime="application/zip")

        st.write("### Regions Detail")
        st.dataframe(regions_df)
        if analysis_type and not filtered_DF_GOI.empty:
            st.write(f"### {GOI} Mapping Detail")
            st.dataframe(filtered_DF_GOI)

st.caption("This program is part of the Disruptomics project. Please cite Balouz V, Buscaglia CA. Methods Mol Biol. 2026;2982:31-45. doi: 10.1007/978-1-0716-4848-3_3. PMID: 41182609.")
