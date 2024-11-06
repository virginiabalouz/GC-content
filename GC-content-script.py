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

st.title("GC Content and Compartment Analysis")
st.write("Upload a FASTA file to analyze the GC content and predict core/disruptive compartments.")

fasta_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])
min_len = st.number_input("Minimum length to analyze:", min_value=0, max_value=1000000, value=50000)  # Control para la cantidad de gráficos
window_size = st.number_input("Window size:", min_value=0, max_value=1000, value=500)
step_size = st.number_input("Step size:", min_value=0, max_value=500, value=300)
smooth_f = st.slider("Smoothing factor (Lowess):", 0, 200, 100)
cutoff_value = st.slider("Core/Disruptive Cutoff:", 0.1, 1.0, 0.51)
num_plots = st.slider("Number of plots to display:", 1, 20, 10)  # Control para la cantidad de gráficos

def f_smoothcal(dataframe, column, f_smooth_points=100):
    points = len(dataframe[column])
    f_val = (f_smooth_points / points)
    return min(f_val, 1)

def calculate_gc_content(sequence, window, step):
    gc_content = []
    sequence_length = len(sequence)
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
    df['Core'] = np.where(df["Smoothed_GC"] < cutoff, 1, 0)
    df['Changes'] = df['Core'].diff().fillna(0).abs()
    hits = df[df['Changes'] == 1].index
    regions = []
    cutoff_status = "disruptive" if df["Smoothed_GC"].iloc[0] > cutoff else "core"

    for i, change in enumerate(hits):
        start = df['Start'].iloc[hits[i-1]] + 1 if i > 0 else 1
        end = df['End'].iloc[change]
        regions.append((start, end, cutoff_status))
        cutoff_status = "core" if cutoff_status == "disruptive" else "disruptive"

    final_start = regions[-1][1] + 1 if regions else 1
    regions.append((final_start, fasta_len, cutoff_status))
    return regions

# Almacena gráficos y regiones para exportación
all_figs = []
all_regions = []

# Análisis cuando se sube el archivo
if fasta_file is not None:
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = SeqIO.parse(fasta_content, "fasta")
    plot_count = 0

    for record in sequences:
        if len(str(record.seq))>=min_len:

            gc_content = calculate_gc_content(str(record.seq), window_size, step_size)
            df_gc_content = pd.DataFrame(gc_content, columns=["Start", "End", "GC_Content"])

            smooth_factor = f_smoothcal(df_gc_content, "GC_Content", smooth_f)
            df_gc_content['Smoothed_GC'] = sm.nonparametric.lowess(df_gc_content["GC_Content"], df_gc_content["Start"], frac=smooth_factor)[:, 1]

            regions = map_regions(df_gc_content, len(record.seq), cutoff_value)
            for start, end, region_type in regions:
                all_regions.append({"Contig": record.id, "Start": start, "End": end, "Region_Type": region_type})

            df_gc_content["Core/Disruptive"] = np.where(df_gc_content["Smoothed_GC"] < cutoff_value, "Core", "Disruptive")
            # Visualización de datos
            if plot_count >= num_plots:
                continue
            # st.write(f"Analyzing Sequence: {record.id}")
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=df_gc_content, x="Start", y="GC_Content", ax=ax, label="GC content").set_title(f"Sequence name: {record.id}, window size:{window_size}, step size:{step_size}, smoothing points:{smooth_f}, cutoff:{cutoff_value}")
            sns.lineplot(data=df_gc_content, x="Start", y="Smoothed_GC", color="red", ax=ax, label="Lowess smoothing")
            ax.axhline(y=cutoff_value, xmin=0, xmax=max(df_gc_content["End"]), color="white")
            ax.fill_between(df_gc_content["Start"], 0, 1, where=df_gc_content["Core/Disruptive"] == "Core", color="green", alpha=0.2, transform=ax.get_xaxis_transform())
            ax.fill_between(df_gc_content["Start"], 0, 1, where=df_gc_content["Core/Disruptive"] == "Disruptive", color="red", alpha=0.2, transform=ax.get_xaxis_transform())
            ax.set_ylabel("GC Content")
            ax.set_ylabel("GC content and lowess smoothing")
            ax.set_xlabel("Sequence position")
            ax.legend()
            st.pyplot(fig)
            all_figs.append((record.id, fig))
            plot_count += 1

    regions_df = pd.DataFrame(all_regions)
    regions_df["Len-region"]=regions_df["End"]-regions_df["Start"]+1
    dis=sum(regions_df[regions_df["Region_Type"]=="disruptive"]["Len-region"])
    core=sum(regions_df[regions_df["Region_Type"]=="core"]["Len-region"])
    tot= sum(regions_df["Len-region"])
    df = pd.DataFrame({'Compartment': ['Core', 'Disruptive'],'%': [core/tot*100, dis/tot*100]})
    st.write("Summary statistics")
    # st.write(f"Disruptive bases: {dis}\nCore bases: {core}\nTotal bases: {tot} ")


    fig, ax = plt.subplots(figsize=(5, 3))
    ax = sns.barplot(df, x='Compartment', y='%',palette=["green", "red"], alpha=0.2)
    for i in ax.containers:
        ax.bar_label(i,)
     # sns.barplot(df, x="Compartment", y="%", palette=["green", "red"], alpha=0.2)
    st.pyplot(fig)


    # Exportar gráficos como SVG
    with tempfile.TemporaryDirectory() as tmpdirname:
        svg_files = []
        for figure in all_figs:
            svg_path = f"{tmpdirname}/plot_{figure[0]}.svg"
            figure[1].savefig(svg_path, format="svg")
            svg_files.append(svg_path)

        zip_path = f"{tmpdirname}/plots_svg.zip"
        with zipfile.ZipFile(zip_path, "w") as zipf:
            for svg_file in svg_files:
                zipf.write(svg_file, arcname=os.path.basename(svg_file))

        with open(zip_path, "rb") as f:
            st.download_button("Download figures", data=f, file_name="plots_svg.zip", mime="application/zip")

    # Exportar datos de regiones como CSV
        csv_path = f"{tmpdirname}/mapped_regions.csv"
        regions_df.to_csv(csv_path, index=False)

        with open(csv_path, "rb") as f:
            st.download_button("Download regions (CSV)", data=f, file_name="mapped_regions.csv", mime="text/csv")
