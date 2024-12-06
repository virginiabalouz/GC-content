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

#functions
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
    cutoff_status = "Disruptive" if df["Smoothed_GC"].iloc[0] >= cutoff else "Core"

    for i, change in enumerate(hits):
        start = df['Start'].iloc[hits[i-1]] + 1 if i > 0 else 1
        end = df['Start'].iloc[change]
        regions.append((start, end, cutoff_status))
        cutoff_status = "Core" if cutoff_status == "Disruptive" else "Disruptive"
    final_start = regions[-1][1] + 1 if regions else 1
    regions.append((final_start, fasta_len, cutoff_status))
    return regions

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc= 'upper right')



### SIDEBAR
with st.sidebar:
    st.subheader("Select analysis")
    analysis_type = st.toggle("GOI-based analysis")
    if analysis_type:
        GOI= st.text_input("Search by gene of interest (keyword)", value="MASP", max_chars=100,label_visibility="visible")
        GOI= GOI.upper()
    st.subheader("General settings")
    min_len = st.number_input("Minimum sequence length to analyze:", min_value=0, max_value=1000000, value=50000)  # Control para la cantidad de gráficos
    window_size = st.number_input("Window size (0-1000):", min_value=0, max_value=1000, value=500)
    step_size = st.number_input("Step size (0-500):", min_value=0, max_value=500, value=300)
    smooth_f = st.slider("Smoothing factor (Lowess):", 0, 200, 100)# Control para la cantidad de gráficos
    cutoff_value = st.slider("Core/Disruptive Cutoff:", 0.1, 1.0, 0.51)
    num_plots = st.slider("Number of GC content plots to display:", 0, 20, 1)  # Control para la cantidad de gráficos
    st.subheader("Graph settings")
    ymin_graph = st.number_input("Y min:", min_value=0.00, step=0.01, max_value=0.50, value=0.20)
    ymax_graph = st.number_input("Y max:", min_value=0.50, step=0.01, max_value=1.00, value=0.80)

###BODY
st.markdown("<h1 style='text-align: center; color: red; line-height: 0.1'>GCanner</h1><h5 style='text-align: center; color: grey;line-height:1 '> An interactive tool for genome-wide GC analysis and isochore identification</h5>", unsafe_allow_html=True)
st.subheader("", divider="red")
st.markdown("<h4 style='text-align: left; color: black;line-height:0.3; margin-top: 0px; '>Upload FASTA file</h5>", unsafe_allow_html=True)
fasta_file = st.file_uploader("", type=["fasta", "fa", "fna"])
st.subheader("", divider="red")

######################FASTA FILE ONLY
# Almacena gráficos y regiones para exportación
all_figs = []
all_regions = []

# Análisis cuando se sube el archivo
if fasta_file is not None and not analysis_type:
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = SeqIO.parse(fasta_content, "fasta")
    plot_count = 0
    st.subheader("Results: Core and Disruptive regions", divider=False)
    filtered_Contigs_to_analyze=[]
    for record in sequences:
        if len(str(record.seq))>=min_len:
            filtered_Contigs_to_analyze.append(record.id)
            gc_content = calculate_gc_content(str(record.seq), window_size, step_size)
            df_gc_content = pd.DataFrame(gc_content, columns=["Start", "End", "GC_Content"])

            smooth_factor = f_smoothcal(df_gc_content, "GC_Content", smooth_f)
            df_gc_content['Smoothed_GC'] = sm.nonparametric.lowess(df_gc_content["GC_Content"], df_gc_content["Start"], frac=smooth_factor)[:, 1]

            regions = map_regions(df_gc_content, len(record.seq), cutoff_value)
            for start, end, region_type in regions:
                all_regions.append({"Sequence name"  : record.id, "Start": start, "End": end, "Region Type": region_type})

            df_gc_content["Core/Disruptive"] = np.where(df_gc_content["Smoothed_GC"] < cutoff_value, "Core", "Disruptive")
            # Visualización de datos
            regions_df = pd.DataFrame(all_regions)
            if plot_count >= num_plots:
                continue
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=df_gc_content, x="Start", y="GC_Content", ax=ax, label="GC-Content",alpha=0.5).set_title(f"Sequence name: {record.id}, window size:{window_size}, step size:{step_size}, smoothing points:{smooth_f}, cutoff:{cutoff_value}")
            sns.lineplot(data=df_gc_content, x="Start", y="Smoothed_GC", color="red", ax=ax, label="Lowess smoothing")
            ax.axhline(y=cutoff_value, xmin=0, xmax=max(df_gc_content["End"]), color="white") #curoff line
            for t in regions:
                if t[-1] == "Disruptive":
                    plt.axvspan(t[0], t[1], facecolor='salmon', alpha=0.2)
                elif t[-1] == "Core":
                    plt.axvspan(t[0], t[1], facecolor='springgreen', alpha=0.2)
            ax.set_ylim(ymin_graph,ymax_graph)
            ax.set_ylabel("GC Content")
            ax.set_ylabel("GC content and Lowess smoothing")
            ax.set_xlabel("Sequence position")
            ax.margins(x=0)
            st.pyplot(fig)
            all_figs.append((record.id, fig))
            plot_count += 1

    if len(all_regions)== 0:
        st.markdown(f" **:red-background[The sequence length is shorter than {min_len} bp]**")
        st.stop()
    regions_df["Region Length"]=regions_df["End"]-regions_df["Start"]+1
    dis = sum(regions_df[regions_df["Region Type"]=="Disruptive"]["Region Length"])
    Core= sum(regions_df[regions_df["Region Type"]=="Core"]["Region Length"])
    tot = sum(regions_df["Region Length"])
    df  = pd.DataFrame({'Compartment': ['Core', 'Disruptive'],'%': [Core/tot*100, dis/tot*100], "bases":[Core,dis]})
    # st.write("Summary statistics")
    #st.write(f"Disruptive bases: {dis}\nCore bases: {Core}\nTotal bases: {tot} ")
    col6,col7 = st.columns(2)
    #plots Summary Core Disruptive
    with col6:
        fig, ax = plt.subplots(figsize=(4, 4))
        sns.barplot(df, x='Compartment', y=round(df['%'],1), hue= 'Compartment',palette=["green", "red"], alpha=0.2, ax=ax)
        for i, p in enumerate(ax.patches):
            h, w, x = p.get_height(), p.get_width(), p.get_x()
            xy = (x + w / 2.,max(df["%"])*1.2)
            text1 = f' {round(df.iloc[i]["bases"]/1000,2)}Kb\n {round(df.iloc[i]["%"],2)}%'
            ax.annotate(text=text1, xy=xy, ha='center', va='top')
        ax.set_ylim(0,110)
        if len(filtered_Contigs_to_analyze)== 1:
            ax.set_title(f"Compartment proportions in {filtered_Contigs_to_analyze[0]} ")
        else:
            ax.set_title(f"Compartment proportions in all sequences\n>{min_len} from the FASTA file")
        all_figs.append((f"Compartment_porportions_in_>{min_len}", fig))
        st.pyplot(fig)
    # Export/Download graphs as SVG
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
            st.download_button("Download figures as SVG", data=f, file_name="plots_svg.zip", mime="application/zip")
    st.dataframe(data=regions_df)
    if analysis_type:
        st.dataframe(data=filtered_DF_GOI)


################### FASTA FILE AND GFF
if fasta_file is not None and analysis_type:
    st.markdown("<h4 style='text-align: left; color: black;line-height:0.1; margin-top: 0px; '>Upload GFF file</h5>", unsafe_allow_html=True)
    GFF_file = st.file_uploader("", type=["gff"])
    st.subheader("", divider="red")
    if GFF_file is not None:
        GFF_file_content = StringIO(GFF_file.getvalue().decode("utf-8"))
        all_GOI = []
        contigs_GOI=[]
        for line in GFF_file_content:
            if GOI in line.upper():
                line=line.strip().split("\t")
                if len(line)>8:
                    Contig=line[0].strip()
                    pbi=int(line[3])
                    pbf=int(line[4])
                    desc=line[8]
                    contigs_GOI.append(Contig)
                    all_GOI.append({"Sequence name"  : Contig, "pbi": pbi, "pbf": pbf, "Description": desc})
        df_GOI = pd.DataFrame(data=all_GOI, columns=["Sequence name"  , "pbi", "pbf", "Description"])
        set_contigs_GOI= set(contigs_GOI)
        c_to_map= set(df_GOI["Sequence name"])
        if len(c_to_map) == 1:
            map=1
        elif len(c_to_map) == 0:
            st.markdown(f" **:red-background[No ***{GOI}*** were found]**")
            st.stop()
        elif GOI == "":
            st.markdown(f" **:red-background[Write a keyword for the GOI-based search]**")
            st.stop()

# Almacena gráficos y regiones para exportación
all_figs = []
all_regions = []

# Análisis cuando se sube el archivo

if fasta_file is not None and analysis_type and GFF_file is not None and len(GOI)>=1:
    fasta_content = StringIO(fasta_file.getvalue().decode("utf-8"))
    sequences = SeqIO.parse(fasta_content, "fasta")
    st.subheader("Results: GOI-based search")
    filtered_Contigs_to_analyze=[]
    plot_count = 0
    for record in sequences:
        if len(str(record.seq))>=min_len and record.id in set_contigs_GOI:
            filtered_Contigs_to_analyze.append(record.id)
            gc_content = calculate_gc_content(str(record.seq), window_size, step_size)
            df_gc_content = pd.DataFrame(gc_content, columns=["Start", "End", "GC_Content"])

            smooth_factor = f_smoothcal(df_gc_content, "GC_Content", smooth_f)
            df_gc_content['Smoothed_GC'] = sm.nonparametric.lowess(df_gc_content["GC_Content"], df_gc_content["Start"], frac=smooth_factor)[:, 1]

            regions = map_regions(df_gc_content, len(record.seq), cutoff_value)
            for start, end, region_type in regions:
                all_regions.append({"Sequence name"  : record.id, "Start": start, "End": end, "Region Type": region_type})

            df_gc_content["Core/Disruptive"] = np.where(df_gc_content["Smoothed_GC"] < cutoff_value, "Core", "Disruptive")

            regions_df = pd.DataFrame(all_regions)
            for hit in all_GOI: #mapeo de GFF en regiones
                f_df = regions_df[regions_df["Sequence name"  ] == hit["Sequence name"  ] ] #me quedo con las regiones de ese Contig solamente
                d_Contig = f_df.to_dict('records')
                for region in d_Contig : # list of dict
                    if int(hit["pbi"]) >= region["Start"] and  int(hit["pbf"]) <= region["End"]:
                        hit["Region Type"]= region["Region Type"]
                        hit["Pbi-Region"]= region["Start"]
                        hit["Pbf-Region"]= region["End"]
                        hit["Contig-len"]=  max(regions_df[regions_df["Sequence name"  ] == hit["Sequence name"  ] ]["End"]  )
                    elif int(hit["pbi"]) >= region["Start"] and not int(hit["pbf"]) <= region["End"]:
                        hit["Region Type"]= "Edge"
                        hit["Pbi-Region"]= region["Start"]
                        hit["Pbf-Region"]= region["End"]
                        hit["Contig-len"]= max(regions_df[regions_df["Sequence name"  ] == hit["Sequence name"  ] ]["End"]  )
            if plot_count >= num_plots :# Visualización de datos
                continue

            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=df_gc_content, x="Start", y="GC_Content", ax=ax, label="GC content",alpha=0.5).set_title(f"Sequence name: {record.id}, window size:{window_size}, step size:{step_size}, smoothing points:{smooth_f}, cutoff:{cutoff_value}")
            sns.lineplot(data=df_gc_content, x="Start", y="Smoothed_GC", color="red", ax=ax, label="Lowess smoothing")
            ax.axhline(y=cutoff_value, xmin=0, xmax=max(df_gc_content["End"]), color="white") #curoff line
            for t in regions:
                if t[-1] == "Disruptive":
                    plt.axvspan(t[0], t[1], facecolor='salmon', alpha=0.2)
                elif t[-1] == "Core":
                    plt.axvspan(t[0], t[1], facecolor='springgreen', alpha=0.2)
            for pb in df_GOI[df_GOI["Sequence name"  ]==record.id]["pbi"]:
                ax.axvspan(xmin=pb, ymax=0.05,xmax=df_GOI[ (df_GOI["pbi"]==pb) & (df_GOI["Sequence name"  ]==record.id) ]["pbf"].iloc[0], color="black", alpha=0.5,label=GOI)
                legend_without_duplicate_labels(ax)
            ax.set_ylim(ymin_graph,ymax_graph)
            ax.set_ylabel("GC Content")
            ax.set_ylabel("GC content and Lowess smoothing")
            ax.set_xlabel("Sequence position")
            ax.margins(x=0)
            st.pyplot(fig)
            plot_count+=1
            all_figs.append((record.id, fig))


if fasta_file is not None and analysis_type and GFF_file is not None and len(GOI)>=1:
    try: #handles empty dataframes due to the lack of GOI within it or non annotated contigs
        PRU= regions_df["End"]
    except:
        regions_df= None
        st.markdown(f" **:red-background[Inconsistencies were found]**")

if fasta_file is not None and analysis_type and GFF_file is not None and regions_df is not None and len(GOI)>=1:

    regions_df["Region Length"]=regions_df["End"]-regions_df["Start"]+1
    dis = sum(regions_df[regions_df["Region Type"]=="Disruptive"]["Region Length"])
    Core= sum(regions_df[regions_df["Region Type"]=="Core"]["Region Length"])
    tot = sum(regions_df["Region Length"])
    df  = pd.DataFrame({'Compartment': ['Core', 'Disruptive'],'%': [Core/tot*100, dis/tot*100], "bases":[Core,dis]})

    col6,col7 = st.columns(2)
    #plots Summary Core Disruptive
    with col6:
        fig, ax = plt.subplots(figsize=(4, 5))
        sns.barplot(df, x='Compartment', y=round(df['%'],1), hue= 'Compartment',palette=["green", "red"], alpha=0.2, ax=ax)
        for i, p in enumerate(ax.patches):
            h, w, x = p.get_height(), p.get_width(), p.get_x()
            xy = (x + w / 2.,max(df["%"])*1.2)
            text1 = f' {round(df.iloc[i]["bases"]/1000,2)}Kb\n {round(df.iloc[i]["%"],2)}%'
            ax.annotate(text=text1, xy=xy, ha='center', va='top')
        ax.set_ylim(0,110)
        if len(filtered_Contigs_to_analyze)== 1:
            ax.set_title(f"Compartment proportions in {filtered_Contigs_to_analyze[0]}\n ")
        else:
            ax.set_title(f"Compartments proportions in all sequences\n>{min_len} from the FASTA file")
        all_figs.append((f"Compartments proportions in >{min_len} sequences", fig))
        st.pyplot(fig)

    df_GOI = pd.DataFrame(data=all_GOI , columns=["Sequence name"  , "pbi", "pbf", "Description", "Region Type", "Pbi-Region", "Pbf-Region" , "Contig-len"])
    filtered_DF_GOI= df_GOI[df_GOI['Contig-len']>min_len][["Sequence name"  , "pbi", "pbf", "Description", "Region Type"]]
    df_trans=filtered_DF_GOI.groupby('Region Type').describe()["pbi"]
    tot_counts= sum(df_trans["count"])
    df_trans["%count"]= df_trans["count"]/tot_counts*100
    # col8,col9= st.columns(2)
    with col7:
        try:
            fig, ax = plt.subplots(figsize=(4, 5))
            sns.barplot(df_trans, x='Region Type', y="%count", hue= 'Region Type',palette={"Disruptive":"red","Core":"green","Edge":"grey"}, alpha=0.2, ax=ax)
            for i, p in enumerate(ax.patches):
                h, w, x = p.get_height(), p.get_width(), p.get_x()
                xy = (x + w / 2.,max(df_trans["%count"])*1.15)
                text1 = f' {round(df_trans.iloc[i]["count"],2)}\n {round(df_trans.iloc[i]["%count"],2)}%'
                ax.annotate(text=text1, xy=xy, ha='center', va='top')
            ax.set_ylim(0,120)
            ax.set_ylabel(f"% {GOI}")
            if map==1:
                c=df_GOI["Sequence name"  ].iloc[0]
                ax.set_title(f"Compartment proportions for\n{GOI} in sequence {c}")
            else:
                ax.set_title(f"Compartment proportions for {GOI}")
            all_figs.append((f"{GOI} counts", fig))
            st.pyplot(fig)
        except:
            if map ==0:
                st.write(f"No {GOI} were found in GFF file, try another keyword")
            # st.write(f"No {GOI} were found in GFF file, try another keyword")

    # Export/Download graphs as SVG
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
            st.download_button("Download figures as SVG", data=f, file_name="plots_svg.zip", mime="application/zip")

    st.dataframe(data=regions_df)
    if analysis_type:
        st.dataframe(data=filtered_DF_GOI)

st.caption("This program is part of the Disruptomics project. Please cite Balouz V. et al. 2025")
