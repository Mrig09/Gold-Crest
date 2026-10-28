import streamlit as st
import pandas as pd
from io import StringIO
import plotly.graph_objects as go
from src.bio_pipeline.pipeline import (
    parse_single_fasta, parse_genbank_features, run_pipeline,
)

st.set_page_config(page_title="Genomic Mutation Explorer", layout="wide")

# Sidebar
st.sidebar.title("Mutation Explorer")
st.sidebar.caption("Reference vs Variant • GenBank Annotation")

mode = st.sidebar.radio("Analysis mode", ["Single variant", "Variant vs Variant"], index=0)

ref_file = st.sidebar.file_uploader("Upload Reference FASTA", type=["fa","fasta","fna"], key="ref")
if mode == "Single variant":
    var_file = st.sidebar.file_uploader("Upload Variant FASTA", type=["fa","fasta","fna"], key="var")
    varA_file = None
    varB_file = None
else:
    var_file = None
    varA_file = st.sidebar.file_uploader("Upload Variant A FASTA", type=["fa","fasta","fna"], key="varA")
    varB_file = st.sidebar.file_uploader("Upload Variant B FASTA", type=["fa","fasta","fna"], key="varB")

gb_file  = st.sidebar.file_uploader("Upload GenBank Annotation", type=["gb","gbk"], key="gb")
fast_mode = st.sidebar.checkbox("Fast alignment (edlib)", value=True, help="Greatly speeds up global alignment for ~30kb genomes.")

col_a, col_b = st.sidebar.columns(2)
run = col_a.button("Analyze", type="primary")
reset = col_b.button("Reset")

if reset:
    st.experimental_rerun()

st.title("Genomic Mutation Explorer")
st.write("Compare a reference and variant sequence, detect mutations, map to genes, classify effects, and visualize results.")

# Placeholders for sections
summary_ph = st.container()
metrics_ph = st.container()
charts_ph = st.container()
table_ph = st.container()
interp_ph = st.expander("Scientific Interpretation", expanded=False)
concl_ph = st.container()

# Step 1/2: File validation and preprocessing
ref_rec = None
var_rec = None
varA_rec = None
varB_rec = None
features = []

if ref_file is not None:
    try:
        ref_text = ref_file.getvalue().decode("utf-8", errors="replace")
        ref_rec = parse_single_fasta(StringIO(ref_text))
    except Exception as e:
        st.error(f"Reference FASTA error: {e}")

if var_file is not None:
    try:
        var_text = var_file.getvalue().decode("utf-8", errors="replace")
        var_rec = parse_single_fasta(StringIO(var_text))
    except Exception as e:
        st.error(f"Variant FASTA error: {e}")

if varA_file is not None:
    try:
        varA_text = varA_file.getvalue().decode("utf-8", errors="replace")
        varA_rec = parse_single_fasta(StringIO(varA_text))
    except Exception as e:
        st.error(f"Variant A FASTA error: {e}")

if varB_file is not None:
    try:
        varB_text = varB_file.getvalue().decode("utf-8", errors="replace")
        varB_rec = parse_single_fasta(StringIO(varB_text))
    except Exception as e:
        st.error(f"Variant B FASTA error: {e}")

if gb_file is not None:
    try:
        gb_text = gb_file.getvalue().decode("utf-8", errors="replace")
        features = parse_genbank_features(StringIO(gb_text))
    except Exception as e:
        st.warning(f"GenBank parse warning: {e}")

with summary_ph:
    st.subheader("File Summary")
    if mode == "Single variant":
        cols = st.columns(3)
        with cols[0]:
            st.write("Reference:")
            if ref_rec:
                st.write(f"ID: {ref_rec.id}")
                st.write(f"Length: {len(ref_rec.seq)} bp")
            else:
                st.write("—")
        with cols[1]:
            st.write("Variant:")
            if var_rec:
                st.write(f"ID: {var_rec.id}")
                st.write(f"Length: {len(var_rec.seq)} bp")
            else:
                st.write("—")
        with cols[2]:
            st.write("GenBank:")
            st.write(f"Features parsed: {len(features)}")
    else:
        cols = st.columns(4)
        with cols[0]:
            st.write("Reference:")
            if ref_rec:
                st.write(f"ID: {ref_rec.id}")
                st.write(f"Length: {len(ref_rec.seq)} bp")
            else:
                st.write("—")
        with cols[1]:
            st.write("Variant A:")
            if varA_rec:
                st.write(f"ID: {varA_rec.id}")
                st.write(f"Length: {len(varA_rec.seq)} bp")
            else:
                st.write("—")
        with cols[2]:
            st.write("Variant B:")
            if varB_rec:
                st.write(f"ID: {varB_rec.id}")
                st.write(f"Length: {len(varB_rec.seq)} bp")
            else:
                st.write("—")
        with cols[3]:
            st.write("GenBank:")
            st.write(f"Features parsed: {len(features)}")

# Run analysis
result = None
resultA = None
resultB = None
if run:
    if not ref_rec:
        st.error("Please upload a Reference FASTA file.")
    else:
        if mode == "Single variant":
            if not var_rec:
                st.error("Please upload a Variant FASTA file.")
            else:
                with st.spinner("Running alignment and mutation analysis..."):
                    result = run_pipeline(ref_rec, var_rec, features, fast=fast_mode)
        else:
            if not varA_rec or not varB_rec:
                st.error("Please upload both Variant A and Variant B FASTA files.")
            else:
                with st.spinner("Running alignment and mutation analysis for Variant A and Variant B..."):
                    resultA = run_pipeline(ref_rec, varA_rec, features, fast=fast_mode)
                    resultB = run_pipeline(ref_rec, varB_rec, features, fast=fast_mode)

# Display results
if mode == "Single variant" and result:
    with metrics_ph:
        st.subheader("Key Stats")
        mcols_top = st.columns(4)
        mcols_top[0].metric("Total Mutations", result["metrics"]["total"])
        mcols_top[1].metric("Substitutions", result["metrics"]["substitutions"])
        mcols_top[2].metric("Insertions", result["metrics"]["insertions"])
        mcols_top[3].metric("Deletions", result["metrics"]["deletions"])
        df = result["mutations_df"]
        sev_counts = df["severity_level"].value_counts() if not df.empty else {}
        high_count = int(sev_counts.get("high", 0)) if hasattr(sev_counts, "get") else 0
        med_count = int(sev_counts.get("medium", 0)) if hasattr(sev_counts, "get") else 0
        low_count = int(sev_counts.get("low", 0)) if hasattr(sev_counts, "get") else 0
        mcols_bottom = st.columns(3)
        mcols_bottom[0].metric("High-impact", high_count)
        mcols_bottom[1].metric("Medium-impact", med_count)
        mcols_bottom[2].metric("Low-impact", low_count)

    with charts_ph:
        st.subheader("Charts")
        c1, c2 = st.columns([2,1])
        if not result["per_gene"].empty:
            c1.bar_chart(result["per_gene"].set_index("gene")["count"])
        else:
            c1.info("No gene mapping available.")
        if not result["per_type"].empty:
            c2.plotly_chart(result["per_type_fig"])
        else:
            c2.info("No mutation types available.")
        if not result["per_severity"].empty:
            st.plotly_chart(result["per_severity_fig"])
        else:
            st.info("No severity information available.")

        st.subheader("Protein Domain Impact Map")
        protein_map = result.get("protein_map", {}) or {}
        df_mut = result["mutations_df"]
        if not protein_map:
            st.info("Domains not available in GenBank — showing gene region only.")
        else:
            protein_names = sorted(protein_map.keys())
            selected_protein = st.selectbox("Choose protein to visualise", protein_names)
            pinfo = protein_map.get(selected_protein, {})
            length_aa = int(pinfo.get("length_aa") or 0)
            domains = pinfo.get("domains", [])

            df_p = df_mut[df_mut["protein_name"] == selected_protein].copy()
            if df_p.empty:
                st.info("No mutations mapped to this protein.")
            else:
                type_symbols = []
                colors = []
                for _, row in df_p.iterrows():
                    eff = row.get("effect")
                    mtype = row.get("type")
                    sev = row.get("severity_level")
                    if eff == "nonsynonymous":
                        sym = "circle"
                    elif eff == "synonymous":
                        sym = "circle-open"
                    elif mtype in ("deletion",) or eff in ("in_frame_indel", "frameshift"):
                        sym = "triangle-down"
                    elif mtype in ("insertion",):
                        sym = "triangle-up"
                    else:
                        sym = "diamond"
                    if sev == "high":
                        col = "red"
                    elif sev == "medium":
                        col = "orange"
                    else:
                        col = "steelblue"
                    type_symbols.append(sym)
                    colors.append(col)

                fig = go.Figure()
                if length_aa > 0:
                    fig.add_trace(go.Scatter(
                        x=[1, length_aa],
                        y=[0, 0],
                        mode="lines",
                        line=dict(color="black", width=6),
                        hoverinfo="skip",
                        showlegend=False,
                    ))

                for d in domains:
                    crit = d.get("criticality", "none")
                    if crit == "critical":
                        fill = "rgba(220,20,60,0.4)"  # subtle crimson for critical regions
                    else:
                        fill = "rgba(100,149,237,0.35)"  # cornflower blue for other domains
                    fig.add_shape(
                        type="rect",
                        x0=d["aa_start"], x1=d["aa_end"],
                        y0=-0.4, y1=0.4,
                        line=dict(color="rgba(0,0,0,0)"),
                        fillcolor=fill,
                    )
                    fig.add_annotation(
                        x=(d["aa_start"] + d["aa_end"]) / 2,
                        y=0.6,
                        text=d.get("label", "domain"),
                        showarrow=False,
                        font=dict(size=10),
                    )

                fig.add_trace(go.Scatter(
                    x=df_p["aa_pos_in_protein"],
                    y=[0.0] * len(df_p),
                    mode="markers",
                    marker=dict(symbol=type_symbols, color=colors, size=10, line=dict(width=1, color="black")),
                    hovertemplate=(
                        "AA position: %{x}<br>"+
                        "Effect: %{customdata[0]}<br>"+
                        "Severity: %{customdata[1]}<br>"+
                        "AA change: %{customdata[2]}<br>"+
                        "Domain: %{customdata[3]}<extra></extra>"
                    ),
                    customdata=df_p[["effect","severity_level","aa_change","domain_label"]].values,
                    name="Mutations",
                ))

                fig.update_layout(
                    xaxis_title="Amino acid position",
                    yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    showlegend=False,
                    margin=dict(l=40, r=20, t=40, b=40),
                )

                if not domains:
                    st.info("Domains not available in GenBank — showing gene region only.")

                st.plotly_chart(fig, use_container_width=True)

    with table_ph:
        st.subheader("Mutation Table")
        df = result["mutations_df"]
        st.dataframe(df, use_container_width=True)
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("Download CSV", data=csv, file_name="mutations.csv", mime="text/csv")

    with interp_ph:
        st.markdown(result["interpretation"])

    with concl_ph:
        st.subheader("Conclusion")
        st.write(result["conclusion"])

elif mode == "Variant vs Variant" and resultA and resultB:
    with metrics_ph:
        st.subheader("Comparison Stats")
        dfA = resultA["mutations_df"]
        dfB = resultB["mutations_df"]
        totalA = int(len(dfA))
        totalB = int(len(dfB))
        scoreA = int(dfA["severity_score"].sum()) if not dfA.empty else 0
        scoreB = int(dfB["severity_score"].sum()) if not dfB.empty else 0
        mcols = st.columns(4)
        mcols[0].metric("Variant A mutations", totalA)
        mcols[1].metric("Variant B mutations", totalB)
        mcols[2].metric("Variant A severity score", scoreA)
        mcols[3].metric("Variant B severity score", scoreB)

    dfA = resultA["mutations_df"].copy()
    dfB = resultB["mutations_df"].copy()
    if not dfA.empty:
        dfA["_key"] = dfA[["ref_pos","type","ref_base","var_base","length","gene","effect"]].astype(str).agg("|".join, axis=1)
    else:
        dfA["_key"] = []
    if not dfB.empty:
        dfB["_key"] = dfB[["ref_pos","type","ref_base","var_base","length","gene","effect"]].astype(str).agg("|".join, axis=1)
    else:
        dfB["_key"] = []

    keysA = set(dfA["_key"]) if not dfA.empty else set()
    keysB = set(dfB["_key"]) if not dfB.empty else set()
    shared_keys = keysA & keysB
    onlyA_keys = keysA - keysB
    onlyB_keys = keysB - keysA

    shared_df = dfA[dfA["_key"].isin(shared_keys)].copy() if shared_keys else dfA.iloc[0:0].copy()
    uniqueA_df = dfA[dfA["_key"].isin(onlyA_keys)].copy() if onlyA_keys else dfA.iloc[0:0].copy()
    uniqueB_df = dfB[dfB["_key"].isin(onlyB_keys)].copy() if onlyB_keys else dfB.iloc[0:0].copy()

    with charts_ph:
        st.subheader("Severity and Gene-wise Comparison")
        c1, c2 = st.columns(2)
        if not resultA["per_severity"].empty:
            c1.plotly_chart(resultA["per_severity_fig"], use_container_width=True)
        else:
            c1.info("No severity data for Variant A.")
        if not resultB["per_severity"].empty:
            c2.plotly_chart(resultB["per_severity_fig"], use_container_width=True)
        else:
            c2.info("No severity data for Variant B.")

        per_geneA = resultA["per_gene"].copy()
        per_geneB = resultB["per_gene"].copy()
        if not per_geneA.empty or not per_geneB.empty:
            per_geneA = per_geneA.rename(columns={"count": "count_A"})
            per_geneB = per_geneB.rename(columns={"count": "count_B"})
            gene_stats = per_geneA.merge(per_geneB, on="gene", how="outer").fillna(0)
            gene_stats = gene_stats.set_index("gene")[["count_A","count_B"]]
            st.bar_chart(gene_stats)
        else:
            st.info("No gene-level mutation data to compare.")

    with table_ph:
        st.subheader("Variant Comparison Tables")
        st.markdown("**Shared mutations (in both Variant A and Variant B)**")
        st.dataframe(shared_df.drop(columns=["_key"]), use_container_width=True)
        st.markdown("**Unique to Variant A**")
        st.dataframe(uniqueA_df.drop(columns=["_key"]), use_container_width=True)
        st.markdown("**Unique to Variant B**")
        st.dataframe(uniqueB_df.drop(columns=["_key"]), use_container_width=True)

    with interp_ph:
        interp_text = []
        if scoreA > scoreB:
            interp_text.append("Variant A shows a higher overall severity score than Variant B.")
        elif scoreB > scoreA:
            interp_text.append("Variant B shows a higher overall severity score than Variant A.")
        else:
            interp_text.append("Variant A and Variant B have similar overall severity scores.")
        interp_text.append(f"Shared mutations: {len(shared_df)}")
        interp_text.append(f"Unique to Variant A: {len(uniqueA_df)}")
        interp_text.append(f"Unique to Variant B: {len(uniqueB_df)}")
        st.markdown("\n".join(f"- {line}" for line in interp_text))

    with concl_ph:
        st.subheader("Conclusion")
        st.write("Comparison mode highlights shared and unique mutations and their severity across Variant A and Variant B against the same reference.")
