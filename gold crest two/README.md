# Bioinformatics Mutation Analysis App

This app compares a reference vs variant sequence, detects mutations, maps them to genes using a GenBank file, classifies effects, and presents visualizations.

## Quickstart

1) Install dependencies

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

2) Run the app

```
streamlit run app.py
```

3) Upload files in the sidebar
- Reference FASTA (single sequence)
- Variant FASTA (single sequence)
- GenBank annotation (.gb/.gbk with features)

4) Explore results
- File summary
- Alignment-derived mutations
- Charts and table
- Interpretation

## Notes
- Works best for viral genomes (e.g., SARS-CoV-2). For other organisms, SARS-CoV-2 specific insights are skipped.
- GenBank should include CDS features for gene/protein mapping.
