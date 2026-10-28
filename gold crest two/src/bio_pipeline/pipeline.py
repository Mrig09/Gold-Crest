from typing import List, Dict, Any, Optional, Tuple
from io import BytesIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio import pairwise2
import pandas as pd
import numpy as np
import plotly.express as px
import re
try:
    import edlib
except Exception:  # edlib optional, fallback to Biopython
    edlib = None
try:
    from Bio.Align import PairwiseAligner
except Exception:
    PairwiseAligner = None

# ---------------------------
# Parsing
# ---------------------------

def parse_single_fasta(handle: BytesIO) -> SeqRecord:
    records = list(SeqIO.parse(handle, "fasta"))
    if len(records) == 0:
        raise ValueError("No sequences found in FASTA.")
    if len(records) > 1:
        raise ValueError("Multiple sequences found. Provide a single-sequence FASTA.")
    return records[0]


def parse_genbank_features(handle: BytesIO) -> List[Dict[str, Any]]:
    features: List[Dict[str, Any]] = []
    try:
        record = next(SeqIO.parse(handle, "genbank"))
    except StopIteration:
        return []

    def _criticality_from_text(text: str) -> str:
        t = text.lower()
        critical_keywords = [
            "active site",
            "catalytic",
            "transmembrane",
            "tm domain",
            "cleavage site",
            "polymerase",
            "receptor binding",
            "binding site",
            "antigenic site",
            "drug-binding",
            "drug binding",
            "resistance",
        ]
        domain_keywords = ["domain", "motif", "repeat", "loop", "region"]
        if any(k in t for k in critical_keywords):
            return "critical"
        if any(k in t for k in domain_keywords):
            return "domain"
        return "none"

    for feat in record.features:
        if feat.type in ("CDS", "gene", "region", "misc_feature", "mat_peptide"):
            try:
                start = int(feat.location.start)
                end = int(feat.location.end) - 1
                strand = int(feat.location.strand or 1)
            except Exception:
                continue
            gene = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", ["unknown"]))[0]
            product = feat.qualifiers.get("product", ["unknown"])[0]
            notes = []
            if "note" in feat.qualifiers:
                notes.extend(feat.qualifiers["note"])
            if "function" in feat.qualifiers:
                notes.extend(feat.qualifiers["function"])
            if "region_name" in feat.qualifiers:
                notes.extend(feat.qualifiers["region_name"])
            note_text = "; ".join(notes)
            crit = _criticality_from_text(" ".join([product, note_text])) if (product or note_text) else "none"
            features.append({
                "type": feat.type,
                "start": start + 1,  # to 1-based
                "end": end + 1,
                "strand": strand,
                "gene": gene,
                "product": product,
                "note": note_text,
                "criticality": crit,
            })
    # sort by start
    features.sort(key=lambda x: (x["start"], x["end"]))
    return features

# ---------------------------
# Alignment
# ---------------------------

def global_align(ref_seq: str, var_seq: str) -> Tuple[str, str]:
    if len(ref_seq) == 0 or len(var_seq) == 0:
        return ref_seq, var_seq
    # Prefer fast C-backed PairwiseAligner if available
    if PairwiseAligner is not None:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        aln = next(iter(aligner.align(ref_seq, var_seq)))
        # Reconstruct aligned strings from coordinates
        t_blocks, q_blocks = aln.aligned
        a_ref_parts = []
        a_var_parts = []
        prev_t = 0
        prev_q = 0
        for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
            # gaps before this block
            if ts > prev_t and qs == prev_q:
                # deletion in variant (gap in var)
                a_ref_parts.append(ref_seq[prev_t:ts])
                a_var_parts.append('-' * (ts - prev_t))
            elif qs > prev_q and ts == prev_t:
                # insertion in variant (gap in ref)
                a_ref_parts.append('-' * (qs - prev_q))
                a_var_parts.append(var_seq[prev_q:qs])
            elif ts > prev_t and qs > prev_q:
                # mismatches without gaps
                n = min(ts - prev_t, qs - prev_q)
                a_ref_parts.append(ref_seq[prev_t:prev_t + n])
                a_var_parts.append(var_seq[prev_q:prev_q + n])
                # handle any residual as gaps if lengths differ
                if (ts - prev_t) > n:
                    a_ref_parts.append(ref_seq[prev_t + n:ts])
                    a_var_parts.append('-' * (ts - prev_t - n))
                if (qs - prev_q) > n:
                    a_ref_parts.append('-' * (qs - prev_q - n))
                    a_var_parts.append(var_seq[prev_q + n:qs])
            # aligned block (match/mismatch run)
            a_ref_parts.append(ref_seq[ts:te])
            a_var_parts.append(var_seq[qs:qe])
            prev_t, prev_q = te, qe
        # tails
        if prev_t < len(ref_seq) and prev_q == len(var_seq):
            a_ref_parts.append(ref_seq[prev_t:])
            a_var_parts.append('-' * (len(ref_seq) - prev_t))
        elif prev_q < len(var_seq) and prev_t == len(ref_seq):
            a_ref_parts.append('-' * (len(var_seq) - prev_q))
            a_var_parts.append(var_seq[prev_q:])
        elif prev_t < len(ref_seq) and prev_q < len(var_seq):
            n = min(len(ref_seq) - prev_t, len(var_seq) - prev_q)
            a_ref_parts.append(ref_seq[prev_t:prev_t + n])
            a_var_parts.append(var_seq[prev_q:prev_q + n])
            if (len(ref_seq) - prev_t) > n:
                a_ref_parts.append(ref_seq[prev_t + n:])
                a_var_parts.append('-' * (len(ref_seq) - prev_t - n))
            if (len(var_seq) - prev_q) > n:
                a_ref_parts.append('-' * (len(var_seq) - prev_q - n))
                a_var_parts.append(var_seq[prev_q + n:])
        a_ref = ''.join(a_ref_parts)
        a_var = ''.join(a_var_parts)
        return a_ref, a_var
    # Fallback to pairwise2 (Python-level, slower)
    alignments = pairwise2.align.globalms(ref_seq, var_seq, 2, -1, -5, -1, one_alignment_only=True)
    a_ref, a_var, *_ = alignments[0]
    return a_ref, a_var


def _expand_cigar_to_alignment(ref_seq: str, var_seq: str, cigar: str) -> Tuple[str, str]:
    # CIGAR tokens like 10=1X5=2I3=... where:
    # '=' or 'M' -> match, 'X' -> mismatch, 'I' -> insertion to target (var), 'D' -> deletion from target
    tokens = re.findall(r"(\d+)([=XIDM])", cigar)
    i = j = 0
    a_ref = []
    a_var = []
    for length_str, op in tokens:
        n = int(length_str)
        if op in ('=', 'M', 'X'):
            a_ref.append(ref_seq[i:i+n])
            a_var.append(var_seq[j:j+n])
            i += n
            j += n
        elif op == 'I':  # insertion to target (var) => gap in ref
            a_ref.append('-' * n)
            a_var.append(var_seq[j:j+n])
            j += n
        elif op == 'D':  # deletion from target (var) => gap in var
            a_ref.append(ref_seq[i:i+n])
            a_var.append('-' * n)
            i += n
        else:
            # Unknown op, fallback (shouldn't happen)
            a_ref.append(ref_seq[i:i+n])
            a_var.append(var_seq[j:j+n])
            i += n
            j += n
    return ''.join(a_ref), ''.join(a_var)


def global_align_fast(ref_seq: str, var_seq: str) -> Tuple[str, str]:
    if edlib is None:
        return global_align(ref_seq, var_seq)
    if len(ref_seq) == 0 or len(var_seq) == 0:
        return ref_seq, var_seq
    res = edlib.align(ref_seq, var_seq, task="path", mode="NW", k=-1, additionalEqualities=None, allowSecondary=False, gapPenalty=1)
    cigar = res.get('cigar')
    if not cigar:
        return global_align(ref_seq, var_seq)
    return _expand_cigar_to_alignment(ref_seq, var_seq, cigar)


def is_fast_available() -> bool:
    return edlib is not None

# ---------------------------
# Mutation detection
# ---------------------------

def detect_point_events(aln_ref: str, aln_var: str) -> List[Dict[str, Any]]:
    events: List[Dict[str, Any]] = []
    rpos = 0  # 1-based genomic positions
    vpos = 0
    for i, (r, v) in enumerate(zip(aln_ref, aln_var)):
        if r != "-":
            rpos += 1
        if v != "-":
            vpos += 1
        if r == v:
            continue
        if r == "-" and v != "-":
            events.append({
                "aln_index": i,
                "type": "insertion",
                "ref_base": "-",
                "var_base": v,
                "ref_pos": rpos if r != "-" else rpos,  # position after previous base
                "var_pos": vpos,
            })
        elif v == "-" and r != "-":
            events.append({
                "aln_index": i,
                "type": "deletion",
                "ref_base": r,
                "var_base": "-",
                "ref_pos": rpos,
                "var_pos": vpos if v != "-" else vpos,
            })
        else:
            events.append({
                "aln_index": i,
                "type": "substitution",
                "ref_base": r,
                "var_base": v,
                "ref_pos": rpos,
                "var_pos": vpos,
            })
    return events


def collapse_indels(events: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if not events:
        return []
    collapsed: List[Dict[str, Any]] = []
    i = 0
    while i < len(events):
        e = events[i]
        if e["type"] == "substitution":
            collapsed.append({**e, "length": 1})
            i += 1
            continue
        # merge consecutive same-type indels with consecutive alignment indices
        j = i + 1
        seq = [e]
        while j < len(events) and events[j]["type"] == e["type"] and events[j]["aln_index"] == events[j-1]["aln_index"] + 1:
            seq.append(events[j])
            j += 1
        # build merged record
        length = len(seq)
        if e["type"] == "insertion":
            var_bases = ''.join(x["var_base"] for x in seq)
            rec = {
                "type": "insertion",
                "ref_pos": seq[0]["ref_pos"],
                "var_pos": seq[0]["var_pos"],
                "ref_base": "-",
                "var_base": var_bases,
                "length": length,
            }
        else:  # deletion
            ref_bases = ''.join(x["ref_base"] for x in seq)
            rec = {
                "type": "deletion",
                "ref_pos": seq[0]["ref_pos"],
                "var_pos": seq[0]["var_pos"],
                "ref_base": ref_bases,
                "var_base": "-",
                "length": length,
            }
        collapsed.append(rec)
        i = j
    return collapsed

# ---------------------------
# Gene mapping and effect classification
# ---------------------------

def map_gene(pos1: int, features: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    for f in features:
        if f["start"] <= pos1 <= f["end"] and f["type"] == "CDS":
            return f
    return None


def assign_severity(m: Dict[str, Any], features: List[Dict[str, Any]]) -> Tuple[str, int]:
    pos = m.get("ref_pos") or m.get("var_pos") or 0
    overlapping = [
        f for f in features
        if f.get("start") is not None and f.get("end") is not None and f["start"] <= pos <= f["end"]
    ]
    cds = None
    for f in overlapping:
        if f.get("type") == "CDS":
            cds = f
            break
    critical_hits = [f for f in overlapping if f.get("criticality") == "critical"]
    domain_hits = [
        f
        for f in overlapping
        if f.get("criticality") == "domain" or f.get("type") in ("region", "misc_feature", "mat_peptide")
    ]

    effect = m.get("effect", "intergenic")
    mtype = m.get("type", "")
    length = int(m.get("length", 1) or 1)

    if effect in ("synonymous", "intergenic"):
        if cds is None:
            return "low", 1

    if mtype in ("insertion", "deletion") and length > 10:
        return "high", 5

    if effect in ("nonsense", "frameshift", "stop_loss"):
        return "high", 5

    if critical_hits:
        return "high", 5

    if effect == "in_frame_indel" and cds is not None:
        if length <= 10:
            if domain_hits:
                return "high", 5
            return "medium", 3
        return "high", 5

    if effect == "nonsynonymous" and cds is not None:
        if domain_hits:
            return "high", 5
        return "medium", 3

    if cds is not None and effect in ("coding_region",):
        return "medium", 3

    return "low", 1


def codon_index_within_cds(genomic_pos1: int, cds_start1: int, strand: int) -> int:
    if strand >= 0:
        return (genomic_pos1 - cds_start1) // 3
    else:
        return (cds_start1 - genomic_pos1) // 3


def get_codon(seq: Seq, genomic_pos1: int, cds: Dict[str, Any]) -> Optional[str]:
    cds_start = cds["start"]
    cds_end = cds["end"]
    strand = cds["strand"]
    if not (cds_start <= genomic_pos1 <= cds_end):
        return None
    offset = (genomic_pos1 - cds_start) if strand >= 0 else (cds_end - genomic_pos1)
    frame_start = (offset // 3) * 3
    if strand >= 0:
        s = cds_start + frame_start
        e = s + 3
        codon = str(seq[s-1:e-1])
    else:
        e = cds_end - frame_start
        s = e - 3
        codon = str(seq[s-1:e-1].reverse_complement())
    if len(codon) == 3:
        return codon
    return None


def translate_codon(codon: str) -> str:
    try:
        aa = str(Seq(codon).translate(to_stop=False))
        return aa[0] if aa else "?"
    except Exception:
        return "?"


def classify_effect(m: Dict[str, Any], ref_seq: Seq, var_seq: Seq, features: List[Dict[str, Any]]) -> Dict[str, Any]:
    gene = map_gene(m.get("ref_pos", m.get("var_pos", 0)), features)
    m["gene"] = gene["gene"] if gene else None
    m["product"] = gene["product"] if gene else None

    # default type
    m["effect"] = "intergenic"
    m["aa_change"] = None

    if not gene or gene["type"] != "CDS":
        return m

    # Substitution effect
    if m["type"] == "substitution":
        ref_codon = get_codon(ref_seq, m["ref_pos"], gene)
        var_codon = get_codon(var_seq, m["var_pos"], gene)
        if ref_codon and var_codon:
            ref_aa = translate_codon(ref_codon)
            var_aa = translate_codon(var_codon)
            m["aa_change"] = f"{ref_aa}->{var_aa}"
            if var_aa == ref_aa:
                m["effect"] = "synonymous"
            elif var_aa == "*":
                m["effect"] = "nonsense"
            elif ref_aa == "*" and var_aa != "*":
                m["effect"] = "stop_loss"
            else:
                m["effect"] = "nonsynonymous"
        else:
            m["effect"] = "coding_region"
        return m

    # Indels
    if m["type"] in ("insertion", "deletion"):
        L = int(m.get("length", 1))
        if L % 3 == 0:
            m["effect"] = "in_frame_indel"
        else:
            m["effect"] = "frameshift"
        return m

    return m


def build_protein_map(features: List[Dict[str, Any]]) -> Dict[str, Any]:
    proteins: Dict[str, Any] = {}
    cds_features = [f for f in features if f.get("type") == "CDS"]
    other_features = [f for f in features if f.get("type") != "CDS"]

    for cds in cds_features:
        gene = cds.get("gene") or "unknown"
        product = cds.get("product") or "unknown"
        label = f"{gene} ({product})"
        if label in proteins:
            continue
        cds_start = cds.get("start")
        cds_end = cds.get("end")
        strand = cds.get("strand", 1)
        if cds_start is None or cds_end is None:
            continue
        length_nt = max(0, int(cds_end) - int(cds_start) + 1)
        length_aa = length_nt // 3 if length_nt > 0 else 0
        doms = []
        for f in other_features:
            fs = f.get("start")
            fe = f.get("end")
            if fs is None or fe is None:
                continue
            if fe < cds_start or fs > cds_end:
                continue
            ov_start = max(fs, cds_start)
            ov_end = min(fe, cds_end)
            if ov_end < ov_start:
                continue
            if strand >= 0:
                offset_start = ov_start - cds_start
                offset_end = ov_end - cds_start
            else:
                offset_start = cds_end - ov_end
                offset_end = cds_end - ov_start
            aa_start = offset_start // 3 + 1
            aa_end = offset_end // 3 + 1
            dom_label = f.get("product") or f.get("note") or f.get("type") or "region"
            doms.append({
                "label": dom_label,
                "aa_start": int(aa_start),
                "aa_end": int(aa_end),
                "criticality": f.get("criticality", "none"),
            })
        proteins[label] = {
            "gene": gene,
            "product": product,
            "cds_start": cds_start,
            "cds_end": cds_end,
            "strand": strand,
            "length_aa": length_aa,
            "domains": doms,
        }
    return proteins


def annotate_protein_context(df: pd.DataFrame, features: List[Dict[str, Any]], ref_seq: Seq) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    proteins = build_protein_map(features)
    if df.empty or not proteins:
        for col in ["protein_name", "aa_pos_in_protein", "domain_label"]:
            if col not in df.columns:
                df[col] = pd.Series([None] * len(df))
        return df, proteins

    label_to_cds = {}
    for f in features:
        if f.get("type") == "CDS":
            gene = f.get("gene") or "unknown"
            product = f.get("product") or "unknown"
            label = f"{gene} ({product})"
            label_to_cds[label] = f

    protein_names = []
    aa_positions = []
    domain_labels = []

    for _, row in df.iterrows():
        gene_name = row.get("gene")
        if not gene_name:
            protein_names.append(None)
            aa_positions.append(None)
            domain_labels.append(None)
            continue
        protein_label = None
        cds = None
        for label, info in proteins.items():
            if info.get("gene") == gene_name:
                protein_label = label
                cds = label_to_cds.get(label)
                break
        if protein_label is None or cds is None:
            protein_names.append(None)
            aa_positions.append(None)
            domain_labels.append(None)
            continue
        cds_start = cds.get("start")
        cds_end = cds.get("end")
        strand = cds.get("strand", 1)
        pos = row.get("ref_pos") or row.get("var_pos")
        if cds_start is None or cds_end is None or not pos:
            protein_names.append(protein_label)
            aa_positions.append(None)
            domain_labels.append(None)
            continue
        if strand >= 0:
            offset = int(pos) - int(cds_start)
        else:
            offset = int(cds_end) - int(pos)
        if offset < 0:
            aa_idx = None
        else:
            aa_idx = offset // 3
        aa_pos = aa_idx + 1 if aa_idx is not None else None
        dom_label = None
        if aa_pos is not None:
            doms = proteins[protein_label].get("domains", [])
            for d in doms:
                if d["aa_start"] <= aa_pos <= d["aa_end"]:
                    dom_label = d["label"]
                    break
        protein_names.append(protein_label)
        aa_positions.append(aa_pos)
        domain_labels.append(dom_label)

    df = df.copy()
    df["protein_name"] = protein_names
    df["aa_pos_in_protein"] = aa_positions
    df["domain_label"] = domain_labels
    return df, proteins

# ---------------------------
# Heuristics (SARS-CoV-2)
# ---------------------------

SARS2_REGIONS = [
    {"gene": "S", "label": "Spike RBD", "start": 23013, "end": 23651},  # approx Wuhan ref positions (AA 319-541)
    {"gene": "S", "label": "Spike NTD", "start": 21563, "end": 22552},
    {"gene": "S", "label": "Furin site", "start": 23603, "end": 23618},
]


def heuristic_annotation(m: Dict[str, Any]) -> Optional[str]:
    g = m.get("gene")
    pos = m.get("ref_pos") or m.get("var_pos")
    if not (g and pos):
        return None
    for r in SARS2_REGIONS:
        if r["gene"] == g and r["start"] <= pos <= r["end"]:
            if r["label"] == "Spike RBD":
                return "Spike RBD region: possible impact on ACE2/antibody binding"
            if r["label"] == "Spike NTD":
                return "Spike NTD: potential antigenic supersite involvement"
            if r["label"] == "Furin site":
                return "Near furin cleavage site: possible effect on activation"
    return None

# ---------------------------
# Pipeline
# ---------------------------

def run_pipeline(ref_rec: SeqRecord, var_rec: SeqRecord, features: List[Dict[str, Any]], fast: bool = False) -> Dict[str, Any]:
    ref_seq = ref_rec.seq.upper()
    var_seq = var_rec.seq.upper()

    if fast:
        a_ref, a_var = global_align_fast(str(ref_seq), str(var_seq))
    else:
        a_ref, a_var = global_align(str(ref_seq), str(var_seq))
    point_events = detect_point_events(a_ref, a_var)
    events = collapse_indels(point_events)

    # Classify, annotate, and assign severity
    annotated: List[Dict[str, Any]] = []
    for m in events:
        m = classify_effect(m, ref_seq, var_seq, features)
        hint = heuristic_annotation(m)
        m["impact_note"] = hint
        sev_level, sev_score = assign_severity(m, features)
        m["severity_level"] = sev_level
        m["severity_score"] = sev_score
        annotated.append(m)

    # DataFrame
    df = pd.DataFrame(annotated)
    df, protein_map = annotate_protein_context(df, features, ref_seq)
    if df.empty:
        df = pd.DataFrame(columns=[
            "ref_pos","var_pos","type","ref_base","var_base","length","gene","product","aa_change","effect","impact_note","severity_level","severity_score"
        ])

    # Metrics
    metrics = {
        "total": int(len(df)),
        "substitutions": int((df["type"] == "substitution").sum()) if not df.empty else 0,
        "insertions": int((df["type"] == "insertion").sum()) if not df.empty else 0,
        "deletions": int((df["type"] == "deletion").sum()) if not df.empty else 0,
    }

    # Per gene
    if not df.empty and "gene" in df:
        per_gene = df.groupby(df["gene"].fillna("intergenic")).size().reset_index(name="count").rename(columns={"gene":"gene"})
    else:
        per_gene = pd.DataFrame(columns=["gene","count"])

    # Per type pie
    if not df.empty:
        per_type = df.groupby("effect").size().reset_index(name="count")
        per_type_fig = px.pie(per_type, names="effect", values="count", title="Mutation Effects")
    else:
        per_type = pd.DataFrame(columns=["effect","count"])
        per_type_fig = px.pie(per_type, names="effect", values="count", title="Mutation Effects")

    # Per severity
    if not df.empty and "severity_level" in df:
        per_severity = df.groupby("severity_level").size().reset_index(name="count")
        per_severity_fig = px.pie(per_severity, names="severity_level", values="count", title="Mutation Severity")
    else:
        per_severity = pd.DataFrame(columns=["severity_level","count"])
        per_severity_fig = px.pie(per_severity, names="severity_level", values="count", title="Mutation Severity")

    # Interpretation (brief)
    interp_lines = []
    interp_lines.append("- Most mutated gene: " + (per_gene.sort_values("count", ascending=False).iloc[0]["gene"] if not per_gene.empty else "n/a"))
    if not per_type.empty:
        nonsyn = int(per_type.loc[per_type["effect"] == "nonsynonymous", "count"].sum())
        syn = int(per_type.loc[per_type["effect"] == "synonymous", "count"].sum())
        interp_lines.append(f"- Nonsynonymous: {nonsyn}, Synonymous: {syn}")
    # Spike highlights
    spike_hits = df[(df["gene"] == "S") & df["impact_note"].notna()] if not df.empty else pd.DataFrame()
    if not spike_hits.empty:
        interp_lines.append("- Spike critical-region mutations detected (RBD/NTD/furin). Potential effects on entry/antibody binding.")

    interpretation = "\n".join(interp_lines) if interp_lines else "No notable patterns detected."

    conclusion = (
        "Results are preliminary. Effects are inferred from simple rules and annotations. "
        "For drug resistance or immune escape, consult literature and structural analyses."
    )

    return {
        "metrics": metrics,
        "mutations_df": df,
        "per_gene": per_gene,
        "per_type": per_type,
        "per_type_fig": per_type_fig,
        "per_severity": per_severity,
        "per_severity_fig": per_severity_fig,
        "protein_map": protein_map,
        "interpretation": interpretation,
        "conclusion": conclusion,
        "alignment": {"ref": a_ref, "var": a_var},
    }
