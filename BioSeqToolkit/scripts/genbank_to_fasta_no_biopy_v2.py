#!/usr/bin/env python3
"""
genbank_to_fasta_no_biopy_v2.py

Improved GenBank/FASTA extractor without Biopython.

Features:
 - Accepts files / directories / globs (batch mode).
 - Detects sequence from ORIGIN or SEQUENCE blocks.
 - Parses FEATURES (CDS, mat_peptide, sig_peptide) and qualifiers.
 - Extracts ranges supporting join(...) and complement(...).
 - Writes:
     - <basename>_full.fasta
     - <basename>_cds_seqs.fasta
     - <basename>_cds_proteins.fasta (if --translate)
     - <basename>_<feature>_seqs.fasta for mat_peptide & sig_peptide if present
 - Translation honors /codon_start and supports --to-stop (stop at first in-frame '*').
 - Output per input stored in timestamped folder; aggregate summary CSV written.
 - Uses only Python standard library (pandas optional).

Usage:
  python genbank_to_fasta_no_biopy_v2.py inputs... [--outdir OUTDIR] [--translate] [--to-stop]

Examples:
  python genbank_to_fasta_no_biopy_v2.py data\seq.gb --translate --to-stop
  python genbank_to_fasta_no_biopy_v2.py "data\*.gb" --outdir results
  python genbank_to_fasta_no_biopy_v2.py data --translate

"""
import argparse
import glob
import os
import re
from pathlib import Path
from datetime import datetime
import csv
from typing import List, Tuple, Dict

# Optional pandas
try:
    import pandas as pd
    PANDAS = True
except Exception:
    PANDAS = False

# ---------------- utilities ----------------
def resolve_inputs(inputs: List[str]) -> List[Path]:
    exts = ('*.fa','*.fasta','*.fna','*.ffn','*.gb','*.gbk','*.gbff','*.embl')
    out = []
    for inp in inputs:
        p = Path(inp)
        if p.is_file():
            out.append(p.resolve())
        elif p.is_dir():
            for e in exts:
                out.extend(sorted(p.glob(e)))
        else:
            for m in glob.glob(inp):
                out.append(Path(m).resolve())
    seen = set(); uniq = []
    for p in out:
        if p not in seen:
            seen.add(p); uniq.append(p)
    return uniq

def safe_basename(p: Path) -> str:
    name = p.stem
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)

def make_outdir(base_out: str, infile: Path) -> Path:
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    folder = Path(base_out) / f"{safe_basename(infile)}_{stamp}"
    folder.mkdir(parents=True, exist_ok=True)
    return folder

# ---------------- detection / reading ----------------
def detect_format_by_content(path: Path) -> str:
    with open(path, 'r', errors='ignore') as fh:
        head = fh.read(4096)
    if re.search(r'^\s*LOCUS\b', head, flags=re.MULTILINE):
        return 'genbank'
    if re.search(r'^\s*ORIGIN\b', head, flags=re.MULTILINE):
        return 'genbank'
    if re.search(r'^\s*>', head, flags=re.MULTILINE):
        return 'fasta'
    ext = path.suffix.lower()
    if ext in ('.gb','.gbk','.gbff'):
        return 'genbank'
    if ext in ('.fa','.fasta','.fna'):
        return 'fasta'
    return 'fasta'

def read_file_lines(path: Path) -> List[str]:
    with open(path, 'r', errors='ignore') as fh:
        return [line.rstrip('\n') for line in fh]

# ---------------- GenBank parsing ----------------
def parse_sections(lines: List[str]) -> Tuple[List[str], List[str]]:
    """Return (annotation_lines, origin_lines_or_sequence_lines)."""
    anot = []; seq_lines = []; in_origin = False; in_seqblock = False
    for line in lines:
        if re.match(r'^\s*ORIGIN', line, flags=re.IGNORECASE):
            in_origin = True; in_seqblock = True; continue
        # Some GenBank-like exports use a SEQUENCE: block -- detect 'SEQUENCE' line
        if re.match(r'^\s*SEQUENCE', line, flags=re.IGNORECASE):
            in_seqblock = True; continue
        if not in_seqblock:
            anot.append(line)
        else:
            # if we hit a terminator (//) or next record marker, stop collecting sequence
            if re.match(r'^\s*//', line):
                in_seqblock = False
                continue
            seq_lines.append(line)
    return anot, seq_lines

def build_sequence(seq_lines: List[str]) -> str:
    s = ''.join(seq_lines)
    s = re.sub(r'[^A-Za-z]', '', s)
    s = s.upper().replace('U','T')
    return s

def extract_features_block(annotation_lines: List[str]) -> List[str]:
    start = None
    for i, line in enumerate(annotation_lines):
        if re.match(r'^\s*FEATURES', line, flags=re.IGNORECASE):
            start = i; break
    if start is None:
        return [l for l in annotation_lines if re.match(r'^\s{5,}\w+', l)]
    return annotation_lines[start+1:]

def parse_features(lines: List[str]) -> List[Dict]:
    features = []; cur = None
    for raw in lines:
        m = re.match(r'^\s{5,}(\S+)\s+(.+)', raw)
        if m:
            if cur: features.append(cur)
            key = m.group(1); loc = m.group(2).strip()
            cur = {'key': key, 'location': loc, 'qualifiers': {}}
            continue
        qm = re.match(r'^\s{21,}(/[-\w]+)(?:=(?:\"([^\"]*)\"|([^\"].*)))?', raw)
        if qm and cur is not None:
            qual_name = qm.group(1).lstrip('/')
            qual_val = qm.group(2) if qm.group(2) is not None else (qm.group(3) or "")
            cur['qualifiers'][qual_name] = cur['qualifiers'].get(qual_name, []) + [qual_val.strip()]
        else:
            cont = re.match(r'^\s{21,}(.+)', raw)
            if cont and cur is not None and cur['qualifiers']:
                last = list(cur['qualifiers'].keys())[-1]
                cur['qualifiers'][last][-1] += ' ' + cont.group(1).strip()
    if cur: features.append(cur)
    return features

# ---------------- Location parsing & extraction ----------------
def reverse_complement(seq: str) -> str:
    trans = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                          'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')
    return seq.translate(trans)[::-1]

def split_top_level_comma(s: str) -> List[str]:
    parts=[]; depth=0; cur=[]
    for ch in s:
        if ch=='(': depth+=1; cur.append(ch)
        elif ch==')': depth-=1; cur.append(ch)
        elif ch==',' and depth==0:
            parts.append(''.join(cur).strip()); cur=[]
        else: cur.append(ch)
    if cur: parts.append(''.join(cur).strip())
    return parts

def extract_by_location(location: str, full_seq: str) -> str:
    location = location.strip()
    def _extract(loc):
        loc = loc.strip()
        mcomp = re.match(r'^(?:complement)\((.+)\)$', loc)
        if mcomp:
            inner = _extract(mcomp.group(1)); return reverse_complement(inner)
        mjoin = re.match(r'^(?:join)\((.+)\)$', loc)
        if mjoin:
            parts = split_top_level_comma(mjoin.group(1))
            return ''.join(_extract(p) for p in parts)
        m = re.search(r'(\d+)\s*\.\.\s*(\d+)', loc)
        if m:
            a=int(m.group(1)); b=int(m.group(2)); return full_seq[a-1:b]
        m2 = re.search(r'(\d+)', loc)
        if m2:
            a=int(m2.group(1)); return full_seq[a-1:a]
        return ''
    return _extract(location)

# ---------------- Translation ----------------
CODON_TABLE = {
'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W'
}

def translate_nt(seq_nt: str, codon_start: int = 1, to_stop: bool = False) -> str:
    if codon_start not in (1,2,3):
        codon_start = 1
    s = seq_nt.upper().replace('U','T')
    # adjust frame by codon_start: GenBank codon_start is 1-based
    offset = (codon_start - 1) % 3
    aa = []
    for i in range(offset, len(s)-2, 3):
        codon = s[i:i+3]
        aa_res = CODON_TABLE.get(codon, 'X')
        if aa_res == '*' and to_stop:
            aa.append('*')
            break
        aa.append(aa_res)
    return ''.join(aa)

# ---------------- QC ----------------
def basic_seq_qc(seq: str, min_len: int = 30) -> Tuple[bool,str]:
    if not seq: return False, "empty"
    sequ = seq.upper()
    valid = sum(1 for c in sequ if c in "ACGTN-")
    if valid == 0: return False, "no canonical bases"
    if len(sequ) < min_len: return False, f"too short ({len(sequ)} nt)"
    return True, f"OK len={len(sequ)} valid={valid}"

# ---------------- Main processing ----------------
def make_header(q, default_id, loc):
    gene = q.get('gene',[None])[0] if q.get('gene') else None
    product = q.get('product',[None])[0] if q.get('product') else None
    prot = q.get('protein_id',[None])[0] if q.get('protein_id') else None
    codon_start = q.get('codon_start',[None])[0] if q.get('codon_start') else None
    parts = [gene or default_id]
    if prot: parts.append(prot)
    if product: parts.append(product)
    if codon_start: parts.append(f"codon_start={codon_start}")
    header_id = "|".join(str(p) for p in parts if p)
    return f"{header_id} loc={loc}"

def process_genbank_file(infile: Path, outdir_base: str, translate: bool=False, to_stop: bool=False) -> Dict:
    lines = read_file_lines(infile)
    annotation_lines, origin_lines = parse_sections(lines)
    full_seq = build_sequence(origin_lines)
    outdir = make_outdir(outdir_base, infile)

    # write full sequence fasta
    header = infile.name
    full_fasta = outdir / (infile.stem + "_full.fasta")
    with open(full_fasta, 'w') as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(full_seq), 60):
            fh.write(full_seq[i:i+60] + "\n")

    feat_block = extract_features_block(annotation_lines)
    features = parse_features(feat_block)

    # feature groups to extract
    feature_keys = ['CDS', 'mat_peptide', 'sig_peptide']
    grouped = {k: [] for k in feature_keys}
    for f in features:
        if f['key'] in feature_keys:
            grouped[f['key']].append(f)

    # write cds nt + aa
    combined_cds_fasta = outdir / (infile.stem + "_cds_seqs.fasta")
    combined_aa_fasta = outdir / (infile.stem + "_cds_proteins.fasta")
    with open(combined_cds_fasta, 'w') as fh_cds, open(combined_aa_fasta, 'w') as fh_aa:
        count = 0
        for idx, cds in enumerate(grouped.get('CDS', []), start=1):
            loc = cds.get('location','')
            seq_nt = extract_by_location(loc, full_seq)
            if not seq_nt: continue
            q = cds.get('qualifiers', {})
            hdr = make_header(q, f"CDS_{idx}", loc)
            fh_cds.write(f">{hdr}\n")
            for i in range(0, len(seq_nt), 60):
                fh_cds.write(seq_nt[i:i+60] + "\n")
            # translation - check codon_start qualifier
            cs_val = None
            if 'codon_start' in q:
                try:
                    cs_val = int(q['codon_start'][0])
                except Exception:
                    cs_val = 1
            aa = translate_nt(seq_nt, codon_start=cs_val or 1, to_stop=to_stop) if translate else ''
            if translate:
                fh_aa.write(f">{hdr}\n")
                for i in range(0, len(aa), 60):
                    fh_aa.write(aa[i:i+60] + "\n")
            count += 1

    # write mat_peptide and sig_peptide sequences if present
    for key in ('mat_peptide','sig_peptide'):
        items = grouped.get(key, [])
        if not items: continue
        outf = outdir / (infile.stem + f"_{key}_seqs.fasta")
        with open(outf, 'w') as fh:
            for idx, item in enumerate(items, start=1):
                loc = item.get('location','')
                seq_nt = extract_by_location(loc, full_seq)
                if not seq_nt: continue
                hdr = make_header(item.get('qualifiers', {}), f"{key}_{idx}", loc)
                fh.write(f">{hdr}\n")
                for i in range(0, len(seq_nt), 60):
                    fh.write(seq_nt[i:i+60] + "\n")

    ok, msg = basic_seq_qc(full_seq)
    return {
        'file': str(infile),
        'format': 'genbank',
        'n_cds_found': len(grouped.get('CDS',[])),
        'full_length': len(full_seq),
        'qc': msg,
        'outdir': str(outdir)
    }

def process_fasta_file(infile: Path, outdir_base: str) -> Dict:
    recs = []
    with open(infile, 'r', errors='ignore') as fh:
        header=None; seq_lines=[]
        for line in fh:
            line=line.rstrip('\n')
            if not line: continue
            if line.startswith('>'):
                if header: recs.append((header, ''.join(seq_lines)))
                header = line[1:].split()[0]; seq_lines=[]
            else:
                seq_lines.append(line.strip())
        if header: recs.append((header, ''.join(seq_lines)))
    outdir = make_outdir(outdir_base, infile)
    full_fasta = outdir / (infile.stem + "_sequences.fasta")
    with open(full_fasta, 'w') as fh:
        for h,s in recs:
            fh.write(f">{h}\n")
            for i in range(0,len(s),60):
                fh.write(s[i:i+60] + "\n")
    concat = ''.join(s for _,s in recs)
    ok, msg = basic_seq_qc(concat)
    return {'file': str(infile), 'format': 'fasta', 'n_seqs': len(recs), 'full_length': len(concat), 'qc': msg, 'outdir': str(outdir)}

# ---------------- CLI ----------------
def main():
    p = argparse.ArgumentParser(description="GenBank/FASTA extractor (no Biopython required) v2")
    p.add_argument('inputs', nargs='+', help='input file(s), directory(ies), or glob(s)')
    p.add_argument('--outdir', '-o', default='results', help='base output directory')
    p.add_argument('--translate', '-t', action='store_true', help='also translate CDS -> proteins for GenBank inputs')
    p.add_argument('--to-stop', action='store_true', help='stop translation at first in-frame stop codon')
    args = p.parse_args()

    inputs = resolve_inputs(args.inputs)
    if not inputs:
        print("No inputs found. Check paths/globs."); return

    aggregate=[]
    for infile in inputs:
        fmt = detect_format_by_content(infile)
        if fmt == 'genbank':
            info = process_genbank_file(infile, args.outdir, translate=args.translate, to_stop=args.to_stop)
        else:
            info = process_fasta_file(infile, args.outdir)
        aggregate.append(info)
        print(f"Processed {infile} -> {info['outdir']} (qc: {info.get('qc')})")

    results_dir = Path(args.outdir); results_dir.mkdir(exist_ok=True)
    agg_path = results_dir / 'aggregate_summary.csv'
    if PANDAS:
        df = pd.DataFrame(aggregate); df.to_csv(agg_path, index=False)
    else:
        keys = set().union(*(a.keys() for a in aggregate))
        with open(agg_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(keys)); w.writeheader()
            for a in aggregate:
                w.writerow(a)
    print(f"Aggregate summary written to {agg_path}")

if __name__ == '__main__':
    main()
