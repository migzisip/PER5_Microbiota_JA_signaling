#!/usr/bin/env python3
"""
Fetch Arabidopsis thaliana (TAIR10) JA biosynthesis/signaling proteins from UniProt.

Outputs:
  - out/ja_uniprot.faa       (FASTA, headers: >ATxG...|GENE|UniProt:ACC|pathway=...)
  - out/ja_uniprot_map.tsv   (AGI, gene, pathway, UniProt acc, reviewed?, URL, protein_name, go_mf_terms)

Requires internet access to UniProt REST.
"""

from pathlib import Path
import sys, time, csv, requests, re

OUTDIR = Path("out"); OUTDIR.mkdir(parents=True, exist_ok=True)
FASTA_OUT = OUTDIR / "ja_uniprot.faa"
MAP_OUT   = OUTDIR / "ja_uniprot_map.tsv"
BUILD_DMND = False  # set True if you want to `diamond makedb` at the end
DMND_PREFIX = "data/ref_ja_uniprot"  # -> data/ref_ja_uniprot.dmnd

# ---------------- Curated JA gene list (TAIR10) ----------------
# Columns: agi, gene, pathway
JA_GENES = [
  # LOX family (biosynthesis)
  ("AT3G45140","LOX2","JA biosynthesis - lipoxygenase"),
  ("AT1G17420","LOX3","JA biosynthesis - lipoxygenase"),
  ("AT1G72520","LOX4","JA biosynthesis - lipoxygenase"),
  ("AT1G67560","LOX6","JA biosynthesis - lipoxygenase"),
  # AOS / AOC / OPR3 / OPCL1
  ("AT5G42650","AOS","JA biosynthesis - allene oxide synthase"),
  ("AT3G25760","AOC1","JA biosynthesis - allene oxide cyclase"),
  ("AT3G25770","AOC2","JA biosynthesis - allene oxide cyclase"),
  ("AT3G25780","AOC3","JA biosynthesis - allene oxide cyclase"),
  ("AT1G13280","AOC4","JA biosynthesis - allene oxide cyclase"),
  ("AT2G06050","OPR3","JA biosynthesis - OPDA reductase"),
  ("AT1G20510","OPCL1","JA biosynthesis - OPC-8:0 CoA ligase"),
  # Peroxisomal beta-oxidation (ACX/MFP2/KAT)
  ("AT4G16760","ACX1","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT5G65110","ACX2","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT1G06290","ACX3","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT3G51840","ACX4","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT2G35690","ACX1.2","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT1G06310","ACX3.2","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT3G06860","MFP2","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT2G33150","KAT2","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT1G04710","KAT1","JA biosynthesis - peroxisomal beta-oxidation"),
  ("AT5G48880","KAT5","JA biosynthesis - peroxisomal beta-oxidation"),
  # Conjugation
  ("AT2G46370","JAR1","JA conjugation - JA?JA-Ile"),
  # Signaling (receptor, repressors, TFs)
  ("AT2G39940","COI1","JA signaling - receptor (F-box)"),
  ("AT1G19180","JAZ1","JA signaling - JAZ repressor"),
  ("AT1G74950","JAZ2","JA signaling - JAZ repressor"),
  ("AT3G17860","JAZ3","JA signaling - JAZ repressor"),
  ("AT1G48500","JAZ4","JA signaling - JAZ repressor"),
  ("AT5G13220","JAZ10","JA signaling - JAZ repressor"),
  ("AT5G20900","JAZ12","JA signaling - JAZ repressor"),
  ("AT1G32640","MYC2","JA signaling - bHLH TF"),
  ("AT5G46760","MYC3","JA signaling - bHLH TF"),
  ("AT4G17880","MYC4","JA signaling - bHLH TF"),
  # JA-Ile catabolism
  ("AT3G48520","CYP94B3","JA-Ile catabolism - CYP94"),
  ("AT2G27690","CYP94C1","JA-Ile catabolism - CYP94"),
  ("AT5G63450","CYP94B1","JA-Ile catabolism - CYP94"),
]

# ---------------- UniProt query helpers ----------------
BASE_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,reviewed,protein_name,gene_primary,gene_names,organism_name,go_f"  # includes protein_name + GO MF
HEADERS_TSV = {"Accept": "text/plain; format=tsv"}  # explicit TSV per UniProt REST

def _tsv_query(query: str, size: int = 5):
    """Run a UniProtKB search returning TSV rows as list[dict]."""
    r = requests.get(BASE_SEARCH, params={"query": query, "fields": FIELDS, "format": "tsv", "size": size}, timeout=30, headers=HEADERS_TSV)
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        vals = line.split("\t")
        row = {h: vals[i] if i < len(vals) else "" for i, h in enumerate(header)}
        rows.append(row)
    return rows

def _parse_go_mf_from_tsv_field(go_f_field: str):
    """
    Parse the UniProt 'Gene Ontology (molecular function)' TSV field into a clean '; '-joined list of MF term names.
    Typical cell content contains repeated patterns like 'GO:xxxx; F:term; evidence...'.
    We extract the 'F:term' fragments.
    """
    if not go_f_field:
        return ""
    # Split by ' | ' or ' ; ' sequences, then capture substrings after 'F:'
    # Keep robust for cases with commas/pipes/semicolons.
    terms = []
    # First split into records using ' | ' as primary separator if present
    records = re.split(r'\s*\|\s*', go_f_field)
    for rec in records:
        # Extract all F:term occurrences in the record
        for m in re.finditer(r'\bF:([^;|]+)', rec):
            term = m.group(1).strip()
            if term and term not in terms:
                terms.append(term)
    # Fallback: if nothing captured, try simple semicolon split and look for 'F:'
    if not terms:
        for part in re.split(r'\s*;\s*', go_f_field):
            m = re.search(r'\bF:([^;|]+)', part)
            if m:
                term = m.group(1).strip()
                if term and term not in terms:
                    terms.append(term)
    return "; ".join(terms)

def _fetch_entry_json(acc: str):
    """Fetch full UniProtKB entry JSON for robust fallbacks."""
    url = f"https://rest.uniprot.org/uniprotkb/{acc}"
    r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

def _protein_name_from_json(j):
    """Recommended protein full name, with safe fallbacks."""
    try:
        return j["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except Exception:
        # Try submittedName or alternativeNames if recommendedName missing
        try:
            return j["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
        except Exception:
            try:
                return j["proteinDescription"]["alternativeNames"][0]["fullName"]["value"]
            except Exception:
                return ""

def _go_mf_from_json(j):
    """Collect GO Molecular Function (aspect 'F') terms from JSON cross-references."""
    terms = []
    try:
        for x in j.get("uniProtKBCrossReferences", []):
            if x.get("database") != "GO":
                continue
            props = {p.get("key"): p.get("value") for p in x.get("properties", []) if "key" in p and "value" in p}
            # Aspect 'F' = molecular function
            aspect = props.get("Aspect") or props.get("GoAspect") or ""
            if aspect == "F":
                # Prefer human-readable 'Term' if present; else leave blank
                t = props.get("Term") or ""
                if t and t not in terms:
                    terms.append(t)
    except Exception:
        pass
    return "; ".join(terms)

def query_uniprot_single(agi: str, gene: str):
    """
    Resolve a single UniProtKB accession for the given Arabidopsis gene.
    Returns: (acc, reviewed_bool, protein_name, go_mf_terms)
    Strategy:
      1) reviewed Swiss-Prot by exact gene name
      2) reviewed by AGI locus name
      3) allow unreviewed
    """
    # 1) reviewed, by gene_exact
    q1 = f'(organism_id:3702) AND (gene_exact:{gene})'
    rows = _tsv_query(q1, size=5)
    for r in rows:
        if r.get("Status","").lower() == "reviewed":
            acc = r.get("Entry","")
            prot = r.get("Protein names","")
            mf   = _parse_go_mf_from_tsv_field(r.get("Gene Ontology (molecular function)",""))
            if (not prot) or (not mf):
                try:
                    j = _fetch_entry_json(acc)
                    prot = prot or _protein_name_from_json(j)
                    mf   = mf   or _go_mf_from_json(j)
                except Exception:
                    pass
            return acc, True, prot, mf

    # 2) reviewed, by AGI locus name in gene names
    q2 = f'(organism_id:3702) AND (gene:{agi})'
    rows = _tsv_query(q2, size=5)
    for r in rows:
        if r.get("Status","").lower() == "reviewed":
            acc = r.get("Entry","")
            prot = r.get("Protein names","")
            mf   = _parse_go_mf_from_tsv_field(r.get("Gene Ontology (molecular function)",""))
            if (not prot) or (not mf):
                try:
                    j = _fetch_entry_json(acc)
                    prot = prot or _protein_name_from_json(j)
                    mf   = mf   or _go_mf_from_json(j)
                except Exception:
                    pass
            return acc, True, prot, mf

    # 3) fallback: allow unreviewed
    q3 = f'(organism_id:3702) AND (gene_exact:{gene} OR gene:{agi})'
    rows = _tsv_query(q3, size=1)
    if rows:
        r = rows[0]
        acc = r.get("Entry","")
        reviewed = r.get("Status","").lower() == "reviewed"
        prot = r.get("Protein names","")
        mf   = _parse_go_mf_from_tsv_field(r.get("Gene Ontology (molecular function)",""))
        if (not prot) or (not mf):
            try:
                j = _fetch_entry_json(acc)
                prot = prot or _protein_name_from_json(j)
                mf   = mf   or _go_mf_from_json(j)
            except Exception:
                pass
        return acc, reviewed, prot, mf

    return None, None, "", ""

def fetch_fasta(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text

# ---------------- Main ----------------
def main():
    mapping_rows = []
    with open(FASTA_OUT, "w") as fo:
        for agi, gene, pathway in JA_GENES:
            try:
                acc, reviewed, protein_name, go_mf = query_uniprot_single(agi, gene)
            except Exception as e:
                sys.stderr.write(f"[{agi} {gene}] UniProt query failed: {e}\n")
                acc, reviewed, protein_name, go_mf = None, None, "", ""

            if not acc:
                sys.stderr.write(f"[{agi} {gene}] No UniProt match found.\n")
                mapping_rows.append([agi, gene, pathway, "", "", "", "", ""])
                continue

            # fetch FASTA
            try:
                fa = fetch_fasta(acc)
            except Exception as e:
                sys.stderr.write(f"[{agi} {gene}] FASTA fetch failed for {acc}: {e}\n")
                mapping_rows.append([agi, gene, pathway, acc, "reviewed" if reviewed else "unreviewed",
                                     f"https://www.uniprot.org/uniprotkb/{acc}", protein_name, go_mf])
                continue

            # rewrite header to our preferred form (unchanged)
            lines = fa.strip().splitlines()
            seq = "\n".join(lines[1:])
            new_header = f">{agi}|{gene}|UniProt:{acc}|pathway={pathway}"
            fo.write(new_header + "\n" + seq + "\n")

            mapping_rows.append([
                agi, gene, pathway, acc,
                "reviewed" if reviewed else "unreviewed",
                f"https://www.uniprot.org/uniprotkb/{acc}",
                protein_name, go_mf
            ])
            print(f"[OK] {agi} {gene} -> {acc} ({'reviewed' if reviewed else 'unreviewed'}) | protein='{protein_name}' | MF='{go_mf}'")
            time.sleep(0.25)  # be polite

    # write map (with new columns)
    with open(MAP_OUT, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["agi","gene","pathway","uniprot_acc","review_status","uniprot_url","protein_name","go_mf_terms"])
        for row in mapping_rows:
            w.writerow(row)

    print("\nWrote:")
    print("  FASTA:", FASTA_OUT)
    print("  MAP  :", MAP_OUT)

    if BUILD_DMND:
        print("\nBuilding DIAMOND DB...")
        import subprocess
        subprocess.check_call(["diamond","makedb","--in",str(FASTA_OUT),"-d",DMND_PREFIX])
        print("  ->", DMND_PREFIX + ".dmnd")

if __name__ == "__main__":
    main()
