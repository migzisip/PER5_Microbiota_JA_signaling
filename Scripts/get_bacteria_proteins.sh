#!/usr/bin/env bash
# safer, more portable flags
set -e
set -u
set -o pipefail 2>/dev/null || true

# --- Config (edit paths if you like) ---
ASSEMBLY_LIST="${1:-assemblies.txt}"
OUTDIR="${OUTDIR:-bacteria_genomes}"
ZIP="${ZIPFILE:-bacteria.zip}"
PROT_FASTA="${PROT_FASTA:-data/ref_bacteria_proteins.faa}"  # combined output

mkdir -p data "$OUTDIR"

# --- Sanity checks ---
if [[ ! -s "$ASSEMBLY_LIST" ]]; then
  echo "ERROR: assembly list not found: $ASSEMBLY_LIST"
  echo "Create a file with one accession per line (e.g., GCA_001427525.1)."
  exit 1
fi

# --- Ensure datasets CLI ---
choose_datasets () {
  if [[ -x "./datasets" ]]; then
    echo "./datasets"; return
  fi
  if command -v datasets >/dev/null 2>&1; then
    echo "datasets"; return
  fi
  echo "[info] downloading NCBI datasets binary (no sudo) ..."
  wget -q https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets -O ./datasets
  chmod +x ./datasets
  echo "./datasets"
}
DATASETS_BIN="$(choose_datasets)"
echo "[ok] using datasets at: $DATASETS_BIN"

# Detect CLI flavor
HELP="$($DATASETS_BIN download genome accession --help 2>/dev/null || true)"
if   echo "$HELP" | grep -q -- "--input-file"; then  SYNTAX="old_input_file"
elif echo "$HELP" | grep -q -- "--inputfile" ; then  SYNTAX="old_inputfile"
else                                                SYNTAX="v2_positional"
fi
echo "[ok] datasets syntax detected: $SYNTAX"

# --- Download assemblies (proteins included) ---
echo "[step] downloading assemblies listed in $ASSEMBLY_LIST ..."
rm -f "$ZIP"
set +e
case "$SYNTAX" in
  v2_positional)
    "$DATASETS_BIN" download genome accession "$ASSEMBLY_LIST" --include protein --filename "$ZIP"
    ;;
  old_input_file)
    "$DATASETS_BIN" download genome accession --input-file "$ASSEMBLY_LIST" --include protein --filename "$ZIP"
    ;;
  old_inputfile)
    "$DATASETS_BIN" download genome accession --inputfile "$ASSEMBLY_LIST" --include protein --filename "$ZIP"
    ;;
esac
RC=$?
set -e
if [[ $RC -ne 0 || ! -s "$ZIP" ]]; then
  echo "ERROR: download failed (no $ZIP created)."
  echo "Try: $DATASETS_BIN download genome accession --help"
  exit 1
fi
echo "[ok] downloaded zip: $ZIP"

# --- Unpack (unzip or python fallback) ---
echo "[step] extracting $ZIP -> $OUTDIR ..."
if command -v unzip >/dev/null 2>&1; then
  unzip -q "$ZIP" -d "$OUTDIR"
else
  python - "$ZIP" "$OUTDIR" <<'PY'
import sys, zipfile
zip_path, outdir = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(outdir)
print("extracted with Python")
PY
fi
echo "[ok] extracted."

# --- Collect protein FASTAs (concatenate into one file) ---
echo "[step] collecting protein FASTAs ..."
> "$PROT_FASTA"  # truncate/create

# Prefer plain .faa first
FOUND_PLAIN=0
while IFS= read -r -d '' f; do
  cat "$f" >> "$PROT_FASTA"
  FOUND_PLAIN=1
done < <(find "$OUTDIR" -type f -name "protein.faa" -size +0 -print0)

# If none found, try .faa.gz
if [[ $FOUND_PLAIN -eq 0 ]]; then
  echo "[info] no plain protein.faa found; trying .faa.gz"
  while IFS= read -r -d '' gz; do
    gunzip -c "$gz" >> "$PROT_FASTA"
  done < <(find "$OUTDIR" -type f -name "protein.faa.gz" -print0)
fi

if [[ ! -s "$PROT_FASTA" ]]; then
  echo "ERROR: could not find any protein FASTA in the downloaded dataset."
  echo "Check contents under $OUTDIR/ncbi_dataset/data/*/ ."
  exit 1
fi

echo "[ok] combined proteins: $PROT_FASTA"
echo "[info] $(grep -c '^>' "$PROT_FASTA") sequences combined."

# --- Optional cleanup of headers & duplicates (if seqkit available) ---
if command -v seqkit >/dev/null 2>&1; then
  echo "[step] seqkit cleanup (rmdup + header trim)"
  seqkit rmdup -s "$PROT_FASTA" | seqkit replace -p '^>(\\S+).*' -r '>\\$1' > "${PROT_FASTA%.faa}.clean.faa"
  PROT_FASTA="${PROT_FASTA%.faa}.clean.faa"
  echo "[ok] cleaned proteins: $PROT_FASTA"
  echo "[info] $(grep -c '^>' "$PROT_FASTA") sequences after cleanup."
fi

echo "Done. FASTAs are under: $OUTDIR/ncbi_dataset/data/*/ (per-assembly) and combined at: $PROT_FASTA"