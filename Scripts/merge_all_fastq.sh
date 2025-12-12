#!/usr/bin/env bash
set -euo pipefail

manifest="merge-manifest.tsv"
root="/mnt/c/Users/user/Desktop/symcom/miseq/extra_1G/FLASH2/flash2_results/usearch/"
pattern="*.extendedFrags.filtered.fastq"   # adjust if your suffix is different

# header (tabs)
printf "sample-id\tabsolute-filepath\tdirection\n" > "$manifest"

i=1
find "$root" -type f -name "$pattern" -print0 \
  | sort -z \
  | while IFS= read -r -d '' f; do
      abs="$(readlink -f "$f")"
      printf "%d\t%s\tforward\n" "$i" "$abs" >> "$manifest"
      i=$((i+1))
    done
