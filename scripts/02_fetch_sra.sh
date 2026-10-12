#!/usr/bin/env bash
# Download one SRA run as a .sra file and verify it against NCBI's own manifest.
#
# Why not `prefetch`? sra-tools resolves accessions through a name-resolver
# service that is not reachable from every network; when it is blocked, prefetch
# and fasterq-dump fail with "Failed to call external services" even for local
# files. We therefore ask the NCBI SDL endpoint for the delivery URL + md5
# ourselves, pull the object with ranged parallel HTTP, and run sra-tools fully
# offline (see 03_reconstruct_locus.sh, which sets NCBI_SETTINGS).
#
# Usage: scripts/02_fetch_sra.sh SRR3948161
set -euo pipefail
cd "$(dirname "$0")/.."
source config/params.sh
ACC="${1:?usage: 02_fetch_sra.sh <SRR accession>}"
mkdir -p "$WORK/sra" "$RESULTS/$ACC"
SRA="$WORK/sra/${ACC}.sra"

read -r SIZE MD5 URL < <(curl -sS "https://trace.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=${ACC}" \
  | python3 -c "import sys,json;f=json.load(sys.stdin)['result'][0]['files'][0];print(f['size'],f['md5'],f['locations'][0]['link'])")
echo "[02] $ACC: $SIZE bytes, md5 $MD5"
printf 'accession\t%s\nsize_bytes\t%s\nmd5_expected\t%s\nurl\t%s\n' "$ACC" "$SIZE" "$MD5" "$URL" \
  > "$RESULTS/$ACC/sra_manifest.txt"

# wc -c, unlike BSD stat -f%z, follows symlinks -- the .sra is often a link
# into shared scratch.
if [[ -s "$SRA" && "$(wc -c < "$SRA" | tr -d ' ')" == "$SIZE" ]]; then
  echo "[02] already downloaded, skipping"
else
  CH=$(( (SIZE + DL_STREAMS - 1) / DL_STREAMS ))
  echo "[02] downloading in $DL_STREAMS ranged streams"
  for i in $(seq 0 $((DL_STREAMS-1))); do
    s=$(( i * CH )); e=$(( s + CH - 1 )); [[ $e -ge $SIZE ]] && e=$(( SIZE - 1 ))
    curl -sS --retry 5 --retry-all-errors -r "$s-$e" -o "$WORK/sra/${ACC}.part$i" "$URL" &
  done
  wait
  cat "$WORK/sra/${ACC}".part* > "$SRA" && rm -f "$WORK/sra/${ACC}".part*
fi

OBS=$(md5 -q "$SRA" 2>/dev/null || md5sum "$SRA" | cut -d' ' -f1)
printf 'md5_observed\t%s\n' "$OBS" >> "$RESULTS/$ACC/sra_manifest.txt"
if [[ "$OBS" != "$MD5" ]]; then echo "[02] ERROR: md5 mismatch ($OBS != $MD5)"; exit 1; fi
echo "[02] md5 verified against the NCBI manifest"
