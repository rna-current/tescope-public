#!/usr/bin/env bash
# Download the two public ENCODE files used by this example (~1 GB)
# into ./data and verify their checksums (as published by ENCODE).
# ENCODE data use policy: free download, analysis and publication
# without restriction (https://www.encodeproject.org/help/citing-encode/).
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p data; cd data
fetch() { [ -s "$2" ] || curl -L --fail -o "$2" "https://www.encodeproject.org/files/$1/@@download/$2"; }
fetch ENCFF601JGK ENCFF601JGK.bigWig     # H3K9me3 fold change over control, GRCh38, K562 (ENCSR000APE)
fetch ENCFF963GZJ ENCFF963GZJ.bed.gz     # replicated peaks, GRCh38
md5sum -c <<'SUMS'
a893ae9a929c36457c6386f4a6475a84  ENCFF601JGK.bigWig
8e71cfdb547cdb5057d99d12d7be3ba6  ENCFF963GZJ.bed.gz
SUMS
