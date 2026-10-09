#!/usr/bin/env bash
# Download the raw public catalogues (about 190 MB) into tools/raw/.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p raw; cd raw
get(){ echo "-> $2"; curl -fL --retry 3 -o "$2" "$1"; }
get https://raw.githubusercontent.com/astronexus/HYG-Database/main/hyg/CURRENT/hygdata_v41.csv hyg.csv
get https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv NGC.csv
get https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_all.csv lv_dwarf_all.csv
get https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_harris.csv lv_gc_harris.csv
get https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_mw_new.csv lv_gc_mw_new.csv
get https://raw.githubusercontent.com/astroML/astroML-data/main/datasets/SDSSspecgalsDR8_1.fit.gz sdss1.fit.gz
get https://raw.githubusercontent.com/astroML/astroML-data/main/datasets/SDSSspecgalsDR8_2.fit.gz sdss2.fit.gz
get https://raw.githubusercontent.com/astroML/astroML-data/main/datasets/nsa_v0_1_2_reduced.npy nsa.npy
echo "Done. Now run: python3 tools/build_data.py"
