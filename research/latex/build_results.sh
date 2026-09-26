#!/usr/bin/env bash
# Build the PDF of the focused results document.  Usage: bash build_results.sh
set -euo pipefail
cd "$(dirname "$0")"
SRC="../reports/Dreikörper-Gitterenergie – Ergebnisse.md"
sed -e '1{/^# /d}' \
    -e 's/^\(##\+\) [0-9]\+\(\.[0-9]\+\)*\.\? /\1 /' \
    -e 's#(\.\./verification/figures/\([^)]*\)\.png)#(../verification/figures/\1.pdf)#g' \
    -e 's#(\.\./extended/figures/\([^)]*\)\.png)#(../extended/figures/\1.pdf)#g' \
    "$SRC" | awk 'BEGIN{prev=""} { if ($0 ~ /^([-*] |[0-9]+\. )/ && prev != "" && prev !~ /^([-*] |[0-9]+\. |  )/ && prev !~ /^\|/) print ""; print; prev=$0 }' > body_results.md
pandoc body_results.md -f markdown-auto_identifiers+pipe_tables -t latex --lua-filter=report.lua \
    --shift-heading-level-by=-1 --no-highlight --wrap=preserve -o body_results.tex
latexmk -lualatex -interaction=nonstopmode -halt-on-error main_results.tex > build_results.log 2>&1 || { tail -40 build_results.log; exit 1; }
cp main_results.pdf "../reports/Dreikörper-Gitterenergie – Ergebnisse.pdf"
echo "PDF written: research/reports/Dreikörper-Gitterenergie – Ergebnisse.pdf"
