#!/usr/bin/env bash
# Build the PDF report from the Markdown source.
# Usage: bash build.sh   (needs pandoc >= 3 and lualatex)
set -euo pipefail
cd "$(dirname "$0")"
SRC="../reports/Neue Anwendungsfelder für GZL.md"
# drop the H1 title (it is on the title page), drop manual section numbers, use PDF figures
sed -e '1{/^# /d}' \
    -e 's/^\(##\+\) [0-9]\+\(\.[0-9]\+\)*\.\? /\1 /' \
    -e 's#(\.\./verification/figures/\([^)]*\)\.png)#(../verification/figures/\1.pdf)#g' \
    "$SRC" | awk 'BEGIN{prev=""} { if ($0 ~ /^([-*] |[0-9]+\. )/ && prev != "" && prev !~ /^([-*] |[0-9]+\. |  )/ && prev !~ /^\|/) print ""; print; prev=$0 }' > body.md
pandoc body.md -f markdown-auto_identifiers+pipe_tables-implicit_figures+implicit_figures \
    -t latex --lua-filter=report.lua --shift-heading-level-by=-1 --no-highlight --wrap=preserve -o body.tex
latexmk -lualatex -interaction=nonstopmode -halt-on-error main.tex > build.log 2>&1 || { tail -40 build.log; exit 1; }
cp main.pdf "../reports/Neue Anwendungsfelder für GZL.pdf"
echo "PDF written: research/reports/Neue Anwendungsfelder für GZL.pdf"
