#!/bin/zsh
setopt nullglob

ROOT=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/smoking-networks/analysis-versions/version001
DIR=results/003
cd $ROOT/$DIR || exit 1

# echo "Making flowchart"
# for f in *.dot; do
#   echo "flowing $f"
#   dot -Tsvg $f -o "${f%.dot}.svg" 
# done

for f in *.pdf; do
  [[ $f == crpd-* ]] && continue
  echo "cropping $f" 
  pdfcrop --quiet "$f" "crpd-$f"
done

for f in *.png *.jpg *.jpeg; do
  [[ $f == crpd-* ]] && continue
  echo "trimming $f"
  magick "$f" -fuzz 2% -trim +repage "crpd-$f" 
done

