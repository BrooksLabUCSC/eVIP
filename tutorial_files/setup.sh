#!/bin/bash

"""
script to download RNF43 kallisto outputs from GEO
and properly format for the tutorial
"""

cd tutorial_files
mkdir RNF43_kallisto_outputs
cd RNF43_kallisto_outputs

#download kallisto files
wget -i ../ftp_links.txt

#unzip kallisto files
gunzip *.tsv.gz

#moving abundance file to named dir
for x in *abundance.tsv; do
  mkdir "${x%.*}" && mv "$x" "${x%.*}"
done

# renaming directories
for x in *_abundance; do
  count=$(echo $x | tr -cd '_' | wc -c)

  # only for GFP samples
  if [ "$count" -eq 3 ]; then
    A="$(cut -d'_' -f2 <<< $x )"
    B="$(cut -d'_' -f3 <<< $x )"
    C="$A""_""$B"
    mv "$x" "$C"

  else
    A="$(cut -d'_' -f2 <<< $x )"
    B="$(cut -d'_' -f3 <<< $x )"
    C="$(cut -d'_' -f4 <<< $x )"
    D="$A""_""$B""_""$C"
    mv "$x" "$D"

  fi

done

#renaming files to abundance.tsv
find . -type f -name "*_abundance.tsv" -not -name "abundance.tsv" \
           -execdir mv -v {} abundance.tsv \;
