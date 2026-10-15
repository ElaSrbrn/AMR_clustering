#!/bin/bash
# Builds clean S1-based read folders for 10x, 20x, 40x and 60x (symlinks) and writes the
# sample sheet for the Hybracter depth series. Runs in seconds on the login node, no SLURM needed.
#
# Usage (from revision_scripts, where s1_ids_all.txt lives):
#   bash build_depth_inputs.sh
#
# Output:
#   ../reads_{10,20,40,60}x_S1/            symlinks to the read files
#   s1_min_chrom.tsv                        isolate, species, minimum chromosome length used by Hybracter
#   hybracter_depth_series.csv              sample_name,fastq,min_chrom_length (no header) (one row per isolate x depth)
#   build_depth_inputs.report               MISSING / MULTIPLE lines to check by hand

set -uo pipefail
NEWCRE=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/new_cre
BASE=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project
IDS=s1_ids_all.txt
[ -f "$IDS" ] || { echo "ERROR: $IDS not found; run from revision_scripts"; exit 1; }

# Minimum chromosome length per isolate (Hybracter -c): a lower bound set below the typical
# chromosome size of the genus/species (species from S1, Pathogenwatch). Check these values.
: > s1_min_chrom.tsv
while IFS=$'\t' read -r id sp; do
  case "$sp" in
    "Escherichia coli"|"Shigella dysenteriae")                       len=4000000 ;;
    "Klebsiella pneumoniae"|"K. pneumoniae"|"Klebsiella variicola"|\
    "Klebsiella quasipneumoniae"|"Klebsiella aerogenes")              len=4500000 ;;
    "Klebsiella michiganensis"|"Raoultella planticola")               len=5000000 ;;
    Citrobacter*|"C. freundii"|Enterobacter*|"E. roggenkampi"|\
    "Hafnia alvei")                                                   len=4000000 ;;
    "Serratia marcescens")                                            len=4300000 ;;
    "Providencia stuartii")                                           len=3700000 ;;
    "Proteus mirabilis")                                              len=3500000 ;;
    *) echo "WARN: no length rule for $id ($sp), using 3500000" >&2;  len=3500000 ;;
  esac
  echo -e "$id\t$sp\t$len" >> s1_min_chrom.tsv
done <<'EOF'
3989	Citrobacter freundii
4704	Citrobacter freundii
4838	Klebsiella pneumoniae
4920	Klebsiella pneumoniae
5043	Proteus mirabilis
5244	Klebsiella pneumoniae
5305	Providencia stuartii
5371	Enterobacter hormaechei
5399	Escherichia coli
5425	Citrobacter freundii
5590	Enterobacter hormaechei
5607	Enterobacter ludwigii
5668	Citrobacter freundii
5860	Klebsiella pneumoniae
5861	Klebsiella aerogenes
5871	Klebsiella variicola
5969	Citrobacter freundii
5970	Citrobacter freundii
5976	Enterobacter hormaechei
6009	Escherichia coli
6014	K. pneumoniae
6035	Klebsiella pneumoniae
6036	Enterobacter hormaechei
6056	Proteus mirabilis
6057	Escherichia coli
6067	Klebsiella michiganensis
6068	Escherichia coli
6233	Escherichia coli
6243	Klebsiella variicola
6244	Klebsiella pneumoniae
6247	Escherichia coli
6287	Klebsiella quasipneumoniae
6312	Escherichia coli
6365	Klebsiella variicola
6378	Klebsiella pneumoniae
6394	Klebsiella pneumoniae
6403	Escherichia coli
6433	Serratia marcescens
6454	Klebsiella variicola
6455	Klebsiella pneumoniae
6466	Escherichia coli
6472	Citrobacter freundii
6482	Enterobacter hormaechei
6500	Klebsiella pneumoniae
6520	Klebsiella pneumoniae
6567	Raoultella planticola
6571	Citrobacter freundii
6601	Citrobacter freundii
6623	Enterobacter hormaechei
6668	Klebsiella aerogenes
6669	Klebsiella pneumoniae
6670	Citrobacter freundii
6701	Klebsiella michiganensis
6711	Citrobacter europaeus
6712	Escherichia coli
6726	Klebsiella pneumoniae
6735	Citrobacter freundii
6737	Klebsiella pneumoniae
6738	Citrobacter freundii
6748	Escherichia coli
6774	Klebsiella pneumoniae
6775	Klebsiella aerogenes
6776	Citrobacter freundii
6778	Klebsiella aerogenes
6779	Klebsiella aerogenes
6795	Citrobacter amalonaticus
6797	Citrobacter freundii
6804	Enterobacter hormaechei
6813	Hafnia alvei
6817	Klebsiella pneumoniae
6838	Escherichia coli
6855	Klebsiella michiganensis
6856	Enterobacter hormaechei
6888	Klebsiella pneumoniae
6900	Citrobacter portucalensis
6901	Citrobacter freundii
6951	Escherichia coli
6970	Citrobacter freundii
6989	Klebsiella pneumoniae
7000	Proteus mirabilis
7067	Escherichia coli
7077	Serratia marcescens
7083	Escherichia coli
7084	E. roggenkampi
7088	Citrobacter koseri
7089	Klebsiella pneumoniae
7090	Citrobacter koseri
7101	Enterobacter bugandensis
7136	Klebsiella pneumoniae
7176	Klebsiella pneumoniae
7181	Klebsiella pneumoniae
7213	Proteus mirabilis
7221	Escherichia coli
7243	Enterobacter hormaechei
7273	Escherichia coli
7285	Citrobacter braakii
7310	Escherichia coli
7316	Klebsiella pneumoniae
7320	Klebsiella pneumoniae
7365	Klebsiella aerogenes
7370	Shigella dysenteriae
7387	Escherichia coli
7418	Citrobacter braakii
7419	Enterobacter hormaechei
7454	Escherichia coli
7475	Klebsiella pneumoniae
7476	Klebsiella pneumoniae
7488	Escherichia coli
7503	C. freundii
7521	Klebsiella pneumoniae
EOF

: > hybracter_depth_series.csv
: > build_depth_inputs.report
shopt -s nullglob
for d in 10 20 40 60; do
  DEST=$BASE/reads_${d}x_S1
  mkdir -p "$DEST"
  n=0
  while read -r id; do
    f=( $NEWCRE/carblot_pos/${d}x/isolate_${id}_*${d}x.fastq* $NEWCRE/carblot_neg/${d}x/isolate_${id}_*${d}x.fastq* )
    if   [ ${#f[@]} -eq 0 ]; then echo "MISSING ${d}x $id" >> build_depth_inputs.report; continue
    elif [ ${#f[@]} -gt 1 ]; then echo "MULTIPLE ${d}x $id: ${f[*]}" >> build_depth_inputs.report; continue
    fi
    ln -sf "${f[0]}" "$DEST/"
    len=$(awk -F'\t' -v i="$id" '$1==i {print $3}' s1_min_chrom.tsv)
    echo "${id}_${d}x,$DEST/$(basename "${f[0]}"),$len" >> hybracter_depth_series.csv
    n=$((n+1))
  done < "$IDS"
  echo "${d}x: $n isolates linked in $DEST"
done
shopt -u nullglob

echo "Rows in hybracter_depth_series.csv: $(wc -l < hybracter_depth_series.csv)"
echo "Problems (expected: MISSING 60x for 5371 5668 6804 7418):"
cat build_depth_inputs.report
