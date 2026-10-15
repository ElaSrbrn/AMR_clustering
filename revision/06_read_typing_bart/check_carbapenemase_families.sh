#!/usr/bin/env bash
# Checks whether carbapenemase alleles share a catalogue gene family with non-carbapenemase alleles
# (e.g. blaOXA-48-like with blaOXA-1), which would let a co-carried allele mask the carbapenemase
# under the rank-1 rule. Output: catalog_family_check.txt
CAT=${1:-$HOME/bart_schemes/ReferenceGeneCatalog.txt}
OUT=catalog_family_check.txt
{
echo "### 1. Catalogue: gene_family and subclass of selected alleles"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
  $h["allele"]~/^bla(OXA-(1|9|10|48|181|232|244)|KPC-2|GES-(1|5|28))$/{print $h["allele"],$h["gene_family"],$h["subclass"]}' "$CAT" | column -t

echo; echo "### 2. Catalogue: gene families containing both carbapenemase and other alleles"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["type"]=="AMR"{
  f=$h["gene_family"]; if($h["subclass"]~/CARBAPENEM/) c[f]++; else o[f]++}
  END{for(f in c) if(o[f]) print f, c[f]" carbapenemase alleles", o[f]" other alleles"}' "$CAT" | sort | column -t

echo; echo "### 3. Samples where a carbapenemase hit (rel_depth >= 0.2) exists but is not rank 1 in its family"
awk -F'\t' 'NR>1{k=$1"\t"$2"\t"$4; if($14==1) top[k]=$5" ["$9"] rel "$13;
  if($9~/CARBAPENEM/){ if($14==1) ok[k]=1; else if($13>=0.2) c[k]=c[k]" "$5"(rank "$14", rel "$13")"}}
  END{for(k in c) if(!ok[k]) print k"\trank 1: "top[k]"\tcarbapenemase:"c[k]}' bart_results_amr.tsv | sort | column -t -s$'\t'

echo; echo "### 4. Lowest rel_depth of rank-1 carbapenemase hits (true-positive side of the threshold window)"
awk -F'\t' 'NR>1 && $9~/CARBAPENEM/ && $14==1{print $13"\t"$1"\t"$2"\t"$5}' bart_results_amr.tsv | sort -g | head -10 | column -t -s$'\t'

echo; echo "### 5. Carbapenemase calls that differ from the gold standard (reads and assembly)"
column -t -s$'\t' summary_carbapenemase_discordant.tsv
} > "$OUT" 2>&1
echo "written: $OUT ($(wc -l < "$OUT") lines)"
