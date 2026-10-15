#!/bin/bash
#SBATCH --job-name=ref_check
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs_ref_check/%x_%A_%a.out
#
# Phase 1 check of the new Autocycler (+Plassembler) references, one isolate per array task.
# Per isolate it writes refcheck/<isolate>/ with:
#   contigs_new.tsv / contigs_old.tsv   contig, length, circular flag (new = with Plassembler, old = without)
#   depth.tsv                           mean read depth per contig (60x reads mapped back to the new reference)
#   amrfinder.tsv                       AMRFinderPlus on the new reference
#   mobtyper.tsv                        MOB-typer on every contig except the longest (chromosome)
# Afterwards, merge everything into one table with:  MODE=merge bash reference_check.sh
#
# Usage (from revision_scripts):
#   mkdir -p logs_ref_check refcheck
#   LIST=isolate_dirs_S1.txt sbatch --array=1 reference_check.sh            # test one isolate
#   LIST=isolate_dirs_S1.txt sbatch --array=2-106%20 reference_check.sh
#   MODE=merge bash reference_check.sh
#
# Set these to the environments that contain the tools (check with: micromamba env list).
ENV_MAP=${ENV_MAP:-autocycler}     # needs minimap2, samtools, seqkit
ENV_AMR=${ENV_AMR:-amrfinder}      # needs amrfinder (database installed)
ENV_MOB=${ENV_MOB:-mobsuite}       # needs mob_typer

set -uo pipefail
LIST=${LIST:-isolate_dirs_S1.txt}
WD=$(pwd)
OUTROOT=$WD/refcheck

if [ "${MODE:-run}" = "merge" ]; then
python3 - "$OUTROOT" <<'PY'
import sys, os, csv, glob
root = sys.argv[1]
def rd(p):
    if not os.path.exists(p) or os.path.getsize(p) == 0: return []
    with open(p) as f: return list(csv.DictReader(f, delimiter='\t'))
def col(row, *names):
    for n in names:
        if n in row: return row[n]
    return ''
out = open(os.path.join(root, 'reference_check_summary.tsv'), 'w')
out.write('\t'.join(['isolate','contig','length','circular','mean_depth','depth_ratio_to_chrom',
                     'carbapenemase_genes','all_amr_genes','rep_type','relaxase','predicted_mobility',
                     'in_old_reference_same_length']) + '\n')
for d in sorted(glob.glob(os.path.join(root, '*'))):
    if not os.path.isdir(d): continue
    iso = os.path.basename(d)
    new = rd(os.path.join(d, 'contigs_new.tsv'))
    old = rd(os.path.join(d, 'contigs_old.tsv'))
    old_lens = [int(r['length']) for r in old]
    dep = {r['#rname']: float(r['meandepth']) for r in rd(os.path.join(d, 'depth.tsv'))}
    chrom = max(new, key=lambda r: int(r['length']))['contig'] if new else None
    cdep = dep.get(chrom, 0) or float('nan')
    amr = {}
    carb = {}
    for r in rd(os.path.join(d, 'amrfinder.tsv')):
        c = col(r, 'Contig id', 'Contig ID')
        g = col(r, 'Element symbol', 'Gene symbol')
        sub = col(r, 'Subclass', 'Element subclass')
        amr.setdefault(c, []).append(g)
        if 'CARBAPENEM' in sub.upper(): carb.setdefault(c, []).append(g)
    mob = {}
    for r in rd(os.path.join(d, 'mobtyper.tsv')):
        sid = col(r, 'sample_id').split()[0]
        sid = sid.split(':')[-1]
        mob[sid] = r
    for r in new:
        c = r['contig']; L = int(r['length'])
        m = mob.get(c, {})
        same = 'yes' if any(abs(L - o) <= max(2, 0.001 * L) for o in old_lens) else 'NO'
        dp = dep.get(c, float('nan'))
        out.write('\t'.join(map(str, [iso, c, L, r['circular'], f'{dp:.1f}', f'{dp/cdep:.2f}' if cdep == cdep else 'NA',
            ','.join(carb.get(c, [])), ','.join(amr.get(c, [])),
            col(m, 'rep_type(s)'), col(m, 'relaxase_type(s)'), col(m, 'predicted_mobility'), same])) + '\n')
out.close()
print('Written:', os.path.join(root, 'reference_check_summary.tsv'))
PY
exit 0
fi

eval "$(micromamba shell hook --shell bash)"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WD/$LIST")
NAME=$(basename "$DIR")
ISO=$(echo "$NAME" | cut -d_ -f2)
OUT=$OUTROOT/$ISO
mkdir -p "$OUT"
NEWREF=$DIR/autocycler_out/consensus_assembly.fasta
OLDREF=$DIR/autocycler_out_noplassembler/consensus_assembly.fasta
READS=$(dirname "$(dirname "$DIR")")/$NAME.fastq
echo "Isolate $ISO  dir $DIR"
[ -s "$NEWREF" ] || { echo -e "STATUS\t$ISO\tFAIL_no_reference"; exit 1; }
[ -s "$READS" ]  || { echo -e "STATUS\t$ISO\tFAIL_no_reads ($READS)"; exit 1; }

# Contig table from FASTA headers (Autocycler writes e.g. ">1 length=5333834 circular=true")
contig_table() {
    awk 'BEGIN{OFS="\t"; print "contig","length","circular"}
         /^>/{ if (n) print n,l,c; n=substr($1,2); l=0; c="NA";
               for(i=2;i<=NF;i++) if ($i ~ /^circular=/) {split($i,a,"="); c=a[2]} ; next }
         { l+=length($0) } END{ if (n) print n,l,c }' "$1"
}
contig_table "$NEWREF" > "$OUT/contigs_new.tsv"
if [ -s "$OLDREF" ]; then contig_table "$OLDREF" > "$OUT/contigs_old.tsv"; else : > "$OUT/contigs_old.tsv"; fi

# Read depth per contig
micromamba activate "$ENV_MAP"
minimap2 -ax map-ont -t 8 "$NEWREF" "$READS" 2> /dev/null | samtools sort -@ 4 -o "$OUT/reads_vs_ref.bam" -
samtools index "$OUT/reads_vs_ref.bam"
samtools coverage "$OUT/reads_vs_ref.bam" > "$OUT/depth.tsv"

# AMR genes
micromamba activate "$ENV_AMR"
amrfinder -n "$NEWREF" --threads 8 -o "$OUT/amrfinder.tsv" || echo "WARN: amrfinder failed"

# Replicon / mobility for every contig except the longest
micromamba activate "$ENV_MOB"
CHROM=$(awk -F'\t' 'NR>1{print $2"\t"$1}' "$OUT/contigs_new.tsv" | sort -k1,1nr | head -1 | cut -f2)
awk -v c="$CHROM" '/^>/{p = (substr($1,2) != c)} p' "$NEWREF" > "$OUT/non_chrom_contigs.fasta"
if [ -s "$OUT/non_chrom_contigs.fasta" ]; then
    mob_typer --multi -i "$OUT/non_chrom_contigs.fasta" -o "$OUT/mobtyper.tsv" -n 8 > "$OUT/mobtyper.log" 2>&1 \
        || echo "WARN: mob_typer failed (see mobtyper.log)"
else
    : > "$OUT/mobtyper.tsv"
fi

echo -e "STATUS\t$ISO\tOK"
