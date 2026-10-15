#!/bin/bash
#SBATCH --job-name=pcs_indepth
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/06_indepth_%A_%a.out
#
# Step 06b: in-depth annotation of one selected plasmid (EMBL setup):
#   Bakta (reference plasmid)                       -> bakta/<plasmid_id>.gbff/.gff3/.tsv
#   ISfinder BLAST (reference and divergent contigs) -> isfinder_ref.blast6, isfinder_test.blast6
#   self-alignment repeats (reference)               -> self.blast6
# Usage (from the pipeline folder):
#   $PY scripts/06_prepare_indepth.py
#   sbatch --array=1-$(wc -l < work/indepth_list.txt)%10 scripts/06_indepth.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"
PID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" work/indepth_list.txt)
D=work/indepth/$PID
echo "Plasmid $PID"

micromamba activate "$ENV_BAKTA"
echo "Bakta version: $(bakta --version)"
bakta --db "$BAKTA_DB" --output "$D/bakta" --prefix "$PID" --threads 8 --tmp-dir "$JOBTMP" --force "$D/ref.fasta" \
  > "$D/bakta.log" 2>&1 || echo "WARN: bakta failed (see $D/bakta.log)"
# test contigs (non-concordant Flye/Hybracter contigs), same settings, for the synteny figures
if [ -s "$D/test_contigs.fasta" ]; then
  bakta --db "$BAKTA_DB" --output "$D/bakta_test" --prefix "${PID}_test" --threads 8 --tmp-dir "$JOBTMP" --force \
    "$D/test_contigs.fasta" > "$D/bakta_test.log" 2>&1 || echo "WARN: bakta (test contigs) failed (see $D/bakta_test.log)"
fi

micromamba activate "$ENV_BLAST"
FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
blastn -query "$D/ref.fasta" -db "$ISFINDER_DB" -outfmt "$FMT" -num_threads 8 \
  > "$D/isfinder_ref.blast6" || echo "WARN: ISfinder BLAST (ref) failed"
if [ -s "$D/test_contigs.fasta" ]; then
  blastn -query "$D/test_contigs.fasta" -db "$ISFINDER_DB" -outfmt "$FMT" -num_threads 8 \
    > "$D/isfinder_test.blast6" || echo "WARN: ISfinder BLAST (test) failed"
fi
blastn -query "$D/ref.fasta" -subject "$D/ref.fasta" -outfmt "$FMT" -evalue 1e-20 > "$D/self.blast6" || echo "WARN: self BLAST failed"
touch "$D/DONE"
echo -e "STATUS\t$PID\tOK"
