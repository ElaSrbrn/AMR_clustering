#!/usr/bin/env bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=64G





THREADS="${THREADS:-16}"
ANNOT_GLOB="${ANNOT_GLOB:-annotation_*}"
PLING_ENV="${PLING_ENV:-pling}"
MASH_CMD="${MASH_CMD:-mash}"
PLING_CMD="${PLING_CMD:-pling}"
MASH_THRESHOLD="${MASH_THRESHOLD:-0.005}"
PLING_CONTAINMENT="${PLING_CONTAINMENT:-0.3}"
LINK_MODE="${LINK_MODE:-symlink}"
SKIP_PLING="${SKIP_PLING:-0}"

EXTRA_ARGS=()
if [[ "${SKIP_PLING}" == "1" ]]; then
  EXTRA_ARGS+=("--skip-pling")
fi

micromamba run -n "${PLING_ENV}" python3 group_plasmids_mash_then_pling.py \
  --annotation-glob "${ANNOT_GLOB}" \
  --mash-cmd "${MASH_CMD}" \
  --pling-cmd "${PLING_CMD}" \
  --threads "${THREADS}" \
  --mash-threshold "${MASH_THRESHOLD}" \
  --pling-containment "${PLING_CONTAINMENT}" \
  --link-mode "${LINK_MODE}" \
  "${EXTRA_ARGS[@]}"
