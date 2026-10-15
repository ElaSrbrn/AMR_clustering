# Shared DCJ-indel colour scale for Figure 4 (panels A and B).
#
# Sourced by Fig4_panels.Rmd, Fig4A_IncN.Rmd and Fig4B_IncLM.Rmd so that the
# same colour means the same distance in both panels. 
# Change the scale here only, never in the individual scripts.

library(readr)
library(circlize)
library(ComplexHeatmap)
library(grid)

# Input tables for both panels. Both are read by every script, because the
# colour scale depends on the pooled distances of both panels.
fig4_files <- c(
  A = "all_IncN_distances_with_species.tsv",
  B = "inclm_plasmids_distances_with_species.tsv"
)

# Upper end of the colour scale: this quantile of the pooled off-diagonal
# distances of both panels. 1 = pooled maximum, no clamping. DCJ-indel distances
# here are small integers (60x Flye: IncN 0-5, IncL/M 0-11) without outliers, so
# clamping at a lower quantile would merge real differences into one colour.
# If set below 1, distances above the cap are drawn in the darkest colour and the
# legend shows a "≥" label.
fig4_cap_quantile <- 1

# DCJ-indel subcommunity threshold used by pling. If set, the white midpoint of
# the scale is placed at this value, so blue means "within the threshold" and
# orange/red means "beyond it". run_mash_then_pling_by_group_v2.slurm calls pling
# without --dcj, so the pling default of 4 applies (as in the Methods).
# Set to NA for an unanchored scale.
fig4_threshold <- 4

fig4_ramp <- c("#2166ac", "#f7f7f7", "#fdae61", "#d73027")

# One distance per unordered plasmid pair, self-comparisons excluded.
fig4_offdiag <- function(path) {
  d <- read_tsv(path, show_col_types = FALSE)
  d <- d[d$plasmid_1_labeled != d$plasmid_2_labeled, ]
  key <- ifelse(d$plasmid_1_labeled < d$plasmid_2_labeled,
                paste(d$plasmid_1_labeled, d$plasmid_2_labeled, sep = "||"),
                paste(d$plasmid_2_labeled, d$plasmid_1_labeled, sep = "||"))
  d$distance[!duplicated(key)]
}

fig4_pooled <- unlist(lapply(fig4_files, fig4_offdiag), use.names = FALSE)
fig4_cap <- as.numeric(quantile(fig4_pooled, fig4_cap_quantile, na.rm = TRUE))

fig4_n_clamped <- sum(fig4_pooled > fig4_cap, na.rm = TRUE)
cat(sprintf(
  "Shared colour scale 0 - %.2f (q = %.2f of %d pooled pairs; %d pairs above the cap)\n",
  fig4_cap, fig4_cap_quantile, length(fig4_pooled), fig4_n_clamped))
for (p in names(fig4_files)) {
  x <- fig4_offdiag(fig4_files[[p]])
  cat(sprintf("  Panel %s: %d pairs, max %.2f, %d above the cap\n",
              p, length(x), max(x, na.rm = TRUE), sum(x > fig4_cap, na.rm = TRUE)))
}

if (is.na(fig4_threshold)) {
  fig4_col_fun <- colorRamp2(seq(0, fig4_cap, length.out = 4), fig4_ramp)
} else {
  if (fig4_threshold <= 0 || fig4_threshold >= fig4_cap) {
    stop("fig4_threshold must lie between 0 and the colour scale cap.")
  }
  fig4_col_fun <- colorRamp2(
    c(0, fig4_threshold, (fig4_threshold + fig4_cap) / 2, fig4_cap),
    fig4_ramp)
}

# Legend breaks: rounded values below the cap, then the cap itself labelled
# with ">=" when distances above it exist.
fig4_at <- pretty(c(0, fig4_cap), n = 4)
fig4_at <- fig4_at[fig4_at < 0.85 * fig4_cap]
fig4_at <- c(fig4_at, fig4_cap)
fig4_labels <- formatC(fig4_at, format = "fg", digits = 3)
if (fig4_n_clamped > 0) {
  fig4_labels[length(fig4_labels)] <- paste0("≥", round(fig4_cap, 1))
} else {
  fig4_labels[length(fig4_labels)] <- format(round(fig4_cap, 1))
}

fig4_legend_param <- function(legend_pt) {
  list(
    title = "DCJ-indel distance",
    at = fig4_at,
    labels = fig4_labels,
    direction = "horizontal",
    legend_width = unit(35, "mm"),
    title_position = "topleft",
    title_gp = gpar(fontsize = legend_pt),
    labels_gp = gpar(fontsize = legend_pt)
  )
}
