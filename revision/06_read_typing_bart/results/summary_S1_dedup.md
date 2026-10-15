## AMR: read-based (bart) vs assembly-based (AMRFinderPlus) at the same depth
Paired isolates only. Reads: top-n rule (n from the gold standard). Allele concordance: detected families whose gold call is exact (EXACT/ALLELE); Gold not exact = gold families without such a call.

| Group | Method | Isolates | Carbapenemase family detected | Carbapenemase allele concordant | All core AMR family detected | All core AMR allele concordant | Gold not exact |
|---|---|---|---|---|---|---|---|
| S1a/10x | reads (bart) | 98 | 103/103 (100.0%) | 102/102 (100.0%) | 1249/1329 (94.0%) | 985/1112 (88.6%) | 193 |
| S1a/10x | assembly | 98 | 99/103 (96.1%) | 91/98 (92.9%) | 1287/1329 (96.8%) | 1070/1103 (97.0%) | 193 |
| S1a/20x | reads (bart) | 99 | 104/104 (100.0%) | 103/103 (100.0%) | 1257/1332 (94.4%) | 989/1119 (88.4%) | 193 |
| S1a/20x | assembly | 99 | 104/104 (100.0%) | 101/103 (98.1%) | 1324/1332 (99.4%) | 1127/1134 (99.4%) | 193 |
| S1b/10x | reads (bart) | 10 | 3/3 (100.0%) | 2/3 (66.7%) | 50/72 (69.4%) | 33/37 (89.2%) | 31 |
| S1b/10x | assembly | 10 | 3/3 (100.0%) | 3/3 (100.0%) | 70/72 (97.2%) | 36/39 (92.3%) | 31 |
| S1b/20x | reads (bart) | 10 | 3/3 (100.0%) | 2/3 (66.7%) | 53/72 (73.6%) | 34/39 (87.2%) | 31 |
| S1b/20x | assembly | 10 | 3/3 (100.0%) | 3/3 (100.0%) | 71/72 (98.6%) | 40/40 (100.0%) | 31 |

Carbapenemase calls that differ from the gold standard: summary_S1_dedup_carbapenemase_discordant.tsv (detail = AMRFinderPlus method for assemblies, read depths for bart)

### Isolate level: all carbapenemase alleles correct
| Group | Isolates with carbapenemase | Reads correct | Assembly correct | Reads only | Assembly only | McNemar exact p |
|---|---|---|---|---|---|---|
| S1a/10x | 98 | 98 | 87 | 11 | 0 | 0.000977 |
| S1a/20x | 99 | 99 | 97 | 2 | 0 | 0.5 |
| S1b/10x | 2 | 1 | 2 | 0 | 1 | 1 |
| S1b/20x | 2 | 1 | 2 | 0 | 1 | 1 |

## False positives: carbapenemases (calls made without the gold standard)
Paired isolates with a gold standard. Family FP = carbapenemase family called but absent from the gold standard. Secondary alleles = hits of rank >= 2 within a gold family, allele not in the gold standard, passing the depth threshold (the keep rule itself reports only rank 1, so these show what a depth filter without the rank restriction would add).

| Group | Method | Isolates | Gold carbapenemase families | Detected | Missed | FP families | Secondary alleles | Sensitivity (isolate) | Specificity (isolate) |
|---|---|---|---|---|---|---|---|---|---|
| S1a/10x | reads, unfiltered | 98 | 103 | 103 | 0 | 0 | 82 | 98/98 (100.0%) | NA |
| S1a/10x | reads, keep >= 0.1 | 98 | 103 | 103 | 0 | 0 | 42 | 98/98 (100.0%) | NA |
| S1a/10x | reads, keep >= 0.2 | 98 | 103 | 103 | 0 | 0 | 7 | 98/98 (100.0%) | NA |
| S1a/10x | reads, keep >= 0.3 | 98 | 103 | 102 | 1 | 0 | 1 | 97/98 (99.0%) | NA |
| S1a/10x | reads, keep >= 0.5 | 98 | 103 | 93 | 10 | 0 | 0 | 89/98 (90.8%) | NA |
| S1a/10x | assembly | 98 | 103 | 99 | 4 | 0 | NA | 94/98 (95.9%) | NA |
| S1a/20x | reads, unfiltered | 99 | 104 | 104 | 0 | 3 | 166 | 99/99 (100.0%) | NA |
| S1a/20x | reads, keep >= 0.1 | 99 | 104 | 104 | 0 | 0 | 20 | 99/99 (100.0%) | NA |
| S1a/20x | reads, keep >= 0.2 | 99 | 104 | 104 | 0 | 0 | 2 | 99/99 (100.0%) | NA |
| S1a/20x | reads, keep >= 0.3 | 99 | 104 | 104 | 0 | 0 | 1 | 99/99 (100.0%) | NA |
| S1a/20x | reads, keep >= 0.5 | 99 | 104 | 96 | 8 | 0 | 0 | 91/99 (91.9%) | NA |
| S1a/20x | assembly | 99 | 104 | 104 | 0 | 0 | NA | 99/99 (100.0%) | NA |
| S1b/10x | reads, unfiltered | 10 | 3 | 3 | 0 | 0 | 3 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/10x | reads, keep >= 0.1 | 10 | 3 | 3 | 0 | 0 | 3 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/10x | reads, keep >= 0.2 | 10 | 3 | 2 | 1 | 0 | 1 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/10x | reads, keep >= 0.3 | 10 | 3 | 2 | 1 | 0 | 0 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/10x | reads, keep >= 0.5 | 10 | 3 | 2 | 1 | 0 | 0 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/10x | assembly | 10 | 3 | 3 | 0 | 0 | NA | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | reads, unfiltered | 10 | 3 | 3 | 0 | 0 | 6 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | reads, keep >= 0.1 | 10 | 3 | 3 | 0 | 0 | 4 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | reads, keep >= 0.2 | 10 | 3 | 3 | 0 | 0 | 1 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | reads, keep >= 0.3 | 10 | 3 | 2 | 1 | 0 | 0 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | reads, keep >= 0.5 | 10 | 3 | 2 | 1 | 0 | 0 | 2/2 (100.0%) | 8/8 (100.0%) |
| S1b/20x | assembly | 10 | 3 | 3 | 0 | 0 | NA | 2/2 (100.0%) | 8/8 (100.0%) |

### Isolate level: all carbapenemase alleles correct, reads with keep >= 0.2 (no gold information)
Correct = every gold carbapenemase family called (rank 1, depth >= threshold, core AMR) with an allele contained in the gold allele set (presence only for gold calls that are not exact).
| Group | Isolates with carbapenemase | Reads correct | Assembly correct | Reads only | Assembly only | McNemar exact p |
|---|---|---|---|---|---|---|
| S1a/10x | 98 | 98 | 87 | 11 | 0 | 0.000977 |
| S1a/20x | 99 | 99 | 97 | 2 | 0 | 0.5 |
| S1b/10x | 2 | 1 | 2 | 0 | 1 | 1 |
| S1b/20x | 2 | 1 | 2 | 0 | 1 | 1 |

Details of all false positives: summary_S1_dedup_carbapenemase_fp.tsv

## MLST
Exact = single ST with all alleles 100% identity and coverage (bart) or a numeric ST (mlst).
Reference ST: mlst on the highest-depth assembly >= 40x with a numeric ST (100x: 66, 60x: 33, 40x: 3 isolates).

| Group | Samples | No scheme | Reads exact ST | Reads ST = gold | Assembly ST = gold |
|---|---|---|---|---|---|
| S1a/10x | 100 | 10 | 69/90 (76.7%) | 65/79 (82.3%) | 66/79 (83.5%) |
| S1a/20x | 100 | 9 | 77/91 (84.6%) | 73/80 (91.2%) | 78/80 (97.5%) |
| S1b/10x | 10 | 0 | 9/10 (90.0%) | 9/9 (100.0%) | 9/9 (100.0%) |
| S1b/20x | 10 | 0 | 9/10 (90.0%) | 9/9 (100.0%) | 9/9 (100.0%) |
