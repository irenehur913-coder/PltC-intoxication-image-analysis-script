#!/usr/bin/env python3
"""
analysis_whole_tissue.py  —  robustness control (whole-tissue, NOT restricted to villi).

Identical pipeline to analysis_villus_restricted.py except the readout is taken over
the whole DAPI+ tissue (villus-interior restriction removed). Provided to show that
the conclusion (B6_Intoxication > all other conditions) does not depend on the villus
restriction. Deterministic: same inputs -> identical outputs.

Usage:
    python analysis_whole_tissue.py [path_to_"2D Images"_folder]

Output: whole_tissue_metrics.csv + per-condition summary and two-sided
Mann-Whitney comparisons (B6_Intoxication vs each other condition).
"""
import sys, numpy as np
from scipy import stats
from pltc_analysis_common import run_folder

DEFAULT_BASE = "2D Images"
ORDER = ['B6_NoToxin', 'B6_Intoxication', 'Tek-Cre_NoToxin', 'Tek-Cre_Intoxication']

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def summarize(rows):
    by = {c: [r for r in rows if r['condition'] == c] for c in ORDER}
    print("\n%-22s %5s %8s %8s %8s %8s" % ("condition", "n", "puncta", "M2", "enrich", "cover"))
    for c in ORDER:
        v = by[c]; g = lambda k: np.mean([r[k] for r in v])
        print("%-22s %5d %8.2f %8.3f %8.3f %8.2f" %
              (c, len(v), g('puncta_density'), g('M2'), g('enrichment'), g('coverage')))
    print("\nB6_Intoxication vs each condition (two-sided Mann-Whitney U):")
    B = by['B6_Intoxication']
    for c in ORDER:
        if c == 'B6_Intoxication':
            continue
        line = "  vs %-20s" % c
        for k in ['puncta_density', 'M2', 'enrichment']:
            p = stats.mannwhitneyu([r[k] for r in B], [r[k] for r in by[c]], alternative='two-sided')[1]
            line += "  %s: p=%.3g %s" % (k.split('_')[0], p, stars(p))
        print(line)

if __name__ == '__main__':
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    rows = run_folder(base, region='whole', out_csv='whole_tissue_metrics.csv')
    summarize(rows)
