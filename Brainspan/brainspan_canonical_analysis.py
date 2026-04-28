#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu

TARGET_GENES: List[str] = [
    'AMIGO3', 'AP1G1', 'APEH', 'APOB', 'ASTN2', 'ASXL1', 'ASXL3', 'BAI2', 'BBX', 'BPTF',
    'BSN', 'C20orf112', 'C6orf106', 'CAMKV', 'CCDC36', 'CDHR4', 'CEP170', 'CPS1',
    'CTD-2330K9.3', 'DCAKD', 'DCC', 'DDX20', 'DNM1', 'ECM1', 'ERBB4', 'FAF1', 'FAM120A',
    'FAM120AOS', 'FAM129A', 'FHL5', 'FOXP2', 'GABRB2', 'GMPPB', 'GNAT1', 'GRK4', 'HEXIM2',
    'INTS3', 'IP6K1', 'KCNC2', 'KLHDC8B', 'LAMB2', 'MAML3', 'MLLT10', 'MON1A', 'MRPS21',
    'MSL2', 'MST1R', 'NCAM1', 'NMT1', 'NRXN1', 'NUMB', 'PABPC4', 'PHF2', 'PPP1R13B', 'PRPF3',
    'RAB5B', 'RABGAP1L', 'RBM5', 'RBM6', 'RERG', 'RNF123', 'RP11-3B7.1', 'RPRD2', 'SDK1',
    'SEMA3F', 'SLC24A3', 'SLC25A13', 'SNRPC', 'SP4', 'SPDEF', 'SPHKAP', 'STAG1', 'TARS2',
    'TM9SF4', 'TRPC7', 'UFL1', 'UHRF1BP1', 'VPS33B'
]

KEY_CASCADE_GENES = [
    'BSN', 'ERBB4', 'DNM1', 'LAMB2', 'KCNC2', 'NCAM1', 'SEMA3F', 'NRXN1',
    'GABRB2', 'FOXP2', 'NUMB', 'MAML3', 'DCC', 'SDK1'
]

PERIOD_LABELS = {
    1: '1_Early_Prenatal',
    2: '2_Early_Mid_Prenatal',
    3: '3_Late_Mid_Prenatal',
    4: '4_Late_Prenatal',
    5: '5_Neonatal_Early_Infancy',
    6: '6_Late_Infancy',
    7: '7_Early_Childhood',
    8: '8_Middle_Childhood',
    9: '9_Adolescence',
    10: '10_Young_Adulthood',
    11: '11_Middle_Adulthood',
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Canonical BrainSpan prenatal peak analysis.')
    p.add_argument('--input', required=True, help='Path to BrainSpan zip or directory containing the three CSV files.')
    p.add_argument('--output-dir', required=True, help='Directory for outputs.')
    p.add_argument('--permutations', type=int, default=10000, help='Number of permutation draws.')
    p.add_argument('--seed', type=int, default=42, help='Random seed for permutation test.')
    return p.parse_args()


def _open_csv_bytes(input_path: Path, name: str):
    if input_path.is_dir():
        return open(input_path / name, 'rb')
    if input_path.suffix.lower() == '.zip':
        zf = zipfile.ZipFile(input_path)
        return zf.open(name, 'r')
    raise ValueError(f'Unsupported input path: {input_path}')


def read_brainspan_tables(input_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], Dict[str, List[str]]]:
    resources = {}
    filenames = ['rows_metadata.csv', 'columns_metadata.csv', 'expression_matrix.csv']
    for fname in filenames:
        resources[fname] = _open_csv_bytes(input_path, fname)

    rows = pd.read_csv(resources['rows_metadata.csv'])
    cols = pd.read_csv(resources['columns_metadata.csv'])
    expr = pd.read_csv(resources['expression_matrix.csv'], header=None)
    expr = expr.set_index(0)
    expr.columns = pd.Index(range(1, expr.shape[1] + 1), name='column_num')

    raw_preview = {}
    if input_path.is_dir():
        expr_path = input_path / 'expression_matrix.csv'
        with open(expr_path, 'r', newline='') as fh:
            for line in fh:
                pieces = line.strip().split(',')
                if pieces and pieces[0] in {'11351', '16767', '8849'}:
                    raw_preview[pieces[0]] = pieces[:6]
    else:
        with zipfile.ZipFile(input_path) as zf:
            with zf.open('expression_matrix.csv', 'r') as fh:
                for raw in io.TextIOWrapper(fh):
                    pieces = raw.strip().split(',')
                    if pieces and pieces[0] in {'11351', '16767', '8849'}:
                        raw_preview[pieces[0]] = pieces[:6]
                    if len(raw_preview) == 3:
                        break

    for f in resources.values():
        try:
            f.close()
        except Exception:
            pass
    return rows, cols, expr, filenames, raw_preview


def assign_period(age: str) -> int:
    parts = str(age).strip().split()
    if len(parts) < 2:
        raise ValueError(f'Could not parse age string: {age!r}')
    value = int(parts[0])
    unit = parts[1].lower()

    if unit == 'pcw':
        if 8 <= value <= 9:
            return 1
        if 10 <= value <= 12:
            return 2
        if 13 <= value <= 18:
            return 3
        if 19 <= value <= 24:
            return 4
        if 25 <= value <= 38:
            return 5
    elif unit == 'mos':
        return 5 if 0 <= value <= 5 else 6
    elif unit == 'yrs':
        if 1 <= value <= 2:
            return 6
        if 3 <= value <= 5:
            return 7
        if 6 <= value <= 11:
            return 8
        if 12 <= value <= 17:
            return 9
        if 18 <= value <= 29:
            return 10
        if value >= 30:
            return 11
    raise ValueError(f'Age value {age!r} did not map to a developmental period')


def build_period_masks(columns: pd.DataFrame) -> Dict[int, np.ndarray]:
    columns = columns.copy()
    columns['period_number'] = columns['age'].map(assign_period)
    masks: Dict[int, np.ndarray] = {}
    for period in range(1, 12):
        col_nums = columns.loc[columns['period_number'] == period, 'column_num'].astype(int).tolist()
        masks[period] = np.array(col_nums, dtype=int)
    return masks


def filter_background_genes(rows: pd.DataFrame) -> pd.DataFrame:
    rows2 = rows.copy()
    rows2['gene_symbol'] = rows2['gene_symbol'].astype(str)
    valid_symbol = rows['gene_symbol'].notna() & (rows2['gene_symbol'].str.strip() != '')
    valid_entrez = rows['entrez_id'].notna() & pd.to_numeric(rows['entrez_id'], errors='coerce').fillna(0).astype(float).ne(0)
    bg = rows.loc[valid_symbol & valid_entrez].copy()
    bg = bg.drop_duplicates(subset=['gene_symbol'], keep='first').copy()
    bg['row_num'] = bg['row_num'].astype(int)
    bg['entrez_id'] = pd.to_numeric(bg['entrez_id'], errors='coerce').astype('Int64')
    return bg


def compute_peak_table(meta: pd.DataFrame, expr: pd.DataFrame, period_masks: Dict[int, np.ndarray]) -> pd.DataFrame:
    meta2 = meta[['gene_symbol', 'row_num']].copy()
    row_nums = meta2['row_num'].astype(int).tolist()
    subexpr = expr.loc[row_nums].to_numpy(dtype=float)
    period_mean_arrays = []
    for period in range(1, 12):
        cols = period_masks[period] - 1  # expr columns are labeled 1..524; numpy arrays are zero-based
        period_mean_arrays.append(subexpr[:, cols].mean(axis=1))
    mean_matrix = np.column_stack(period_mean_arrays)
    peak_idx = mean_matrix.argmax(axis=1)
    peak_period_numbers = peak_idx + 1
    out = meta2.copy()
    out['peak_period_number'] = peak_period_numbers
    out['peak_period_label'] = [PERIOD_LABELS[p] for p in peak_period_numbers]
    out['prenatal_peak'] = peak_period_numbers <= 4
    for period in range(1, 12):
        out[f'mean_rpkm_period_{period}'] = mean_matrix[:, period - 1]
    return out


def mann_whitney_two_sided(prenatal: np.ndarray, postnatal: np.ndarray) -> float:
    try:
        return float(mannwhitneyu(prenatal, postnatal, alternative='two-sided').pvalue)
    except ValueError:
        return float('nan')


def compute_target_detailed(target_peak: pd.DataFrame, expr: pd.DataFrame, period_masks: Dict[int, np.ndarray]) -> pd.DataFrame:
    prenatal_cols = np.concatenate([period_masks[p] for p in [1, 2, 3, 4]])
    postnatal_cols = np.concatenate([period_masks[p] for p in [6, 7, 8, 9, 10, 11]])
    records = []
    for _, row in target_peak.iterrows():
        values = expr.loc[int(row['row_num'])].astype(float)
        prenatal_vals = values.loc[prenatal_cols].to_numpy(dtype=float)
        postnatal_vals = values.loc[postnatal_cols].to_numpy(dtype=float)
        mean_pre = float(np.mean(prenatal_vals))
        mean_post = float(np.mean(postnatal_vals))
        log2fc = float(np.log2((mean_pre + 0.1) / (mean_post + 0.1)))
        if log2fc > 0.5:
            cls = 'Prenatal-biased'
        elif log2fc < -0.5:
            cls = 'Postnatal-biased'
        else:
            cls = 'Stable'
        records.append({
            'gene_symbol': row['gene_symbol'],
            'row_num': int(row['row_num']),
            'peak_period_number': int(row['peak_period_number']),
            'peak_period_label': row['peak_period_label'],
            'prenatal_peak': bool(row['prenatal_peak']),
            'log2FC': log2fc,
            'classification': cls,
            'mann_whitney_p': mann_whitney_two_sided(prenatal_vals, postnatal_vals),
            'mean_rpkm_prenatal': mean_pre,
            'mean_rpkm_postnatal': mean_post,
        })
    return pd.DataFrame.from_records(records)


def compute_zscored_by_period(target_peak: pd.DataFrame, expr: pd.DataFrame, period_masks: Dict[int, np.ndarray]) -> pd.DataFrame:
    records = []
    for _, row in target_peak.iterrows():
        values = expr.loc[int(row['row_num'])].to_numpy(dtype=float)
        sd = float(np.std(values, ddof=0))
        if sd == 0:
            z = np.zeros_like(values)
        else:
            z = (values - np.mean(values)) / sd
        rec = {'gene_symbol': row['gene_symbol']}
        for period in range(1, 12):
            cols = period_masks[period] - 1  # z array is zero-based positional now
            rec[PERIOD_LABELS[period]] = float(np.mean(z[cols]))
        records.append(rec)
    return pd.DataFrame.from_records(records)


def permutation_test(all_prenatal: np.ndarray, observed: int, draw_size: int, n_perm: int, seed: int) -> Tuple[np.ndarray, float, float, float]:
    np.random.seed(seed)
    n = len(all_prenatal)
    null_counts = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        idx = np.random.choice(n, size=draw_size, replace=False)
        null_counts[i] = int(all_prenatal[idx].sum())
    p = float(np.mean(null_counts >= observed))
    percentile = float(100.0 * np.mean(null_counts <= observed))
    return null_counts, p, float(np.mean(null_counts)), float(np.std(null_counts, ddof=0)), percentile


def write_histogram(null_counts: np.ndarray, observed: int, p_value: float, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))
    plt.hist(null_counts, bins=range(int(null_counts.min()), int(null_counts.max()) + 2), edgecolor='black')
    plt.axvline(observed, color='red', linewidth=2)
    plt.title('Permutation null: prenatal peak count in 78 random brain-expressed genes')
    plt.xlabel('Number of genes peaking prenatally (out of 78)')
    plt.ylabel('Frequency')
    plt.text(0.98, 0.95, f'Permutation P = {p_value:.4g}', transform=plt.gca().transAxes,
             ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_console_summary(
    verification: List[str],
    sample_counts: Dict[int, int],
    n_background: int,
    n_background_prenatal: int,
    background_rate: float,
    n_targets_present: int,
    n_targets_bg: int,
    observed_all78: int,
    observed_bg: int,
    binom_all78: float,
    binom_bg: float,
    null_mean: float,
    null_sd: float,
    observed_perm: int,
    perm_p: float,
    percentile: float,
    detailed: pd.DataFrame,
) -> str:
    lines = []
    lines.append('=== ROW MAPPING VERIFICATION ===')
    lines.extend(verification)
    lines.append('')
    lines.append('=== SAMPLE COUNTS PER PERIOD ===')
    for period in range(1, 12):
        lines.append(f"{PERIOD_LABELS[period]}: {sample_counts[period]} samples")
    lines.append('')
    lines.append('=== BACKGROUND ===')
    lines.append(f'N_background: {n_background}')
    lines.append(f'N_background_prenatal: {n_background_prenatal}')
    lines.append(f'Background prenatal rate: {background_rate*100:.2f}%')
    lines.append('')
    lines.append('=== TARGET GENES ===')
    lines.append(f'N target genes present in BrainSpan: {n_targets_present}/78')
    lines.append(f'N target genes in protein-coding background: {n_targets_bg}/78')
    lines.append(f'N target genes peaking prenatally (all 78): {observed_all78}/78 = {observed_all78/78*100:.2f}%')
    lines.append(f'N target genes peaking prenatally (77 protein-coding background genes): {observed_bg}/{n_targets_bg} = {observed_bg/n_targets_bg*100:.2f}%')
    lines.append('')
    lines.append('=== BINOMIAL TEST (uniform null) ===')
    lines.append(f'Binomial P (one-sided, null=4/11) using all 78 targets: {binom_all78:.10g}')
    lines.append(f'Binomial P (one-sided, null=4/11) using 77 protein-coding background targets: {binom_bg:.10g}')
    lines.append('')
    lines.append('=== PERMUTATION TEST (brain-gene background null) ===')
    lines.append(f'Null mean: {null_mean:.4f}')
    lines.append(f'Null SD: {null_sd:.4f}')
    lines.append(f'Observed: {observed_perm}')
    lines.append(f'Permutation P: {perm_p:.4g}')
    lines.append(f'Observed percentile (<=): {percentile:.1f}')
    lines.append('')
    lines.append('=== KEY CASCADE GENES ===')
    d = detailed.set_index('gene_symbol')
    for gene in KEY_CASCADE_GENES:
        r = d.loc[gene]
        lines.append(f"{gene}: peak = {r['peak_period_label']}, log2FC = {r['log2FC']:.3f}, class = {r['classification']}")
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows, cols, expr, _, raw_preview = read_brainspan_tables(input_path)

    # verification
    verification = []
    for gene in ['BSN', 'DCC', 'GABRB2']:
        row_num = int(rows.loc[rows['gene_symbol'] == gene, 'row_num'].iloc[0])
        first5 = expr.loc[row_num].astype(float).iloc[:5].tolist()
        verification.append(f"{gene}: row_num={row_num}, first 5 RPKM values from expr.loc[row_num]={[round(x, 6) for x in first5]}")
        preview = raw_preview.get(str(row_num))
        if preview is not None:
            verification.append(f"{gene}: raw CSV line starts with {preview}")

    period_masks = build_period_masks(cols)
    sample_counts = {p: len(period_masks[p]) for p in range(1, 12)}

    background = filter_background_genes(rows)
    all_peak = compute_peak_table(background[['gene_symbol', 'row_num']], expr, period_masks)
    all_peak.to_csv(outdir / 'all_genes_peak_periods.csv', index=False)

    n_background = len(all_peak)
    n_background_prenatal = int(all_peak['prenatal_peak'].sum())
    background_rate = n_background_prenatal / n_background

    target_present = rows.loc[rows['gene_symbol'].isin(TARGET_GENES), ['gene_symbol', 'row_num']].copy()
    target_present['row_num'] = target_present['row_num'].astype(int)
    target_peak = compute_peak_table(target_present, expr, period_masks)
    target_peak['target_order'] = target_peak['gene_symbol'].map({g: i for i, g in enumerate(TARGET_GENES)})
    target_peak = target_peak.sort_values('target_order').drop(columns='target_order')
    target_peak.to_csv(outdir / 'somatic78_canonical_peak_periods.csv', index=False)

    target_detailed = compute_target_detailed(target_peak, expr, period_masks)
    target_detailed['target_order'] = target_detailed['gene_symbol'].map({g: i for i, g in enumerate(TARGET_GENES)})
    target_detailed = target_detailed.sort_values('target_order').drop(columns='target_order')
    target_detailed.to_csv(outdir / 'somatic78_canonical_detailed.csv', index=False)

    z_by_period = compute_zscored_by_period(target_peak, expr, period_masks)
    z_by_period['target_order'] = z_by_period['gene_symbol'].map({g: i for i, g in enumerate(TARGET_GENES)})
    z_by_period = z_by_period.sort_values('target_order').drop(columns='target_order')
    z_by_period.to_csv(outdir / 'somatic78_zscored_by_period.csv', index=False)

    background_symbols = set(background['gene_symbol'])
    target_bg = target_peak[target_peak['gene_symbol'].isin(background_symbols)].copy()
    observed_all78 = int(target_peak['prenatal_peak'].sum())
    observed_bg = int(target_bg['prenatal_peak'].sum())
    n_targets_present = len(target_peak)
    n_targets_bg = len(target_bg)

    binom_all78 = float(binomtest(observed_all78, n=78, p=4/11, alternative='greater').pvalue)
    binom_bg = float(binomtest(observed_bg, n=n_targets_bg, p=4/11, alternative='greater').pvalue)

    all_prenatal = all_peak['prenatal_peak'].to_numpy(dtype=bool)
    null_counts, perm_p, null_mean, null_sd, percentile = permutation_test(
        all_prenatal, observed_bg, n_targets_bg, args.permutations, args.seed
    )
    pd.DataFrame({'prenatal_count': null_counts}).to_csv(outdir / 'permutation_null_distribution.csv', index=False)
    write_histogram(null_counts, observed_bg, perm_p, outdir / 'brainspan_canonical_permutation_histogram.png')

    summary = pd.DataFrame([
        {'metric': 'N_background', 'value': n_background},
        {'metric': 'N_background_prenatal', 'value': n_background_prenatal},
        {'metric': 'background_prenatal_rate', 'value': background_rate},
        {'metric': 'N_targets_present', 'value': n_targets_present},
        {'metric': 'N_targets_in_background', 'value': n_targets_bg},
        {'metric': 'observed_prenatal_all78', 'value': observed_all78},
        {'metric': 'observed_prenatal_background_targets', 'value': observed_bg},
        {'metric': 'binomial_p_all78', 'value': binom_all78},
        {'metric': 'binomial_p_background_targets', 'value': binom_bg},
        {'metric': 'null_mean', 'value': null_mean},
        {'metric': 'null_sd', 'value': null_sd},
        {'metric': 'permutation_p', 'value': perm_p},
        {'metric': 'observed_percentile', 'value': percentile},
    ])
    summary.to_csv(outdir / 'brainspan_canonical_summary.csv', index=False)

    console = make_console_summary(
        verification, sample_counts, n_background, n_background_prenatal, background_rate,
        n_targets_present, n_targets_bg, observed_all78, observed_bg, binom_all78, binom_bg,
        null_mean, null_sd, observed_bg, perm_p, percentile, target_detailed
    )
    (outdir / 'brainspan_canonical_console_summary.txt').write_text(console)
    print(console)

    metadata = {
        'input': str(input_path),
        'output_dir': str(outdir),
        'permutations': args.permutations,
        'seed': args.seed,
        'target_genes': TARGET_GENES,
    }
    (outdir / 'run_metadata.json').write_text(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
