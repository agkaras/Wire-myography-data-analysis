"""
Wire Myography Analysis
=======================
Author : Agnieszka Karas, PhD
Contact: agaakaras@gmail.com | linkedin.com/in/agnieszka-karas

Edit INPUT_FILE, OUTPUT_DIR and the settings block below, then run:

    python myograph.py
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 10,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 1.1, 'lines.linewidth': 1.6, 'figure.dpi': 130,
})
# =============================================================================
#Load file: insert full path of the txt file for analysis and desired output directory
# =============================================================================
INPUT_FILE = '//Users/agnieszkakaras/Wire-myography-data-analysis/sample_data.txt'
OUTPUT_DIR = '/Users/agnieszkakaras/output'

CHANNEL_COLS = [1, 2, 3, 4, 5, 6, 7, 8]

df = pd.read_csv(
     INPUT_FILE, sep='\t', decimal=',', skiprows=9,
    header=None, names=['t', 1, 2, 3, 4, 5, 6, 7, 8, 'comment'],
    dtype={c: float for c in ['t'] + CHANNEL_COLS},
    keep_default_na=False)

df['comment'] = df['comment'].astype(str).str.strip()
for c in CHANNEL_COLS:
    df[c] = pd.to_numeric(df[c], errors='coerce')


def preview_comments(df): 
    return df[df['comment'].str.startswith('#*')][['t', 'comment']].reset_index(drop=True)
preview_comments(df)
# =============================================================================
# Settings — edit here
# =============================================================================

# Choose Phe protocol - multiple doses or one!
PHE_MODE = 'dose_response'   # 'dose_response' or 'single_dose'

### Edit labels to match comments in Labchart ###
# Cumulative dose labels and concentrations 
PHE_LABELS = ['Phe0,01', 'Phe0,03', 'Phe0,1', 'Phe0,3', 'Phe1', 'Phe3']
PHE_CONC   = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]   # uM

ACH_LABELS = ['Ach0,001', 'Ach0,01', 'Ach0,1', 'Ach1', 'Ach10']
ACH_CONC   = [0.001, 0.01, 0.1, 1.0, 10.0]       # uM

SNP_LABELS = ['SNP0,001', 'SNP0,01', 'SNP0,1', 'SNP1']
SNP_CONC   = [0.001, 0.01, 0.1, 1.0]              # uM

# Default protocol labels — edit if your file uses different annotations
L_KCL60       = 'KCl60'    # KCl 60 mM addition
L_KCL60_END   = 'P2'       # rinse after KCl60
L_PHE_END     = 'PP'       # rinse after Phe: dose-response or maximal contraction
L_PHE_SINGLE  = 'Phe3uM'  # single-dose Phe label (maximal contraction)
L_SUBPHE      = 'subPhe'   # submaximal Phe precontraction before ACh
L_SUBPHE_END  = 'Ach0,001' # first ACh dose (= end of Phe precontraction window)
L_SUBPHE2     = '2subPhe'  # submaximal Phe precontraction before SNP
L_SUBPHE2_END = 'SNP0,001' # first SNP dose (= end of Phe precontraction window)
L_ACH_END     = 'P3'       # rinse after ACh
L_SNP_END     = 'K'        # end of recording

CH_COLORS    = plt.cm.tab10(np.linspace(0, 0.8, 8))


# ── Functions ──────────────────────────────────────────────────────────────

def _find_row(df, label):
    search = label.strip().lstrip('#').lstrip('*').strip()
    norm   = df['comment'].str.lstrip('#').str.lstrip('*').str.strip()
    hits   = np.where(norm == search)[0]
    if len(hits) == 0:
        avail = df.loc[df['comment'].str.startswith('#*'), 'comment'].unique()
        raise ValueError(f"Label '{label}' not found.\nAvailable: {list(avail)}")
    return int(hits[0])


def compute_contraction(df, start_label, end_label):
    s = _find_row(df, start_label)
    e = _find_row(df, end_label)
    return (df.iloc[s:e][CHANNEL_COLS].max() - df.iloc[s - 1][CHANNEL_COLS]).rename('mN')


def compute_cumulative_contraction(df, dose_labels, end_label):
    baseline = df.iloc[_find_row(df, dose_labels[0]) - 1][CHANNEL_COLS]
    records  = []
    for i, label in enumerate(dose_labels):
        s = _find_row(df, label)
        e = (_find_row(df, dose_labels[i + 1]) if i < len(dose_labels) - 1
             else _find_row(df, end_label))
        records.append(df.iloc[s:e][CHANNEL_COLS].max() - baseline)
    return pd.DataFrame(records, index=dose_labels)


def compute_cumulative_relaxation(df, dose_labels, end_label):
    baseline = df.iloc[_find_row(df, dose_labels[0]) - 1][CHANNEL_COLS]
    records  = []
    for i, label in enumerate(dose_labels):
        s = _find_row(df, label)
        e = (_find_row(df, dose_labels[i + 1]) if i < len(dose_labels) - 1
             else _find_row(df, end_label))
        records.append(baseline - df.iloc[s:e][CHANNEL_COLS].min())
    return pd.DataFrame(records, index=dose_labels)
# =============================================================================
# Curve fitting functions
# =============================================================================

def hill_equation(x, bottom, top, ec50, n):
    """4-parameter logistic (Hill) equation for dose-response fitting."""
    return bottom + (top - bottom) / (1 + (ec50 / x) ** n)


def fit_dose_response(conc, responses_pct):
    """
    Fit Hill equation to dose-response data for each channel.

    Parameters
    ----------
    conc : list of float
        Concentrations (same units as PHE_CONC / ACH_CONC / SNP_CONC).
    responses_pct : pd.DataFrame
        Rows = doses, columns = channels (1–8), values = % response.

    Returns
    -------
    fit_params : dict  {channel: {'bottom', 'top', 'ec50', 'n'} or None}
    fit_curves : dict  {channel: (x_fine, y_fine) or None}
    """
    x = np.array(conc, dtype=float)
    x_fine = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)

    fit_params = {}
    fit_curves = {}

    for ch in CHANNEL_COLS:
        y = responses_pct[ch].values.astype(float)

        # Skip channels with all-NaN or flat response
        if np.all(np.isnan(y)) or np.ptp(y[~np.isnan(y)]) < 1:
            fit_params[ch] = None
            fit_curves[ch] = None
            continue

        # Initial guesses: bottom ~ min, top ~ max, EC50 ~ geometric midpoint, n = 1
        p0 = [np.nanmin(y), np.nanmax(y), np.sqrt(x.min() * x.max()), 1.0]
        bounds = (
            [-10,   0,   x.min() * 0.01, 0.1],
            [50,  150,  x.max() * 100,  10.0],
        )

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(
                    hill_equation, x, y,
                    p0=p0, bounds=bounds,
                    maxfev=5000
                )
            fit_params[ch] = {
                'bottom': popt[0],
                'top':    popt[1],
                'ec50':   popt[2],
                'n':      popt[3],
            }
            fit_curves[ch] = (x_fine, hill_equation(x_fine, *popt))
        except Exception:
            fit_params[ch] = None
            fit_curves[ch] = None

    return fit_params, fit_curves


def summarise_ec50(fit_params, label=""):
    """Print EC50 summary table for all channels."""
    print(f"\n{'─'*50}")
    print(f"  EC50 summary — {label}")
    print(f"{'─'*50}")
    print(f"  {'Channel':<10} {'EC50 (uM)':<14} {'Emax (%)':<12} {'Hill n'}")
    print(f"  {'─'*7:<10} {'─'*9:<14} {'─'*8:<12} {'─'*6}")
    for ch in CHANNEL_COLS:
        p = fit_params.get(ch)
        if p:
            print(f"  Ch {ch:<7} {p['ec50']:<14.4f} {p['top']:<12.1f} {p['n']:.2f}")
        else:
            print(f"  Ch {ch:<7} {'fit failed':<14} {'—':<12} {'—'}")
    print(f"{'─'*50}\n")

# ── Main analysis ─────────────────────────────────────────────────────────────

def analyse_experiment(df, phe_mode, phe_labels, ach_labels, snp_labels, l):

    rows_mn  = []
    rows_pct = []

    def add_mn(label, series):
        row = {'metric': label}
        for ch in CHANNEL_COLS:
            row[f'ch{ch}'] = round(float(series[ch]), 4)
        rows_mn.append(row)

    def add_pct(label, series):
        row = {'metric': label}
        for ch in CHANNEL_COLS:
            row[f'ch{ch}'] = round(float(series[ch]), 4)
        rows_pct.append(row)

    # KCl 60 mM
    kcl60 = compute_contraction(df, l['KCl60'], l['KCl60_end'])
    add_mn('KCl60_mN', kcl60)

    # Phenylephrine - dose response or single dose induced contraction
    if phe_mode == 'dose_response':
        phe_abs = compute_cumulative_contraction(df, phe_labels, l['phe_end'])
        phe_pct = phe_abs.div(kcl60, axis=1) * 100
        for label in phe_labels:
            add_mn( f'Phe_{label}_mN',    phe_abs.loc[label])
            add_pct(f'Phe_{label}_%_KCl', phe_pct.loc[label])
    else:  # single_dose
        phe_max = compute_contraction(df, l['phe_single'], l['phe_end'])
        add_mn( 'Phe_3uM_mN',    phe_max)
        add_pct('Phe_3uM_%_KCl', phe_max / kcl60 * 100)
        phe_pct = None

    # Submaximal Phe (before ACh and SNP)
    phe_sub_Ach = compute_contraction(df, l['subPhe'],  l['subPhe_end'])
    phe_sub_SNP = compute_contraction(df, l['subPhe2'], l['subPhe2_end'])
    add_mn('subPhe_ACh_mN', phe_sub_Ach)
    add_mn('subPhe_SNP_mN', phe_sub_SNP)

    # ACh relaxation
    ach_abs = compute_cumulative_relaxation(df, ach_labels, l['ach_end'])
    ach_pct = ach_abs.div(phe_sub_Ach, axis=1) * 100
    for label in ach_labels:
        add_mn( f'ACh_{label}_mN', ach_abs.loc[label])
        add_pct(f'ACh_{label}_%',  ach_pct.loc[label])

    # SNP relaxation
    snp_abs = compute_cumulative_relaxation(df, snp_labels, l['snp_end'])
    snp_pct = snp_abs.div(phe_sub_SNP, axis=1) * 100
    for label in snp_labels:
        add_mn( f'SNP_{label}_mN', snp_abs.loc[label])
        add_pct(f'SNP_{label}_%',  snp_pct.loc[label])

    # mN rows first, then % rows
    result_df = pd.DataFrame(rows_mn + rows_pct).set_index('metric')

    # Store subtables for plotting
    result_df.attrs['kcl60']      = kcl60
    result_df.attrs['ach_%']      = ach_pct
    result_df.attrs['snp_%']      = snp_pct
    result_df.attrs['phe_%']      = phe_pct
    result_df.attrs['phe_labels'] = phe_labels if phe_mode == 'dose_response' else None
    # Fit dose-response curves
    phe_fit_params, phe_fit_curves = fit_dose_response(PHE_CONC, phe_pct) if phe_pct is not None else ({}, {})
    ach_fit_params, ach_fit_curves = fit_dose_response(ACH_CONC, ach_pct)
    snp_fit_params, snp_fit_curves = fit_dose_response(SNP_CONC, snp_pct)
    
    result_df.attrs['phe_fit'] = (phe_fit_params, phe_fit_curves)
    result_df.attrs['ach_fit'] = (ach_fit_params, ach_fit_curves)
    result_df.attrs['snp_fit'] = (snp_fit_params, snp_fit_curves)
    
    return result_df


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_summary(results, filename):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax_phe, ax_ach, ax_snp, ax_kcl = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    def _ylim(values, margin=0.05):
        mn, mx = np.nanmin(values), np.nanmax(values)
        span = mx - mn if mx != mn else max(abs(mx), 1)
        return mn - span * margin, mx + span * margin

    # ── Phe panel ────────────────────────────────────────────────────────────
    if results.attrs['phe_%'] is not None:
        phe_pct = results.attrs['phe_%']
        for i, ch in enumerate(range(1, 9)):
            ax_phe.plot(PHE_CONC, phe_pct[ch].values, 'o-',
                        color=CH_COLORS[i], label=f'Ch {ch}',
                        markerfacecolor='white', markeredgewidth=1.5, markersize=5)
        ax_phe.set_xscale('log')
        ax_phe.set_xlabel('Phenylephrine (uM)')
        ax_phe.set_ylabel('Contraction (% KCl 60 mM)')
        ax_phe.set_title('Phenylephrine dose-response')
        vals = [v for ch in range(1,9) for v in phe_pct[ch].values]
        ax_phe.set_ylim(*_ylim(phe_pct.values))
        ax_phe.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        ax_phe.legend(fontsize=8, ncol=2, frameon=False)
        
        phe_fit_params, phe_fit_curves = results.attrs.get('phe_fit', ({}, {}))
        for i, ch in enumerate(range(1, 9)):
            curve = phe_fit_curves.get(ch)
            if curve:
                ax_phe.plot(curve[0], curve[1], '-', color=CH_COLORS[i], alpha=0.4, linewidth=1)
    else:
        phe_max_pct = results.loc['Phe_3uM_%_KCl']
        bars = ax_phe.bar([f'Ch {ch}' for ch in range(1, 9)], phe_max_pct.values,
                          color=CH_COLORS, edgecolor='white', width=0.65)
        for bar, val in zip(bars, phe_max_pct.values):
            ax_phe.text(bar.get_x() + bar.get_width() / 2, val + 0.5, f'{val:.1f}',
                        ha='center', va='bottom', fontsize=8)
        ax_phe.set_ylabel('Contraction (% KCl 60 mM)')
        ax_phe.set_title('Phenylephrine 3 uM — max contraction')
        ax_phe.legend(frameon=False)
        plt.setp(ax_phe.get_xticklabels(), rotation=30, ha='right')

    # ── ACh relaxation ────────────────────────────────────────────────────────
    ach_pct = results.attrs['ach_%']
    for i, ch in enumerate(range(1, 9)):
        ax_ach.plot(ACH_CONC, ach_pct[ch].values, 'o-',
                    color=CH_COLORS[i], label=f'Ch {ch}',
                    markerfacecolor='white', markeredgewidth=1.5, markersize=5)
    ax_ach.set_xscale('log')
    ax_ach.set_xlabel('Acetylcholine (uM)')
    ax_ach.set_ylabel('Relaxation (% Phe pre-contraction)')
    ax_ach.set_title('Endothelium-dependent relaxation (ACh)')
    vals = [v for ch in range(1,9) for v in ach_pct[ch].values]
    ax_ach.set_ylim(*_ylim(ach_pct.values))
    ax_ach.invert_yaxis()
    ax_ach.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
    ax_ach.legend(fontsize=8, ncol=2, frameon=False, loc='lower left')
    ach_fit_params, ach_fit_curves = results.attrs.get('ach_fit', ({}, {}))
    for i, ch in enumerate(range(1, 9)):
        curve = ach_fit_curves.get(ch)
        if curve:
            ax_ach.plot(curve[0], curve[1], '-', color=CH_COLORS[i], alpha=0.4, linewidth=1)
    # ── SNP relaxation ────────────────────────────────────────────────────────
    snp_pct = results.attrs['snp_%']
    for i, ch in enumerate(range(1, 9)):
        ax_snp.plot(SNP_CONC, snp_pct[ch].values, 'o-',
                    color=CH_COLORS[i], label=f'Ch {ch}',
                    markerfacecolor='white', markeredgewidth=1.5, markersize=5)
    ax_snp.set_xscale('log')
    ax_snp.set_xlabel('SNP (uM)')
    ax_snp.set_ylabel('Relaxation (% Phe pre-contraction)')
    ax_snp.set_title('Endothelium-independent relaxation (SNP)')
    vals = [v for ch in range(1,9) for v in snp_pct[ch].values]
    ax_snp.set_ylim(*_ylim(snp_pct.values))
    ax_snp.invert_yaxis()
    ax_snp.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
    ax_snp.legend(fontsize=8, ncol=2, frameon=False, loc='lower left')
    snp_fit_params, snp_fit_curves = results.attrs.get('snp_fit', ({}, {}))
    for i, ch in enumerate(range(1, 9)):
        curve = snp_fit_curves.get(ch)
        if curve:
            ax_snp.plot(curve[0], curve[1], '-', color=CH_COLORS[i], alpha=0.4, linewidth=1)

    # ── KCl 60 mM barplot ────────────────────────────────────────────────────
    kcl = results.attrs['kcl60']
    bars = ax_kcl.bar([f'Ch {ch}' for ch in range(1, 9)], kcl.values,
                      color=CH_COLORS, edgecolor='white', width=0.65)
    for bar, val in zip(bars, kcl.values):
        ax_kcl.text(bar.get_x() + bar.get_width() / 2, val + 0.04, f'{val:.2f}',
                    ha='center', va='bottom', fontsize=8)
    ax_kcl.set_ylabel('Contraction (mN)')
    ax_kcl.set_title('KCl 60 mM — viability check')
    ax_kcl.legend(frameon=False)
    plt.setp(ax_kcl.get_xticklabels(), rotation=30, ha='right')

    fig.suptitle(f'Wire Myography — {filename}  [{PHE_MODE}]',
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig


# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    stem = Path(INPUT_FILE).stem

    print()

    LABELS = {
        'KCl60':       L_KCL60,
        'KCl60_end':   L_KCL60_END,
        'phe_end':     L_PHE_END,
        'phe_single':  L_PHE_SINGLE,
        'subPhe':      L_SUBPHE,
        'subPhe_end':  L_SUBPHE_END,
        'subPhe2':     L_SUBPHE2,
        'subPhe2_end': L_SUBPHE2_END,
        'ach_end':     L_ACH_END,
        'snp_end':     L_SNP_END,
    }

    print('Running analysis ...')
    results = analyse_experiment(df, PHE_MODE, PHE_LABELS, ACH_LABELS, SNP_LABELS, LABELS)
    phe_fp, _ = results.attrs.get('phe_fit', ({}, {}))
    ach_fp, _ = results.attrs.get('ach_fit', ({}, {}))
    snp_fp, _ = results.attrs.get('snp_fit', ({}, {}))
    summarise_ec50(phe_fp, label="Phenylephrine")
    summarise_ec50(ach_fp, label="Acetylcholine")
    summarise_ec50(snp_fp, label="SNP")
    print(f'Results table: {results.shape[0]} metrics x {results.shape[1]} channels\n')
    print(results.to_string())
    print()

    print('Generating plot ...')
    fig = plot_summary(results, Path(INPUT_FILE).name)
    out_png = Path(OUTPUT_DIR) / f'{stem}.png'
    fig.savefig(out_png, bbox_inches='tight', dpi=150)
    plt.show()
    print(f'Plot saved: {out_png}')

    print('Exporting Excel ...')
    out_xlsx = Path(OUTPUT_DIR) / f'{stem}.xlsx'
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        results.to_excel(writer, sheet_name='results')
    print(f'Excel saved: {out_xlsx}')
