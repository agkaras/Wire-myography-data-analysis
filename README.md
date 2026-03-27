[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/agkaras/Wire-myography-data-analysis/blob/main/myograph_colab.ipynb)
Python analysis pipeline for data measured with wire myograph 620M (Danish Myo Technology) exported from LabChart (AD Instruments). Calculates vascular contraction and relaxation responses from `.txt` files and produces summary figures and an Excel results table. Analysis condition: every event (reagent addition, rinse) has to be marked as a comment during recording in LabChart.

## What it does

Parses time-series force recordings from 8-channel myograph experiments. Contraction and relaxation values are calculated from comment-anchored time windows:

- **KCl 60 mM contraction** — smooth muscle viability reference
- **Phenylephrine (Phe)** — α1-adrenergic contractility, either as a full cumulative dose-response curve (0.01–3 μM) or a single maximal dose (3 μM)
- **Acetylcholine (ACh)** — endothelium-dependent relaxation, expressed as % of submaximal Phe pre-contraction
- **SNP** — endothelium-independent relaxation (NO donor), expressed as % of second submaximal Phe pre-contraction

All results land in a single DataFrame (and one Excel sheet): absolute values in mN first, then all percentage values — so columns are easy to copy into GraphPad or R for statistics.

## Files

```
myograph_colab.ipynb    — main analysis notebook, runs in Google Colab
sample_data.txt    — synthetic LabChart export for testing
```

## Usage

Open in Google Colab: **File → Open notebook → GitHub**, paste the repo URL, select `myograph_colab.ipynb`. Or download the .ipynb file and upload in Colab.

Run cells in order. The only cell you need to edit is **Step 3**:

```python
# Protocol mode
PHE_MODE = 'dose_response'   # or 'single_dose'

# Cumulative dose labels — match exactly what was typed in LabChart
PHE_LABELS = ['Phe0,01', 'Phe0,03', 'Phe0,1', 'Phe0,3', 'Phe1', 'Phe3']
ACH_LABELS = ['Ach0,001', 'Ach0,01', 'Ach0,1', 'Ach1', 'Ach10']
SNP_LABELS = ['SNP0,001', 'SNP0,01', 'SNP0,1', 'SNP1']

# Concentrations — must be in the same order as the labels above
PHE_CONC = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]   # uM
ACH_CONC = [0.001, 0.01, 0.1, 1.0, 10.0]        # uM
SNP_CONC = [0.001, 0.01, 0.1, 1.0]              # uM

# Fixed protocol markers — edit if your LabChart annotations differ
L_KCL60       = 'KCl60'    # KCl 60 mM addition
L_KCL60_END   = 'P2'       # rinse after KCl60
L_PHE_END     = 'PP'       # plateau marker after Phe (both modes)
L_PHE_SINGLE  = 'Phe3uM'  # single-dose Phe label (single_dose mode only)
L_SUBPHE      = 'subPhe'   # submaximal Phe before ACh
L_SUBPHE_END  = 'Ach0,001' # first ACh dose (= end of subPhe window)
L_SUBPHE2     = '2subPhe'  # submaximal Phe before SNP
L_SUBPHE2_END = 'SNP0,001' # first SNP dose (= end of 2subPhe window)
L_ACH_END     = 'P3'       # rinse after ACh
L_SNP_END     = 'K'        # end of recording
```

If a label doesn't match, the notebook raises a `ValueError` listing all labels actually found in the file — use **Step 2** (comment preview) to check before running the analysis.

**Output** (Step 6, downloaded automatically):
- `<filename>_results.xlsx` — single sheet, metrics as rows, channels as columns
- `myograph_results.png` — 2×2 summary figure: Phe, ACh, SNP curves + KCl barplot

## Input file format

Standard LabChart tab-separated `.txt` export, decimal comma, 9 header lines. Each channel (1–8) is one aortic ring. Protocol events must be annotated with the `#*` prefix in LabChart (comment option) before export.

## Dependencies

All standard — no install needed in Colab. `openpyxl` is pip-installed automatically in the first cell.
