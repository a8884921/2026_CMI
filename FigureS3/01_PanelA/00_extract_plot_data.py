import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from publication_data import numeric_id, numeric_ids, read_csv
HERE = Path(__file__).resolve().parent
INPUTS = Path(os.environ.get('SCI_INPUT_DIR', str(HERE / 'inputs')))
SOURCE = INPUTS / 'source'
PLOT_DATA = INPUTS / 'plot_data'
CELL_TYPES = ['CD4T', 'CD8T', 'DC', 'Endothelial', 'Macrophage', 'Mast', 'NK', 'Normal_B', 'Stromal', 'Mal_1', 'Mal_2', 'Mal_3']
MALIGNANT_TYPES = {'Mal_1', 'Mal_2', 'Mal_3'}
DISEASE_ORDER = ['DLBCL', 'PCNSL', 'MZL', 'BL']
WEIGHT_COLUMNS = ['RCTD_weight_' + cell_type for cell_type in CELL_TYPES]
UNASSIGNED = 'Unassigned'
files = sorted(SOURCE.glob('*_spot_master.csv'), key=lambda file: int(numeric_id(file.name.removesuffix('_spot_master.csv'))))
if len(files) != 42:
    raise ValueError('The spatial atlas requires 42 anonymized spot tables.')
parts = []
sample_rows = []
seen = set()
next_spot_id = 1
for file in files:
    sample_id = numeric_id(file.name.removesuffix('_spot_master.csv'))
    if sample_id in seen:
        raise ValueError('Each sample must have exactly one input table.')
    seen.add(sample_id)
    table = read_csv(file, usecols=['sample_id', 'disease', 'x', 'y'] + WEIGHT_COLUMNS, low_memory=False)
    if table.empty:
        raise ValueError('Spot tables cannot be empty.')
    if table['sample_id'].isna().any() or set(table['sample_id']) != {sample_id}:
        raise ValueError('The numeric input filename and sample_id column must agree.')
    if table['disease'].nunique(dropna=False) != 1 or table['disease'].iloc[0] not in DISEASE_ORDER:
        raise ValueError('Each spot table must have one recognized disease group.')
    weights = table[WEIGHT_COLUMNS].to_numpy(dtype='float64')
    all_missing = np.isnan(weights).all(axis=1)
    filled = np.where(np.isnan(weights), -np.inf, weights)
    top1 = np.array(CELL_TYPES, dtype=object)[filled.argmax(axis=1)]
    top1[all_missing | (filled.max(axis=1) == -np.inf)] = UNASSIGNED
    table = table.assign(top1_celltype=top1)
    table['spot_id'] = np.arange(next_spot_id, next_spot_id + len(table)).astype(str)
    next_spot_id += len(table)
    n = len(table)
    n_malignant = int(np.isin(top1, list(MALIGNANT_TYPES)).sum())
    disease = str(table['disease'].iloc[0])
    parts.append(table[['spot_id', 'sample_id', 'disease', 'x', 'y', 'top1_celltype']])
    sample_rows.append(dict(sample_id=sample_id, disease=disease, n_spot=n, n_malignant=n_malignant, mal_frac=n_malignant / n, n_unassigned_spots=int(all_missing.sum())))
spots = pd.concat(parts, ignore_index=True)
samples = pd.DataFrame(sample_rows)
samples['disease_order'] = samples['disease'].map({disease: i for i, disease in enumerate(DISEASE_ORDER)})
samples['sample_order'] = samples['sample_id'].map(int)
samples = samples.sort_values(['disease_order', 'mal_frac', 'sample_order'], ascending=[True, False, True]).reset_index(drop=True)
samples['panel_order'] = np.arange(1, len(samples) + 1)
panel_order = dict(zip(samples['sample_id'], samples['panel_order']))
spots['panel_order'] = spots['sample_id'].map(panel_order)
spots['spot_order'] = spots['spot_id'].map(int)
spots = spots.sort_values(['panel_order', 'spot_order'], kind='stable').drop(columns='spot_order').reset_index(drop=True)
samples = samples.drop(columns=['disease_order', 'sample_order'])
total_spots = len(spots)
class_counts = spots['top1_celltype'].value_counts().to_dict()
malignant_spots = int(sum((class_counts.get(cell_type, 0) for cell_type in MALIGNANT_TYPES)))
types = CELL_TYPES + ([UNASSIGNED] if UNASSIGNED in class_counts else [])
counts = pd.DataFrame([dict(cell_type=cell_type, N=int(class_counts.get(cell_type, 0)), frac=class_counts.get(cell_type, 0) / total_spots, is_malignant=cell_type in MALIGNANT_TYPES) for cell_type in types])
checks = {'numbers': {'total_spots': total_spots, 'malignant_spots': malignant_spots, 'malignant_fraction': malignant_spots / total_spots, 'n_samples': len(samples), 'class_counts': {cell_type: int(class_counts.get(cell_type, 0)) for cell_type in types}}}
PLOT_DATA.mkdir(parents=True, exist_ok=True)
spots.to_csv(PLOT_DATA / 'panelA_r01_spatial_atlas.csv', index=False)
samples.to_csv(PLOT_DATA / 'panelA_r02_sample_statistics.csv', index=False)
counts.to_csv(PLOT_DATA / 'panelA_r03_cell_type_counts.csv', index=False)
(PLOT_DATA / 'panelA_r05_checks.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
print('Spatial atlas data extracted for ' + str(len(samples)) + ' samples.')
