from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import os
import numpy as np
import pandas as pd
_HERE = os.path.dirname(os.path.abspath(__file__))

def _first_dir(*cands):
    for c in cands:
        if c and os.path.isdir(c):
            return c
    return next((c for c in cands if c))
INDEX = os.environ.get('SCI_INPUT_DIR', os.path.join(_HERE, 'inputs'))
F5 = _first_dir(os.environ.get('FIG5_ROOT'), os.path.join(INDEX, 'figure5_data'))
S5 = _first_dir(os.environ.get('FIGS3_ROOT'), os.path.join(os.path.dirname(os.path.dirname(_HERE)), 'FigureS3'))
ASSETS = os.environ.get('FIG5_ASSETS', os.path.join(INDEX, 'assets'))

def p(*a):
    return os.path.join(F5, *a)
_spots = None

def spots():
    global _spots
    if _spots is None:
        _spots = read_parquet(p('shared', 'c01_r01_spot_analysis.parquet'))
    return _spots

def territory_outline(state, sample):
    v = read_csv(p('08_panelH_territory', 'inputs', 'source', 'c02.4_r01_territory_vertices.csv'))
    v = v[(v.state == state) & (v.sample_id == sample)].sort_values('vertex_order')
    return (v.x.values, v.y.values)

def territory_rings():
    return read_parquet(p('08_panelH_territory', 'inputs', 'source', 'c02.4_r02_spot_territory_distances.parquet'))

def territory_summary():
    return read_csv(p('08_panelH_territory', 'inputs', 'source', 'c02.4_r03_territory_summary.csv'))

def composition():
    d = read_csv(p('01_panelA_state_composition', 'results', 'panelA_state_composition_values.csv'))
    return d[d['row_type'] == 'bar'].sort_values('plot_order').reset_index(drop=True)

def moran():
    d = read_csv(p('02_panelB_moran_i', 'results', 'panelB_moran_i_values.csv'))
    return d

def self_z():
    d = read_csv(p('03_panelC_self_enrichment', 'results', 'panelC_self_enrichment_values.csv'))
    return d

def nhood3x3():
    d = read_csv(p('04_panelD_neighborhood_enrichment', 'inputs', 'source', 'c03_r08_figure3_plot_data.csv'))
    d = d[d.facet == 'ALL (pooled)']
    assert len(d) == 9
    return d

def cooccur():
    d = read_csv(p('05_panelE_cooccurrence', 'inputs', 'source', 'c07.1_r02_supplement_pooled_plot_data.csv'))
    return d[d.caliber == 'top1']

def tme_composition():
    return read_csv(p('09_panelI_tme_neighborhood', 'inputs', 'plot_data', 'c03_r06_figure1A_composition_plot_data.csv'))

def tme_access():
    d = read_csv(p('10_panelJ_abundance_delta_z', 'inputs', 'source', 'c05_r03_subject_paired_tests.csv'))
    return d[(d.disease == 'DLBCL') & (d.caliber == 'top1')].set_index('tme_cell')

def importance():
    return read_csv(p('11_panelK_joint_model_importance', 'inputs', 'source', 'c07.3_r05_joint_model_importance.csv'))

def roc():
    return read_csv(p('12_panelL_roc_comparison', 'inputs', 'plot_data', 'c07.3_r02_roc_data.csv'))

def roc_perf():
    return read_csv(p('12_panelL_roc_comparison', 'inputs', 'source', 'c07.3_r01_performance.csv'))

def visiumhd_values():
    return read_csv(p('13_panelM_visium_hd_transformation', 'results', 'panelM_visium_hd_transformation_values.csv'))

def geomx_pairs():
    return read_csv(p('14_panelN_geomx_pairs', 'inputs', 'plot_data', 'c04_r01_subject_pairs.csv'))

def geomx_tests():
    return read_csv(p('14_panelN_geomx_pairs', 'inputs', 'source', 'c04_r02_paired_tests.csv')).set_index('indicator')

def geomx_models():
    return read_csv(p('15_panelO_geomx_effects', 'inputs', 'plot_data', 'c05_r01_model_effects.csv'))

def cosmx_deciles():
    return read_csv(p('17_panelQ_cosmx_neighborhood_trends', 'inputs', 'plot_data', 'c06_r04_panelB_plot_data.csv'))

def cosmx_models():
    return read_csv(p('17_panelQ_cosmx_neighborhood_trends', 'inputs', 'plot_data', 'c05.1_r01_CMI_adjusted_effects.csv'))

def cosmx_cells(cols=None):
    return read_csv(p('17_panelQ_cosmx_neighborhood_trends', 'inputs', 'source', 'c04.2_r01_tumor_B_neighborhood_access.csv.gz'), usecols=cols)

def dose_base():
    return read_csv(os.path.join(S5, '02_PurityDoseResponse', 'inputs', 'source', 'c01_r01_analysis_table.csv'))

def dose_points():
    return read_csv(os.path.join(S5, '02_PurityDoseResponse', 'inputs', 'source', 'c02_r09_figure1A_scatter_data.csv'))

def dose_stats():
    return read_csv(os.path.join(S5, '02_PurityDoseResponse', 'inputs', 'source', 'c02_r10_scatter_slopes_explained_variance.csv'))

def atlas():
    return read_csv(os.path.join(S5, '01_PanelA', 'inputs', 'plot_data', 'panelA_r01_spatial_atlas.csv'))

def atlas_samples():
    return read_csv(os.path.join(S5, '01_PanelA', 'inputs', 'plot_data', 'panelA_r02_sample_statistics.csv'))
_terr_all = None

def territory_sizes_all():
    global _terr_all
    if _terr_all is not None:
        return _terr_all
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import NearestNeighbors
    S = spots()
    rows = []
    for sid, d in S.groupby('sample_id'):
        xy = d[['x', 'y']].values
        n = len(d)
        _, idx = NearestNeighbors(n_neighbors=7).fit(xy).kneighbors(xy)
        r = np.repeat(np.arange(n), 6)
        c = idx[:, 1:].ravel()
        A = coo_matrix((np.ones(len(r)), (r, c)), shape=(n, n)).tocsr()
        A = (A + A.T > 0).astype(int)
        lab = d.CMI_state_top1.values
        for st in ['Mal_1', 'Mal_2', 'Mal_3']:
            m = np.where(lab == st)[0]
            if len(m) == 0:
                continue
            _, labs = connected_components(A[m][:, m], directed=False)
            rows.append(dict(sample_id=sid, disease=d.disease.iloc[0], state=st, n_state=len(m), n_territory=int(np.bincount(labs).max())))
    _terr_all = pd.DataFrame(rows)
    _terr_all['share'] = _terr_all.n_territory / _terr_all.n_state
    return _terr_all
