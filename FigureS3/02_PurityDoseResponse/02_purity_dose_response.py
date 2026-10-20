import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
STEP = HERE
PROJECT_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'shared'))
try:
    import dose_response_common as CO
except ModuleNotFoundError as exc:
    if exc.name != 'dose_response_common':
        raise
    raise SystemExit('The required shared module dose_response_common.py was not included in the supplied archive.') from None
CO.draw_panel(plot_data_dir=os.path.join(os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs')), 'plot_data'), results_dir=os.environ.get('SCI_OUTPUT_DIR', os.path.join(HERE, 'results')), figure_name='02_Panel_purity_dose_response', operationalization='Local purity (malignant load)', panel_letter='A', PROJECT_ROOT=PROJECT_ROOT)
