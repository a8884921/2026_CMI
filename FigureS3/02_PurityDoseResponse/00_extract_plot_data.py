import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
STEP = HERE
PROJECT_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'shared'))
try:
    import dose_response_extract as EX
except ModuleNotFoundError as exc:
    if exc.name != 'dose_response_extract':
        raise
    raise SystemExit('The required shared module dose_response_extract.py was not included in the supplied archive.') from None
upstream_dir = os.environ.get('DOSE_SOURCE_DIR', os.path.join(os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs')), 'source'))
EX.extract_panel(raw_data_dir=os.path.join(os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs')), 'source'), plot_data_dir=os.path.join(os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs')), 'plot_data'), operationalization='Local purity (malignant load)', panel_letter='A', upstream_dir=upstream_dir)
