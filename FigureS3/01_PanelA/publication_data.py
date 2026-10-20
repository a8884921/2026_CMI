import os
import re
from numbers import Integral, Real

import pandas as pd

ID_COLUMNS = frozenset(('sample_id', 'patient_id', 'subject_id', 'case_id', 'spot_id', 'cell_id', 'donor_id', 'specimen_id', 'roi_id', 'aoi_id', 'fov_id'))
PRIVATE_COLUMNS = frozenset(('name', 'full_name', 'first_name', 'last_name', 'patient_name', 'subject_name', 'sample_name', 'donor_name', 'patient_number', 'subject_number', 'hospital_id', 'hospital_number', 'medical_record_number', 'record_number', 'original_id', 'original_patient_id', 'original_sample_id', 'serum_sample_id', 'urine_sample_id', 'serum_id', 'urine_id', 'barcode', 'raw_barcode', 'original_barcode', 'address', 'email', 'phone', 'date_of_birth', 'dob', 'notes', 'note', 'remarks', 'remark', 'comments', 'comment', 'source', 'source_id', 'dataset_source', 'source_path', 'local_path', 'file_path', 'source_file', 'filename', 'file_name', 'sha256', 'md5', 'size_bytes'))
_SAMPLE_SUBJECTS = {}


def normalized_field(value):
    return re.sub(r'[^a-z0-9]', '', str(value).lower())


_ID_FIELDS = {normalized_field(column): column for column in ID_COLUMNS}
_PRIVATE_FIELDS = frozenset(normalized_field(column) for column in PRIVATE_COLUMNS)


def numeric_id(value):
    if isinstance(value, bool) or pd.isna(value):
        raise ValueError('An anonymized positive integer identifier is required.')
    if isinstance(value, Integral):
        value = str(int(value))
    elif isinstance(value, Real) and float(value).is_integer():
        value = str(int(value))
    else:
        value = str(value)
    if re.fullmatch(r'[1-9][0-9]*', value) is None:
        raise ValueError('An anonymized positive integer identifier is required.')
    return value


def numeric_ids(values, allow_missing=False):
    values = pd.Series(values, copy=True)
    if not allow_missing and values.isna().any():
        raise ValueError('Anonymized identifiers cannot be missing.')
    return values.map(lambda value: pd.NA if pd.isna(value) else numeric_id(value)).astype('string')


def configured_id(name, value=None):
    selected = os.environ.get(name) if value is None else value
    if selected is None:
        raise ValueError('Set ' + name + ' to the anonymized positive integer identifier of the intended sample.')
    return numeric_id(selected)


def prepare_table(table):
    table = table.drop(columns=[column for column in table.columns if normalized_field(column) in _PRIVATE_FIELDS]).copy()
    table = table.rename(columns={column: _ID_FIELDS[normalized_field(column)] for column in table.columns if normalized_field(column) in _ID_FIELDS})
    if table.columns.duplicated().any():
        raise ValueError('The publication schema cannot contain duplicate columns.')
    if any(not str(column).isascii() for column in table.columns):
        raise ValueError('Input column names must use the English publication schema.')
    for column in table.columns:
        if str(column).lower() in ID_COLUMNS:
            table[column] = numeric_ids(table[column], allow_missing=True).array
    if normalized_field(table.index.name) in _ID_FIELDS:
        table.index = pd.Index(numeric_ids(table.index), name=_ID_FIELDS[normalized_field(table.index.name)])
    if 'sample_id' in table:
        for column in ('patient_id', 'subject_id'):
            if column in table:
                for sample_id, subject_id in table[['sample_id', column]].dropna().itertuples(index=False, name=None):
                    previous = _SAMPLE_SUBJECTS.get(sample_id)
                    if previous is not None and previous != subject_id:
                        raise ValueError('One sample cannot refer to different subjects across input tables.')
                    _SAMPLE_SUBJECTS[sample_id] = subject_id
    for column in table:
        if normalized_field(column) in ('cohort', 'cohortid', 'dataset', 'datasetid', 'cohortlabel', 'datasetlabel'):
            values = table[column].dropna().astype(str)
            if not values.str.fullmatch(r'(?:Cohort_[A-Z]+(?:_(?:DLBCL|FL|BL|MCL))?|Training|Test|All)').all():
                raise ValueError('Dataset fields must use neutral publication aliases.')
    return table


def read_csv(path, *args, **kwargs):
    dtype = kwargs.get('dtype')
    if dtype is None:
        kwargs['dtype'] = {column: 'string' for column in ID_COLUMNS}
    elif isinstance(dtype, dict):
        kwargs['dtype'] = dict(dtype, **{column: 'string' for column in ID_COLUMNS})
    try:
        table = pd.read_csv(path, *args, **kwargs)
    except FileNotFoundError:
        raise FileNotFoundError('A required anonymized input file is unavailable.') from None
    return prepare_table(table)


def read_parquet(path, *args, **kwargs):
    try:
        table = pd.read_parquet(path, *args, **kwargs)
    except FileNotFoundError:
        raise FileNotFoundError('A required anonymized input file is unavailable.') from None
    return prepare_table(table)


def validate_cohorts(table, column='group'):
    if column not in table:
        return table
    values = table[column].dropna().astype(str)
    if not values.str.fullmatch(r'(?:Cohort_[A-Z]+(?:_(?:DLBCL|FL|BL|MCL))?|Training|Test|All)').all():
        raise ValueError('Cohort identifiers must use the neutral publication aliases.')
    table = table.copy()
    if 'label' in table:
        table['label'] = table[column].map(lambda value: '_'.join(str(value).split('_')[:2]) if str(value).startswith('Cohort_') else str(value))
    return table
