import pandas as pd
import numpy as np
import re

def preprocess_data(file_path):
    """
    Preprocesses the longitudinal proteomics dataset to generate:
    - working_df_all: all patients with corrected statuses and slopes
    - working_df: patients who transition from Normal to Abnormal
    - df_normal_only: patients who stay Normal throughout
    - df_abnormal_only: patients who stay Abnormal throughout

    Parameters
    ----------
    file_path : str
        Path to the input CSV file containing the raw data with onset information.

    Returns
    -------
    dict
        A dictionary with the following keys:
        - 'working_df_all': DataFrame with all patients, corrected Status_NTA
        - 'working_df': Status change patients cleaned to remove initially abnormal cases
        - 'df_normal_only': Normal-only patients
        - 'df_abnormal_only': Abnormal-only patients
    """
    df = pd.read_csv(file_path)

    # Clean column names
    df.columns = df.columns.map(
        lambda col: re.sub(r'[\s.\-|]+', '_', col) if any(char in col for char in " .-|") else col
    )

    protein_list = df.columns[14:-3].tolist()
    filtered_data = df[
        ['SUBID','PROCEDURE_DATE', 'BASELINE_AGE', 'DECAGE', 'PROCEDURE_AGE',
         'SEX','CATEGORY','Status_H', 'TRANSITION_AGE', 'ONSET_AGE', 'Status'] + protein_list
    ].copy()

    filtered_data['DECAGE_Mutated'] = filtered_data['DECAGE']

    # Identify Normal-only patients and mutate DECAGE
    normal_patients = (
        filtered_data.groupby('SUBID')
        .filter(lambda x: set(x['Status_H']) == {'Normal'})
        ['SUBID']
        .unique()
    )
    for subid in normal_patients:
        last_procedure_age = filtered_data.loc[filtered_data['SUBID'] == subid, 'PROCEDURE_AGE'].max()
        filtered_data.loc[filtered_data['SUBID'] == subid, 'DECAGE_Mutated'] = last_procedure_age

    filtered_data['ONSET_AGE'] = filtered_data['DECAGE_Mutated'].fillna(filtered_data['TRANSITION_AGE'])
    filtered_data['Status_NTA'] = filtered_data['Status_H']

    # Initial status corrections
    filtered_data.loc[
        (filtered_data['ONSET_AGE'] > filtered_data['PROCEDURE_AGE']) &
        (filtered_data['Status_NTA'] == 'Abnormal'),
        'Status_NTA'
    ] = 'Normal'

    filtered_data.loc[
        (filtered_data['ONSET_AGE'] < filtered_data['PROCEDURE_AGE']) &
        (filtered_data['Status_NTA'] == 'Normal'),
        'Status_NTA'
    ] = 'Abnormal'

    filtered_data = filtered_data.sort_values(by=['SUBID', 'PROCEDURE_DATE'])

    def correct_status(df_group):
        abnormal_flag = False
        for idx, row in df_group.iterrows():
            if row['Status_NTA'] == 'Abnormal':
                abnormal_flag = True
            if abnormal_flag:
                df_group.at[idx, 'Status_NTA'] = 'Abnormal'
        return df_group

    filtered_data = filtered_data.groupby('SUBID', group_keys=False).apply(correct_status)
    working_df_all = filtered_data.copy()

    def identify_status_change(df):
        df_sorted = df.sort_values(by=["SUBID", "PROCEDURE_DATE"])
        return df_sorted.groupby("SUBID").filter(
            lambda x: (x["Status_NTA"] == "Normal").any() and
                      (x["Status_NTA"] != "Normal").any()
        )

    changed_status_df = identify_status_change(filtered_data)

    def remove_initial_abnormal(df):
        df_sorted = df.sort_values(by=["SUBID", "PROCEDURE_DATE"])
        return df_sorted.groupby("SUBID").filter(lambda group: group.iloc[0]["Status_NTA"] == "Normal")

    changed_status_df_cleaned = remove_initial_abnormal(changed_status_df)
    working_df = changed_status_df_cleaned.copy()

    # Create piecewise terms
    for df_ in [working_df, working_df_all]:
        df_['piecewise_age'] = df_['PROCEDURE_AGE'] - df_['ONSET_AGE']
        df_['years_since_onset'] = df_['PROCEDURE_AGE'] - df_['ONSET_AGE']
        df_.dropna(subset=['piecewise_age', 'years_since_onset'], inplace=True)

    working_df_all.loc[
        (working_df_all['ONSET_AGE'] > working_df_all['PROCEDURE_AGE']) &
        (working_df_all['Status_NTA'] == 'Abnormal'),
        'Status_NTA'
    ] = 'Normal'

    working_df_all.loc[
        (working_df_all['ONSET_AGE'] < working_df_all['PROCEDURE_AGE']) &
        (working_df_all['Status_NTA'] == 'Normal'),
        'Status_NTA'
    ] = 'Abnormal'

    # Subsetting Normal-only and Abnormal-only
    normal_patients = working_df_all.groupby('SUBID')['Status_NTA'].apply(lambda x: (x == 'Normal').all())
    abnormal_patients = working_df_all.groupby('SUBID')['Status_NTA'].apply(lambda x: (x != 'Normal').all())

    df_normal_only = working_df_all[working_df_all['SUBID'].isin(normal_patients[normal_patients].index)]
    df_abnormal_only = working_df_all[working_df_all['SUBID'].isin(abnormal_patients[abnormal_patients].index)]

    return {
        "working_df_all": working_df_all,
        "working_df": working_df,
        "df_normal_only": df_normal_only,
        "df_abnormal_only": df_abnormal_only
    }
