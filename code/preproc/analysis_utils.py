import pandas as pd
import os
import numpy as np

def get_mri_use_list(participants_tsv_path):
    """
    Takes in a path to the omni-filled BIDS participants.tsv and returns lists of people to use and people to exclude (MRI data).

    For now we are excluding based on 'censor_runs-all' in participats.tsv

    Input:
        participants_tsv_path: full path (including filename) to a bids-y participants.tsv that has omni-expected columns (see below)
        example: participants_tsv_path = '/Volumes/data/HMM/HMM_MRI/omni_BIDS/derivatives/fmriprep/participants.tsv'

    Output:
        use_peeps_bids: list of bids subject numbers to use
        use_peeps_omni: list of PRISMAF style subject numbers to use (as used in the main dataframe)

        exclude_peeps_bids: list of bids subject numbers to EXCLUDE 
        exclude_peeps_omni: list of PRISMAF style subject numbers to exclude (as used in the main dataframe)
    """

    if not os.path.exists(participants_tsv_path):
        raise ValueError(
            participants_tsv_path + " doesn" "t seem to be a valid path and/or filename"
        )

    else:
        participants_df = pd.read_csv(participants_tsv_path, sep="\t")

    # check for expected censoring columns -- these are omni-specific and don't exist in base bids participants.tsv
    censor_cols = [
        "censor_run-all",
        "censor_run-01",
        "censor_run-02",
        "censor_run-03",
        "censor_run-04",
        "censor_run-05",
        "NOTES",
    ]

    for c in censor_cols:
        if c not in participants_df.columns:
            raise ValueError("participants_df doesn" "t have expected column " + c)
    # ----------------

    use_peeps_bids = participants_df.loc[
        participants_df["censor_run-all"] == 0, "participant_id"
    ].values.tolist()
    exclude_peeps_bids = participants_df.loc[
        participants_df["censor_run-all"] == 1, "participant_id"
    ].values.tolist()

    use_peeps_omni = [p.replace("sub-", "PRISMAF") for p in use_peeps_bids]
    exclude_peeps_omni = [p.replace("sub-", "PRISMAF") for p in exclude_peeps_bids]

    return use_peeps_bids, use_peeps_omni, exclude_peeps_bids, exclude_peeps_omni


def flatten_cols(index, prefix, col_names, sep="_"):
    new_cols = [prefix + sep + str(t) for t in col_names]
    return list(index) + new_cols

def block_z_rel(df):
    x = np.array(df['min_study_block'])
    y = np.array(df['fisher_z'])
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    r = np.corrcoef(x, y)[0, 1]
    last_first = y[-1] - y[0]
    max_min = max(y) - min(y)
    rel_df = pd.DataFrame({'slope': [m], 'corr': [r], 'last_first': [last_first], 'max_min': [max_min]})
    return rel_df

def lower_triangular_corr(corr_mat):
    mask = np.tril(np.ones_like(corr_mat, dtype=np.bool), k=-1)
    return corr_mat.where(mask)

def lower_triangular_corr_filter(df, filter_cols):
    corr_mat = df.drop(columns=filter_cols)
    mask = np.tril(np.ones_like(corr_mat, dtype=np.bool), k=-1)
    corr_mat_lt = corr_mat.where(mask)
    corr_df = pd.concat([df[filter_cols], corr_mat_lt], axis=1)
    return corr_df

#[au.flatten_col(col) for col in df.columns]
def flatten_col(multi_col):
    if multi_col[1] != '':
        if type(multi_col) == tuple:
            new_col = '_'.join([str(col) for col in multi_col])
    else:
        new_col = str(multi_col[0])
    return new_col

def get_pattern_corr_wc_df(pattern_df):
    # get df in position to compute correlations
    pattern_df_pivot = pd.pivot_table(pattern_df, values='vox_val', columns=['word', 'run'], index=['sub', 'roi', 'voxnum'])
    pattern_df_pivot.columns = pattern_df_pivot.columns.to_flat_index()
    pattern_df_pivot.columns = [t[0] + '_' + str(t[1]) for t in pattern_df_pivot.columns]
    pattern_df_pivot = pattern_df_pivot.reset_index()
    pattern_df_pivot
    
    #compute correlation
    pattern_corrs = pattern_df_pivot.drop(columns='voxnum').groupby(['sub', 'roi']).corr()

    #get only lower triangle so each correlation only shows up once
    pattern_corrs_lt = pattern_corrs.groupby(['sub', 'roi']).apply(lower_triangular_corr)
    pattern_corrs_lt_ri = pattern_corrs_lt.reset_index().rename(columns={'level_2': 'word_run'})

    # melt into tidy frame
    pattern_corrs_melt = pattern_corrs_lt_ri.melt(id_vars=['sub', 'roi', 'word_run'], value_vars=pattern_df_pivot.columns[3:], var_name='word_run_2', value_name='corr')
    pattern_corrs_melt = pattern_corrs_melt.dropna()

    #split word and run into multiple columns
    word_run_1 = pattern_corrs_melt['word_run'].str.split(pat="_", expand=True)
    pattern_corrs_melt['word_1'] = word_run_1[0]
    pattern_corrs_melt['run_1'] = word_run_1[1]
    pattern_corrs_melt['run_1'] = pattern_corrs_melt['run_1'].astype(int)
    word_run_2 = pattern_corrs_melt['word_run_2'].str.split(pat="_", expand=True)
    pattern_corrs_melt['word_2'] = word_run_2[0]
    pattern_corrs_melt['run_2'] = word_run_2[1]
    pattern_corrs_melt['run_2'] = pattern_corrs_melt['run_2'].astype(int)

    #split subject
#     subject = pattern_corrs_melt['sub'].str.split(pat="-", expand=True)
#     pattern_corrs_melt['subject'] = subject[1].astype(int)

    #create type of correlation
    pattern_corrs_melt['within_cross'] = 'Cross'
    pattern_corrs_melt.loc[pattern_corrs_melt['word_1'] == pattern_corrs_melt['word_2'], 'within_cross'] = 'Within'
    return pattern_corrs_melt

def get_pattern_corr_df(pattern_df, values='vox_val', 
                        columns=['word', 'run'], 
                        group_index=['sub', 'roi'], 
                        corr_index=['voxnum'],
                       lower_triangular=True):

    # get df in position to compute correlations
    if len(columns) == 1:
        pattern_df_pivot = pattern_df.pivot_table(values=[values], columns=columns, 
                                      index=(group_index+corr_index)).reset_index() #maybe whether values is in a list depends on # of columns?
    else: 
        pattern_df_pivot = pattern_df.pivot_table(values=values, columns=columns, 
                                      index=(group_index+corr_index)).reset_index() #maybe whether values is in a list depends on # of columns?
    pattern_df_pivot.columns = [flatten_col(col) for col in pattern_df_pivot.columns]
    
    #compute correlation
    pattern_corrs = pattern_df_pivot.drop(columns=corr_index).groupby(group_index).corr()
    
    #get only lower triangle so each correlation only shows up once
    if lower_triangular:
        pattern_corrs = pattern_corrs.groupby(group_index).apply(lower_triangular_corr)
    colname = '_'.join(columns)
    colname_1 = colname + '_1'
    if len(group_index) == 1:
        old_colname = 'level_1'
    else:
        old_colname = 'level_2'
    pattern_corrs_ri = pattern_corrs.reset_index().rename(columns={old_colname: colname_1})

    # melt into tidy frame
    colname_2 = colname + '_2'
    id_vars = group_index + [colname_1]
    value_vars = list(pattern_corrs_ri.columns[len(id_vars):])
    pattern_corrs_melt = pattern_corrs_ri.melt(id_vars=id_vars, value_vars=value_vars, 
                                                  var_name=colname_2, value_name='corr')
    pattern_corrs_melt = pattern_corrs_melt.dropna()
    
    #split word and run into multiple columns
    if len(columns) == 2:
        temp_col_1 = pattern_corrs_melt[colname_1].str.split(pat="_", expand=True)
        pattern_corrs_melt[columns[0] + '_1'] = temp_col_1[0]
        pattern_corrs_melt[columns[1] + '_1'] = temp_col_1[1]
        pattern_corrs_melt[columns[1] + '_1'] = pattern_corrs_melt[columns[1] + '_1'].astype(int)
        temp_col_2 = pattern_corrs_melt[colname_2].str.split(pat="_", expand=True)
        pattern_corrs_melt[columns[0] + '_2'] = temp_col_2[0]
        pattern_corrs_melt[columns[1] + '_2'] = temp_col_2[1]
        pattern_corrs_melt[columns[1] + '_2'] = pattern_corrs_melt[columns[1] + '_2'].astype(int)
    else:
        pattern_corrs_melt[colname_1] = pattern_corrs_melt[colname_1].str.replace(values+'_', '').astype(int)
        pattern_corrs_melt[colname_2] = pattern_corrs_melt[colname_2].str.replace(values+'_', '').astype(int)
    return pattern_corrs_melt

def create_long_activity_df(rois, base_dir, var_group, roi_col='roi', sub_col='sub_id', study_block_col='study_block', 
    item_col='lith_word', value_col='value', lat_rois=['bilat']):
    activity_dfs = []
    for roi in rois:
        for lat in lat_rois:
            if len(lat_rois) > 1:
                roi_name = roi + '_' + lat
            else:
                roi_name = roi
            print(roi_name)
            pattern_file = base_dir + var_group + '_' + str(roi_name) + '_df.csv'
            pattern_df = pd.read_csv(pattern_file)
            if len(lat_rois) > 1:
                pattern_df[roi_col] = roi_name
            activity_df = pattern_df.groupby([sub_col, study_block_col, roi_col, item_col]).agg({value_col: 'mean'}).reset_index()
            activity_dfs.append(activity_df)

    activity_df_long = pd.concat(activity_dfs)

    if roi_col != 'roi':
        activity_df_long.rename(columns={roi_col: 'roi'}, inplace=True)
    activity_df_long['var'] = 'activity'
    activity_df_long['var_group'] = var_group

    return activity_df_long

def create_lt_pattern_sim_df(rois, base_dir, var_group, roi_col='roi', sub_col='sub_id', study_block_col='study_block', 
    item_col='lith_word', value_col='value', id_col='xyz', lat_rois=['bilat']):
    pattern_corr_dfs = []
    for roi in rois:
        for lat in lat_rois:
            if len(lat_rois) > 1:
                roi_name = roi + '_' + lat
            else:
                roi_name = roi
            print(roi_name)
            pattern_file = base_dir + var_group + '_' + str(roi_name) + '_df.csv'
            pattern_df = pd.read_csv(pattern_file)
            if len(lat_rois) > 1:
                pattern_df[roi_col] = roi_name
            pattern_df[id_col] = pattern_df['x'].astype(str) + '_' + pattern_df['y'].astype(str) + '_' + pattern_df['z'].astype(str)
            pattern_corr_df = get_pattern_corr_df(pattern_df, 
                                                     columns=[item_col, study_block_col],
                                                     group_index=[sub_col, roi_col],
                                                    values=value_col,
                                                    corr_index=[id_col])
            pattern_corr_df = pattern_corr_df[pattern_corr_df[study_block_col + '_1'] != pattern_corr_df[study_block_col + '_2']]
            pattern_corr_df = pattern_corr_df.drop(columns=[(item_col + '_' + study_block_col + '_1'), (item_col + '_' + study_block_col + '_2')])
            pattern_corr_dfs.append(pattern_corr_df)

    # n entries is (((n_words * n_reps) ^ 2 - n_reps * (n_words ^ 2)) * n_rois * n_subs) / 2 because subtracting off same rep and lower triagular
    pattern_corr_df_long_lt = pd.concat(pattern_corr_dfs)

    pattern_corr_df_long_lt['within_cross'] = 'Cross'
    pattern_corr_df_long_lt.loc[pattern_corr_df_long_lt[(item_col + '_1')] == pattern_corr_df_long_lt[(item_col + '_2')], 'within_cross'] = 'Within'
    pattern_corr_df_long_lt['fisher_z'] = np.arctanh(pattern_corr_df_long_lt['corr'])
    if roi_col != 'roi':
        pattern_corr_df_long_lt.rename(columns={roi_col: 'roi'}, inplace=True)

    return pattern_corr_df_long_lt

def rename_pattern_cols(pattern_corr_df_longer, roi_col='roi', sub_col='sub_id', study_block_col='study_block', 
    item_col='lith_word', use_roi_col=True):
    pattern_corr_df_longer[('min_' + study_block_col)] = pattern_corr_df_longer[[(study_block_col + '_1'), (study_block_col + '_2')]].min(axis=1)
    pattern_corr_df_longer[('max_' + study_block_col)] = pattern_corr_df_longer[[(study_block_col + '_1'), (study_block_col + '_2')]].max(axis=1)

    pattern_corr_df_longer[('target_' + study_block_col)] = pattern_corr_df_longer[(study_block_col + '_1')]
    pattern_corr_df_longer.loc[pattern_corr_df_longer['item_id'] == 2, ('target_' + study_block_col)] = pattern_corr_df_longer[(study_block_col + '_2')]
    pattern_corr_df_longer[('opp_' + study_block_col)] = pattern_corr_df_longer[(study_block_col + '_1')]
    pattern_corr_df_longer.loc[pattern_corr_df_longer['item_id'] == 1, ('target_' + study_block_col)] = pattern_corr_df_longer[(study_block_col + '_2')]

    pattern_corr_df_longer[study_block_col] = pattern_corr_df_longer[('target_' + study_block_col)].astype(str) + '_' + pattern_corr_df_longer[('opp_' + study_block_col)].astype(str)
    pattern_corr_df_longer.loc[(pattern_corr_df_longer['within_cross'] == 'Within'), study_block_col] = pattern_corr_df_longer[('min_' + study_block_col)].astype(str) + '_' + pattern_corr_df_longer[('max_' + study_block_col)].astype(str)
    
    if use_roi_col: 
        groupby_cols = [sub_col, roi_col, item_col, 'within_cross', study_block_col]
    else:
        groupby_cols = [sub_col, item_col, 'within_cross', study_block_col]

    pattern_corr_df_long = pattern_corr_df_longer.groupby(groupby_cols).agg(
        {'fisher_z': 'mean'}).reset_index()
    
    return pattern_corr_df_long

def create_long_pattern_sim_df(pattern_corr_df_long_lt, var_group, roi_col='roi', sub_col='sub_id', study_block_col='study_block', 
    item_col='lith_word', value_col='value', use_chunks=False, use_roi_col=True):
    print('start')
    if use_chunks:
        print('use_chunks')
        pattern_corr_dfs = []
        grouped = pattern_corr_df_long_lt.groupby(roi_col)
        for roi, roi_pattern_corr_df_long_lt in grouped:
            print(roi)
            pattern_corr_df_longer = pd.wide_to_long(roi_pattern_corr_df_long_lt, 
                        stubnames=item_col, 
                        i=['within_cross', sub_col, 'fisher_z'], 
                        j='item_id', 
                        sep='_').reset_index()

            pattern_corr_df = rename_pattern_cols(pattern_corr_df_longer, roi_col=roi_col, sub_col=sub_col, 
                                                   study_block_col=study_block_col, item_col=item_col)
            pattern_corr_dfs.append(pattern_corr_df)
        pattern_corr_df_long = pd.concat(pattern_corr_dfs)
    else:
        if use_roi_col:
            i = [roi_col,'within_cross', sub_col, 'fisher_z']
        else:
            i = ['within_cross', sub_col, 'fisher_z']
        pattern_corr_df_longer = pd.wide_to_long(pattern_corr_df_long_lt, 
                    stubnames=item_col, 
                    i=i, 
                    j='item_id', 
                    sep='_').reset_index()

        pattern_corr_df_long = rename_pattern_cols(pattern_corr_df_longer, roi_col=roi_col, sub_col=sub_col, 
                                                   study_block_col=study_block_col, item_col=item_col, use_roi_col=use_roi_col)
   
    pattern_corr_df_long.rename(columns={'fisher_z': 'value'}, inplace=True)
    pattern_corr_df_long['var'] = 'pattern_sim_z_' + pattern_corr_df_long['within_cross'].astype(str)
    pattern_corr_df_long.drop(columns=['within_cross'], inplace=True)
    pattern_corr_df_long['var_group'] = var_group
    return pattern_corr_df_long

def create_long_ISC_df(rois, base_dir, var_group, roi_col='roi', sub_col='sub_id', study_block_col='study_block', 
    item_col='lith_word', value_col='value', id_col='xyz', lat_rois=['bilat']):
    ISC_dfs = []
    for roi in rois:
        for lat in lat_rois:
            if len(lat_rois) > 1:
                roi_name = roi + '_' + lat
            else:
                roi_name = roi
            pattern_file = base_dir + var_group + '_' + str(roi_name) + '_df.csv'
            pattern_df = pd.read_csv(pattern_file)
            if len(lat_rois) > 1:
                pattern_df[roi_col] = roi_name
            print(roi)
            pattern_df_pivot = pd.pivot_table(pattern_df, 
                                              values=value_col, 
                                              index=[item_col, study_block_col, roi_col, 'x', 'y', 'z'],
                                              columns=[sub_col]).reset_index()
            pattern_df_pivot.drop(columns=['x', 'y', 'z'], inplace=True)
            ISC_df = pattern_df_pivot.groupby(
                [item_col, study_block_col, roi_col]).corr().reset_index()
            ISC_melt = ISC_df.melt(
                id_vars=[item_col, roi_col, study_block_col, sub_col],
                value_vars=pattern_df_pivot.columns[4:], 
                var_name=(sub_col + '_2'),
                value_name='ISC').dropna()
            ISC_dfs.append(ISC_melt)
    ISC_df_long = pd.concat(ISC_dfs)
    ISC_df_long['fisher_z'] = np.arctanh(ISC_df_long['ISC'])
    ISC_df_long.query(sub_col + ' != ' + (sub_col + '_2'), inplace=True)
    ISC_df_long = ISC_df_long.groupby([item_col, roi_col, sub_col, study_block_col]).agg({'fisher_z': 'mean'}).reset_index()
    ISC_df_long.rename(columns={'fisher_z': 'value', roi_col: 'roi'}, inplace=True)
    ISC_df_long['var'] = 'ISC_z'
    ISC_df_long['var_group'] = var_group
    return ISC_df_long
