import pickle
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

work_dir = 'CIMA/scenicplus/'
tmp_dir = 'CIMA/Temp/scenicplus/'

if not os.path.exists(work_dir):
    os.makedirs(work_dir)
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

cistopic_obj = pickle.load(open(work_dir+'pycisTopic/CIMA_scATAC_PeakRecalling_Metacell_cisTopic_obj.pkl', 'rb'))


# Binarize Topics
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)


# Calculating DARs
# Calculating DARs for each cell_type_l3
markers_dict_cell_type_l3 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l3',
                                           var_features=variable_regions, n_cpu=95,**{'_temp_dir':tmp_dir,'_system_config':{
        # Allow spilling until the local disk is 99% utilized.
        # This only affects spilling to the local file system.
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})


# Calculating DARs for age_group
cistopic_obj.cell_data['age_group'] = cistopic_obj.cell_data['age'].apply(lambda x: 'young' if x < 50 else 'old')

def age_sex_group(row):
    if row['age_group'] == 'young':
        return 'YM' if row['sex'] == 'Male' else 'YF'
    else:
        return 'OM' if row['sex'] == 'Male' else 'OF'
    
cistopic_obj.cell_data['age_sex_group'] = cistopic_obj.cell_data.apply(age_sex_group, axis=1)

def contrasts_group(cell_types, age_sex_groups):
    contrasts = []
    for cell_type in cell_types:
        for group in age_sex_groups:
            inner_list = []
            for age_sex_group in group:
                inner_list.append([f'{cell_type}_{age_sex_group}'])
            contrasts.append(inner_list)
    return contrasts

# Calculating DARs for age_group in each cell_type_l3
cistopic_obj.cell_data['cell_type_l3_age_sex_group'] = cistopic_obj.cell_data['cell_type_l3'].astype(str) + '_' + cistopic_obj.cell_data['age_sex_group'].astype(str)
cell_types = list(cistopic_obj.cell_data['cell_type_l3'].unique())  # get all cell type l3

# Young Female vs. Old Female
age_sex_groups = [['YF','OF'], ['OF','YF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YF_v_OF_l3 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l3_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Young Male vs. Old Male
age_sex_groups = [['YM','OM'], ['OM','YM']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YM_v_OM_l3 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l3_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Old Female vs. Old Male
age_sex_groups = [['OF','OM'], ['OM','OF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_OF_v_OM_l3 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l3_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Young Female vs. Young Male
age_sex_groups = [['YF','YM'], ['YM','YF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YF_v_YM_l3 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l3_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})


# Calculating DARs for each cell_type_l4
markers_dict_cell_type_l4 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Calculating DARs for age_group in each cell_type_l4
cistopic_obj.cell_data['cell_type_l4_age_sex_group'] = cistopic_obj.cell_data['cell_type_l4'].astype(str) + '_' + cistopic_obj.cell_data['age_sex_group'].astype(str)
cell_types = list(cistopic_obj.cell_data['cell_type_l4'].unique())  # get all cell type l4
# remove cell_type_l4 which metacell number less than 50
cell_types.remove('pre-Switched_Bm_IFIT3')
cell_types.remove('NKT_IFNG')
cell_types.remove('Unswitched_Bm_IL6')

# Young Female vs. Old Female
age_sex_groups = [['YF','OF'], ['OF','YF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YF_v_OF_l4 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l4_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Young Male vs. Old Male
age_sex_groups = [['YM','OM'], ['OM','YM']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YM_v_OM_l4 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l4_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Old Female vs. Old Male
age_sex_groups = [['OF','OM'], ['OM','OF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_OF_v_OM_l4 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l4_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Young Female vs. Young Male
age_sex_groups = [['YF','YM'], ['YM','YF']]
contrasts = contrasts_group(cell_types, age_sex_groups)
markers_dict_YF_v_YM_l4 = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type_l4_age_sex_group', var_features=variable_regions,n_cpu=96, contrasts = contrasts,
                                         **{'_temp_dir':tmp_dir, '_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})


# Calculating DARs for cell type l4 in cell type l1
# Calculating DARs for cell type l4 in B cells
B_cell = CistopicObject.subset(cistopic_obj,cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='B'].index,copy=True)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='B'].index, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict_B_l4 = find_diff_features(B_cell, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Calculating DARs for cell type l4 in CD4 T cells
CD4_cell = CistopicObject.subset(cistopic_obj,cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='CD4T'].index,copy=True)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='CD4T'].index, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict_CD4_l4 = find_diff_features(CD4_cell, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Calculating DARs for cell type l4 in CD8 T cells
CD8_cell = CistopicObject.subset(cistopic_obj,cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='CD8T'].index,copy=True)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='CD8T'].index, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict_CD8_l4 = find_diff_features(CD8_cell, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Calculating DARs for cell type l4 in NK cells
NK_cell = CistopicObject.subset(cistopic_obj,cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='NK'].index,copy=True)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='NK'].index, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict_NK_l4 = find_diff_features(NK_cell, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})

# Calculating DARs for cell type l4 in Myeloid cells
Myeloid_cell = CistopicObject.subset(cistopic_obj,cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='Myeloid'].index,copy=True)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=cistopic_obj.cell_data[cistopic_obj.cell_data['cell_type_l1']=='Myeloid'].index, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict_Myeloid_l4 = find_diff_features(Myeloid_cell, imputed_acc_obj, variable='cell_type_l4',
                                           var_features=variable_regions, n_cpu=96,**{'_temp_dir':tmp_dir,'_system_config':{
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"})


# save DARs results
if not os.path.exists(os.path.join(work_dir, 'candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'candidate_enhancers'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))

pickle.dump(markers_dict_cell_type_l3, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_cell_type_l3.pkl'), 'wb'))
pickle.dump(markers_dict_YF_v_OF_l3, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_OF_l3.pkl'), 'wb'))
pickle.dump(markers_dict_YM_v_OM_l3, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YM_v_OM_l3.pkl'), 'wb'))
pickle.dump(markers_dict_OF_v_OM_l3, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_OF_v_OM_l3.pkl'), 'wb'))
pickle.dump(markers_dict_YF_v_YM_l3, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_YM_l3.pkl'), 'wb'))

pickle.dump(markers_dict_cell_type_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_cell_type_l4.pkl'), 'wb'))
pickle.dump(markers_dict_YF_v_OF_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_OF_l4.pkl'), 'wb'))
pickle.dump(markers_dict_YM_v_OM_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YM_v_OM_l4.pkl'), 'wb'))
pickle.dump(markers_dict_OF_v_OM_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_OF_v_OM_l4.pkl'), 'wb'))
pickle.dump(markers_dict_YF_v_YM_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_YM_l4.pkl'), 'wb'))

pickle.dump(markers_dict_B_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_B_l4.pkl'), 'wb'))
pickle.dump(markers_dict_CD4_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_CD4_l4.pkl'), 'wb'))
pickle.dump(markers_dict_CD8_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_CD8_l4.pkl'), 'wb'))
pickle.dump(markers_dict_NK_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_NK_l4.pkl'), 'wb'))
pickle.dump(markers_dict_Myeloid_l4, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_Myeloid_l4.pkl'), 'wb'))