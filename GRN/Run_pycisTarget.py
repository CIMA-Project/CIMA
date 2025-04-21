import pickle
import os
import ray

ray._private.ray_constants.DEFAULT_OBJECT_STORE_MEMORY_PROPORTION = 0.80
ray._private.ray_constants.DEFAULT_OBJECT_STORE_MAX_MEMORY_BYTES = 1500 * 10^9

work_dir = 'CIMA/scenicplus/'

def dict_slice(adict, start, end):
    keys = adict.keys()
    dict_slice = {}
    for k in list(keys)[start:end]:
        dict_slice[k] = adict[k]
    return dict_slice

region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict_cell_type_l3 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_cell_type_l3.pkl'), 'rb'))
markers_dict_YF_v_OF_l3 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_OF_l3.pkl'), 'rb'))
markers_dict_YM_v_OM_l3 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YM_v_OM_l3.pkl'), 'rb'))
markers_dict_OF_v_OM_l3 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_OF_v_OM_l3.pkl'), 'rb'))
markers_dict_YF_v_YM_l3 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_YM_l3.pkl'), 'rb'))
markers_dict_cell_type_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_cell_type_l4.pkl'), 'rb'))
markers_dict_YF_v_OF_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_OF_l4.pkl'), 'rb'))
markers_dict_YM_v_OM_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YM_v_OM_l4.pkl'), 'rb'))
markers_dict_OF_v_OM_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_OF_v_OM_l4.pkl'), 'rb'))
markers_dict_YF_v_YM_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_YF_v_YM_l4.pkl'), 'rb'))
markers_dict_B_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_B_l4.pkl'), 'rb'))
markers_dict_CD4_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_CD4_l4.pkl'), 'rb'))
markers_dict_CD8_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_CD8_l4.pkl'), 'rb'))
markers_dict_NK_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_NK_l4.pkl'), 'rb'))
markers_dict_Myeloid_l4 = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_Myeloid_l4.pkl'), 'rb'))

region_bin_topics_otsu_part1 = dict_slice(region_bin_topics_otsu,0,45)
region_bin_topics_otsu_part2 = dict_slice(region_bin_topics_otsu,45,60)
region_bin_topics_otsu_part3 = dict_slice(region_bin_topics_otsu,60,80)
region_bin_topics_otsu_part4 = dict_slice(region_bin_topics_otsu,80,100)
region_bin_topics_otsu_part5 = dict_slice(region_bin_topics_otsu,100,130)

pickle.dump(region_bin_topics_otsu_part1, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu_part1.pkl'), 'wb'))
pickle.dump(region_bin_topics_otsu_part2, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu_part2.pkl'), 'wb'))
pickle.dump(region_bin_topics_otsu_part3, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu_part3.pkl'), 'wb'))
pickle.dump(region_bin_topics_otsu_part4, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu_part4.pkl'), 'wb'))
pickle.dump(region_bin_topics_otsu_part5, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu_part5.pkl'), 'wb'))


import pyranges as pr
from pycistarget.utils import region_names_to_coordinates

datasets = {
    'topics_otsu_part1': region_bin_topics_otsu_part1,
    'topics_otsu_part2': region_bin_topics_otsu_part2,
    'topics_otsu_part3': region_bin_topics_otsu_part3,
    'topics_otsu_part4': region_bin_topics_otsu_part4,
    'topics_otsu_part5': region_bin_topics_otsu_part5,
    'DARs_cell_type_l3': markers_dict_cell_type_l3,
    'DARs_YF_v_OF_l3': markers_dict_YF_v_OF_l3,
    'DARs_YM_v_OM_l3': markers_dict_YM_v_OM_l3,
    'DARs_OF_v_OM_l3': markers_dict_OF_v_OM_l3,
    'DARs_YF_v_YM_l3': markers_dict_YF_v_YM_l3,
    'DARs_cell_type_l4': markers_dict_cell_type_l4,
    'DARs_YF_v_OF_l4': markers_dict_YF_v_OF_l4,
    'DARs_YM_v_OM_l4': markers_dict_YM_v_OM_l4,
    'DARs_OF_v_OM_l4': markers_dict_OF_v_OM_l4,
    'DARs_YF_v_YM_l4': markers_dict_YF_v_YM_l4,
    'DARs_B_l4': markers_dict_B_l4,
    'DARs_CD4_l4': markers_dict_CD4_l4,
    'DARs_CD8_l4': markers_dict_CD8_l4,
    'DARs_NK_l4': markers_dict_NK_l4,
    'DARs_Myeloid_l4': markers_dict_Myeloid_l4,
}


region_sets = {key: {} for key in datasets.keys()}


for output_key, data_dict in datasets.items():
    for key, df in data_dict.items():
        if df.shape[0] > 10:  
            regions = df.index[df.index.str.startswith('chr')] 
            region_sets[output_key][key] = pr.PyRanges(region_names_to_coordinates(regions))


for key, regions in region_sets.items():
    print(f'{key}: {list(regions.keys())}')

db_fpath = '/CIMA/scenicplus/'
rankings_db = os.path.join(db_fpath, 'C4_PBMC_1kb_bg_with_mask.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'C4_PBMC_1kb_bg_with_mask.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(db_fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs_custom'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 96,save_partial=True,
    _temp_dir = os.path.join('/CIMA/Temp/'),
    annotation_version = 'v10nr_clust',
    **{
         '_system_config':{
        # Allow spilling until the local disk is 99% utilized.
        # This only affects spilling to the local file system.
        "local_fs_capacity_threshold": 0.99},"dashboard_host":"0.0.0.0"
    }
)