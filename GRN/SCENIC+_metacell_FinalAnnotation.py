import anndata as ad
import os
import pickle
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

work_dir = 'CIMA/scenicplus/'
tmp_dir = 'CIMA/Temp/scenicplus/'

if not os.path.exists(work_dir):
    os.makedirs(work_dir)
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

atac = ad.read_h5ad('CIMA/scATAC/CIMA_scATAC_PeakRecalling_Metacell.h5ad')
atac.obs.index = atac.obs['new_barcode']

from pycisTopic.cistopic_class import *
cistopic_obj = create_cistopic_object(fragment_matrix=atac.X.T.tocsr(),
                                     cell_names = atac.obs_names.to_list(),
                                      region_names= atac.var_names.to_list(),
                                      tag_cells = False,
                                      split_pattern = ' ')
cistopic_obj.add_cell_data(atac.obs)

if not os.path.exists(os.path.join(work_dir,'pycisTopic')):
    os.makedirs(os.path.join(work_dir,'pycisTopic'))
pickle.dump(cistopic_obj,open(work_dir+'pycisTopic/CIMA_scATAC_PeakRecalling_Metacell_cisTopic_obj.pkl', 'wb'))

# run_cgs_models_mallet
cistopic_obj = pickle.load(open('CIMA/scenicplus/CIMA_scATAC_PeakRecalling_Metacell_cisTopic_obj.pkl','rb'))

from pycisTopic.cistopic_class import *
path_to_mallet_binary="/home/zhengyuhui/mambaforge/envs/scenicplus/bin/mallet"
# Run models
models=run_cgs_models_mallet(path_to_mallet_binary,
                    cistopic_obj,
                    n_topics=[5, 20, 40, 60, 80, 100, 130, 150, 180],
                    n_cpu=96,
                    n_iter=150,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path=tmp_dir, #Use SCRATCH if many models or big data set
                    save_path=work_dir + 'pycisTopic/TopicModel_iter150/',
                    reuse_corpus=True)

topic_model_dir = os.listdir(work_dir + 'pycisTopic/TopicModel_iter150/')
models = []
for file in topic_model_dir:
    topic = pickle.load(open(work_dir + 'pycisTopic/TopicModel_iter150/' + file, 'rb'))
    models.append(topic)
    print(file)

from pycisTopic.lda_models import *
model = evaluate_models(models,
                       select_model=130,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)

cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,open(work_dir+'pycisTopic/CIMA_scATAC_PeakRecalling_Metacell_cisTopic_obj.pkl', 'wb'))

