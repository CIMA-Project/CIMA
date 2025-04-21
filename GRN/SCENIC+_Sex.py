import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

AUC_df = pd.read_parquet('scenicplus_custom_FinalAnnotation_otsu_AUC_with_Metadata.parquet')

cell_types = AUC_df['ACC_cell_type_l4'].unique()


transcription_factors = AUC_df.columns.to_list()[0:-8] 

filtered_factors = [factor for factor in transcription_factors if '_+_' in factor]


all_p_values = []


for cell_type in tqdm(cell_types):

    cell_type_data = AUC_df[AUC_df['ACC_cell_type_l4'] == cell_type]
    

    for factor in filtered_factors:

        male_data = cell_type_data[cell_type_data['ACC_sex'] == 'Male'][factor]
        female_data = cell_type_data[cell_type_data['ACC_sex'] == 'Female'][factor]
        

        _, p_value = mannwhitneyu(male_data, female_data, alternative='two-sided', nan_policy='omit')
        

        all_p_values.append((cell_type, factor, p_value))


p_values_df = pd.DataFrame(all_p_values, columns=['Cell_Type', 'Factor', 'P_Value'])


p_values_df['Corrected_P_Value'] = multipletests(p_values_df['P_Value'], method='fdr_bh')[1]


for index, row in p_values_df.iterrows():
    if row['Corrected_P_Value'] < 0.05:
        cell_type = row['Cell_Type']
        factor = row['Factor']
        cell_type_data = AUC_df[AUC_df['ACC_cell_type_l4'] == cell_type]
        
        plt.figure(figsize=(8, 8))
        sns.boxplot(x='ACC_sex', y=factor, data=cell_type_data)
        sns.stripplot(x='ACC_sex', y=factor, data=cell_type_data, color='black', alpha=0.5)
        

        plt.title(f'Sex Difference in Transcription Factor Activity for {cell_type} - {factor}\n(p-value={row["Corrected_P_Value"]:.3e})')
        plt.xlabel('Sex')
        plt.ylabel('Transcription Factor Activity')
        

        if not os.path.exists(f'Sex_Difference/{factor}'):
            os.makedirs(f'Sex_Difference/{factor}')
        filename = f'Sex_Difference/{factor}/{cell_type}_{factor}.pdf'
        plt.savefig(filename, format='pdf')
        

        plt.close()