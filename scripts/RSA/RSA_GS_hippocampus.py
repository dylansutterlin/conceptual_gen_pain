#%%
import numpy as np
import pandas as pd
import RSA_utils as RSA_utils
import os
from nilearn.plotting import plot_stat_map, view_img
#%%
SAVE_DIR = "/home/dsutterlin/projects/genPain/results/imaging/RSA"
os.makedirs(SAVE_DIR, exist_ok=True)


#%%
# Code for RSA of Generalization stimuli in the hippocampus

gen_stim_order = ['GCA1', 'GCA2', 'GCA3', 'GCV1', 'GCV2', 'GCV3',
                  'GPA1', 'GPA2', 'GPA3', 'GPV1', 'GPV2', 'GPV3',
                  'GWA1', 'GWA2', 'GWA3', 'GWV1', 'GWV2', 'GWV3']

# =============================================================================
# STIMULUS METADATA
# =============================================================================
stimulus_mapping = {
    'GCA1': {'exemplar': 'dog',   'category': 'animal'},
    'GPA1': {'exemplar': 'dog',   'category': 'animal'},
    'GWA1': {'exemplar': 'dog',   'category': 'animal'},

    'GCA2': {'exemplar': 'horse', 'category': 'animal'},
    'GPA2': {'exemplar': 'horse', 'category': 'animal'},
    'GWA2': {'exemplar': 'horse', 'category': 'animal'},

    'GCA3': {'exemplar': 'cow',   'category': 'animal'},
    'GPA3': {'exemplar': 'cow',   'category': 'animal'},
    'GWA3': {'exemplar': 'cow',   'category': 'animal'},

    'GCV1': {'exemplar': 'car',    'category': 'vehicle'},
    'GPV1': {'exemplar': 'car',    'category': 'vehicle'},
    'GWV1': {'exemplar': 'car',    'category': 'vehicle'},

    'GCV2': {'exemplar': 'train',  'category': 'vehicle'},
    'GPV2': {'exemplar': 'train',  'category': 'vehicle'},
    'GWV2': {'exemplar': 'train',  'category': 'vehicle'},

    'GCV3': {'exemplar': 'truck',  'category': 'vehicle'},
    'GPV3': {'exemplar': 'truck',  'category': 'vehicle'},
    'GWV3': {'exemplar': 'truck',  'category': 'vehicle'},
}

gen_stim_full_name = {

    "GCA1": "Cartoon-Dog",
    "GCA2": "Cartoon-Horse",
    "GCA3": "Cartoon-Cow",

    "GCV1": "Cartoon-Car",
    "GCV2": "Cartoon-Train",
    "GCV3": "Cartoon-Truck",

    "GPA1": "Picture-Dog",
    "GPA2": "Picture-Horse",
    "GPA3": "Picture-Cow",

    "GPV1": "Picture-Car",
    "GPV2": "Picture-Train",
    "GPV3": "Picture-Truck",

    "GWA1": "Word-Dog",
    "GWA2": "Word-Horse",
    "GWA3": "Word-Cow",

    "GWV1": "Word-Car",
    "GWV2": "Word-Train",
    "GWV3": "Word-Truck",
}

# =============================================================================
# DATA DIRECTORIES
# =============================================================================
#------------------------------
# Extract Trial by trial activation maps in HPC
import os

# Path and file
project_dir = '/home/dsutterlin/projects/genPain/'

data_dir = os.path.join(project_dir, 'DATA/Imaging/test_genBetas/model3_generalization_ST/')

output_dir = os.path.join(project_dir, 'results/imaging/RSA_GS_hippocampus/')
os.makedirs(output_dir, exist_ok=True)

subjects = sorted([
    folder_name
    for folder_name in os.listdir(data_dir)
    if folder_name.startswith('SCEBL')
    and not folder_name.endswith('.zip')
])

# Hard-coded beta filenames for AllPains trials (4 runs × 18 trials = 72)
beta_names = []
# run1: 0002..0019
beta_names += [f"beta_{i:04d}.img" for i in range(2, 20)]
# run2: 0052..0069
beta_names += [f"beta_{i:04d}.img" for i in range(52, 70)]
# run3: 0103..0120
beta_names += [f"beta_{i:04d}.img" for i in range(103, 121)]
# run4: 0158..0175
beta_names += [f"beta_{i:04d}.img" for i in range(158, 176)]

assert len(beta_names) == 72

func_files = {sub : [os.path.join(data_dir, sub, beta_file) for beta_file in beta_names] for sub in subjects}

residual_files = {sub : os.path.join(data_dir, sub, "ResMS.img") for sub in subjects}

behavioral_cond = pd.read_csv('/home/dsutterlin/projects/genPain/results/behavioral/SCEBLmri_Gendata_TxT_N36.csv')
gs_cond_column = behavioral_cond[['sub_TxTgrp', 'cscat']]
# convert to dict, based on unique values in sub_TxTgrp. Output is sub: [list where gs_cond_column['cscat]] and sub_TxTgrp is sub_i
gs_cond_dict = {subjects[i]: gs_cond_column[gs_cond_column['sub_TxTgrp'] == sub_id]['cscat'].to_numpy() for i, sub_id in enumerate(gs_cond_column['sub_TxTgrp'].unique())}

# %%

func_labels = {}
for sub in func_files.keys():
    # assign gen stim as 4 times 18 stimuli
    func_labels[sub] = np.tile(np.array(gen_stim_order), 4)


# %%
# compute mean image per generalization stimulus per subject
import nibabel as nib
import numpy as np
from nilearn.image import mean_img

# def beta_to_t_imgs(beta_img, resms_img):

#     resms_data = resms_img.get_fdata()
#     beta_data = beta_img.get_fdata()

#     t_data = beta_data / np.sqrt(
#         resms_data + 1e-8
#     )

#     t_img = nib.Nifti1Image(
#         t_data,
#         beta_img.affine
#     )
#     return t_img


sub_mean_cond_imgs = {}

for sub in func_files.keys():
    sub_mean_cond_imgs[sub] = {}
    print(f"Processing subject {sub}...")
    img_paths = func_files[sub]
    labels = func_labels[sub]

    # Load all beta images for this subject
    imgs = [nib.load(img_path) for img_path in img_paths]
    print(f"Loaded {len(imgs)} images for subject {sub}.")
    
    # load residual
    residual_img = nib.load(residual_files[sub])

    # Compute mean image for each generalization stimulus
    for cond in gen_stim_order:
        cond_indices = np.where(labels == cond)[0]  # indices of trials with this condition
        cond_data = [imgs[i] for i in cond_indices]  # shape: (num_trials_cond, x, y, z)
        mean_cond_img = mean_img(cond_data)  # average across trials

       # t_mean_map = beta_to_t_imgs(mean_cond_img, residual_img)
        sub_mean_cond_imgs[sub][cond] = mean_cond_img
    

#%%
# =============================================================================
# BUILD STIMULUS DATAFRAME
# =============================================================================

stim_info = []

for stim in gen_stim_order:

    stim_info.append({
        "stim": stim,
        "exemplar": stimulus_mapping[stim]["exemplar"],
        "category": stimulus_mapping[stim]["category"]
    })

stim_df = pd.DataFrame(stim_info)

n = len(stim_df)

# =============================================================================
# EXEMPLAR MODEL
# =============================================================================

exemplar_model = np.zeros((n, n))

for i in range(n):

    for j in range(n):

        if stim_df.loc[i, "exemplar"] == stim_df.loc[j, "exemplar"]:

            exemplar_model[i, j] = 1

        else:

            exemplar_model[i, j] = 0

# =============================================================================
# CATEGORY MODEL
# =============================================================================

category_model = np.zeros((n, n))

for i in range(n):

    for j in range(n):

        if stim_df.loc[i, "category"] == stim_df.loc[j, "category"]:

            category_model[i, j] = 1

        else:

            category_model[i, j] = 0


# =============================================================================
# DATAFRAMES
# =============================================================================

exemplar_model_df = pd.DataFrame(
    exemplar_model,
    index=gen_stim_order,
    columns=gen_stim_order
)

category_model_df = pd.DataFrame(
    category_model,
    index=gen_stim_order,
    columns=gen_stim_order
)

#%%
# =============================================================================
# CS HIGH / LOW DEFINITIONS
# =============================================================================

cs_high_raw = [
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LA1dog_sm.bmp'
]

cs_low_raw = [
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV3truck_sm.bmp',
    'Stimuli/LA3cow_sm.bmp',
    'Stimuli/LV1car_sm.bmp',
    'Stimuli/LA2horse_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LA1dog_sm.bmp',
    'Stimuli/LV2train_sm.bmp',
    'Stimuli/LV1car_sm.bmp'
]

# =============================================================================
# CLEAN LABELS
# =============================================================================

def clean_learning_stim(stim):

    stim = stim.replace("Stimuli/", "")
    stim = stim.replace("_sm.bmp", "")

    if "dog" in stim:
        return "dog"

    elif "horse" in stim:
        return "horse"

    elif "cow" in stim:
        return "cow"

    elif "car" in stim:
        return "car"

    elif "train" in stim:
        return "train"

    elif "truck" in stim:
        return "truck"

# =============================================================================
# SUBJECT -> CSHIGH / CSLOW
# =============================================================================

subject_cs_dict = {}

for sub, csh, csl in zip(subjects, cs_high_raw, cs_low_raw):

    subject_cs_dict[sub] = {
        "CShigh": clean_learning_stim(csh),
        "CSlow": clean_learning_stim(csl)
    }

# =============================================================================
# PERSONALIZED GENERALIZATION MODEL 
# =============================================================================

custom_exemplar_model = {}

for sub in subjects:

    cshigh_exemplar = subject_cs_dict[sub]["CShigh"] #e.g. car
    cslow_exemplar = subject_cs_dict[sub]["CSlow"] #e.g. dog

    model = np.zeros((len(gen_stim_order), len(gen_stim_order)))

    for i, stim_i_code in enumerate(gen_stim_order):

        stim_i = stimulus_mapping[stim_i_code]['exemplar'] #e.g.GCA1 --> dog

        for j, stim_j_code in enumerate(gen_stim_order):
            
            stim_j = stimulus_mapping[stim_j_code]['exemplar'] 
           
            # stim_i_lower = stim_i.lower()
            # stim_j_lower = stim_j.lower()

            # --------------------------------------------------
            # Same CS HIGH exemplar
            # --------------------------------------------------

            same_cshigh = (
                cshigh_exemplar in stim_i
                and
                cshigh_exemplar in stim_j
            )

            # --------------------------------------------------
            # Same CS LOW exemplar
            # --------------------------------------------------

            same_cslow = (
                cslow_exemplar in stim_i
                and
                cslow_exemplar in stim_j
            )

            if same_cshigh or same_cslow:
                model[i, j] = 1

            else:
                model[i, j] = 0

    custom_exemplar_model[sub] = pd.DataFrame(
        model,
        index=gen_stim_order,
        columns=gen_stim_order
    )

#%%
# GRaded personalized models

animal_exemplars = ["dog", "horse", "cow"]
vehicle_exemplars = ["car", "train", "truck"]

def get_category(exemplar):

    if exemplar in animal_exemplars:
        return "animal"

    elif exemplar in vehicle_exemplars:
        return "vehicle"

personalized_graded_model = {}

for sub in subjects:

    cshigh_exemplar = subject_cs_dict[sub]["CShigh"]
    cslow_exemplar  = subject_cs_dict[sub]["CSlow"]

    cshigh_category = get_category(cshigh_exemplar)
    cslow_category  = get_category(cslow_exemplar)

    model = np.zeros((len(gen_stim_order), len(gen_stim_order)))

    # -------------------------------------------------------------------------
    # Loop over GS pairs
    # -------------------------------------------------------------------------

    for i, stim_i_code in enumerate(gen_stim_order):

        stim_i_exemplar = stimulus_mapping[stim_i_code]["exemplar"]
        stim_i_category = stimulus_mapping[stim_i_code]["category"]

        for j, stim_j_code in enumerate(gen_stim_order):

            stim_j_exemplar = stimulus_mapping[stim_j_code]["exemplar"]
            stim_j_category = stimulus_mapping[stim_j_code]["category"]

            similarity = 0

            # ==============================================================
            # 1) SAME LEARNED EXEMPLAR
            # ==============================================================

            same_cshigh = (
                stim_i_exemplar == cshigh_exemplar
                and
                stim_j_exemplar == cshigh_exemplar
            )

            same_cslow = (
                stim_i_exemplar == cslow_exemplar
                and
                stim_j_exemplar == cslow_exemplar
            )

            if same_cshigh or same_cslow:

                similarity = 1

            # ==============================================================
            # 2) SAME CATEGORY AS LEARNED EXEMPLAR
            # ==============================================================

            else:

                # animal generalization
                same_high_category = (
                    stim_i_category == cshigh_category
                    and
                    stim_j_category == cshigh_category
                )

                # vehicle generalization
                same_low_category = (
                    stim_i_category == cslow_category
                    and
                    stim_j_category == cslow_category
                )

                if same_high_category or same_low_category:

                    similarity = 0.5

            model[i, j] = similarity

    personalized_graded_model[sub] = pd.DataFrame(
        model,
        index=gen_stim_order,
        columns=gen_stim_order
    )

#%%
# high_category_model only where similarity is 1 for all the GS category for the high similus

animal_exemplars = ["dog", "horse", "cow"]
vehicle_exemplars = ["car", "train", "truck"]


high_category_model = {}

for sub in subjects:

    cshigh_exemplar = subject_cs_dict[sub]["CShigh"]
    cshigh_category = get_category(cshigh_exemplar)

    model = np.zeros((len(gen_stim_order), len(gen_stim_order)))

    # -------------------------------------------------------------------------
    # Loop over GS pairs
    # -------------------------------------------------------------------------

    for i, stim_i_code in enumerate(gen_stim_order):

        stim_i_category = stimulus_mapping[stim_i_code]["category"]

        for j, stim_j_code in enumerate(gen_stim_order):

            stim_j_category = stimulus_mapping[stim_j_code]["category"]

            similarity = 0


            same_cshigh = (
                stim_i_category == cshigh_category
                and
                stim_j_category == cshigh_category
            )

            if same_cshigh:

                similarity = 1

            model[i, j] = similarity

    high_category_model[sub] = pd.DataFrame(
        model,
        index=gen_stim_order,
        columns=gen_stim_order
    )
# %%

gen_names_plotting = stimulus_mapping.keys()  # same order as in the models

RSA_utils.plot_similarity_matrix(
    exemplar_model_df,
    title="Exemplar Model",
    show_labels=True,
    rename_labels=gen_stim_full_name,

)

RSA_utils.plot_similarity_matrix(
    category_model_df,
    title="Category Model",
    show_labels=True,
    rename_labels=gen_stim_full_name,
)

RSA_utils.plot_similarity_matrix(
    category_model_df,
    title="Category Model",
    show_labels=False,
    reorder_columns=gen_names_plotting,
    rename_labels=gen_stim_full_name
)

RSA_utils.plot_similarity_matrix(
    exemplar_model_df,
    title="Exemplar Model",
    show_labels=False,
    rename_labels=gen_stim_full_name,
    reorder_columns=gen_names_plotting

)
#%%
RSA_utils.plot_similarity_matrix(
    custom_exemplar_model[subjects[0]],
    title=f"Custom Exemplar Model - {subjects[0]}",
    show_labels=False,
    reorder_columns=gen_names_plotting,
    rename_labels=gen_stim_full_name,
)

RSA_utils.plot_similarity_matrix(high_category_model[subjects[0]],
                       title=f"High Category Model - {subjects[0]}",
                          show_labels=False,
                          reorder_columns=gen_names_plotting,
                            rename_labels=gen_stim_full_name
)

# test matrix rotation for pattern separation
# suppose maximal pattern separation for most distinct stimuli.
pattern_separation_model = pd.DataFrame(np.rot90(np.array(category_model_df), k=1), index=gen_stim_order, columns=gen_stim_order)

RSA_utils.plot_similarity_matrix(pattern_separation_model,
                       title=f"Pattern Separation Model",
                          show_labels=False,
                          reorder_columns=gen_names_plotting,
                            rename_labels=gen_stim_full_name
                          )
# %%
# Example plot

example_sub = subjects[0]

RSA_utils.plot_similarity_matrix(
    personalized_graded_model[subjects[0]],
    title="Personalized graded",
    reorder_columns=gen_names_plotting,
    rename_labels=gen_stim_full_name,
    show_labels=False
)

# %%
#==============================================================================
# HIPPOCAMPUS MASK and ROI
#==============================================================================
# Load hippocampus mask and extract mean activation for each condition and subject

#extract hippocampus activity patterns from path a
from nilearn.maskers import NiftiMapsMasker, NiftiMasker
from nilearn import plotting, image

template_img = sub_mean_cond_imgs[subjects[0]]['GCA1']  # example image to get shape and affine

import importlib
import RSA_utils

importlib.reload(RSA_utils)
hip_mask = RSA_utils.extract_hpc_mask("sphere", template_img)
view_img(hip_mask, title="Hippocampus Mask", threshold=0.5)

#%%
from scipy.stats import zscore

# Extract activation pattern in the HPC per condition and subject
masker = NiftiMasker(mask_img=hip_mask, resampling_target='data')
hpc_patterns = {}

masker.fit(template_img)

for sub in sub_mean_cond_imgs.keys():
    hpc_patterns[sub] = {}
    for cond in gen_stim_order:
        mean_img_cond = sub_mean_cond_imgs[sub][cond]
        hpc_pattern = masker.transform(mean_img_cond)[0]  # shape: (n_regions,)
        #print("Applying z-scoring to pattern")
        #hpc_pattern = zscore(hpc_pattern)
        hpc_patterns[sub][cond] = hpc_pattern

# %%
#---------------------
# compute RSA matrix per subject and compare to models
from scipy.spatial.distance import euclidean
import numpy as np
import pandas as pd
from importlib import reload
reload(RSA_utils)

subject_hpc_similarity = {}

for sub in hpc_patterns.keys():

    print(f"Computing HPC similarity matrix: {sub}")

    subject_hpc_similarity[sub] = RSA_utils.compute_neural_similarity(
        hpc_patterns[sub],
        gen_stim_order, # model order
        metric="pearson"
    )
# Mean similarity matrix across subjects
mean_hpc_similarity = np.mean(
    [subject_hpc_similarity[sub].values for sub in subject_hpc_similarity.keys()],
    axis=0
)
RSA_utils.plot_similarity_matrix(mean_hpc_similarity, title=f"HPC Similarity Matrix - {subjects[0]}", show_labels=True)


#%%
from importlib import reload
reload(RSA_utils)
# RUN RSA and permutation test for each subject and model
rsa_results = []

comparator_method = 'spearmanr'  # or 'pearson', 'spearman', etc.

for sub in subject_hpc_similarity.keys():

    neural_mat = subject_hpc_similarity[sub]

    # EXEMPLAR
    rho_ex, p_ex, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        exemplar_model_df,
        n_permutations=10000,
        comparator_method = comparator_method
    )

    # CATEGORY
    rho_cat, p_cat, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        category_model_df,
        n_permutations=10000,
        comparator_method = comparator_method
    )

    # Personalized exemplar models
    custom_exemplar_model_sub = custom_exemplar_model[sub]
    rho_personalized, p_personalized, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        custom_exemplar_model_sub,
        n_permutations=1000,
        comparator_method = comparator_method
    )
    
    personalized_graded_model_sub = personalized_graded_model[sub]
    rho_personalized_graded, p_personalized_graded, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        personalized_graded_model_sub,
        n_permutations=1000,
        comparator_method = comparator_method
    )

    high_category_model_sub = high_category_model[sub]
    rho_high_category, p_high_category, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        high_category_model_sub,
        n_permutations=1000,
        comparator_method = comparator_method
    )

    # Pattern separation model
    pattern_separation_model = pd.DataFrame(np.rot90(np.array(category_model_df), k=1), index=gen_stim_order, columns=gen_stim_order)
    rho_pattern_separation, p_pattern_separation, _ = RSA_utils.rsa_permutation_test(
        neural_mat,
        pattern_separation_model,
        n_permutations=1000,
        comparator_method = comparator_method
    )

    rsa_results.append({
        "subject": sub,
        "rho_exemplar": rho_ex,
        "p_exemplar": p_ex,
        "rho_category": rho_cat,
        "p_category": p_cat,
        "rho_personalized": rho_personalized,
        "p_personalized": p_personalized,
        "rho_personalized_graded": rho_personalized_graded,
        "p_personalized_graded": p_personalized_graded,
        "rho_high_category": rho_high_category,
        "p_high_category": p_high_category,
        "rho_pattern_separation": rho_pattern_separation,
        "p_pattern_separation": p_pattern_separation
    })

rsa_results_df = pd.DataFrame(rsa_results)

#%%
# plot example subject results
# Mean similarity matrix across subjects
reload(RSA_utils)


RSA_utils.plot_similarity_matrix(subject_hpc_similarity['SCEBL220_1527'], title=f"HPC Similarity Matrix - SCEBL220_1527 (rho = {rsa_results_df[rsa_results_df['subject'] == 'SCEBL220_1527']['rho_category'].values[0]:.2f})", show_labels=True,
                       reorder_columns=gen_names_plotting, rename_labels=gen_stim_full_name)
RSA_utils.plot_similarity_matrix(subject_hpc_similarity['SCEBL217_1505'], title=f"HPC Similarity Matrix - SCEBL217_1505 (rho = {rsa_results_df[rsa_results_df['subject'] == 'SCEBL217_1505']['rho_category'].values[0]:.2f})", show_labels=True,
                                                           reorder_columns=gen_names_plotting,
                            rename_labels=gen_stim_full_name )
RSA_utils.plot_similarity_matrix(subject_hpc_similarity['SCEBL209_1495'], title=f"HPC Similarity Matrix - SCEBL209_1495 (rho = {rsa_results_df[rsa_results_df['subject'] == 'SCEBL209_1495']['rho_category'].values[0]:.2f})", show_labels=True,
                                                           reorder_columns=gen_names_plotting,
                            rename_labels=gen_stim_full_name,
                            rm_diagonal=True,
                            add_contour = True)
RSA_utils.plot_similarity_matrix(category_model_df, title=f"Category Model", show_labels=True,
                       reorder_columns=gen_names_plotting, rename_labels=gen_stim_full_name) 

# %%

# behavioral expectations
#load behavioral data
csv_path = os.path.join(project_dir, 'results/behavioral/gen_exp_learning_betas.csv')
# Load the data
data = pd.read_csv(csv_path)
gen_betas = data['gen_betas']
learn_betas = data['exp_learning_beta']

# %%
#----------------------------
# PERMUTATION TEST FOR BEHAVIOR-RSA CORRELATION
from scipy.stats import pearsonr, spearmanr, kendalltau
import numpy as np
import matplotlib.pyplot as plt
import RSA_utils
reload(RSA_utils)

def behavioral_permutation_test(
    x,
    y,
    n_permutations=10000,
    random_state=22,
    comparator_method="pearson"
):
    """
    Permutation test for behavioral-RSA correlations.
    """

    rng = np.random.default_rng(random_state)

    x = np.asarray(x)
    y = np.asarray(y)

    if comparator_method == "pearson":
        observed_r, _ = pearsonr(x, y)
    elif comparator_method == "kendalltau":
        observed_r, _ = kendalltau(x, y)
    elif comparator_method == "spearman":
        observed_r, _ = spearmanr(x, y)
    elif comparator_method == "regress":
        slope, intercept, observed_r, p_value, std_err = stats.linregress(x, y)
    else:
        raise ValueError("Invalid comparator method. Use 'pearson' or 'kendalltau'.")

    permuted_r = np.zeros(n_permutations)

    for i in range(n_permutations):

        shuffled_y = rng.permutation(y)

        if comparator_method == "pearson":
            permuted_r[i], _ = pearsonr(x, shuffled_y)
        elif comparator_method == "kendalltau":
            permuted_r[i], _ = kendalltau(
                x,
                shuffled_y
            )
        elif comparator_method == "spearman":
            permuted_r[i], _ = spearmanr(
                x,
                shuffled_y
            )
        elif comparator_method == "regress":
            slope, intercept, permuted_r[i], p_value, std_err = stats.linregress(
                x,
                shuffled_y
            )

    p_value = (
        np.sum(
            np.abs(permuted_r)
            >=
            np.abs(observed_r)
        ) + 1
    ) / (n_permutations + 1)

    return observed_r, p_value, permuted_r

def plot_scatter(x, y, xlabel, ylabel, title):
    
    plt.scatter(x, y)
    
    # Add regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.array([x.min(), x.max()])
    line_y = slope * line_x + intercept
    plt.plot(line_x, line_y, 'r-', label=f'R² = {r_value**2:.3f}')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_null_dist(null_dist, observed_r, title, save_as=None):
    plt.figure(figsize=(9, 6))
    plt.hist(null_dist, bins=30, alpha=0.7, color = "#a37171")

    plt.axvline(observed_r, color='red', linestyle='dashed', linewidth=4)
    plt.xlabel('Correlation (rho)', fontsize=30, labelpad=5)
    plt.ylabel('Frequency', fontsize=30, labelpad=5)
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.tick_params(axis='both', which='minor', labelsize=30)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(30)
    plt.title(title, fontsize=30, pad=20)
    plt.legend()
    plt.tight_layout()
    
    if save_as is not None:
        plt.savefig(save_as, dpi=1000)
    plt.show()



model_names = [
    "rho_exemplar",
    "rho_category",
    "rho_personalized",
    "rho_personalized_graded",
    "rho_high_category",
    "rho_pattern_separation"
]

print("\nRSA ~ learning betas")
for model in model_names:
   
   r, p_perm, null_dist = behavioral_permutation_test(
       learn_betas,
       rsa_results_df[model],
       comparator_method="spearman"
   )

   print(
       f"{model}: "
       f"r={r:.3f}, "
       f"perm_p={p_perm:.4f}"
   )
#    if model in ['rho_high_category']:
#        RSA_utils.jointplot(learn_betas, rsa_results_df[model],
#                            x_label=r'Pain learning ($\beta$)',
#                              y_label=f'RSA Rho', title=f'{model} RSA vs Learning Beta',
#                              r_p_input=(r, p_perm),
#                              save_as =os.path.join(SAVE_DIR, f"scatter_{model}_rsa_vs_learning_beta.png"))
#        plot_null_dist(null_dist, r, title=f'Null Distribution of {model}',
#                       save_as=os.path.join(SAVE_DIR, f"null_dist_{model}_rsa_vs_learning_beta.png"))

# print("RSA ~ learning betas")


for model in model_names:

    r, p_perm, null_dist = behavioral_permutation_test(
        gen_betas,
        rsa_results_df[model],
        comparator_method="spearman"
    )

    print(
        f"{model}: "
        f"r={r:.8f}, "
        f"perm_p={p_perm:.4f}"
    )

    if model in ['rho_exemplar', 'rho_category']:


        RSA_utils.jointplot(gen_betas, rsa_results_df[model],
                            x_label=r'Pain generalization ($\beta$)',
                              y_label=f'RSA Rho', title=f'{model} RSA vs Generalization Beta',
                              r_p_input=(r, p_perm),
                              save_as =os.path.join(SAVE_DIR, f"scatter_{model}_rsa_vs_gen_beta.png"))
        plot_null_dist(null_dist, r, title=f'Null Distribution of {model}',
                       save_as=os.path.join(SAVE_DIR, f"null_dist_{model}_rsa_vs_gen_beta.png"))

# %%
# one-sample t-test to see if RSA correlations are significantly greater than 0 across subjects

from scipy.stats import ttest_1samp, ttest_rel
import numpy as np
from scipy import stats

models = [
    "rho_exemplar",
    "rho_category",
    "rho_personalized",
    "rho_personalized_graded",
    "rho_high_category",
    "rho_pattern_separation"
]

for model in models:

    t, p = ttest_1samp(
        rsa_results_df[model],
        popmean=0
    )

    print(
        f"{model}: "
        f"mean={rsa_results_df[model].mean():.3f}, "
        f"t={t:.3f}, "
        f"p={p:.4f}"
    )

# compare exemplar vs category models
t_exemplar_category, p_exemplar_category = ttest_rel(
    rsa_results_df['rho_exemplar'],
    rsa_results_df['rho_category']
)
print(
    f"Exemplar vs Category: "
    f"t={t_exemplar_category:.3f}, "
    f"p={p_exemplar_category:.4f}"
        )

# comapre category vs high pain
t_category_high, p_category_high = ttest_rel(
    rsa_results_df['rho_category'],
    rsa_results_df['rho_high_category']
)
print(
    f"Category vs High Category: "
    f"t={t_category_high:.3f}, "
    f"p={p_category_high:.4f}"
        )

reload(RSA_utils)

RSA_utils.plot_rsa_model_fits(
    rsa_results_df,
    models=models,
    title="RSA Model Fits",

)
# %%
rsa_results_df.describe()
# %%
# test high leverage point influence
current_df = rsa_results_df["rho_category"]
current_model = "rho_category_outlier_test"
#filter highleverage points
std = current_df.std()
print(f"Standard deviation: {std:.3f}")

mask = np.abs(current_df - current_df.mean()) < 3 * std
filtered_df = current_df[mask]
filtered_gen_betas = gen_betas[mask]

# plot
import matplotlib.pyplot as plt
plt.scatter(filtered_gen_betas, filtered_df)
plt.xlabel('Behavioral Generalization Beta')
plt.ylabel('RSA High Category Rho')
plt.title(f'High Category RSA vs Generalization Beta\nr={corr_personalized_graded:.3f}, p={p_personalized_graded:.4f}')
plt.tight_layout()
plt.show()

r, p_perm, null_dist = behavioral_permutation_test(
        filtered_gen_betas,
        filtered_df
    )

print(
    f"{current_model}: "
    f"r={r:.3f}, "
    f"perm_p={p_perm:.4f}"
)
# %%
# PUBLISHED OUTPUT 
from nilearn.datasets import load_mni152_template
high_res_template = load_mni152_template(resolution=1)

plotting.plot_stat_map(hip_mask, bg_img=high_res_template,black_bg=False, title="Hippocampus Mask", threshold=0.5, display_mode='ortho', cut_coords=(-32, -22, -12),cmap='coolwarm', colorbar=False)


