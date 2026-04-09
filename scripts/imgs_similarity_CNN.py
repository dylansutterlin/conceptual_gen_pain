#%%
from transformers import AutoFeatureExtractor, ResNetModel
from PIL import Image
import torch
import os
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_resnet_embeddings(image_paths):
    extractor = AutoFeatureExtractor.from_pretrained("microsoft/resnet-50")
    model = ResNetModel.from_pretrained("microsoft/resnet-50")
    model.eval()

    embeddings = {}
    for img_path in image_paths:
        img = Image.open(img_path).convert("RGB")
        inputs = extractor(images=img, return_tensors="pt")
        with torch.no_grad():
            outputs = model(**inputs)
        # Use global pooled output
        embeddings[os.path.basename(img_path)] = outputs.pooler_output.squeeze().numpy()
    return embeddings

#%%
from transformers import CLIPProcessor, CLIPModel

def get_clip_embeddings(image_paths):
    processor = CLIPProcessor.from_pretrained("openai/clip-vit-base-patch32")
    model = CLIPModel.from_pretrained("openai/clip-vit-base-patch32")
    model.eval()

    embeddings = {}
    for img_path in image_paths:
        img = Image.open(img_path).convert("RGB")
        inputs = processor(images=img, return_tensors="pt")
        with torch.no_grad():
            outputs = model.get_image_features(**inputs)
        embeddings[os.path.basename(img_path)] = outputs.squeeze().numpy()
    return embeddings

#%%
from transformers import CLIPProcessor, CLIPModel

def get_openclip_embeddings(image_paths):
    processor = CLIPProcessor.from_pretrained("laion/CLIP-ViT-B-32-laion2B-s34B-b79K")
    model = CLIPModel.from_pretrained("laion/CLIP-ViT-B-32-laion2B-s34B-b79K")
    model.eval()


    embeddings = {}
    for img_path in image_paths:
        img = Image.open(img_path).convert("RGB")
        inputs = processor(images=img, return_tensors="pt")
        with torch.no_grad():
            outputs = model.get_image_features(**inputs)
        embeddings[os.path.basename(img_path)] = outputs.squeeze().numpy()
    return embeddings


#%%


def compute_similarity_matrix(embeddings):
    names = list(embeddings.keys())
    emb_matrix = np.stack([embeddings[name] for name in names])
    sim_matrix = cosine_similarity(emb_matrix)
    return pd.DataFrame(sim_matrix, index=names, columns=names)



def compute_cross_similarity(embeddings_A, embeddings_B):
    names_A = list(embeddings_A.keys())
    names_B = list(embeddings_B.keys())

    matrix_A = np.stack([embeddings_A[k] for k in names_A])
    matrix_B = np.stack([embeddings_B[k] for k in names_B])

    sim_matrix = cosine_similarity(matrix_A, matrix_B)

    return pd.DataFrame(sim_matrix, index=names_A, columns=names_B)


#%%
#Run

import glob

base_path = '/home/dsutterlin/projects/genPain/docs/all_img_stimuli'
image_paths_learning = glob.glob(os.path.join(base_path, "Learning_Stim", "*.png"))
img_path_generalization = glob.glob(os.path.join(base_path, "Generalization_Stim", "*.png"))

#learning 
resnet_embeddings_learning = get_resnet_embeddings(image_paths_learning)
clip_embeddings_learning = get_openclip_embeddings(image_paths_learning)

#generalization
resnet_embeddings_generalization = get_resnet_embeddings(img_path_generalization)
clip_embeddings_generalization = get_openclip_embeddings(img_path_generalization)

resnet_sim_df = compute_cross_similarity(resnet_embeddings_learning, resnet_embeddings_generalization)
clip_sim_df = compute_cross_similarity(clip_embeddings_learning, clip_embeddings_generalization)

#%%

gen_names_ordered = ['GCA1_dog.png', 'GCA2_horse.png', 'GCA3_cow.png', 'GCV1_car.png', 'GCV2_train.png', 'GCV3_truck.png', 
'GPA1_dog.png', 'GPA2_horse.png', 'GPA3_cow.png', 'GPV1_car.png', 'GPV2_train.png', 'GPV3_truck.png',
'GWA1_dog.png', 'GWA2_horse.png', 'GWA3_cow.png', 'GWV1_car.png', 'GWV2_train.png', 'GWV3_truck.png']

learn_names_ordered = [
'LA1_dog.png', 'LA2_horse.png', 'LA3_cow.png', 'LV1_car.png', 'LV2_train.png', 'LV3_truck.png']

resnet_sim_df = resnet_sim_df.loc[learn_names_ordered, gen_names_ordered]
clip_sim_df = clip_sim_df.loc[learn_names_ordered, gen_names_ordered]

#%%
cmap = plt.get_cmap('OrRd')
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
fig, ax = plt.subplots(figsize=(6, 8))  # adjust size as needed
fig.subplots_adjust(left=0.5, right=0.6)  # narrow figure to show only bar
cb = mpl.colorbar.ColorbarBase(
    ax,
    cmap=cmap,
    norm=norm,
    orientation='vertical'
)
cb.set_ticks([])          # remove tick marks
cb.set_ticklabels([])     # remove labels
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap


def plot_similarity_heatmap(sim_df_raw, title, figsize=(8, 24),cmap = 'viridis', save_path=None):
    
    reordered_cols = [0, 1, 2, 6, 7, 8, 12, 13, 14, 3, 4, 5, 9, 10, 11, 15, 16, 17]
    sim_df = sim_df_raw.iloc[:, reordered_cols].T # fliped!! 

    plt.figure(figsize=figsize)
    sns.heatmap(sim_df, cmap='OrRd', cbar=False, xticklabels=False, yticklabels=False)
    # plt.title(title, fontsize=14)
    # plt.xlabel("Generalization Images")
    # plt.ylabel("Learning Images")
    # plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=1000)
    plt.show()

#%% 
# PLOT SIMILARITY MODELS + compute model fit (from empirical to other models)

save_matrices_to = '/home/dsutterlin/projects/genPain/results/behavioral/similarity_CS_GS_matrices'
def save_similarity_matrix(matrix, filename):
    filepath = os.path.join(save_matrices_to, filename)
    matrix.to_csv(filepath, index=True)

#load empirical similarity matrix
empirical_similarity = pd.read_csv(os.path.join(save_matrices_to, "empiricalSimmat.csv"),
                                    header=None)
empirical_similarity.index = learn_names_ordered
empirical_similarity.columns = gen_names_ordered

# Boulder semantic similarity matrix (6 × 6)
A_boulder = np.array([
    [1,     0.483, 0.402, 0.31,  0.216, 0.247],
    [0.483, 1,     0.494, 0.249, 0.25,  0.226],
    [0.402, 0.494, 1,     0.176, 0.174, 0.232],
    [0.31,  0.249, 0.176, 1,     0.34,  0.674],
    [0.216, 0.25,  0.174, 0.34,  1,     0.38],
    [0.247, 0.226, 0.232, 0.674, 0.38,  1]
])
# Expand to 6 × 18
boulder_similarity = np.tile(A_boulder, (1, 3))
save_similarity_matrix(pd.DataFrame(boulder_similarity), "boulder_similarity.csv")

A_soft = np.array([
    [1, 1, 1, -1, -1, -1],
    [1, 1, 1, -1, -1, -1],
    [1, 1, 1, -1, -1, -1],
    [-1, -1, -1, 1, 1, 1],
    [-1, -1, -1, 1, 1, 1],
    [-1, -1, -1, 1, 1, 1]
])
categ_similarity = np.tile(A_soft, (1, 3))  # 6 × 18
save_similarity_matrix(pd.DataFrame(categ_similarity), "categ_similarity.csv")

A_strict = np.array([
    [ 1, -1, -1, -1, -1, -1],
    [-1,  1, -1, -1, -1, -1],
    [-1, -1,  1, -1, -1, -1],
    [-1, -1, -1,  1, -1, -1],
    [-1, -1, -1, -1,  1, -1],
    [-1, -1, -1, -1, -1,  1]
])
# Also expand if needed (not specified, but for symmetry):
examplar_similarity = np.tile(A_strict, (1, 3))  # 6 × 18
save_similarity_matrix(pd.DataFrame(A_strict), "examplar_similarity.csv")

save_sim_fig = '/home/dsutterlin/projects/genPain/results/py_figures/similarity'

plot_similarity_heatmap(pd.DataFrame(empirical_similarity),
                         "Empirical Similarity Matrix",
                           figsize= (10, 24),
                           save_path = os.path.join(save_sim_fig, 'sim_empirical.png'))

plot_similarity_heatmap(pd.DataFrame(boulder_similarity), "Boulder Semantic Similarity",
                        save_path = os.path.join(save_sim_fig, 'sim_boulder.png'))
plot_similarity_heatmap(pd.DataFrame(categ_similarity), "Categorical Similarity",
                        save_path = os.path.join(save_sim_fig, 'sim_categ.png'))
plot_similarity_heatmap(pd.DataFrame(examplar_similarity), "Exemplar Similarity",
                        save_path = os.path.join(save_sim_fig, 'sim_exemplar.png'))


plot_similarity_heatmap(resnet_sim_df, "Visual Similarity (ResNet-50)",
                        save_path = os.path.join(save_sim_fig, 'sim_resnet50.png'))
# plot_similarity_heatmap(clip_sim_df, "Conceptual Similarity (CLIP)")


#%%
from scipy.stats import spearmanr

def compute_similarity_fit(matrix_a, matrix_b):
    a = np.asarray(matrix_a).flatten()
    b = np.asarray(matrix_b).flatten()
    assert a.shape == b.shape, "Matrix dimensions must match"
    rho, pval = spearmanr(a, b)
    return rho, pval

rho_resnet, p_resnet = compute_similarity_fit(empirical_similarity, resnet_sim_df.values)
rho_openclip, p_openclip = compute_similarity_fit(empirical_similarity, clip_sim_df.values)
rho_examplar, p_examplar = compute_similarity_fit(empirical_similarity, examplar_similarity)
rho_categ, p_categ = compute_similarity_fit(empirical_similarity, categ_similarity)
rho_boulder, p_boulder = compute_similarity_fit(empirical_similarity, boulder_similarity)

# Print results
print(f"Spearman correlation (ResNet-50): {rho_resnet:.3f}, p-value: {p_resnet:.3e}")
print(f"Spearman correlation (OpenCLIP): {rho_openclip:.3f}, p-value: {p_openclip:.3e}")
print(f"Spearman correlation (Exemplar): {rho_examplar:.3f}, p-value: {p_examplar:.3e}")
print(f"Spearman correlation (Categorical): {rho_categ:.3f}, p-value: {p_categ:.3e}")
print(f"Spearman correlation (Boulder): {rho_boulder:.3f}, p-value: {p_boulder:.3e}")

# %%
