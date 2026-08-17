from nilearn import image, plotting
from atlasreader.atlasreader import process_img
import nibabel as nib
import os
import pandas as pd
import tqdm


def extract_and_save_clusters(
    stat_img,
    output_dir,
    voxel_thresh=3.1,
    cluster_extent=20,
    atlas_mask=None,
    save_nifti=True,
    save_png=True
):
    """
    Extract connected clusters from statistical map.

    Parameters
    ----------
    stat_img : nib image
        Statistical image.

    output_dir : str

    voxel_thresh : float

    cluster_extent : int
        Minimum cluster size (voxels).

    direction : str
        both / positive / negative

    atlas_mask : image or None
        Optional anatomical mask.

    Returns
    -------
    cluster_dict : dict
        cluster_id -> cluster image
    """

    os.makedirs(output_dir, exist_ok=True)

    # --------------------------------------------------
    # constrain to anatomical mask if desired
    # --------------------------------------------------

    if atlas_mask is not None:

        stat_img = image.math_img(
            "img1 * (img2 > 0)",
            img1=stat_img,
            img2=atlas_mask
        )

    # --------------------------------------------------
    # atlasreader cluster extraction
    # --------------------------------------------------

    cluster_4d = process_img(
        stat_img,
        voxel_thresh=voxel_thresh,
        cluster_extent=cluster_extent,
    )

    cluster_dict = {}

    # --------------------------------------------------
    # save clusters
    # --------------------------------------------------

    for i, cluster_img in enumerate(
        tqdm.tqdm(
            image.iter_img(cluster_4d),
            desc="Saving clusters"
        )
    ):

        cluster_id = i + 1

        cluster_dict[cluster_id] = cluster_img

        if save_nifti:

            nii_path = os.path.join(
                output_dir,
                f"cluster_{cluster_id}.nii.gz"
            )

            nib.save(
                cluster_img,
                nii_path
            )

        if save_png:

            fig_path = os.path.join(
                output_dir,
                f"cluster_{cluster_id}.png"
            )

            display = plotting.plot_stat_map(
                cluster_img,
                threshold=0,
                display_mode="ortho",
                colorbar=True,
                title=f"Cluster {cluster_id}"
            )

            display.savefig(
                fig_path,
                dpi=300
            )

            display.close()

    return cluster_dict


#EXTRACT ROI MASK
from nilearn import datasets, image, plotting
from nilearn.maskers import nifti_spheres_masker
from nilearn.masking import _unmask_3d
from nibabel import Nifti1Image
import nibabel as nib
import numpy as np


def extract_hpc_mask(type, template_img):


    import nilearn
    import nibabel as nib
    from nilearn import datasets, image, plotting

    pruned_ab_map = nib.load("/home/dsutterlin/projects/genPain/results/py_figures/mediation_derivatives/ab_pruned_FDR_01unc.nii.gz")
    nilearn.plotting.plot_stat_map(pruned_ab_map, title='pruned ab map', display_mode='z', cut_coords=1)

    save_derivatives_rsa = "/home/dsutterlin/projects/genPain/results/imaging/RSA_GS_hippocampus"

    #cluster 1 = HPC
    cluster_1_path = os.path.join(save_derivatives_rsa, "cluster_1.nii.gz")

    if os.path.exists(cluster_1_path):
        cluster_1_mask = nib.load(cluster_1_path)
    else:
        # run cluster detection
        cluster_masks = extract_and_save_clusters(
            stat_img=pruned_ab_map,
            atlas_mask=None,
            output_dir=save_derivatives_rsa,
            voxel_thresh=0,
            cluster_extent=10,
        )
        cluster_1_mask = cluster_masks.get(1)  # get cluster 1 mask


    binary_cluster_1_mask = image.binarize_img(cluster_1_mask   
    )
    mediation_hcp_mask = image.resample_to_img(binary_cluster_1_mask, template_img, interpolation='nearest')
    plotting.plot_roi(mediation_hcp_mask, title='HPC mask from cluster 1')


        # atlas_img = nib.load(os.path.join(project_dir, 'atlas/Tian2020_schaeffer200_subcortical16', 'combined_schaefer200_tian16.nii.gz'))
        # atlas_df_path = os.path.join(project_dir, 'atlas/Tian2020_schaeffer200_subcortical16', 'Schaefer2018_200Parcels_17Networks_order_Tian_Subcortex_S1_label.txt')
        # with open(atlas_df_path, 'r') as f:
        #     atlas_labels = [line.strip() for line in f][::2] #!! getting only labels, only 16 !!
        # atlas_df = pd.DataFrame({'label': atlas_labels, 'index': range(len(atlas_labels))})

        # atlas_masker = NiftiMapsMasker(maps_img=atlas_img, resampling_target='maps')
        # plotting.view_img(atlas_img, title='Schaefer 200 atlas', cmap='Paired', colorbar=True, symmetric_cmap=False)
        # print('shape', atlas_img.shape)

        # # get region 0 and 8 where hippocampus
        # #assign 0 
        # atlas_data = atlas_img.get_fdata()

        # single_region_masks = {}
        # for i in [200, 208]:  # hippocampus regions in the atlas
        #     atlas_np = np.zeros_like(atlas_data, dtype=float)
        #     atlas_np[atlas_data == i+1] = 1.0
        #     roi_mask = nib.Nifti1Image(atlas_np, atlas_img.affine, atlas_img.header)  
        #     single_region_masks[i] = image.resample_to_img(roi_mask, template_img, interpolation='nearest')

        # plotting.view_img(single_region_masks[208], title='Hippocampus Atlas', cmap='gray', colorbar=True)

    # GET ANATOMICAL HPC from atlas
    from nilearn import datasets, image, plotting
    atlas = datasets.fetch_atlas_juelich(
        "maxprob-thr25-1mm"
    )

    atlas_img = atlas.maps
    labels = atlas.labels

    for i, label in enumerate(labels):

        print(i, label)

    atlas_data = atlas_img.get_fdata()

    atlas_hpc_mask = np.isin(
        atlas_data,
        [9, 10, 11, 12] # manually ID with labels
    )

    # Keep only left hemisphere (x < 0)
    x_coords = np.arange(atlas_hpc_mask.shape[0])
    atlas_hpc_mask[x_coords >= atlas_hpc_mask.shape[0] // 2] = False
    
    atlas_hpc_mask = nib.Nifti1Image(
        atlas_hpc_mask.astype(float),
        atlas_img.affine
    )

    atlas_hpc_mask = image.resample_to_img(
        atlas_hpc_mask,
        template_img,
        interpolation="nearest"
    )


    # SPHERE MASK from coordinate

    # ---------------------------------------------------
    # Resample left hippocampus mask to brain mask space
    # ---------------------------------------------------

    brain_mask = datasets.load_mni152_brain_mask()

    resampled_hpc_mask = image.resample_to_img(
        atlas_hpc_mask,
        brain_mask,
        interpolation="nearest"
    )

    _, A = nifti_spheres_masker._apply_mask_and_get_affinity(
        seeds=[(-32, -22, -12)],
        niimg=None,
        radius=10,
        allow_overlap=False,
        mask_img=resampled_hpc_mask
    )

    # IMPORTANT:
    # use SAME mask used to generate A
    sphere_mask = _unmask_3d(
        X=A.toarray().flatten(),
        mask=resampled_hpc_mask.get_fdata().astype(bool)
    )

    sphere_mask_img = Nifti1Image(
        sphere_mask.astype(np.int8),
        resampled_hpc_mask.affine
    )

    sphere_mask_img = image.resample_to_img(
        sphere_mask_img,
        template_img,
        interpolation="nearest"
    )

    nib.save(
        sphere_mask_img,
        "left_hpc_constrained_sphere.nii.gz"
    )

    plotting.plot_roi(
        sphere_mask_img,
        title="Constrained Left HPC Sphere"
    )

    plotting.show()

    # CHOOSE HPC MASK


    if type == "atlas":
        hip_mask = atlas_hpc_mask
    elif type == "sphere":
        hip_mask = sphere_mask_img
    elif type == "cluster":
        hip_mask = mediation_hcp_mask
    elif type == "intersection":
        hip_mask = image.math_img(
            'img1 * img2',
            img1=mediation_hcp_mask,
            img2=atlas_hpc_mask
        )

    plotting.view_img(hip_mask, title=f"HPC mask - {type}")

    return hip_mask

#=======================================================

from scipy.stats import spearmanr
from scipy.stats import pearsonr
import numpy as np
import pandas as pd


def compute_neural_similarity(
    pattern_dict,
    condition_order,
    metric="pearson"
):
    """
    Compute condition × condition similarity matrix
    from voxel activation patterns.

    Parameters
    ----------
    pattern_dict : dict
        cond -> voxel pattern

    condition_order : list
        ordering of conditions

    metric : str
        "pearson" or "euclidean_similarity"
    """

    n_cond = len(condition_order)

    sim_matrix = np.zeros((n_cond, n_cond))

    for i, cond_i in enumerate(condition_order):

        pattern_i = pattern_dict[cond_i]

        for j, cond_j in enumerate(condition_order):

            pattern_j = pattern_dict[cond_j]

            # -----------------------------------------
            # Pearson correlation
            # -----------------------------------------
            if metric == "pearson":

                sim = np.corrcoef(
                    pattern_i,
                    pattern_j
                )[0, 1]

            # -----------------------------------------
            # Euclidean similarity
            # -----------------------------------------
            elif metric == "euclidean_sim":

                dist = euclidean(
                    pattern_i,
                    pattern_j
                )

                sim = 1 / (1 + dist)
            
            elif metric == "dot_product":

                sim = np.dot(pattern_i, pattern_j)

            else:

                raise ValueError(
                    "metric must be 'pearson' or "
                    "'euclidean_sim' or 'dot_product'"
                )

            sim_matrix[i, j] = sim

    return pd.DataFrame(
        sim_matrix,
        index=condition_order,
        columns=condition_order
    )

from scipy.stats import kendalltau

def rsa_similarity(neural_matrix, model_matrix, comparator_method="spearmanr"):

    """
    RSA between neural similarity matrix
    and theoretical similarity model.
    """

    neural_matrix = np.asarray(neural_matrix)
    model_matrix = np.asarray(model_matrix)

    # upper triangle only
    iu = np.triu_indices(neural_matrix.shape[0], k=1)

    neural_vec = neural_matrix[iu]
    model_vec = model_matrix[iu]


    if comparator_method == "spearmanr":
        rho, pval = spearmanr(neural_vec, model_vec)

    elif comparator_method == "kendalltau":
        rho, pval = kendalltau(neural_vec, model_vec, alternative = 'greater') # one side test for pos.

    return rho, pval

def rsa_permutation_test(
    neural_matrix,
    model_matrix,
    n_permutations=10000,
    random_state=33,
    comparator_method="spearmanr"
):

    rng = np.random.default_rng(random_state)

    neural_matrix = np.asarray(neural_matrix)
    model_matrix = np.asarray(model_matrix)

    # observed statistic
    observed_rho, _ = rsa_similarity(
        neural_matrix,
        model_matrix
    )

    n = neural_matrix.shape[0]

    permuted_rhos = np.zeros(n_permutations)

    for p in range(n_permutations):

        perm_idx = rng.permutation(n)

        permuted_model = model_matrix[perm_idx][:, perm_idx]

        perm_rho, _ = rsa_similarity(
            neural_matrix,
            permuted_model,
            comparator_method = comparator_method
  
        )

        permuted_rhos[p] = perm_rho

    # one-sided
    p_value = (
        np.sum(permuted_rhos >= observed_rho) + 1
    ) / (n_permutations + 1)

    return observed_rho, p_value, permuted_rhos


from scipy.stats import spearmanr
import numpy as np


def rsa_similarity_rectangle(neural_matrix, model_matrix):

    """
    RSA between rectangular neural and model matrices.

    Expected shape:
        GS × CS
        (18 × 2)
    """

    neural_matrix = np.asarray(neural_matrix)
    model_matrix = np.asarray(model_matrix)

    if neural_matrix.shape != model_matrix.shape:

        raise ValueError(
            f"Shape mismatch: "
            f"{neural_matrix.shape} vs "
            f"{model_matrix.shape}"
        )

    neural_vector = neural_matrix.flatten()

    model_vector = model_matrix.flatten()

    rho, pval = spearmanr(
        neural_vector,
        model_vector
    )

    return rho, pval

def rsa_permutation_test_rectangle(
    neural_matrix,
    model_matrix,
    n_permutations=10000,
    random_state=None
):

    rng = np.random.default_rng(random_state)

    neural_matrix = np.asarray(neural_matrix)
    model_matrix = np.asarray(model_matrix)

    observed_rho, _ = rsa_similarity_rectangle(
        neural_matrix,
        model_matrix
    )

    n_rows = neural_matrix.shape[0]

    permuted_rhos = np.zeros(
        n_permutations
    )

    for permutation in range(n_permutations):

        shuffled_rows = rng.permutation(
            n_rows
        )

        permuted_model = model_matrix[
            shuffled_rows,
            :
        ]

        perm_rho, _ = rsa_similarity_rectangle(
            neural_matrix,
            permuted_model
        )

        permuted_rhos[
            permutation
        ] = perm_rho

    p_value = (

        np.sum(
            permuted_rhos >= observed_rho
        ) + 1

    ) / (

        n_permutations + 1
    )

    return (
        observed_rho,
        p_value,
        permuted_rhos
    )


# PLOT SIMILARITY MATRIX

def plot_similarity_matrix(
    sim_df,
    title="",
    figsize=(11,9),
    cmap="OrRd",
    show_labels=False,
    reorder_columns=None,
    rename_labels=None,
    save_path=None,
    vmin=0,
    vmax=1,
    rm_diagonal=False,
    add_contour=False,
    remove_tick_lines = True
):
    """
    Plot similarity matrix for RSA models.

    Parameters
    ----------
    sim_df : pd.DataFrame or np.ndarray
        Similarity matrix.

    title : str
        Figure title.

    figsize : tuple
        Figure size.

    cmap : str
        Colormap.

    show_labels : bool
        Whether to display stimulus labels.

    reorder_columns : list or None
        Optional reordered labels for rows/columns.

    save_path : str or None
        Path to save figure.

    vmin, vmax : float
        Color scaling limits.

    add_contour : bool
        Whether to add dark contour lines around heatmap cells and border.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    # --------------------------------------------------
    # Reorder matrix if requested
    # --------------------------------------------------

    if reorder_columns is not None:

        if isinstance(sim_df, pd.DataFrame):

            sim_df = sim_df.loc[
                reorder_columns,
                reorder_columns
            ]

        else:
            raise ValueError(
                "reorder_columns requires sim_df to be a DataFrame"
            )
        
    # RENAME LABELS IF REQUESTED
    if rename_labels is not None:

        if isinstance(sim_df, pd.DataFrame):

            sim_df = sim_df.copy()

            sim_df.index = [
                rename_labels.get(x, x)
                for x in sim_df.index
            ]

            sim_df.columns = [
                rename_labels.get(x, x)
                for x in sim_df.columns
            ]
    # --------------------------------------------------
    # Plot
    # --------------------------------------------------

    if rm_diagonal:
        sim_df = sim_df.copy()
        np.fill_diagonal(sim_df.values, np.nan)

    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        sim_df,
        cmap=cmap,
        square=True,
        cbar=True,
        xticklabels=show_labels,
        yticklabels=show_labels,
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        linecolor='white'
    )

    if add_contour:
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_color('black')

    # --------------------------------------------------
    # Labels formatting
    # --------------------------------------------------

    if show_labels:

        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=90,
            fontsize=22
        )

        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=0,
            fontsize=22
        )

    else:

        ax.set_xticks([])
        ax.set_yticks([])

    # Remove tick lines 
    if remove_tick_lines:
        ax.tick_params(left=False, bottom=False)

    # --------------------------------------------------
    # Colorbar formatting
    # --------------------------------------------------

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=30)

    # --------------------------------------------------
    # Cosmetics
    # --------------------------------------------------

    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.title(title, fontsize=25, pad=20)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=1000, bbox_inches='tight')

    plt.show()

# import seaborn as sns # bug fix with import!!
import matplotlib.pyplot as plt
from scipy.stats import linregress
import seaborn as sns

def jointplot(x, y, x_label="X", y_label="Y", 
              title=None, default_color=[0.2, 0.5, 1],
              density_color_x='#fd6262', density_color_y="#747272", alpha=1,
              r_p_input=None,
              save_as=None
              ):
    """
    Create a jointplot with regression line for X and Y with KDE marginals.

    Parameters:
    -----------
    x : array-like
        Predictor variable.
    y : array-like
        Outcome variable.
    x_label : str
        Label for the x-axis.
    y_label : str
        Label for the y-axis.
    title : str or None
        Title of the plot.
    default_color : list
        RGB color used for plotting points.
    density_color_x : str
        Color code for the X-axis marginal KDE distribution.
    density_color_y : str
        Color code for the Y-axis marginal KDE distribution.

    Returns:
    --------
    g : seaborn.axisgrid.JointGrid
        The seaborn jointplot object.
    """
    # Handle NaNs
    valid_mask = ~np.isnan(x) & ~np.isnan(y)

    x_valid = x[valid_mask]
    y_valid = y[valid_mask]

    # Compute correlation
    if r_p_input is not None:
        (r_value, p_value) = r_p_input
    else:
        slope, intercept, r_value, p_value, _ = linregress(x_valid, y_valid)
    corr_label = f'r = {r_value:.2f}, p = {p_value:.3f}'

    # Determine consistent limits
    x_margin = (x_valid.max() - x_valid.min()) * 0.08
    y_margin = (y_valid.max() - y_valid.min()) * 0.08
    xlim = (x_valid.min() - x_margin, x_valid.max() + x_margin)
    ylim = (y_valid.min() - y_margin, y_valid.max() + y_margin + 0.05) # small buffer for extreme values

    # Create JointGrid with fixed axis limits
    g = sns.JointGrid(x=x_valid, y=y_valid, height=12, xlim=xlim, ylim=ylim) # height to increase box size


    # Scatter plot
    # g.ax_joint.scatter(x_valid, y_valid, alpha=0.7, s=50, edgecolor='black', color=default_color)
    g.ax_joint.scatter(x_valid, y_valid, alpha=0.7, s=150, edgecolor='black', linewidth=1.5, color=default_color)

    g.ax_joint.tick_params(axis='both', which='major', labelsize=40, width=2.5, length=8)
    
    for spine in g.ax_joint.spines.values():
        spine.set_linewidth(3.5)


    # Add regression line
    sns.regplot(x=x_valid, y=y_valid, scatter=False, ax=g.ax_joint,
                line_kws={'color': 'black', 'linewidth': 5})

    # KDE marginal distributions with separate colors
    
    sns.kdeplot(x=x_valid, ax=g.ax_marg_x, fill=True, color=density_color_x, alpha=alpha)
    sns.kdeplot(y=y_valid, ax=g.ax_marg_y, fill=True, color=density_color_y, alpha=alpha)

    # Add correlation text with larger box
    g.ax_joint.text(0.05, 0.95, corr_label, transform=g.ax_joint.transAxes,
                    fontsize=45, verticalalignment='top', 
                    bbox=dict(boxstyle="round,pad=0.48", alpha=0.5, facecolor = 'lightgray'))

    # Add labels and title
    g.ax_joint.set_xlabel(x_label, fontsize=45, labelpad=24)
    g.ax_joint.set_ylabel(y_label, fontsize=45)
    if title:
        g.fig.suptitle(title, fontsize=45, y=1.10)

    if save_as is not None:
        g.fig.savefig(save_as, dpi=1000, bbox_inches='tight')
    
    plt.tight_layout()

    plt.show()

    return g

    # VIOLIN
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ttest_1samp, ttest_rel


def plot_rsa_model_fits(
    rsa_results_df,
    models,
    model_labels=None,
  
    figsize=None,
    ylabel="RSA model fit (ρ)",
    title=None,
    ylim=None,
    show_points=True,
    show_mean=True,
    show_zero_line=True,
    violin_alpha=0.7,
    point_size=80,
    point_alpha=0.7,
    linewidth=2.5,
    font_scale=1.0,
    save_path=None,
    dpi=1000
):
    """
    Plot subject-level RSA model fits as side-by-side violin plots.

    Parameters
    ----------
    rsa_results_df : pandas.DataFrame
        DataFrame containing one column per model and one row per subject.

    models : list of str
        Column names corresponding to the RSA models to plot.

    model_labels : list of str or None
        Labels displayed on the x-axis. If None, the column names are used.

    figsize : tuple or None
        Figure size. If None, determined automatically from the number of models.

    ylabel : str
        Y-axis label.

    title : str or None
        Optional figure title.

    ylim : tuple or None
        Y-axis limits, e.g. (-0.2, 0.6).

    show_points : bool
        Overlay individual subject values.

    show_mean : bool
        Overlay the mean as a black horizontal marker.

    show_zero_line : bool
        Add a horizontal line at y=0.

    violin_alpha : float
        Transparency of violin distributions.

    point_size : float
        Size of individual subject points.

    point_alpha : float
        Transparency of individual points.

    linewidth : float
        Width of plot borders.

    font_scale : float
        Scaling factor for font sizes.

    save_path : str or None
        Path to save the figure.

    dpi : int
        Figure resolution.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    stats_df : pandas.DataFrame
        Summary statistics for each model.
    """

    # --------------------------------------------------
    # Validate inputs
    # --------------------------------------------------

    missing_models = [
        model for model in models
        if model not in rsa_results_df.columns
    ]

    if missing_models:
        raise ValueError(
            f"The following models are not columns in rsa_results_df: "
            f"{missing_models}"
        )

    if model_labels is None:
        model_labels = models

    if len(model_labels) != len(models):
        raise ValueError(
            "`model_labels` must have the same length as `models`."
        )

    # --------------------------------------------------
    # Prepare long-format dataframe
    # --------------------------------------------------

    plot_df = rsa_results_df[models].copy()

    plot_df.columns = model_labels

    plot_df = plot_df.melt(
        var_name="Model",
        value_name="RSA_fit"
    )

    # Remove NaNs
    plot_df = plot_df.dropna(subset=["RSA_fit"])

    # Preserve requested order
    plot_df["Model"] = pd.Categorical(
        plot_df["Model"],
        categories=model_labels,
        ordered=True
    )

    # --------------------------------------------------
    # Figure size
    # --------------------------------------------------

    if figsize is None:
        figsize = (
            max(10, 2.2 * len(models)),
            9
        )

    fig, ax = plt.subplots(figsize=figsize)

    # --------------------------------------------------
    # Violin plots
    # --------------------------------------------------

    sns.violinplot(
        data=plot_df,
        x="Model",
        y="RSA_fit",
        order=model_labels,
        inner=None,
        cut=0,
        linewidth=linewidth,
        alpha=violin_alpha,
        ax=ax
    )

    # --------------------------------------------------
    # Individual subject values
    # --------------------------------------------------

    if show_points:

        sns.stripplot(
            data=plot_df,
            x="Model",
            y="RSA_fit",
            order=model_labels,
            color="black",
            size=np.sqrt(point_size),
            alpha=point_alpha,
            jitter=0.12,
            ax=ax
        )

    # --------------------------------------------------
    # Means
    # --------------------------------------------------

    if show_mean:

        means = (
            plot_df
            .groupby("Model", observed=False)["RSA_fit"]
            .mean()
            .reindex(model_labels)
        )

        for x, mean in enumerate(means):

            ax.plot(
                x,
                mean,
                marker="o",
                markersize=12,
                markeredgewidth=2,
                markeredgecolor="black",
                markerfacecolor="white",
                zorder=10
            )

    # --------------------------------------------------
    # Zero line
    # --------------------------------------------------

    if show_zero_line:

        ax.axhline(
            0,
            linestyle="--",
            linewidth=2,
            zorder=0
        )

    # --------------------------------------------------
    # Axis formatting
    # --------------------------------------------------

    ax.set_xlabel(
        "",
        fontsize=30 * font_scale
    )

    ax.set_ylabel(
        ylabel,
        fontsize=30 * font_scale
    )

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=25 * font_scale,
        width=2.5,
        length=8
    )

    for spine in ax.spines.values():

        spine.set_linewidth(linewidth)

    if title is not None:

        ax.set_title(
            title,
            fontsize=32 * font_scale,
            pad=20
        )

    if ylim is not None:

        ax.set_ylim(ylim)

    # --------------------------------------------------
    # Summary statistics
    # --------------------------------------------------

    summary_stats = []

    for model, label in zip(models, model_labels):

        values = rsa_results_df[model].dropna().to_numpy()

        t_stat, p_value = ttest_1samp(
            values,
            popmean=0
        )

        summary_stats.append({

            "model": model,

            "label": label,

            "n": len(values),

            "mean": np.mean(values),

            "std": np.std(values, ddof=1),

            "sem": np.std(values, ddof=1)
                   / np.sqrt(len(values)),

            "t": t_stat,

            "p": p_value

        })

    stats_df = pd.DataFrame(summary_stats)

    # --------------------------------------------------
    # Layout and saving
    # --------------------------------------------------

    plt.tight_layout()

    if save_path is not None:

        fig.savefig(
            save_path,
            dpi=dpi,
            bbox_inches="tight"
        )

    plt.show()

    return fig, ax, stats_df