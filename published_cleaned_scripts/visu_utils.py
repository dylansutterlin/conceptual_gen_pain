import os as os
import pandas as pd
import numpy as np
import glob as glob
from os.path import join
import nibabel as nb
import nilearn.image as image
from nilearn import plotting
from nilearn.plotting import plot_stat_map
from nilearn.glm import threshold_stats_img
from nilearn.image import math_img
from nilearn import datasets
from nilearn.plotting import plot_stat_map
import matplotlib.pyplot as plt
from nilearn.image import resample_to_img


def plot_slices(img, axis, coords, **kwargs):

    plot_stat_map(img, display_mode=axis, cut_coords=coords,cmap='hot', **kwargs)


def view_results(img, title = None, cmap='cold_hot'): 

    bg_img =  datasets.load_mni152_template(resolution=1)
    new_img = resample_to_img(img, bg_img, interpolation='nearest')
    
    plots = []
    if title is None:
        title = "Effect Size Map"

    stats = plotting.plot_stat_map(
        new_img,
        bg_img=bg_img,
        threshold=0,  # Only show significant voxels
        display_mode='z',  # Axial slices
        cut_coords=None,
        colorbar=True,
        cmap=cmap,
        title=title,
    )

    plots.append(stats)
    # Display the masked effect image interactively
    view = plotting.view_img(new_img, threshold='95%', cmap='cold_hot')
    plots.append(view)

    return plots


def save_slices(img, coords_to_plot, save_to, img_id,cmap='hot', **kwargs):
    ticksfontsize = 14

    bg_img = datasets.load_mni152_template(resolution=1)
    new_img = resample_to_img(img, bg_img, interpolation='nearest')

    for axis, coord_list in coords_to_plot.items():
        for c in coord_list:
            fig, ax = plt.subplots(figsize=(2, 2))  # matched style
            disp = plot_stat_map(
                new_img,  
                bg_img=bg_img,
                cmap=cmap,
                colorbar=False,
                dim=0,
                black_bg=False,
                display_mode=axis,
                axes=ax,
                vmax=None,
                cut_coords=(c,),
                alpha=1,
                annotate=False,
                interpolation='nearest'
            )
            # disp.annotate(size=ticksfontsize, left_right=False, xy=(0.5, -0.08))

            # disp.annotate(size=ticksfontsize, left_right=False, xy=(0.5, -0.08))
            if axis == 'z': # bigger brain on z cuts
                x_pos = 0.5
                y_pos = -0.029
            elif axis == 'y':
                x_pos = 0.5
                y_pos = 0.1
            elif axis == 'x':
                x_pos = 0.65
                y_pos = 0.145

            ax.text(
                x_pos, y_pos, f'{axis} = {c}',
                transform=ax.transAxes,
                fontsize=ticksfontsize,
                ha='center', va='top',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7, edgecolor='none')
            )

            fig.savefig(
                os.path.join(save_to, f'{img_id}_sig_{axis}{c}.png'),
                transparent=True, bbox_inches='tight', dpi=600
            )
            plt.close(fig)

# ===================================
# color bars
# ===================================

# Adapted from Nilearn (https://github.com/nilearn/nilearn),
# Copyright (c) 2007-2024 The Nilearn developers,
# BSD-3-Clause license

from numbers import Number
import numpy as np



def get_cbar_ticks(vmin, vmax, offset, n_ticks=5):
    """Extracted from Nilearn, adapted for standalone use."""
    if vmin == vmax:
        return np.linspace(vmin, vmax, 1)

    if vmax == 0:
        vmax += np.finfo(np.float32).eps

    ticks = np.linspace(vmin, vmax, n_ticks)
    if offset is not None and offset / vmax > 0.12:
        diff = [abs(abs(tick) - offset) for tick in ticks]
        if diff.count(min(diff)) == 4:
            idx_closest = np.sort(np.argpartition(diff, 4)[:4])
            idx_closest = np.isin(ticks, np.sort(ticks[idx_closest])[1:3])
        else:
            idx_closest = np.sort(np.argpartition(diff, 2)[:2])
            if 0 in ticks[idx_closest]:
                idx_closest = np.sort(np.argpartition(diff, 3)[:3])
                idx_closest = idx_closest[[0, 2]]
        ticks[idx_closest] = [-offset, offset]
    if len(ticks) > 0 and ticks[0] < vmin:
        ticks[0] = vmin
    return ticks


def get_colorbar_and_data_ranges(
    stat_map_data,
    vmin=None,
    vmax=None,
    symmetric_cbar=True,
    force_min_stat_map_value=None,
):
    """Extracted from Nilearn, adapted for standalone use."""
    if (not isinstance(vmin, Number)) or (not np.isfinite(vmin)):
        vmin = None
    if (not isinstance(vmax, Number)) or (not np.isfinite(vmax)):
        vmax = None

    if hasattr(stat_map_data, "_mask"):
        stat_map_data = np.asarray(
            stat_map_data[np.logical_not(stat_map_data._mask)]
        )

    if force_min_stat_map_value is None:
        stat_map_min = np.nanmin(stat_map_data)
    else:
        stat_map_min = force_min_stat_map_value
    stat_map_max = np.nanmax(stat_map_data)

    if symmetric_cbar == "auto":
        if vmin is None or vmax is None:
            min_value = stat_map_min if vmin is None else max(vmin, stat_map_min)
            max_value = stat_map_max if vmax is None else min(stat_map_max, vmax)
            symmetric_cbar = min_value < 0 < max_value
        else:
            symmetric_cbar = np.isclose(vmin, -vmax)

    if symmetric_cbar:
        if vmin is None and vmax is None:
            vmax = max(-stat_map_min, stat_map_max)
            vmin = -vmax
        elif vmin is None:
            vmin = -vmax
        elif vmax is None:
            vmax = -vmin
        elif not np.isclose(vmin, -vmax):
            raise ValueError(
                "vmin must be equal to -vmax unless symmetric_cbar is False."
            )
        cbar_vmin = vmin
        cbar_vmax = vmax
    else:
        negative_range = stat_map_max <= 0
        positive_range = stat_map_min >= 0
        if positive_range:
            cbar_vmin = 0 if vmin is None else vmin
            cbar_vmax = vmax
        elif negative_range:
            cbar_vmax = 0 if vmax is None else vmax
            cbar_vmin = vmin
        else:
            cbar_vmin = vmin
            cbar_vmax = vmax

    if vmin is None:
        vmin = stat_map_min
    if vmax is None:
        vmax = stat_map_max

    return cbar_vmin, cbar_vmax, float(vmin), float(vmax)

import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize


def save_colorbar(img, cmap, outpath, symmetric_cbar=False, offset=None, n_ticks=5, transparent=True):
    img_data = img.get_fdata()

    import numpy as np # bug fix, somehow numpy is not imported in the scritp

    # Automatically set vmax and vmin based on actual data range
    data_max = np.nanmax(img_data)
    data_min = np.nanmin(img_data)

    # Important adjustment here:
    if symmetric_cbar:
        vmax = max(abs(data_min), abs(data_max))
        vmin = -vmax
    else:
        vmax = data_max
        vmin = data_min  # <-- set to actual data min, not 0 explicitly

    cbar_vmin, cbar_vmax, _, _ = get_colorbar_and_data_ranges(
        img_data,
        symmetric_cbar=symmetric_cbar,
        vmax=vmax,
        vmin=vmin
    )

    # Guard against None values
    if cbar_vmin is None:
        cbar_vmin = data_min
    if cbar_vmax is None:
        cbar_vmax = data_max

    # Generate ticks safely
    ticks = np.linspace(cbar_vmin, cbar_vmax, n_ticks)

    # Plotting the standalone colorbar
    fig, ax = plt.subplots(figsize=(0.4, 4))
    norm = Normalize(vmin=cbar_vmin, vmax=cbar_vmax)
    cb = ColorbarBase(ax, cmap=cmap, norm=norm, ticks=ticks, orientation='vertical')
    cb.ax.set_yticklabels([f'{t:.2f}' for t in ticks], fontsize=16)

    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.yaxis.set_label_position('left')
    cb.outline.set_visible(False)

    fig.savefig(outpath, transparent=transparent, bbox_inches='tight', dpi=300)
    plt.close(fig)

    return fig


# ===================================
# Behavioral interaction plot
# ===================================

# import seaborn as sns # bug fix with import!!
import matplotlib.pyplot as plt
from scipy.stats import linregress


def jointplot(x, y, x_label="X", y_label="Y", 
              title=None, default_color=[0.2, 0.5, 1],
              density_color_x='#4c80ff', density_color_y="#fd6262", alpha=1):
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
    slope, intercept, r_value, p_value, _ = linregress(x_valid, y_valid)
    corr_label = f'r = {r_value:.2f}, p = {p_value:.3f}'

    # Determine consistent limits
    x_margin = (x_valid.max() - x_valid.min()) * 0.05
    y_margin = (y_valid.max() - y_valid.min()) * 0.05
    xlim = (x_valid.min() - x_margin, x_valid.max() + x_margin)
    ylim = (y_valid.min() - y_margin, y_valid.max() + y_margin)

    # Create JointGrid with fixed axis limits
    g = sns.JointGrid(x=x_valid, y=y_valid, height=10, xlim=xlim, ylim=ylim) # height to increase box size
    # g.fig.set_size_inches(10, 11)  # width, height in inches


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
                    fontsize=35, verticalalignment='top', 
                    bbox=dict(boxstyle="round,pad=0.48", alpha=0.5, facecolor = 'lightgray'))

    # Add labels and title
    g.ax_joint.set_xlabel(x_label, fontsize=45, labelpad=24)
    g.ax_joint.set_ylabel(y_label, fontsize=45)
    if title:
        g.fig.suptitle(title, fontsize=36, y=0.97)


    plt.tight_layout()

    plt.show()

    return g


import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from matplotlib.ticker import MaxNLocator

def plot_pain_violin(low_pain, high_pain, value_col='value', title=' ', ylabel=' ', save_as=None, x_label = 'Generalization Category',ticks_labels = ['Low Pain', 'High Pain'],
                     y_label = 'Pain Rating', palette=None):

    # Perform independent t-test
    t_stat, p_value = ttest_ind(low_pain, high_pain, nan_policy='omit')

    # Prepare DataFrame for plotting
    data_plot = pd.DataFrame({
        'Pain Rating': np.concatenate([low_pain, high_pain]),
        'Condition': [ticks_labels[0]] * len(low_pain) + [ticks_labels[1]] * len(high_pain)
    })

    # Set white background style
    if palette is None:
        sns.set_style('white')
        palette = ['#4c80ff', '#fd6262']

    fig = plt.figure(figsize=(12, 12))
    ax = sns.violinplot(
        x='Condition', y='Pain Rating', data=data_plot,
        hue='Condition', palette=palette, dodge=False,
        scale='width', inner=None, linewidth=2, cut=2
    )

    # Add lines connecting each subject's low/high pain
    for lp, hp in zip(low_pain, high_pain):
        ax.plot([0, 1], [lp, hp], color='gray', alpha=0.7, zorder=0)

    # Set y-axis starting at 0
    ax.set_ylim(top=100)

    # Clip violins: left half for Low Pain, right half for High Pain
    n_cond = data_plot['Condition'].nunique()
    for i, violin in enumerate(ax.collections[:n_cond]):
        path = violin.get_paths()[0]
        bbox = path.get_extents()
        x0, y0, width, height = bbox.bounds
        if i == 0:
            clip_rect = plt.Rectangle((x0, y0), width / 2, height, transform=ax.transData)
        else:
            clip_rect = plt.Rectangle((x0 + width / 2, y0), width / 2, height, transform=ax.transData)
        violin.set_clip_path(clip_rect)

    # Add boxplot outline over violins
    sns.boxplot(
        x='Condition', y='Pain Rating', data=data_plot,
        showcaps=False,
        boxprops={'facecolor': 'none', 'zorder': 3, 'linewidth': 3},
        whiskerprops={'linewidth': 3, 'color': 'k'},
        showfliers=False,
        showmeans=True,
        meanline=True,
        meanprops={'linewidth': 3, 'color': 'black'},
        medianprops={'visible': False, 'linewidth': 0, 'color': (0, 0, 0, 0)},
        width=0.3, saturation=1, ax=ax
    )

    # Overlay stripplot, shift dots: right for Low Pain, left for High Pain
    old_collections = len(ax.collections)
    sns.stripplot(
        x='Condition', y='Pain Rating', data=data_plot,
        hue='Condition', palette=palette, dodge=False,
        size=12, alpha=0.7, ax=ax
    )
    for i, dots in enumerate(ax.collections[old_collections:]):
        offsets = dots.get_offsets()
        shift = 0.12 if i == 0 else -0.12
        dots.set_offsets(offsets + np.array([shift, 0]))

    # Restore x limits
    xlim = ax.get_xlim()
    ax.set_xlim(xlim)

        # Thicken both x and y-axis ticks
    ax.tick_params(axis='both', width=3.5, length=8)

    # Format
    ax.yaxis.set_major_locator(MaxNLocator(nbins='auto', integer=True))

    #thicker square
        # Make the plot frame (spines) thicker
    for spine in ax.spines.values():
        spine.set_linewidth(3.5)  

    plt.title(title, fontsize=40, y=1.02)

    plt.xlabel(x_label, fontsize=55, labelpad=16)
    plt.ylabel(y_label, fontsize=55)
    plt.xticks(fontsize=55)
    plt.yticks(fontsize=45)
    plt.tight_layout()

    if save_as is not None:
        plt.savefig(save_as, dpi=1000)


    plt.show()

    return fig