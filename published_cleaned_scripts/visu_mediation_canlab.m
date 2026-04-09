    %% 

%addpath(genpath('E:\dsutterlin\projects\genPain\Toolboxes'));
resultDir = '/home/dsutterlin/projects/genPain/results/imaging/mediation/'
%modDir = 'E:\dsutterlin\projects\genPain\results\imaging\mediation';
projectDir = '/home/dsutterlin/projects/genPain/'
behavDir = fullfile(projectDir,'results','behavioral');

load(fullfile(behavDir, 'finalPreproc_simGen_data.mat'));
gendata = SCEBLmri_gendata;
load(fullfile(behavDir, 'SCEBLmri_learndata_FINAL_N36.mat'))
learndata = SCEBLmri_data;

% Learning betas
data.learnBetas = gendata.SClearn_beta; % already scaled
csvwrite('scaled_learnBetas.csv', data.learnBetas);
learn_beta = data.learnBetas;
cd resultDir;

%% Cue-M-Pain model
mask = which('gray_matter_mask.img');
%SETUP = mediation_brain_corrected_threshold('fdr');
SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

%% plotting variables
col_pos = [1.0 0.93 0; 1 0.667 0];
col_neg = [.074 .6235 1; 0 0 .95];
fdrSize = 3;
uncSize1 = 10;
uncSize2 = 10;
unc1 = .01;
unc2 = .05;
a_cols = {[.9 .4 0]; [1 .6 0]};
ab_cols = {[.7 .2 .3]; [1 .4 .55]};
dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('\n\n%s\n%s\n%s\n%s\n%s\n\n', dashes, dashes, str, dashes, dashes);

%% Results variables and tables
%% Final plots and stats

[clpos_gencue_A_fdr, clneg_gencue_A_fdr, clpos_data_gencue_A_fdr, clneg_data_gencue_A_fdr clpos2 clneg2] = mediation_brain_results('a',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1],'prune', 'mask', mask, 'tables', 'names', 'save');
[clpos, tmp, tableApos] = table(cluster2region(clpos_gencue_A_fdr{1}),'nolegend', 'nosort');
[tmp, clneg, tableAneg] = table(cluster2region(clneg_gencue_A_fdr{1}),'nolegend');

% Manually get Zmax in cl fdr structure
T_A_pos_ext = extract_maxZ_from_cluster_any(clpos_gencue_A_fdr{1}, [], 'pos');
T_A_neg_ext = extract_maxZ_from_cluster_any(clneg_gencue_A_fdr{1}, [], 'neg');
assert(height(tableApos) == height(T_A_pos_ext), 'A pos: row count mismatch');
assert(height(tableAneg) == height(T_A_neg_ext), 'A neg: row count mismatch');
tableApos{:,4} = T_A_pos_ext.MaxZ;   % true Z (two-tailed)
tableAneg{:,4} = T_A_neg_ext.MinZ;  
tableApos{:, 4} = round(tableApos{:, 4}, 2); %round
tableAneg{:, 4} = round(tableAneg{:, 4}, 2);
writetable(tableApos, 'results/AposTable_fdr.csv');
writetable(tableAneg, 'results/AnegTable_fdr.csv');

montage(regionA, 'regioncenters', 'colormap');

% Path A Plots
% ------
disp_A = fmridisplay; % [44 4 -38; -22 -24 -20; -40 -10 -6; 28 -44 0; 0 28 4; -34 -34 -18]
%disp_A = montage(disp_A, 'coronal', 'wh_slice', [ 0 -34 0; 0 -32 0; 0 4 0], 'onerow', 'spacing', 2);
disp_A = montage(disp_A, 'coronal', 'wh_slice', [ 0 -4 0 ; 0 -10 0; 0 -32 0], 'onerow', 'spacing', 2); % AMY, Insula Coronal
disp_A = montage(disp_A, 'saggital', 'wh_slice', [-40 0 0;-26 0 0; -22 0 0; 0 0 0; 26 0 0], 'onerow');
disp_A = montage(disp_A, 'axial', 'wh_slice', [0 0 76]);

%annotation('textbox', [0., 0.8, 0.8, .5], 'String', 'Path a : Generalization cues effect', 'FontSize', 30,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'EdgeColor', 'none');
%annotation('textbox', [0.0, 0.8, 0.1, 0], 'String', 'Left', 'EdgeColor', 'none','VerticalAlignment', 'middle', 'FontSize', 18);
%annotation('textbox', [0.01, 0.75, 0.3, 0], 'String', 'Saggital : Left (1), Right (2)', 'EdgeColor', 'none','VerticalAlignment', 'middle', 'FontSize', 12);
genpain_disp = addblobs(disp_A, clpos_gencue_A_fdr{2}, 'color', col_pos(2, :));  % path a
genpain_disp = addblobs(disp_A, clpos_gencue_A_fdr{1}, 'color', col_pos(1, :));  % path a 
genpain_disp = addblobs(disp_A, clneg_gencue_A_fdr{2}, 'color', col_neg(1, :));  % path a neg
genpain_disp = addblobs(disp_A, clneg_gencue_A_fdr{1}, 'color', col_neg(2, :));  % path a neg

% Legend
h = zeros(4, 1);
marker_size = 20;  % Adjust marker size as needed
h(1) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_pos(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(2) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_pos(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(3) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_neg(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(4) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_neg(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
legend(h, 'q < 0.05', 'p < 0.01', 'q < 0.05', 'p < 0.01','location', 'northwestoutside');
set(gca, 'FontSize', 22, 'Position', [01 0.15 .2 .4]); %(left, r, width height)

filename = 'aX-M_clusters'; 
saveas(gcf, fullfile(pwd, 'results', [filename, '.png'])); 
%% path A combined view 

disp_A = fmridisplay;
disp_A = montage(disp_A, 'axial', ...
    'wh_slice', [0 0 -66; 0 0 -10; 0 0 0; 0 0 10; 0 0 20; 0 0 30; 0 0 40; 0 0 50; 0 0 60], ...
    'onerow');
axh = axes('Position', [0 0.25 .3 .7]);  % Position vector: [left bottom width height]
disp_A = montage(disp_A, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
disp_A = removeblobs(disp_A);

% Positive clusters
disp_A = addblobs(disp_A, clpos_gencue_A_fdr{2}, 'color', col_pos(2, :));  % Second positive cluster
disp_A = addblobs(disp_A, clpos_gencue_A_fdr{1}, 'color', col_pos(1, :));  % First positive cluster
% Negative clusters
disp_A = addblobs(disp_A, clneg_gencue_A_fdr{2}, 'color', col_neg(2, :));  % Second negative cluster
disp_A = addblobs(disp_A, clneg_gencue_A_fdr{1}, 'color', col_neg(1, :));  % First negative cluster

drawnow;


%% wedge plot
% ----------
yeoLabels = tableApos.modal_label_descriptions;
unique_yeo = unique(yeoLabels);
radius = histcounts(categorical(yeoLabels));
create_figure('wedge_plot');
hh = tor_wedge_plot(radius, unique_yeo, 'outer_circle_radius', 12, 'colors', {[1 .7 0]}, 'nofigure');

%% Make scatter/violin of each sig clusters
% Load the fmri_data object
XM_indiv = fmri_data('X-M_indiv_effect.img');
all_regA = cluster2region(clpos_gencue_A_fdr{1});


mean_reg_cell = cell(length(all_regA), 1);

for i = 1:length(all_regA)
    reg = all_regA(i);  % Get the current region
    roi = region2fmri_data(reg,XM_indiv);
    mask = roi.dat ~= 0;
    indices = find(mask);
    % XM_indiv.dat is (n_voxels x n_subjects)
    data_in_roi = XM_indiv.dat(indices, :);
    mean_values = mean(data_in_roi, 1);
    mean_reg_cell{i} = mean_values;
end

patha_indices = [16, 21, 7, 20,9, 18, 19, 3];
names_A = {'L post. insula', 'S1', 'Amygdala', 'HCP', 'L PHG', 'subgenual/mPFC', 'R auditory c.', 'R temporal pole'};
selected_data = vertcat(mean_reg_cell{patha_indices})';

create_figure('Mean Betas Across Runs');
barplot_columns(selected_data,'doline','names', names_A, 'nofigure');
ylabel('Path a coefficients', 'FontSize',14);
xlabel('Significant regions (FDR < .05)')
title('GS category effect on pain-related brain activity', 'FontSize',18);

%%  convert mean ROI coeff in path a to csv to do scatter plot 

patha_indices = [16, 21, 7, 20,9, 18, 19, 3];
names_A = {'L post. insula', 'S1', 'Amygdala', 'HCP', 'L PHG', 'subgenual/mPFC', 'R auditory c.', 'R temporal pole'};
mean_reg_mat = NaN(length(mean_reg_cell{1}), length(patha_indices));
for i = 1:length(patha_indices)
    mean_reg_mat(:, i) = mean_reg_cell{patha_indices(i)};
end
T = array2table(mean_reg_mat, 'VariableNames', names_A);
writetable(T, 'results/pathA_mean_coefficients_16-21-7_20-18-19-3.csv');
disp('CSV file saved: mean_coefficients_selected.csv');


roi = clpos_data_gencue_A_fdr{1}(1)
region2fmri_data((cluster2region(clpos_data_gencue_A_fdr{1})), fmri_data('X-M_effect.img'))
[paths, stats] = mediation(gdata.SCcscat, gdata.pain, roi.dat, 'plots', 'verbose', 'boot', 'bootsamples', 10000);

ylabel(y_name_axis);
title('Within-person Generalization effects per run (n=36) (N.S. all)', 'FontSize',14);
drawnow, snapnow

save_path_betas = fullfile(genFX_save, 'genFX_4runs_noModerator_violins.png');
exportgraphics(gcf, save_path_betas, 'Resolution', resolution); % Save with high resolution

% Violin plot for the first region
figure;
violinplot(mean_reg_cell{1});
xlabel('Subjects');
ylabel('Mean Effect Size');
title('Region 1 Mean Effect Across Subjects');

figure;
num_regions = length(all_regA);
for i = 1:num_regions
    subplot(1, num_regions, i);
    violinplot(mean_reg_cell{i});
    xlabel('Subjects');
    ylabel('Mean Effect Size');
    title(['Region ' num2str(i)]);
end

% Scatter plot ROI-indiv FX 
% -------------------------
mask = fmri_data('mask.img')
rdlPFC = cluster2region(clneg_gencue_A_fdr{1}); % no idx bc only one cl
rdlPFC2 = cluster2region(clneg_gencue_A_fdr{2});
dlPFC = region2fmri_data([rdlPFC rdlPFC2], mask);
orthviews(dlPFC);
%
dlPFCFx = cell(1, size(XM_indiv.dat, 2));
meanVal = [];
for i=1:size(XM_indiv.dat, 2) % loop over col (imgs)
    
    voxelIDs = find(dlPFC.dat ~= 0); % get voxel in mask (!= 0)
    dlPFCFx{i} = XM_indiv.dat(voxelIDs, i); %at colum i get values at indices od voxels IDs
    meanVal =  [meanVal mean(XM_indiv.dat(voxelIDs, i))];
end

% Comapre with extracting individually for each subj
% Get directly the 20 vxls in dlPFC from clneg_data..
dlPFCFx = cell(1, size(XM_indiv.dat, 2));
meanFX = [];
for i=1:size(XM_indiv.dat, 2) % loop over col (imgs)
    rdlPFC20 = cluster2region(clneg_data_gencue_A_fdr{i});
    dlPFC20 = region2fmri_data(rdlPFC20, mask);
    voxelIDs = find(dlPFC20.dat ~= 0); % get voxel in mask (!= 0)
    dlPFCFx{i} = XM_indiv.dat(voxelIDs, i); %at colum i get values at indices od voxels IDs
    meanFX = [meanFX mean(XM_indiv.dat(voxelIDs, i))];
end

% Scatter plot dlPFC_fx with scale(learn_beta) on x axis
figure;
scatter(learn_beta, meanFX');
xlabel('Scaled Learn Beta');
ylabel('Mean dlPFC_fx');
title('Path A : Scatter Plot of dlPFC_fx vs Scaled Learn Beta');
grid on;
lsline;
hold off;
[r, p] = corr(learn_beta, meanVal');
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');


%% Path A load moderation mask dlPFC
% Relation btwn dlPFC FX and learning
%-----------------------------------
% %mask = fmri_data('mask.img')
% XM_indivmod = fmri_data('X-M_indiv_effect.img')
% dlPFC_mod_mask = fmri_data('dlPFC_mask_neg_fdr_pathA.nii');
% 
% % Get effect according to mask
% dlPFCFXmodA = cell(1, size(XM_indivmod.dat, 2));
% meanFXmodA = [];
% for i=1:size(XM_indivmod.dat, 2) % loop over col (imgs)
% 
%     voxelIDs = find(dlPFC_mod_mask.dat ~= 0); % get voxel in mask (!= 0)
%     dlPFCFXmodA{i} = XM_indivmod.dat(voxelIDs, i); %at colum i get values at indices od voxels IDs
%     meanFXmodA = [meanFXmodA mean(XM_indivmod.dat(voxelIDs, i))];
% 
% end
% 
% % Final mask dlPFC path A, use only fdr sig cl_gen...{1}
% %------
% dlpfc_idx = 12;
% dlPFCFXmodA = cell(1, size(XM_indivmod.dat, 2));
% meanFXmodA = [];
% 
% for i=1:size(XM_indivmod.dat, 2) % loop over col (imgs)
% 
%     rdlPFCmodA = cluster2region(clneg_gencue_modA_fdr{1}(dlpfc_idx));
%     dlPFCmodA = region2fmri_data(rdlPFCmodA, mask_out);
% 
%     voxelIDs = find(dlPFCmodA.dat ~= 0); 
%     dlPFCFXmodA{i} = XM_indivmod.dat(voxelIDs, i); %at colum i get values at indices od voxels IDs
%     meanFXmodA = [meanFXmodA mean(XM_indivmod.dat(voxelIDs, i))];
% 
% end
% 
% % Scatter plot dlPFC_fx with scale(learn_beta) on x axis
% figure;
% scatter(learn_beta, meanFXmodA');
% xlabel('Scaled Learn Beta');
% ylabel('Mean dlPFC_fx');
% title('Path A mod: Scatter Plot of dlPFC_fx vs Scaled Learn Beta');
% grid on;
% lsline;
% hold off;
% [r, p] = corr(learn_beta, meanFXmodA');
% annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
% annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');


%% Path B
% -------
[clpos_gencue_B_fdr, clneg_gencue_B_fdr, clpos_data_gencue_B_fdr, clneg_data_gencue_B_fdr] = mediation_brain_results('b', 'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1],'prune', 'mask', mask, 'names', 'save');
[clpos, tmp, tableBpos] = table(cluster2region(clpos_gencue_B_fdr{1}),'nolegend');
[tmp, clneg, tableBneg] = table(cluster2region(clneg_gencue_B_fdr{1}),'nolegend');

T_B_pos_ext = extract_maxZ_from_cluster_any(clpos_gencue_B_fdr{1}, [], 'pos');
T_B_neg_ext = extract_maxZ_from_cluster_any(clneg_gencue_B_fdr{1}, [], 'neg');
tableBpos{:,4} = round(T_B_pos_ext.MaxZ,2);
tableBneg{:,4} = round(T_B_neg_ext.MinZ,2);
writetable(tableBpos, 'results/BposTable_fdr.csv');
writetable(tableBneg, 'results/BnegTable_fdr.csv');


%% Path ab
%**cluster size for fdr = 10 (vs 3 with other path) : removes single salt and pepper
fdrSize = 3;
uncSize1 = 10;
uncSize2 = 10;

[clpos_gencue_AB_fdr, clneg_gencue_AB_fdr, clpos_data_gencue_AB_fdr, clneg_data_gencue_AB_fdr] = mediation_brain_results('ab', 'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1],'prune', 'mask', mask, 'names', 'save');
[clpos, tmp, tableABpos] = table(cluster2region(clpos_gencue_AB_fdr{1}),'nolegend');
[clpos, tmp, tableAABpos] = table(cluster2region(clpos_gencue_AB_fdr{1}),'nolegend');
[tmp, clneg, tableABneg] = table(cluster2region(clneg_gencue_AB_fdr{1}),'nolegend');
% Manually get extreme Z in cl fdr structure
T_AB_pos_ext = extract_maxZ_from_cluster_any(clpos_gencue_AB_fdr{1}, [], 'pos');  % Z from ln(1/p)
T_AB_neg_ext = extract_maxZ_from_cluster_any(clneg_gencue_AB_fdr{1}, [], 'neg');
tableABpos{:,4} = round(T_AB_pos_ext.MaxZ,2);
tableABneg{:,4} = round(T_AB_neg_ext.MinZ,2);
writetable(tableABpos, 'results/ABposTable_fdr.csv');
writetable(tableABneg, 'results/ABnegTable_fdr.csv');


xyz = cat(1, clpos_gencue_AB_fdr{1}.mm_center);
% no prune option for plotting only
[clpos_gencue_AB_fdr_30vox, tmp, clpos_data_gencue_AB_fdr, clneg_data_gencue_AB_fdr] = mediation_brain_results('ab', 'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize 30], 'mask', mask, 'names', 'save');

disp_AB = fmridisplay; % [44 4 -38; -22 -24 -20; -40 -10 -6; 28 -44 0; 0 28 4; -34 -34 -18]
disp_AB = montage(disp_AB, 'saggital', 'wh_slice', [-64 0 0;-50 0 0; -18 0 0; -10 0 0; 62 0 0], 'onerow', 'spacing', 1); % L hippoc?
%disp_AB = montage(disp_AB, 'axial', 'wh_slice', [0 0 -6;], 'onerow');
disp_AB = montage(disp_AB, 'coronal', 'wh_slice', [0 8 0; 0 8 0; 0 -30 0; 0 0 0], 'onerow'); % 8 AMY to add to full plot
%dip_AB = montage(disp_AB, 'coronal', 'wh_slice', [0 58 0; 0 8 0; 0 -30 0; 0 -70 0], 'onerow'); % 8 AMY to add to full plot

genpain_disp = addblobs(disp_AB, clpos_gencue_AB_fdr_30vox{2}, 'color', col_pos(2, :));  % path a
genpain_disp = addblobs(disp_AB, clpos_gencue_AB_fdr{1}, 'color', col_pos(1, :));  % path a 


h = zeros(4, 1);
marker_size = 20;  % Adjust marker size as needed
h(1) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_pos(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(2) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_pos(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(3) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_neg(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
h(4) = plot(NaN, NaN, 's', 'MarkerFaceColor', col_neg(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', marker_size);
legend(h, 'q < 0.05', 'p < 0.01', 'q < 0.05', 'p < 0.01','location', 'westoutside');
set(gca, 'FontSize', 18, 'Position', [1 0.30 .1 .2]); %(left, r, width height)

%% path ab apply mvpa signature

indiv_ab_fx = fmri_data('X-M-Y_indiv_effect.img')
[nps_values, image_names, data_objects]  = apply_nps(indiv_ab_fx)
[siips_values, image_names, data_objects, siipspos_exp_by_region, siipsneg_exp_by_region, siipspos, siipsneg] = apply_siips(indiv_ab_fx)
% save

%% Moderation
%% plot variables

cd ../cue_B_pain_L2M
mask = which('gray_matter_mask.img');
%SETUP = mediation_brain_corrected_threshold('fdr');
SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

%% Moderator path A noL2M (just for qualitative comparison : supp mat,
[clpos_gencue_A_fdr, clneg_gencue_A_fdr, clpos_data_gencue_A_fdr, clneg_data_gencue_A_fdr clpos2 clneg2] = mediation_brain_results('a',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1],'prune', 'mask', mask, 'tables', 'names', 'save');
[clpos, tmp, tableApos] = table(cluster2region(clpos_gencue_A_fdr{1}),'nolegend', 'nosort');
[tmp, clneg, tableAneg] = table(cluster2region(clneg_gencue_A_fdr{1}),'nolegend');
writetable(tableApos, 'results/AposTable_fdr.csv');
writetable(tableAneg, 'results/AnegTable_fdr.csv');


%% Moderator of path A L2M
[clpos_gencue_modA_fdr, clneg_gencue_modA_fdr, clpos_data_gencue_modA_fdr, clneg_data_gencue_modA_fdr clpos2 clneg2] = mediation_brain_results('al2mod',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1],'prune', 'mask', mask, 'tables', 'names', 'save');
[clpos, tmp, tablemodApos] = table(cluster2region(clpos_gencue_modA_fdr{1}),'nolegend');
[tmp, clneg, tablemodAneg] = table(cluster2region(clneg_gencue_modA_fdr{1}),'nolegend');
writetable(tablemodApos, 'results/modA_l2m_posTable_fdr.csv');
writetable(tablemodAneg, 'results/modA_l2m_negTable_fdr.csv');

cust_size = 55
[tmp1, clneg_gencue_modA_fdr_55vox, tmp2, tmp3, clpos2 clneg2] = mediation_brain_results('al2mod',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [cust_size uncSize1],'prune', 'mask', mask, 'tables', 'names', 'save');


%% path a l2m : load _positive_ cue_B_pain effect and get data from l2mod regions
% --> scatter plot for interaction in dlPDC and HCP bilat.

XM_indiv_a_nomod = fmri_data('../cue_B_pain/X-M_indiv_effect.img');
%load('../cue_B_pain/cl_X-M_pvals_003_01_k10_5_prune.mat')
% use mod a l2mod : clpos_gencue_modA_fdr 
all_regA_l2m = cluster2region(clpos_gencue_modA_fdr{1});

%positive cluster : HCP bilat
mean_reg_cell = cell(length(all_regA_l2m), 1);
for i = 1:length(all_regA_l2m)
    reg = all_regA_l2m(i);  % Get the current region
    roi = region2fmri_data(reg,XM_indiv_a_nomod);
    mask = roi.dat ~= 0;
    indices = find(mask);
    % XM_indiv.dat is (n_voxels x n_subjects)
    data_in_roi = XM_indiv_a_nomod.dat(indices, :);
    mean_values = mean(data_in_roi, 1);
    mean_reg_cell{i} = mean_values;
end

patha_indices = [8, 7] % L and R HCP
names_A_nomod = {'L hippocampus', 'R hippocampus'}
patha_nomod_pos = vertcat(mean_reg_cell{patha_indices})';

output_path = '../cue_B_pain/results/mean_pos_pathA_activity_in_Al2mod_rois_-7-8.csv';
T = array2table(patha_nomod_pos, 'VariableNames', names_A_nomod);
writetable(T, output_path);

%% Moderation b L2M path
clear all
[clpos_gencue_l2modB_fdr, clneg_gencue_l2modB_fdr, clpos_data_gencue_l2modB_fdr, clneg_data_gencue_l2modB_fdr clpos2 clneg2] = mediation_brain_results('bl2mod',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1], 'prune','mask', mask, 'tables', 'save');
[tmp, clneg, tablel2modBpos] = table(cluster2region(clpos_gencue_l2modB_fdr{1}),'nolegend');
[tmp, clneg, tablel2modBneg] = table(cluster2region(clneg_gencue_l2modB_fdr{1}),'nolegend');
writetable(tablel2modBpos, 'results/Bl2mod_posTable_fdr.csv');
writetable(tablel2modBneg, 'results/Bl2mod_negTable_fdr.csv');

%% saving for interaction b 
MY_indiv_a_nomod = fmri_data('../cue_B_pain/M-Y_indiv_effect.img'); %! from no l2m res
%load('../cue_B_pain/cl_X-M_pvals_003_01_k10_5_prune.mat')
% use mod a l2mod : clpos_gencue_modA_fdr 
all_regB_l2m = cluster2region(clpos_gencue_l2modB_fdr{1}); %!! from l2mod results

%positive cluster : HCP bilat
mean_reg_cell = cell(length(all_regB_l2m), 1);
for i = 1:length(all_regB_l2m)
    reg = all_regB_l2m(i);  % Get the current region
    roi = region2fmri_data(reg,MY_indiv_a_nomod);
    mask = roi.dat ~= 0;
    indices = find(mask);
    % XM_indiv.dat is (n_voxels x n_subjects)
    data_in_roi = MY_indiv_a_nomod.dat(indices, :);
    mean_values = mean(data_in_roi, 1);
    mean_reg_cell{i} = mean_values;
end

patha_indices = [8, 9, 16]
names_A_nomod = {'R OFC', 'R dlPFC', 'R inf. parietal cortex'}
patha_nomod_neg = vertcat(mean_reg_cell{patha_indices})';

output_path = '../cue_B_pain/results/mean_neg_pathB_activity_in_Bl2mod_rois_8_9_16.csv';
T = array2table(patha_nomod_neg, 'VariableNames', names_A_nomod);
writetable(T, output_path);


%% moderation path ab (no l2m) for comparison with cue-B-pain path ab model
fdrSize = 10;
uncSize1 = 5;
uncSize2 = 10;
[clpos_gencue_modAB_fdr, clneg_gencue_modAB_fdr, clpos_data_gencue_modAB_fdr, clneg_data_gencue_modAB_fdr clpos2 clneg2] = mediation_brain_results('ab',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1], 'prune','mask', mask, 'tables', 'save');
[clpos, tmp, tablemodABpos] = table(cluster2region(clpos_gencue_modAB_fdr{1}),'nolegend');
[tmp, clneg, tablemodABneg] = table(cluster2region(clneg_gencue_modAB_fdr{1}),'nolegend');
writetable(tablemodABpos, 'results/ABposTable_fdr.csv');
writetable(tablemodABneg, 'results/ABnegTable_fdr.csv');

%% ModerationPath ab : used in paper!!---------------
[clpos_gencue_modAB_fdr, clneg_gencue_modAB_fdr, clpos_data_gencue_modAB_fdr, clneg_data_gencue_modAB_fdr clpos2 clneg2] = mediation_brain_results('abl2mod',  'thresh', [SETUP.fdr_p_thresh unc1 ], 'size', [fdrSize uncSize1], 'prune','mask', mask, 'tables', 'save');
[clpos, tmp, tablemodABpos] = table(cluster2region(clpos_gencue_modAB_fdr{1}),'nolegend');
[tmp, clneg, tablemodABneg] = table(cluster2region(clneg_gencue_modAB_fdr{1}),'nolegend');
% Manually get extreme Z in cl fdr structure
T_ABmod_pos_ext = extract_maxZ_from_cluster_any(clpos_gencue_modAB_fdr{1}, [], 'pos');  % Z from ln(1/p)
tablemodABpos{:,4} = round(T_ABmod_pos_ext.MaxZ,2);
writetable(tablemodABpos, 'results/modABpos_lm2_Table_fdr.csv');

%% moderation ab ROI from abl2mod, coeff from cue_B_pain ab img
XMY_indiv_ab_nomod = fmri_data('../cue_B_pain/X-M-Y_indiv_effect.img'); %! from no l2m res
%load('../cue_B_pain/cl_X-M_pvals_003_01_k10_5_prune.mat')
% use mod a l2mod : clpos_gencue_modA_fdr 
all_regAB_l2m = cluster2region(clpos_gencue_modAB_fdr{1}); %!! from l2mod results

%positive cluster : HCP bilat
mean_reg_cell = cell(length(all_regAB_l2m), 1);
for i = 1:length(all_regAB_l2m)
    reg = all_regAB_l2m(i);  % Get the current region
    roi = region2fmri_data(reg,XMY_indiv_ab_nomod);
    mask = roi.dat ~= 0;
    indices = find(mask);
    % XM_indiv.dat is (n_voxels x n_subjects)
    data_in_roi = XMY_indiv_ab_nomod.dat(indices, :);
    mean_values = mean(data_in_roi, 1);
    mean_reg_cell{i} = mean_values;
end

pathab_indices = [2, 3]
names_A_nomod = {'L hippocampus', 'L caudate'}
pathab_mod_pos = vertcat(mean_reg_cell{pathab_indices})';

output_path = '../cue_B_pain/results/mean_pos_pathAB_activity_in_ABl2mod_rois_2-3.csv';
T = array2table(pathab_mod_pos, 'VariableNames', names_A_nomod);
writetable(T, output_path);

% HCP moderation scatter :  cluster in relation to learning
% -------------------------------------
XYM_indiv = fmri_data('X-M-Y_indiv_effect.img');
mask = fmri_data('mask.img');
hpc_idx = 1; %see idx in cl struct for ROI of interest
hcpmod = cluster2region(clpos_gencue_modAB_fdr{1}(hpc_idx)); % check if adding prune {2} changes
hcpmod = region2fmri_data(hcpmod, mask);
orthviews(hcpmod); % good!

% Extract effect in mask 
hcpFXmodAB = cell(1, size(XYM_indiv.dat, 2));
hcpmeanFXmodAB = [];
for i=1:size(XYM_indiv.dat, 2) % loop over col (imgs)/subj
  
    voxelIDs = find(hcpmod.dat ~= 0); % get voxel in mask (!= 0)
    hcpFx{i} = XYM_indiv.dat(voxelIDs, i); %at colum i get values at indices od voxels IDs
    hcpmeanFXmodAB = [hcpmeanFXmodAB mean(XYM_indiv.dat(voxelIDs, i))];

end

% Scatter plot HCP modAB with scale(learn_beta) on x axis
figure;
scatter(scale(learn_beta), hcpmeanFXmodAB');
xlabel('Scaled Learn Beta', 'FontSize',16);
ylabel('HCP ab coefficients', 'FontSize',16);
title('Interaction of explicit learning in HCP', 'FontSize',20);
grid on;
lsline;
hold off;
[r, p] = corr(learn_beta, hcpmeanFXmodAB');
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');


