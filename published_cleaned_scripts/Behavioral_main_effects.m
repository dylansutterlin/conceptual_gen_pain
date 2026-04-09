%glmfit_multilevel from canlab toolbox for main effects of learning and
%gen. task


clc
close all
clear all

projectDir = '/home/dsutterlin/projects/genPain/';
dataDir = fullfile(projectDir, 'DATA/Behavioral/')
saveDir = fullfile(projectDir, 'results/behavioral')
addpath(genpath(fullfile(projectDir, 'Toolboxes/')));

cd (saveDir);
load finalPreproc_simGen_data.mat SCEBLmri_gendata % if first time running pipeline
gdata = SCEBLmri_gendata;

load finalPreproc_simGen_data.mat % else
%load simrate.mat
load SCEBLmri_learndata_FINAL_N36.mat SCEBLmri_data
%load SCEBLmri_data_gen36.mat
gdata_df = readtable('SCEBLmri_Gendata_TxT_N36.csv');
ldata_df = readtable('SCEBLmri_Learndata_TxT_N36.csv');

% loading and rename learn and gen data for regression
ldata = SCEBLmri_data;

subjects = filenames(fullfile(dataDir,'SCEBL_MRI2*'));
effectcolors = {[.5 .5 .5] [.7 .4 .4] [.4 .4 .7]};
u = 1:36; %sub id
%[Y, X] = deal(cell(1, length(u)));

%% Add rescaled variables

for cell=1:size(gdata.pain, 2)

    currentCell = gdata.pain{cell}; % pain
    gdata.SCpain{cell} = scale(currentCell);

    currentCell = gdata.cscat{cell}; % Category
    gdata.SCcscat{cell} = scale(currentCell);

    currentCell = gdata.empiricalSimr{cell}; % empiricalSimr
    gdata.SCempiricalSimr{cell} = scale(currentCell);

    currentCell = gdata.boulderLSASimr{cell}; % boulderLSASimr
    gdata.SCboulderLSASimr{cell} = scale(currentCell);

    currentCell = gdata.mean_modelSimr{cell}; % mean_modelSimr
    gdata.SCmean_modelSimr{cell} = scale(currentCell);

end
clear cell
%% LEARN
%% EXP ~ CUES
% Extract Learners form ! Learners (learning betas)
%=========
Y_name = 'exp';
X_names = {'cues'};
X= ldata.cues;
X2_names = {'2nd-level Intercept (average within-person effects)'};
cueExp = glmfit_multilevel(ldata.exp,ldata.cues, [], 'names', X_names, 'beta_names', X2_names, 'weighted')

learn_beta = cueExp.first_level.beta(2,:)';
learners = find(learn_beta > median(learn_beta));
nonlearners = find(learn_beta <= median(learn_beta));


% plot learner vs non-learners slopes
create_figure('Learners cue effect on pain');
line_plot_multisubject(X(learners), ldata.pain(learners));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

create_figure('Non-Learners cue effect on pain');
line_plot_multisubject(X(nonlearners), ldata.pain(nonlearners));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

% Add scaled learn betas to gdata struct !!
gdata.learnBetas = learn_beta;
gdata.SClearn_beta = scale(learn_beta);

%% SAVING GEN data with scaled var, used in future analyses

save('finalPreproc_simGen_data.mat', 'gdata')

dark_blue = [0.2, 0.5, 1];
%% PAIN ~ CUES of 49 deg (medium) trials
%======================================
%Apply filter on trials to keep only trials with medium temp.
for sub = 1:numel(subjects)
    mask = ldata.temp{sub} == 49;
    ldata.cues49deg{sub} = ldata.cues{sub}(mask);
    ldata.pain49deg{sub} = ldata.pain{sub}(mask);
    ldata.exp49deg{sub} = ldata.exp{sub}(mask);
end
Y_name = 'pain';
X_var_names = {'Cues from medium pain trials'};
X2_names = {'2nd-level Intercept (average within-person effects)'};
X2 = learn_beta
cuePain49 = glmfit_multilevel(ldata.pain49deg, ldata.cues49deg, [], 'names', X_var_names, 'beta_names', X2_names, 'weighted')

beta_means = cuePain49.first_level.beta; % Means for the bar plot
pain_modulation_beta = cuePain49.first_level.beta(2,:)'
beta_std = std(cuePain49.first_level.beta, 0, 2); % Standard deviation along the correct dimension

% Multilevel effect plot
create_figure('2nd level stats');
barplot_columns(cuePain49.first_level.beta', 'names', cuePain49.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow

% Colors gradient as a function of learning
min_beta = min(learn_beta);
max_beta = max(learn_beta);
color_scale = colormap(parula); % Using parula colormap for a smooth transition from light yellow to orange
normalized_beta = (learn_beta - min_beta) / (max_beta - min_beta);
color_indices = ceil(normalized_beta * (size(color_scale, 1) - 1)) + 1;
rgb_colors = cell(size(color_indices));
for i = 1:length(color_indices)
    rgb_colors{i} = color_scale(color_indices(i), :);
end

% Multiline with bar plots for mean CS -----------
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain, 'colors', rgb_colors, 'exclude_low_range_Y', 'LineWidth', 2); % Increase line width
set(gca, 'XTick', [-0.411, 0.411], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});
hold on;

conditions = [-1, 1]; % Unique conditions
mean_pain49 = zeros(1, length(conditions)); % Preallocate mean pain array
filtered_49 = (ldata_df.temp == 49);

% Calculate mean pain scores for each condition
for i = 1:length(conditions)
    condition_indices = (ldata_df.cues == conditions(i)) & filtered_49;
    mean_pain49(i) = nanmean(ldata_df.pain(condition_indices));
    std_pain49(i) = nanstd(ldata_df.pain(condition_indices)); 
    mean_exp49(i) = nanmean(ldata_df.exp(condition_indices));
end

% Add bar plots of the mean pain scores
bar_handle = bar(conditions, mean_pain49);
bar_handle.FaceColor = dark_blue; % Light blue color
bar_handle.FaceAlpha = 0.8; % Transparency
ylim([0 80]);
hold on;

x_bar_positions = bar_handle.XData; % Get the x positions of the bars
errorbar(x_bar_positions, mean_pain49, std_pain49, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Add error bars
hold on;


xlabel('Conditioned cue conditions', 'FontSize', 20);
ylabel('Pain ratings', 'FontSize', 20);
title('Conditioned cues (CS) effect on pain', 'FontSize', 17);

hLegend = legend('Individual slopes', 'FontSize', 17);
hold on;

p = patch(NaN, NaN, dark_blue, 'FaceAlpha', 0.5);
set(get(get(p, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'on');
colormap(parula);
cbar = colorbar('Location', 'eastoutside');
cbar.Ticks = linspace(0, 1, 10);
rounded_tick_labels = round(linspace(min_beta, max_beta, 10)); % Round the tick labels
cbar.TickLabels = num2cell(rounded_tick_labels);
%cbar.Label.String = 'Explicit learning (betas)';    
hold off;

%legend(bar_handle, 'Mean pain (49 deg.)'); % Use the bar handle to ensure legend color matches


% Add a horizontal line and asterisk to indicate significance
y = max(mean_pain49) + 5; % Y position of the horizontal line and asterisk
shortening_fraction = 0.288; % Fraction to shorten the line on each side
x1 = conditions(1) + shortening_fraction * (conditions(2) - conditions(1)); % Adjusted start point
x2 = conditions(2) - shortening_fraction * (conditions(2) - conditions(1)); % Adjusted end point
line([x1, x2], [y, y], 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off'); % Horizontal line
text(mean([x1, x2]), y, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18, 'FontWeight', 'bold'); % Asterisk

% Increase figure resolution for export
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector graphics
print(gcf, 'figure_name', '-dpdf', '-r300'); % Save as PDF with high resolution

%% EXP ~ CUES from 49deg (medium pain) trials
%===========================================
Y_name = 'exp';
X_names = {'cues'};
X2_names = {'2nd-level Intercept (average within-person effects)'};

cuePain = glmfit_multilevel(ldata.(Y_name), ldata.cues, [], 'names', X_names, 'beta_names', X2_names, 'weighted')


create_figure('2nd level stats');
barplot_columns(cuePain.first_level.beta', 'names', cuePain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

% bar plot 
clf;
bar_handle = bar(conditions, mean_exp49);
bar_handle.FaceColor = dark_blue; % Light blue color
bar_handle.FaceAlpha = 0.8; % Transparency
ylim([0 80]);
hold on;

y_offset = 5; % Amount to offset the line and asterisk above the bar
y = max(mean_exp49) + y_offset; % Set y to be slightly above the tallest bar
hLegend = legend('Pain expecation');
line([x1, x2], [y, y], 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off'); % Horizontal line
text(mean([x1, x2]), y, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18, 'FontWeight', 'bold'); % Asterisk

% Increase figure resolution for export
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector graphics
set(gca, 'XTick', conditions); % Set the x-axis ticks to match the conditions
set(gca, 'XTickLabel', {'CS low', 'CS high'}, 'FontSize', 18); % Change the labels of the conditions

print(gcf, 'figure_name', '-dpdf', '-r300'); % Save as PDF with high resolution

hold off;

%% PAIN ~ TEMP
Y_name = 'pain';
X_var_names = {'temp'};
X2_names = {'2nd-level Intercept (average within-person effects)'};

cuePain = glmfit_multilevel(ldata.pain, ldata.temp, [], 'names', X_var_names, 'beta_names', X2_names, 'weighted', 'plots')
create_figure('2nd level stats');
barplot_columns(cuePain.first_level.beta', 'names', cuePain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

%%
%=====================
% GENERALIZATION TASK
%====================

plot_correlation_matrix(gdata_df, 'names', gdata_df.Properties.VariableNames)

%% PAIN ~ CUES (cscat) : N.S.
%============
Y_name = 'pain';
X_names = {'cscat'};
X2_names = {'2nd-level Intercept (average within-person effects)'};

genPain = glmfit_multilevel(gdata.pain, gdata.SCcscat, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot');
gen_CSbeta = genPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(genPain.first_level.beta', 'names', genPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Within-person cue effect on pain ratings');
drawnow, snapnow

create_figure('1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

% GENERALIZATION EFFECTS ~ Exp : Sig. relationship

% Create a scatter plot
figure;
scatter(learn_beta, gen_CSbeta, "filled");
xlabel('Explicit Learning (beta)', 'FontSize', 20);
ylabel('Generalization (beta)', 'FontSize', 20);
%title('Generalization effect ~ explicit learning', 'FontSize', 16);
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, gen_CSbeta);
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r =  %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = %.2f', p), 'FontSize', 12, 'EdgeColor', 'none');

% Learners only
Y_name = 'pain';
X_names = {'cscat'}
X2_names = {'2nd-level Intercept (average within-person effects)'};

genPain = glmfit_multilevel(gdata.pain(learners_idx), gdata.cscat(learners_idx), [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
gen_beta = genPain.first_level.beta(2,:)'; 
create_figure('2nd level stats');
barplot_columns(genPain.first_level.beta', 'names', genPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);

% ---------------
%% Pain ~ category with second level expectations : Sig.

X_names = {'cscat'};
X2 = [learn_beta];
 
genPain = glmfit_multilevel(gdata.pain, gdata.SCcscat, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')
gen_CSbetaMod = genPain.first_level.beta(2,:)';

create_figure('2nd level stats');
barplot_columns(genPain.Y_star, 'names', {'Intercept', 'cscat'}, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow


%% SAVE GENDATA WITH NEW BETAS
gdata.genBetas_CS = gen_CSbeta;
gdata.genBetas_CSmod = gen_CSbetaMod

save('finalPreproc_simGen_data.mat', 'gdata')


%% Multiline with bar plots for Pain ~ category model
% Supports second level results


% Colors gradient as a function of learning
min_beta = min(learn_beta);
max_beta = max(learn_beta);
color_scale = colormap(parula); % Using parula colormap for a smooth transition from light yellow to orange
normalized_beta = (learn_beta - min_beta) / (max_beta - min_beta);
color_indices = ceil(normalized_beta * (size(color_scale, 1) - 1)) + 1;
rgb_colors = cell(size(color_indices));
for i = 1:length(color_indices)
    rgb_colors{i} = color_scale(color_indices(i), :);
end


min_beta = min(learn_beta);
max_beta = max(learn_beta);


create_figure('gen1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain,'colors', rgb_colors, 'exclude_low_range_Y');
set(gca, 'XTick', [-.4 .4], 'XLim', [-.6 .6], 'XTickLabel', {'Related to CS low' 'Related to CS high'},'FontSize', 18);
hold on;
conditions = [-1, 1]; % Unique conditions
mean_painGen = zeros(1, length(conditions)); % Preallocate mean pain array
%filtered_49 = (gdata_df.temp == 49);
% Calculate mean pain scores for each condition
for i = 1:length(conditions)
    condition_indices = (gdata_df.cscat == conditions(i));
    mean_painGen(i) = nanmean(gdata_df.pain(condition_indices));
end
% Add bar plots of the mean pain scores
bar_handle = bar(conditions, mean_painGen);
bar_handle.FaceColor = dark_blue; % Light blue color
bar_handle.FaceAlpha = 0.7; % Transparency
ylim([0 80]);
xlabel('Generalization conditions', 'FontSize', 20);
ylabel('Pain ratings', 'FontSize', 20);
%title('Generalization cues (GC) effect on pain', 'FontSize', 22);
hLegend = legend('Individual slopes', 'FontSize', 17);
hold on;
p = patch(NaN, NaN, dark_blue, 'FaceAlpha', 0.5);
set(get(get(p, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'on');

colormap(parula);
cbar = colorbar('Location', 'eastoutside');
cbar.Ticks = linspace(0, 1, 10);
rounded_tick_labels = round(linspace(min_beta, max_beta, 10)); % Round the tick labels
cbar.TickLabels = num2cell(rounded_tick_labels);
%cbar.Label.String = 'Explicit learning (betas)';

line([x1, x2], [y, y], 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off'); % Horizontal line
text(mean([x1, x2]), y, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16); % Asterisk

hold off;

% Increase figure resolution for export
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector graphics
print(gcf, 'figure_name', '-dpdf', '-r300'); % Save as PDF with high resolution


%%


%% PAIN ~ Similarity : N.S.
%==============
Y_name = 'pain';
X_names = {'Empirical similarity'};

simPain = glmfit_multilevel(gdata.pain, gdata.SCempiricalSimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
sim_betas = simPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(simPain.first_level.beta', 'names', simPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Within-person similarity on pain ratings');
drawnow, snapnow

figure; % Sig correlation : high leverage points?
scatter(learn_beta, sim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning ~ similarity betas');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, sim_betas);
text(0.5, max(sim_betas) - 0.1, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, max(sim_betas) - 0.11, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% Remove high leverage points and plot correlation
out_thresh = (mean(learn_beta) + 3 *std(learn_beta)); % 3 SD > mean
disp(['Thresh used to remove betas > 3 S.D. is : ' num2str(out_thresh)])
inBetas = find(learn_beta <= out_thresh);
filt_learnBetas = learn_beta(inBetas);
filt_simBetas = sim_betas(inBetas);

% Then plot corr
figure; 
scatter(filt_learnBetas, filt_simBetas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning ~ similarity betas');
hold on;
lsline;
hold off;
[r, p] = corr(filt_learnBetas, filt_simBetas);
text(0.5, max(filt_simBetas) - 0.1, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, max(filt_simBetas) - 0.11, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

%% Pain ~ Empirical similarity with second level expectation : N.S. !! (weird? since EmpSim ~ category is high so should obtain sim results)
X_names = {'Empirical similarity'};
X2 = [gdata.learnBetas];
 
genPain = glmfit_multilevel(gdata.SCpain, gdata.SCempiricalSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')
gen_beta = genPain.first_level.beta(2,:)';

%% Pain ~ Mean model similarity : N.S but sig. interaction/association with learn betas
Y_name = 'pain';
X_names = {'mean_modelSimr'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain, gdata.SCmean_modelSimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
meanSim_betas = simPain.first_level.beta(2,:)'; 

figure; % Correl to learning betas : Sig correl
scatter(learn_beta, meanSim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning ~ similarity betas');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, meanSim_betas);
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');

%% Pain ~ Mean model similarity with second level expectation : N.S. !! 
X_names = {'Mean model similarity'};
X2 = [learn_beta];

genPain = glmfit_multilevel(gdata.pain, gdata.SCmean_modelSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')

%!! Sig. when scaling Y : SCpain vs pain
X_names = {'Mean model similarity'};
genPain = glmfit_multilevel(gdata.SCpain, gdata.SCmean_modelSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')


%% Pain ~ boulder LSA similarity : N.S.
Y_name = 'pain';
X_names = {'boulderLSASimr'} %,'modality'};
genPain = glmfit_multilevel(gdata.pain, gdata.SCboulderLSASimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
LSASim_betas = genPain.first_level.beta(2,:)'; 

figure; % Correl to learning betas : Sig correl
scatter(learn_beta, LSASim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning ~ LSA similarity betas');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, LSASim_betas);
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');

%% Pain ~ boulder LSA similarity with second level learning betas : N.S. 
X_names = {'Boulder LSA conceptual similarity'};
X2 = [learn_beta];

genPain = glmfit_multilevel(gdata.pain, gdata.SCboulderLSASimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')

% !! Scaling Y (pain) : N.S. trending results p = .07
genPain = glmfit_multilevel(gdata.SCpain, gdata.SCboulderLSASimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')

% Learners only!
%==============
Y_name = 'pain';
X_names = {'empirical_sim'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain(learners_idx), gdata.empiricalSimr(learners_idx), [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
sim_betas = simPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(simPain.first_level.beta', 'names', simPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);


%% PAIN ~ CUES (cscat) + simR from learners only
%=======================================
S
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});



% PAIN ~ Similarity and second level expectation : Not sig!
X_names = {'empiricalSimr'}; %,'modality'};
X2_names = {'Learning beta'};
X2 = learn_beta
genPain = glmfit_multilevel(gdata.pain, gdata.boulderSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plots')
gen_beta = genPain.first_level.beta(2,:)'; 

%% PAIN ~  CSCAT + simmilarity (trials)
%=====================
    Y_name = 'pain';
    X_names = {'cscat', 'mean_modelSimr'}; %theotrialSim
    for sub = 1:numel(subjects);
        v1 = gdata.(X_names{1}){sub};
        v2 = gdata.(X_names{2}){sub};
        X_cat_trials{sub} = horzcat(v1, v2);
    end
    X2_names = {'2d level'};
    X2 = []
    genPain = glmfit_multilevel(gdata.pain, X_cat_trials, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
gen_beta = genPain.first_level.beta(2,:)'; 
sim_betas = genPain.first_level.beta(3,:)';
% scatter plot cue betas
figure;
scatter(learn_beta, gen_beta);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and generalization betas, controlling for similarity');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, gen_beta);
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% scatter plot sim betas
figure;
scatter(learn_beta, sim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and similarity betas, controlling for CS effect');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, sim_betas);
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% Second level expectatation
X2 = [learn_beta]
genPain = glmfit_multilevel(gdata.pain, X_cat_trials, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')


plot_correlation_matrix(genPain.Y_star, 'names', genPain.inputOptions.names);
drawnow, snapnow
figure; plotmatrix(genPain.Y_star);
drawnow, snapnow
% GENERALIZATION EFFECTS ~ Exp
% Create a scatter plot
figure;
scatter(learn_beta, gen_beta);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Scatter Plot');
% Add regression line
hold on;
lsline;
hold off;
% Calculate correlation coefficient and p-value
[r, p] = corr(learn_beta, gen_beta)

% Display correlation coefficient and p-value on the plot
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

%% Second level learning
% Second level expectation on Pain~sim
X_names = {'empiricalSimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.empiricalSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')

X_names = {'mean_modelSimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.mean_modelSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')

X_names = {'boulderLSASimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.boulderLSASimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')


%% Correl between var
colnames = {'blocktrial', 'condition','cscat', 'pain','empiricalSimr', 'mean_modelSimr', 'boulderLSASimr'}
plot_correlation_matrix(gdata_df(:,colnames), 'names', colnames);
drawnow, snapnow

%% Mediation models
M = gdata.empiricalSimr

[paths, stats] = mediation(M,gdata.pain,gdata.cscat, 'plots','boot', 'verbose', 'bootsamples', 10000);
mediation_path_diagram(stats)

[paths, stats] = mediation(gdata.cscat(learners),gdata.pain(learners),M(learners) , 'plots','boot', 'verbose', 'bootsamples', 10000);
mediation_path_diagram(stats)

[r,p] = corr(gdata_df.cscat, gdata_df.empiricalSimr)
[r,p] = corr(gdata_df.cscat, gdata_df.mean_modelSimr)
[r,p] = corr(gdata_df.cscat, gdata_df.boulderLSASimr)
%correl CS and sim model?
%mediation

% /home/dylan.sutterlinguindon/genPain/behav_scripts/LKscripts/analyze_main_SCEBL_fMRI.m
% /home/dylan.sutterlinguindon/genPain/results 
% /crnldata/socialhealth/projects/2024_PainGen/Behavioral/DATA

%% Export results in table
ldata.pain49deg, ldata.cues49deg


result_table = array2table([nps{1},siips{1}], 'VariableNames', {'nps_ab', 'siips_ab'});

%% Export variables to CSV for python scripts (modif DSG 20 nov.)

learn_beta           = learn_beta(:);
pain_modulation_beta = pain_modulation_beta(:);
gen_CSbeta           = gen_CSbeta(:);

% Combine
data = [learn_beta, pain_modulation_beta, gen_CSbeta, ];
headers = {'exp_learning_betas', 'learned_pain_betas', 'gen_betas'};
T = array2table(data, 'VariableNames', headers);

writetable(T, 'final-table_within_person_effects_learn_gen.csv');
