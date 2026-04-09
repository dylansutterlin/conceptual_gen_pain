
clc;
close all;
clear all;


projectDir = '/home/dsutterlin/projects/genPain/';
toolboxdir = fullfile(projectDir, 'Toolboxes');
dataDir = fullfile(projectDir, 'DATA/Imaging/test_genBetas/model3_generalization_ST');
behavDir = fullfile(projectDir,'results','behavioral');
saveDir = fullfile(projectDir, 'results/imaging/mediation/testData');

%% cRNL cluster paths
projectDir = '/crnldata/socialhealth/projects/2024_PainGen/';
toolboxdir = fullfile(projectDir, 'Toolboxes');
dataDir = '/crnldata/socialhealth/data/2024_SCEBL_GEN/Imaging_Data/analyses/model3_generalization_ST';
behavDir = fullfile(projectDir,'results','behavioral');
saveDir = fullfile(projectDir, '/results/imaging/mediation/');

%% 
mkdir(saveDir);


% Add Toolboxes to paths
addpath(genpath(fullfile(toolboxdir, 'RobustToolbox')));
addpath(genpath(fullfile(toolboxdir, 'Neuroimaging_Pattern_Masks')));
addpath(genpath(fullfile(toolboxdir, 'CanlabCore')));
addpath(genpath(fullfile(toolboxdir, 'MediationToolbox')));
addpath(genpath(fullfile(toolboxdir, 'spm12')));

% Except SPM/external for compatibility with Canlab
% toolbox/mediation (?)
spmFolders = dir(fullfile(toolboxdir, 'spm12'));
for i = 1:numel(spmFolders)
    folderName = spmFolders(i).name;
    if strcmp(folderName, 'external')
        path2rm = strsplit(genpath(fullfile(toolboxdir, 'spm12', folderName)), ':');
        for i=1:numel(path2rm)
            rmpath(char(path2rm(i)));
        end
    else
        continue
    end
end

%% --------------
% Load varialbes
% --------------
load(fullfile(behavDir, 'finalPreproc_simGen_data.mat'));
gendata = gdata;


load(fullfile(behavDir, 'SCEBLmri_learndata_FINAL_N36.mat'))
learndata = SCEBLmri_data;

load(fullfile(projectDir, '/results/imaging/STmodel3_cues_vifsQC/subjectsVifs.mat'));
vifStruct = subjects;

% Functionnal images
allFiles = filenames(fullfile(dataDir, 'SCEBL*'));
subjectsDir = [];
% remove .zip folder from files to extract func files
for i = 1:numel(allFiles)
    filename = allFiles(i);
    if ~endsWith(filename, '.zip')
        % Add the filename to the subjects cell array
        subjectsDir = [subjectsDir; filename];
    end
end

% Extract subjects names
subjects = cell(1, numel(subjectsDir));
for p = 1:numel(subjectsDir)
    [~, subname] = fileparts(subjectsDir(p));
    subjects{p} = subname;
end 

% Initialize an empty cell array one for each subject
func_imgs = cell(1, numel(subjects)); % Initialize func_imgs as a cell array

% get func files in a cell/sub
for sub = 1:numel(subjects)
    contrastFiles = {}; % Initialize contrastFiles as an empty cell array
    subFiles = fullfile(dataDir, char(subjects{sub}));
    
    SPM = load(fullfile(subFiles, 'SPM.mat')); % Load SPM.mat for the subject
    table = struct2table(SPM.SPM.Vbeta);

    % Loop through each contrast
    for cont = 1:numel(table.descrip)
        % Check if the description of the contrast contains 'allPains_trials0'
        if contains(table.descrip{cont}, 'AllPains_trial0')
            % add the image path to contrastFiles
            add_img = table.fname{cont};
            img_path = fullfile(dataDir, subjects{sub}, add_img);
            contrastFiles = [contrastFiles, img_path]; % Append image path to contrastFiles
        end
    end
    %func_imgs{sub} = contrastFiles;
    func_imgs{sub} = char(contrastFiles'); % Assign contrastFiles to the corresponding subject in func_imgs
end

%% Prep varialbes, scale and compute 2d level moderator
% X
cscat = gendata.cscat;
empSim = gendata.empiricalSimr;
meanSim = gendata.mean_modelSimr;
lsaSim = gendata.boulderLSASimr;
pain = gendata.pain; % Y
cov = gendata.blocktrial;
data = struct();
data.subjects = subjects;
data.func = func_imgs; % M
data.pain = pain;

% Mean center vairables
for i=1:numel(func_imgs) % mean center each cells
    x = cscat{i};
    data.cscat{i} = (x-mean(x))/std(x);

    x = empSim{i};
    data.empSim{i} = (x-mean(x))/std(x);

    x = meanSim{i};
    data.meanSim{i} = (x-mean(x))/std(x);

    x = lsaSim{i};
    data.lsaSim{i} = (x-mean(x))/std(x);

    x = cov{i};
    data.cov{i} = (x-mean(x))/std(x);

    %y=pain{i};
    %data.pain{i} = (y-mean(y))/std(y);
end 

% L2M : Expectation (during learning) as moderator (expectation ~ learning cues)
cue_exp_eff = glmfit_multilevel(learndata.exp, learndata.cues, [], 'weighted');
learn_beta = cue_exp_eff.first_level.beta(2,:)';
l2m = cell(1, numel(learn_beta));
for m = 1:numel(learn_beta)
    l2m{m} = learn_beta(m);
end
% mean center
data.learnBetas = (learn_beta-mean(learn_beta))/std(learn_beta);

%% Remove outlier from X M and Y based on VIF   

for sub = 1:length(data.subjects)
    mask = vifStruct.outvifs{sub}; % Vector of 0 (keep) or 1 (rm)
    keep_mask = mask == 0; % keep ones and leaves the 0; (inverse of mask)

    data.cscat{sub} = data.cscat{sub}(keep_mask);
    data.empSim{sub} = data.empSim{sub}(keep_mask);
    data.meanSim{sub} = data.meanSim{sub}(keep_mask);
    data.lsaSim{sub} = data.lsaSim{sub}(keep_mask);
    data.pain{sub} = data.pain{sub}(keep_mask);
    data.cov{sub} = data.cov{sub}(keep_mask);
    data.func_vifsCleaned{sub} = data.func{sub}(keep_mask, :);
    
end

%% Apply signatures on ALL data (data.func)

% fetch relevant signatures
%[pine] = load_image_set('pines')
NPS_folder_path = fullfile(projectDir, '/NPS_wager_lab/NPS_share')
NPS_path = filenames(fullfile(NPS_folder_path, 'weights_NSF_grouppred_cvpcr.img'))

%add apply_NPS func to path
addpath(fullfile(NPS_folder_path))
NPS_sig = fmri_data(NPS_path);
%orthviews(NPS_sig)

%load func data as fmri obj 1 cell/subj
subj_func = cell(1, length(data.func)); 
for i = 1:length(data.subjects)
    % Load the fMRI data for each subject using fmri_data and store it in subj_func
    subj_func{i} = fmri_data(data.func{i});
end

% apply all sugnature
[nps_dot, nps_names]  = apply_nps(data.func)
[siips_dot, nps_names, data_objects, siipspos_exp_by_region, siipsneg_exp_by_region, siipspos, siipsneg] = apply_siips(data.func)

dot_path = fullfile(projectDir, 'results/behavioral/dot_betaMaps')
save(fullfile(dot_path, "nps_dot.mat"), "nps_dot");
save(fullfile('siips_dot.mat'), 'siips_dot')

%% checkpoint on gdata, save func files, patterns and func paths

gdata.func_painBetas = data.func;
gdata.vifs_mask = vifStruct.outvifs;
gdata.nps_dot = nps_dot;
gdata.siips_dot = siips_dot;

gdata.vardescrip.func_painBetas = 'fMRI paths for pain betas per subject';
gdata.vardescrip.vifs2_5_mask = 'VIFs (variance inflation factors) for each trial map, keep 0, exclude 1s';
gdata.vardescrip.nps_dot = 'Dot-product of NPS (Neurologic Pain Signature) weights with subject data';
gdata.vardescrip.siips_dot = 'Dot-product of SIIPS (Stimulus-Independent Interoception-Pain Signature) weights with subject data';

save_behav = fullfile(projectDir, 'results/behavioral')
save(fullfile(save_behav, 'finalPreproc_simGen_data.mat'), 'gdata')
%% Mediation with dot results
nboot = 1000;

[paths, stats_cueNpsPain] = mediation(data.cscat,data.pain,nps_dot, 'plots','boot', 'verbose', 'bootsamples', nboot);
mediation_path_diagram(stats_cueNpsPain)

[paths, stats_simNpsPain] = mediation(gdata.empiricalSimr,data.pain,nps_dot, 'plots','boot', 'verbose', 'bootsamples', nboot);
mediation_path_diagram(stats_simNpsPain)

%siips
[paths, stats_cueSiipsPain] = mediation(data.cscat,data.pain,siips_dot, 'plots','boot', 'verbose', 'bootsamples', nboot);
mediation_path_diagram(stats_cueSiipsPain)


% Only learners
[paths, stats] = mediation(data.cscat(learners),data.pain(learners),nps_dot(learners) , 'plots','boot', 'verbose', 'bootsamples', nboot);
mediation_path_diagram(stats)

[paths, stats] = mediation(data.cscat(learners),data.pain(learners),nps_dot(learners) , 'plots','boot', 'verbose', 'bootsamples', nboot);
mediation_path_diagram(stats)

barplot_columns(nps_dot(learners))
barplot_columns(siips_dot(learners))

%% Mediation
SETUP.mask = which('gray_matter_mask.nii');
SETUP.preprocX = 0;
SETUP.preprocY = 0;
SETUP.preprocM = 0;

% cd(saveDir);
% Data used has VIFS > 2.5 removed!!
% data.func = data.func_vifsCleaned
% % Bootstrap with second level moderator
% new_folder = fullfile(saveDir, 'cue_B_pain');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.cscat, data.pain, data.func, SETUP, 'covs', data.cov, 'nopreproc','boot','bootsample',10000);
% % Bootstrap with second level moderator
% new_folder = fullfile(saveDir, 'cue_B_pain_L2Mlearn');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.cscat, data.pain, data.func, SETUP, 'covs', data.cov, 'L2M',data.learnBetas, 'nopreproc','boot','bootsample',10000);
% 
% % Empirical Similarity
% new_folder = fullfile(saveDir, 'empsim_B_pain');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.empSim, data.pain, data.func, SETUP, 'covs', data.cov, 'nopreproc','boot','bootsample',10000);
% % simil with second level moderator
% new_folder = fullfile(saveDir, 'empSim_B_pain_L2Mlearn');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.empSim, data.pain, data.func, SETUP, 'covs', data.cov, 'L2M',data.learnBetas, 'nopreproc','boot','bootsample',10000);
% 
% % theoretical Similarity (meanSim)
% new_folder = fullfile(saveDir, 'meanSim_B_pain');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.meanSim, data.pain, data.func, SETUP, 'covs', data.cov, 'nopreproc','boot','bootsample',10000);
% % meanSim with second level moderator
% new_folder = fullfile(saveDir, 'meanSim_B_pain_L2Mlearn');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.meanSim, data.pain, data.func, SETUP, 'covs', data.cov, 'L2M',data.learnBetas, 'nopreproc','boot','bootsample',10000);
% 
% % Boulder LSA Similarity (lsaSim)
% new_folder = fullfile(saveDir, 'lsaSim_B_pain');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.lsaSim, data.pain, data.func, SETUP, 'covs', data.cov, 'nopreproc','boot','bootsample',10000);
% % lsaSim with second level moderator
% new_folder = fullfile(saveDir, 'lsaSim_B_pain_L2Mlearn');
% mkdir(new_folder);
% cd(new_folder)
% mediation_brain_multilevel(data.lsaSim, data.pain, data.func, SETUP, 'covs', data.cov, 'L2M',data.learnBetas, 'nopreproc','boot','bootsample',10000);
% 
% 
% 
% disp('changed dir to : ')
% disp(pwd);
% %mediation_brain_multilevel(centmeanSim, Y, M, SETUP, 'covs', centCov, 'nopreproc' ,'boot', 'bootsample',10000, 'L2M', centLearn_beta);
% 
disp('Everything done!!');

