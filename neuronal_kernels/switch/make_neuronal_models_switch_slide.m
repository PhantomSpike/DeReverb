function make_neuronal_models_switch_slide(model,h_max_ms,dt_ms,normz,save_dir,kfold,n_cores)
%% Input params

% Define STRF window fit params
total_size_s = 10;
window_size_s = 2;
stride_s = 0.5;
n_steps = floor((total_size_s - window_size_s)/stride_s) + 1;

%NPSP selection criteria
NPSP_th = 40;
coch_type = 'specpower';

if ~exist('model','var') || isempty(model)
    model = 'ridge';
end

if ~exist('normz','var') || isempty(normz)
    normz = 'none';
end

if ~exist('h_max_ms','var') || isempty(h_max_ms)
    h_max_ms = 200;
end

if ~exist('dt_ms','var') || isempty(dt_ms)
    dt_ms = 10;
end

if ~exist('kfold','var') || isempty(n_cores)
    kfold = true;
end

if ~exist('n_cores','var') || isempty(n_cores)
    n_cores = 20;
end

params.h_max_ms = h_max_ms;
params.model = model;
params.normz = normz;
params.kfold = kfold;

%% Define the params for all possible models
coch_file1 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch/Switch_Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file2 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch/Switch_Noah_Derekah/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/For_analysis/All_data'; %Directory with all the clusters

%% Select the clusters to fit
load(fullfile(cluster_dir,'info.mat'),'info'); %Load the info file with the meta data
ix = info.NPSP<NPSP_th; %Find only the neurons that are below the NPSP th
%Get the necessary fields
cluster_ids = info.cluster_id(ix);
pen_names = info.pen_name(ix);
animal_names = info.animal_name(ix);
qualities = info.quality(ix);
NPSPs = info.NPSP(ix);
n_clusters = sum(ix);

%% Fit the model
delete(gcp('nocreate'));
parpool('local',n_cores);

coch1 = load(coch_file1,'coch');
coch2 = load(coch_file2,'coch');

kernels = cell(n_clusters,1);
sprintf('== Fitting %s kernels ==\n',model);tic;

% Define the periods of stimuli

%sl - small-to-large transition
%ls - large-to-small transition

% Define the indices for sl
sl_1_small = [8:stride_s:18-window_size_s];
sl_2_small = [24:stride_s:34-window_size_s];
sl_1_large = [16:stride_s:26-window_size_s];
sl_2_large = [32:stride_s:42-window_size_s];

% Define the indices for ls
ls_1_small = [16:stride_s:26-window_size_s];
ls_2_small = [32:stride_s:42-window_size_s];
ls_1_large = [8:stride_s:18-window_size_s];
ls_2_large = [24:stride_s:34-window_size_s];

for s = 1:n_steps
    
    %Initialize
    sl = []; ls = [];
    
    %Small start (R1)
    sl(1).small = [sl_1_small(s), sl_1_small(s) + window_size_s]; sl(2).small = [sl_2_small(s), sl_2_small(s) + window_size_s];
    ls(1).small = [ls_1_small(s), ls_1_small(s) + window_size_s]; ls(2).small = [ls_2_small(s), ls_2_small(s) + window_size_s];
    
    %Large start (R2)
    sl(1).large = [sl_1_large(s), sl_1_large(s) + window_size_s]; sl(2).large = [sl_2_large(s), sl_2_large(s) + window_size_s];
    ls(1).large = [ls_1_large(s), ls_1_large(s) + window_size_s]; ls(2).large = [ls_2_large(s), ls_2_large(s) + window_size_s];
    
    params.sl = sl;
    params.ls = ls;
    
    parfor c = 1:n_clusters
        fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
        temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
        if ismember(temp.data.params.animal_name,["Cork","Kilkenny","Derry"])
            kernels{c}{s} =  gen_neuronal_kernel_switch_slide(coch1,temp.data, params);
        elseif ismember(temp.data.params.animal_name,["Noah","Derekah"])
            kernels{c}{s} =  gen_neuronal_kernel_switch_slide(coch2,temp.data, params);
        else
            error('Unrecognized animal');
        end
        kernels{c}{s}.params = temp.data.params;
        kernels{c}{s}.params.sl = sl;
        kernels{c}{s}.params.ls = ls;
    end
    
end

fprintf('== Done! Fitting all kernels took %0.fs ==\n',toc);

sprintf('== Saving the results ==\n');tic;

norm_dir = fullfile(save_dir,normz);
if ~exist(norm_dir, 'dir')
    mkdir(norm_dir);
end

model_dir = fullfile(norm_dir,model);
if ~exist(model_dir, 'dir')
    mkdir(model_dir);
end

bin_dir = fullfile(model_dir,[num2str(dt_ms,'%1.0f'),'ms']);
if ~exist(bin_dir, 'dir')
    mkdir(bin_dir);
end

hmax_dir = fullfile(bin_dir,[num2str(h_max_ms,'%1.0f'),'ms']);
if ~exist(hmax_dir, 'dir')
    mkdir(hmax_dir);
end

window_dir = fullfile(hmax_dir,['window_',num2str(window_size_s,'%1.1f'),'_s']);
if ~exist(window_dir, 'dir')
    mkdir(window_dir);
end

save_dir_full = fullfile(window_dir,['stride_',num2str(stride_s,'%1.3f'),'_s']);
if ~exist(save_dir_full, 'dir')
    mkdir(save_dir_full);
end

for c = 1:n_clusters
    kernel = kernels{c};
    save(fullfile(save_dir_full,strjoin({kernel{1}.params.animal_name,kernel{1}.params.pen_name,num2str(kernel{1}.params.cluster_id)},'_')),'kernel');
end


clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.quality = qualities;
info.pen_name = pen_names;
save(fullfile(save_dir_full,'info'),'info');
fprintf('== Done! Saving took %0.fs ==\n',toc);