function make_lnp_model(model,h_max_ms,dt_ms,normz,save_dir,kfold,n_cores)
%% Input params
%NPSP selection criteria
NPSP_th = 40;
coch_type = 'specpower';

if ~exist('model','var') || isempty(model)
    model = 'lasso';
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
%% Define the params for all possible models
coch_file1 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Ronnie_PLP/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file2 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file3 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Noah_Derekah/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Directory with all the clusters

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
coch3 = load(coch_file3,'coch');

kernels = cell(n_clusters,1);

sprintf('== Fitting %s kernels ==\n',model);tic;

parfor c = 1:n_clusters
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
    temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    if ismember(temp.data.params.animal_name,["Ronnie","PLP"])
        kernels{c} =  gen_lnp_kernel(coch1,temp.data,h_max_ms,model,normz,kfold);
    elseif ismember(temp.data.params.animal_name,["Cork","Kilkenny","Derry"])
        kernels{c} =  gen_lnp_kernel(coch2,temp.data,h_max_ms,model,normz,kfold);
    elseif ismember(temp.data.params.animal_name,["Noah","Derekah"])
        kernels{c} =  gen_lnp_kernel2(coch3,temp.data,h_max_ms,model,normz,kfold);
    else
        error('Unrecognized animal');
    end
    kernels{c}.params = temp.data.params;
end

fprintf('== Done! Fitting all kernels took %0.fs ==\n',toc);

sprintf('== Saving the results ==\n');tic;

lnp_dir = fullfile(save_dir,'LNP_model_noneuron_norm');
if ~exist(lnp_dir, 'dir')
    mkdir(lnp_dir);
end

norm_dir = fullfile(lnp_dir,normz);
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

save_dir_full = fullfile(bin_dir,[num2str(h_max_ms,'%1.0f'),'ms']);
if ~exist(save_dir_full, 'dir')
    mkdir(save_dir_full);
end


for c = 1:n_clusters
    kernel = kernels{c};
    save(fullfile(save_dir_full,strjoin({kernel.params.animal_name,kernel.params.pen_name,num2str(kernel.params.cluster_id)},'_')),'kernel');
end

clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.quality = qualities;
info.pen_name = pen_names;
save(fullfile(save_dir_full,'info'),'info');
fprintf('== Done! Saving took %0.fs ==\n',toc);