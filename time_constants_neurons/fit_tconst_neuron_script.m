kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/sep';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Dexp_fits';

%% Criteria
NPSP_th = 40;
n_cores = 20;
%% Load the data
load(fullfile(kernel_dir,'info'),'info');

%% Select the data
ix = info.NPSP<NPSP_th;
NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);

%Sort in increasing NPSP
[~,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
qualities = qualities(ix_select);

%% Fit the time constants
n_clust = length(ix_select);
kh_fits = cell(n_clust,1);
problematic_cluster = [];

delete(gcp('nocreate'));
parpool('local',n_cores);

fprintf('== Getting dexp fits for all neurons ==\n');tic;
parfor c = 1:n_clust
    fprintf('== Cluster %0.f/%0.f ==\n',c,n_clust);
    try 
    kernel_name = fullfile(kernel_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_'));
    kh_fits{c} = fit_tconst_neuron(kernel_name);
    catch
        kh_fits{c} = NaN;
        problematic_cluster = [problematic_cluster;c] 
    end
        
end
fprintf('== Done! This took %0.fs ==\n',toc);
%% Save the results
sprintf('== Saving the results ==\n');tic;
[~,model] = fileparts(kernel_dir);
save_dir_full = fullfile(save_dir,model);
if ~exist(save_dir_full, 'dir')
    mkdir(save_dir_full);
end

for c = 1:n_clust
    kh_fit = kh_fits{c};
    save(fullfile(save_dir_full,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'kh_fit');
end

clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.quality = qualities;
info.pen_name = pen_names;
save(fullfile(save_dir_full,'info'),'info');
fprintf('== Done! Saving took %0.fs ==\n',toc);

%% Plot the betas

%% Plot the taus