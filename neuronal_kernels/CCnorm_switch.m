%% Define params
coch_type = 'specpower';
model = 'ridge';
normz = 'perfreq_noneuro';
h_max_ms = 200;
dt_ms = 10;
kfold = true;
coch_file1 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch/Switch_Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file2 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch/Switch_Noah_Derekah/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/Kernel_fits/kfolds/0.25_4_4.25_8s_version/perfreq_noneuro/ridge/10ms/200ms'; %Directory with all the clusters
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/For_analysis/All_data';
%% Select the kernels to compute CCnorm for
load(fullfile(kernel_dir,'info.mat'),'info'); %Load the info file with the meta data

%Get the necessary fields
cluster_ids = info.cluster_id;
pen_names = info.pen_name;
animal_names = info.animal_name;
qualities = info.quality;
NPSPs = info.NPSP;
n_clusters = length(cluster_ids);

%% Load the cochleagrams
coch1 = load(coch_file1,'coch');
coch2 = load(coch_file2,'coch');

%% Load the kernels and compute CCnorm
CCnorm = cell(n_clusters,1);
kernels = cell(n_clusters,1);

sprintf('== Computing CC_norm ==\n');tic;
parfor c = 1:n_clusters
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
    temp_psth = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    temp_ker = load(fullfile(kernel_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'kernel');
    if ismember(temp_psth.data.params.animal_name,["Cork","Kilkenny","Derry"])
        CCnorm{c} =  get_CCnorm(coch1, temp_psth.data, temp_ker.kernel, normz);
    elseif ismember(temp_psth.data.params.animal_name,["Noah","Derekah"])
        CCnorm{c} =  get_CCnorm(coch2, temp_psth.data, temp_ker.kernel, normz);
    else
        error('Unrecognized animal');
    end
    temp_ker.kernel.params.CCnorm = CCnorm{c};
    kernels{c} = temp_ker.kernel;
end
fprintf('== Done! Computing CC_norm for all kernels took %0.fs ==\n',toc);

%% Save the results
sprintf('== Saving the results ==\n');tic;
info.CCnorm = CCnorm;
mean_cc = reach(cell2mat(info.CCnorm),'mean');
info.CCnorm_mean = mean_cc;
save(fullfile(kernel_dir,'info.mat'),'info'); 

for c = 1:n_clusters
    kernel = kernels{c};
    save(fullfile(kernel_dir, strjoin({kernel.params.animal_name,kernel.params.pen_name,num2str(kernel.params.cluster_id)},'_')),'kernel');
end
fprintf('== Done! Saving the results took %0.fs ==\n',toc);