%% Define params
coch_type = 'specpower';
model = 'ridge';
normz = 'perfreq_noneuro';
h_max_ms = 200;
dt_ms = 10;
kfold = true;
coch_file1 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Ronnie_PLP/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file2 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
coch_file3 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Noah_Derekah/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms'; %Directory with all the clusters
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data';

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
coch3 = load(coch_file2,'coch');

%% Load the kernels and compute CCnorm
CCnorm_small = cell(n_clusters,1);
CCnorm_big = cell(n_clusters,1);
kernels = cell(n_clusters,1);

sprintf('== Computing CC_norm ==\n');tic;
parfor c = 1:n_clusters
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
    temp_psth = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    temp_ker = load(fullfile(kernel_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'kernel');
    animal_name = animal_names{c};
    
    if ismember(temp_psth.data.params.animal_name,["Ronnie", "PLP"])
        [CCnorm_small{c}, CCnorm_big{c}] =  get_CC_norm_gain(coch1, temp_psth.data, temp_ker.kernel, normz, animal_name);
        
    elseif ismember(temp_psth.data.params.animal_name,["Cork", "Kilkenny", "Derry"])
        [CCnorm_small{c}, CCnorm_big{c}] =  get_CC_norm_gain(coch2, temp_psth.data, temp_ker.kernel, normz, animal_name);
        
    elseif ismember(temp_psth.data.params.animal_name,["Noah", "Derekah"])
        [CCnorm_small{c}, CCnorm_big{c}] =  get_CC_norm_gain(coch3, temp_psth.data, temp_ker.kernel, normz, animal_name);
        
    else
        error('Unrecognized animal');
    end
    temp_ker.kernel.params.CCnorm_small = CCnorm_small{c};
    temp_ker.kernel.params.CCnorm_big = CCnorm_big{c};
    kernels{c} = temp_ker.kernel;
    
end
fprintf('== Done! Computing CC_norm for all kernels took %0.fs ==\n',toc);

%% Save the results
sprintf('== Saving the results ==\n');tic;
info.CCnorm.small = CCnorm_small;
info.CCnorm.big = CCnorm_big;
mean_cc_small = reach(cell2mat(info.CCnorm.small),'CCnorm_mean');
mean_cc_big = reach(cell2mat(info.CCnorm.big),'CCnorm_mean');
hasNaN_cc_small = reach(cell2mat(info.CCnorm.small),'hasNaN');
hasNaN_cc_big = reach(cell2mat(info.CCnorm.big),'hasNaN');

info.CCnorm_mean.small = mean_cc_small;
info.CCnorm_mean.big = mean_cc_big;
info.CCnorm_hasNaN.small = hasNaN_cc_small;
info.CCnorm_hasNaN.big = hasNaN_cc_big;
save(fullfile(kernel_dir,'info.mat'),'info'); 

for c = 1:n_clusters
    kernel = kernels{c};
    save(fullfile(kernel_dir, strjoin({kernel.params.animal_name,kernel.params.pen_name,num2str(kernel.params.cluster_id)},'_')),'kernel');
end
fprintf('== Done! Saving the results took %0.fs ==\n',toc);