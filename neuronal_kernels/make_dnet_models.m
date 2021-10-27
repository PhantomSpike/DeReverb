n_batch = 4;
lambda_name = 'lambda_ind';
coch_type = 'spechill';
h_max_ms = 30;
dt_ms = 5;
n_hu = 20;
n_pass = 20;
r_type{1} = 'big';
r_type{2} = 'small';
r_type{3} = 'anech';
n_rooms = length(r_type);
%% Input params
%NPSP selection criteria
NPSP_th = 40;
model_type = 'dnet';

if ~exist('n_cores','var') || isempty(n_cores)
    n_cores = 20;
end

cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Directory with all the clusters
dnet_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/DNet_fits';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits';
%% Select the clusters to fit
load(fullfile(dnet_dir,'noise_ratios'),'NR','neuronlist');
%% Select the data
%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NR,'ascend');
n_clust = length(ix_select);

for k = 1:n_clust
    ix = ix_select(k);
    temp = split(neuronlist{ix},'_');
    animal_names{k} = temp{1};
    pen_names{k} = temp{2};
    cluster_ids(k) = str2num(temp{3}(1:end-4));
end

%Get the lambdas for every room
for r = 1:n_rooms
    room = r_type{r};
    load(fullfile(dnet_dir,['run_dn_mlp_Reverb_Datacond_',room,'_30ms']),'model');
    for k = 1:n_clust
        ix = ix_select(k);
        lambda(k).(room) = model(ix).net_str_optim{4};
    end
end

delete(gcp('nocreate'));
parpool('local',n_cores);
%% Fit the models
params.h_max_ms = h_max_ms;
params.dt_ms = dt_ms;
params.n_hu = n_hu;
params.n_pass = n_pass;
params.model = model_type;

coch_file1 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Ronnie_PLP/',coch_type),[num2str(dt_ms),'ms/coch_all_conditions_',coch_type]);
coch_file2 = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Derry_Kilkenny_Cork/',coch_type),[num2str(dt_ms),'ms/coch_all_conditions_',coch_type]);
coch1 = load(coch_file1,'coch');
coch2 = load(coch_file2,'coch');

params.n_batch = n_batch;
params.lambda = lambda;
kernels = cell(n_clust,1);

fprintf('== Fitting DNet model ==\n');tic;
parfor k = 1:n_clust
    fprintf('== Processing cluster %0.f/%0.f ==\n',k,n_clust);
    y_temp = load(fullfile(cluster_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_')),'data');
    if ismember(y_temp.data.params.animal_name,["Ronnie","PLP"])
        kernels{k} =  run_dnet(coch1,y_temp.data,params,k);
    elseif ismember(y_temp.data.params.animal_name,["Cork","Kilkenny","Derry"])
        kernels{k} =  run_dnet(coch2,y_temp.data,params,k);
    else
        error('Unrecognized animal');
    end
    kernels{k}.params = y_temp.data.params;
end
fprintf('== Done! Fitting all kernels took %0.fs ==\n',toc);

sprintf('== Saving the results ==\n');tic;

model_dir = fullfile(save_dir,model_type);
if ~exist(model_dir, 'dir')
    mkdir(model_dir);
end

coch_dir = fullfile(model_dir,coch_type);
if ~exist(coch_dir, 'dir')
    mkdir(coch_dir);
end

lambda_dir = fullfile(coch_dir,lambda_name);
if ~exist(lambda_dir, 'dir')
    mkdir(lambda_dir);
end

save_dir_full = fullfile(lambda_dir,['batch_',num2str(n_batch)]);
if ~exist(save_dir_full, 'dir')
    mkdir(save_dir_full);
end

for k = 1:n_clust
    kernel = kernels{k};
    save(fullfile(save_dir_full,strjoin({kernel.params.animal_name,kernel.params.pen_name,num2str(kernel.params.cluster_id)},'_')),'kernel');
end

clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.pen_name = pen_names;
save(fullfile(save_dir_full,'info'),'info');
fprintf('== Done! Saving took %0.fs ==\n',toc);
