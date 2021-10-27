%% Input variables

kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/Neuronal_data_orig/perfreq/ridge/10ms/200ms';
ei_info_path = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/E_I_together/clust_info.mat';
plot_kernels = 0;

get_neurons = 'good';
NPSP_th = 40;
freq_up_bound = 10000;
freq_down_bound = 700;

%% Params
chop_ms = 190; %The analysis window for the COM
hist_type = 'bar'; %The type of the histogram
calc_method = 'average'; %How to estimate impotant values
%calc_method -- 'raw', 'average'
%               'raw' - Use the whole receptive field 
%               'average' - Take an average across frequencies first

bf_neurons = 'shared'; %How to get the BF
% bf_method -- 'ind', 'shared'
%              'ind' - Treat the BFs separately 
%              'shared' - Take an log weighted mean BF between the two

bf_method = 'max'; %The method use to extract the bf:
% bf_method -- 'max', 'window_max', 'window_mean'
%              'max' - Take the max across all history
%              'window_max' - Take the max in a given time window
%              'window_mean' - Take the mean in a given time window
bf_window_ms = 50; %The window for calcualting bf if this method is used
exc_small_color = [0.9804    0.4196    0.6431];
exc_big_color = [0.8 0 0];
inh_small_color = [0.0745 0.6235 1.0000];
inh_big_color = 'b';
font_type = 'Liberation Sans';

if ~exist('NPSP_th','var') || isempty(NPSP_th)
    NPSP_th = 40;
end

if ~exist('get_neurons','var') || isempty(get_neurons)
    get_neurons = 'all';
end

% [save_dir,~] = fileparts(kernel_dir);
% plot_dir = fullfile(save_dir,'Plots');
% if ~exist(plot_dir,'dir')
%     mkdir(plot_dir);
% end
% 
% save_dir_full = fullfile(plot_dir,get_neurons);
% if ~exist(save_dir_full,'dir')
%     mkdir(save_dir_full);
% end
% 
% if ~strcmp(bf_method,'max')
%     save_dir_bf = fullfile(save_dir_full,[bf_neurons,'_',bf_method,'_',num2str(bf_window_ms),'ms']);
% else
%     save_dir_bf = fullfile(save_dir_full,[bf_neurons,'_',bf_method]);
% end
% 
% if ~exist(save_dir_bf,'dir')
%     mkdir(save_dir_bf);
% end

% r_type{1} = 'anech';
r_type{1} = 'small';
r_type{2} = 'big';
n_rooms = length(r_type);

%% Load the data
load(fullfile(kernel_dir,'info'),'info');
temp_files = dir([kernel_dir,'/*.mat']);
temp = load(fullfile(temp_files(1).folder,temp_files(1).name));
model = temp.kernel.model;
load(ei_info_path);

%% Select the data
switch get_neurons
    case 'all'
        ix_qualia = ones(length(info.cluster_id),1); %Get all the neurons
    case 'good'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'good'),info.quality,'UniformOutput',false)); %Find all the good units
    case 'mua'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'mua'),info.quality,'UniformOutput',false)); %Find all the mua units
end
ix_npsp = info.NPSP<NPSP_th; %Find all the neurons below certain NPSP
ix = ix_qualia & ix_npsp; %Find the intersection of the two

NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);

%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
qualities = qualities(ix_select);
n_clust = length(cluster_ids);

%Asign E and I labels to the neurons based on the clustering
inh_flag = clust_info.inh_flag;

i = 0; e = 0;
for c = 1:n_clust
    ix_curr = cell2mat(cellfun(@(x) strcmp(x,animal_names{c}),clust_info.animal_name,'UniformOutput',false)) & ...
        cell2mat(cellfun(@(x) strcmp(x,pen_names{c}),clust_info.pen_name,'UniformOutput',false)) & ...
        cluster_ids(c)==clust_info.cluster_id;
    label_curr = clust_info.kmeans_label(ix_curr);
    if label_curr ~= inh_flag
        e = e+1;
        exc.animal_name{e} = animal_names{c};
        exc.pen_name{e} = pen_names{c};
        exc.cluster_id(e) = cluster_ids(c);
    elseif label_curr == inh_flag
        i = i+1;
        inh.animal_name{i} = animal_names{c};
        inh.pen_name{i} = pen_names{c};
        inh.cluster_id(i) = cluster_ids(c);
    end
end
n_exc = e;
n_inh = i;

%% First load all the selected kernels
fprintf('== Loading the data ==\n');tic;
for k = 1:n_exc
    c_name = fullfile(kernel_dir,strjoin({exc.animal_name{k},exc.pen_name{k},num2str(exc.cluster_id(k))},'_'));
    load(c_name,'kernel');
    kernels.exc{k,1}.kernel = kernel;
    k_fhn_exc(:,:,k).small = fliplr(kernel.small.main{end}.k_fh);
    k_fhn_exc(:,:,k).big = fliplr(kernel.big.main{end}.k_fh);
end

for k = 1:n_inh
    c_name = fullfile(kernel_dir,strjoin({inh.animal_name{k},inh.pen_name{k},num2str(inh.cluster_id(k))},'_'));
    load(c_name,'kernel');
    kernels.inh{k,1}.kernel = kernel;
    k_fhn_inh(:,:,k).small = fliplr(kernel.small.main{end}.k_fh);
    k_fhn_inh(:,:,k).big = fliplr(kernel.big.main{end}.k_fh);
end
fprintf('== Done! This took %0.fs ==\n',toc);

%Get the params out
freqs = kernels.exc{1}.kernel.freqs;
n_h = kernels.exc{1}.kernel.n_h;
dt_ms = round(kernels.exc{1}.kernel.dt_ms);
h_max_ms = kernels.exc{1}.kernel.h_max_ms;
%% Plot the kernels separately
[base_dir,~] = fileparts(ei_info_path);
strf_dir = fullfile(base_dir,'STRFs');
exc_dir = fullfile(strf_dir,'Excitatory');
inh_dir = fullfile(strf_dir,'Inhibitory');
if plot_kernels
    ker_per_plot = 5;
    
    
    if ~exist(strf_dir,'dir')
        mkdir(strf_dir);
    end
    
    %Plot the excitatory ones
    if ~exist(exc_dir,'dir')
        mkdir(exc_dir);
    end
    
    plot_neuronal_kernels_ei(kernels.exc,exc_dir,ker_per_plot)
    
    %Plot the inhibitory ones
    if ~exist(inh_dir,'dir')
        mkdir(inh_dir);
    end
    plot_neuronal_kernels_ei(kernels.inh,inh_dir,ker_per_plot)
end

%% Plot the mean kernels
if plot_kernels
    
    %Plot the excitatory
    params.save_dir = exc_dir;
    params.type = 'raw'; params.name = 'neurons';
    params.freqs = freqs; params.n_h = n_h; params.dt_ms = dt_ms; params.h_max_ms = h_max_ms;
    params.font_type = 'Liberation Sans';
    params.exc_small_color = exc_small_color;  params.exc_big_color = exc_big_color;
    params.inh_small_color = inh_small_color;  params.inh_big_color = inh_big_color;
    plot_mean_strf(k_fhn_exc,params)
    
    %Plot the inhibitory
    params.save_dir = inh_dir;
    plot_mean_strf(k_fhn_inh,params)
end
%% First compute the bf, com, pt, ph for every cluster
freqs = fliplr(kernels.exc{1}.kernel.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_h = kernel.n_h;
dt_ms = round(kernel.dt_ms);
chop_ix = round(chop_ms/dt_ms);
bf_window_ix = round(bf_window_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;
n_type{1} = 'exc';
n_type{2} = 'inh';
n_cl.exc = n_exc;
n_cl.inh = n_inh;
n_n_types = 2;
fprintf('== Calcuating com and bf ==\n');tic;
for n = 1:n_n_types
    for k = 1:n_cl.(n_type{n})
        sprintf('== Cluster %0.f/%0.f ==\n',k, n_cl.(n_type{n}));
        for r = 1:n_rooms
            room = r_type{r};
            switch model
                
                case {'ridge','lasso','elastic'}
                    k_fh = fliplr(kernels.(n_type{n}){k}.kernel.(room).main{end}.k_fh);
                    k_fh = k_fh(:,1:chop_ix);
                    k_fhn.(n_type{n})(:,:,k).(room) = k_fh;
                    %Get all +ve and -ve values separately
                    k_fh_pos = abs(max(k_fh,0));
                    k_fh_neg = abs(min(k_fh,0));
                    %Get the mean across freqeuncies
                    k_h_pos = mean(k_fh_pos);
                    k_h_neg = mean(k_fh_neg);
                    %Select which method to use
                    switch calc_method
                        case 'raw'
                            k_h_pos_temp = k_fh_pos;
                            k_h_neg_temp = k_fh_neg;
                        case 'average'
                            k_h_pos_temp = k_h_pos;
                            k_h_neg_temp = k_h_neg;
                    end
                    
                    %Peak Height (PH) values
                    ph_pos.(n_type{n})(k).(room) = max(k_h_pos_temp(:));
                    ph_neg.(n_type{n})(k).(room) = max(k_h_neg_temp(:));
                    %Total Excitation and Inhibition
                    total_pos.(n_type{n})(k).(room) = sum(k_h_pos_temp(:));
                    total_neg.(n_type{n})(k).(room) = sum(k_h_neg_temp(:));
                    %IE ratio
                    ie_ratio.(n_type{n})(k).(room) = sum(k_h_neg_temp(:))/sum(k_h_pos_temp(:));
                    %Peak Time (PT) values
                    [~,max_ix_exc] = max(k_h_pos_temp(:));
                    [~, max_ix_col_exc] = ind2sub(size(k_h_pos_temp),max_ix_exc);
                    pt_pos.(n_type{n})(k).(room) = h(max_ix_col_exc);
                    
                    [~,max_ix_inh] = max(k_h_neg_temp(:));
                    [~, max_ix_col_inh] = ind2sub(size(k_h_neg_temp),max_ix_inh);
                    pt_neg.(n_type{n})(k).(room) = h(max_ix_col_inh);
                    
                    switch bf_method
                        case 'max'
                            k_f = max(k_fh_pos,[],2); %Take the max across history steps
                        case 'window_max'
                            k_f = max(k_fh_pos(:,1:bf_window_ix),[],2); %Take the max in a specified window
                        case 'window_mean'
                            k_f = mean(k_fh_pos(:,1:bf_window_ix),2); %Take the mean in a specified window
                    end
                    
                    [~,ix] = max(k_f); %Find the index of the max frequency
                    bf.(n_type{n})(k).(room) = freqs(ix); %Find the corresponding frequency
            end
            
            %Get the center of mass (COM) values
            k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
            k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
            com_neg.(n_type{n})(k).(room) = (k_h_neg*h); %Compute a weighted sum of all values
            com_pos.(n_type{n})(k).(room) = (k_h_pos*h); %Compute a weighted sum of all values
        end
        bf_mean.(n_type{n})(k) = 2^(log2(bf.(n_type{n})(k).small*bf.(n_type{n})(k).big)/2);
        [~,ix_bf] = min(abs(bf_mean.(n_type{n})(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
        bf_closest.(n_type{n})(k) = freqs(ix_bf);
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Stats
for n = 1:n_n_types
    %COM measures
    com_small_pos.(n_type{n}) = reach(com_pos.(n_type{n}),'small');
    com_big_pos.(n_type{n}) = reach(com_pos.(n_type{n}),'big');
    com_small_neg.(n_type{n}) = reach(com_neg.(n_type{n}),'small');
    com_big_neg.(n_type{n}) = reach(com_neg.(n_type{n}),'big');
    
    %IE ratio
    ie_ratio_small.(n_type{n}) = reach(ie_ratio.(n_type{n}),'small');
    ie_ratio_big.(n_type{n}) = reach(ie_ratio.(n_type{n}),'big');
    
    %Total Inhibiton and excitation
    total_pos_small.(n_type{n}) = reach(total_pos.(n_type{n}),'small');
    total_pos_big.(n_type{n}) = reach(total_pos.(n_type{n}),'big');
    total_neg_small.(n_type{n}) = reach(total_neg.(n_type{n}),'small');
    total_neg_big.(n_type{n}) = reach(total_neg.(n_type{n}),'big');

    %Inhibition Max value
    ph_pos_small.(n_type{n}) = reach(ph_pos.(n_type{n}),'small');
    ph_pos_big.(n_type{n}) = reach(ph_pos.(n_type{n}),'big');
    ph_neg_small.(n_type{n}) = reach(ph_neg.(n_type{n}),'small');
    ph_neg_big.(n_type{n}) = reach(ph_neg.(n_type{n}),'big');

    %Time of Max Inhibition
    pt_pos_small.(n_type{n}) = reach(pt_pos.(n_type{n}),'small');
    pt_pos_big.(n_type{n}) = reach(pt_pos.(n_type{n}),'big');
    pt_neg_small.(n_type{n}) = reach(pt_neg.(n_type{n}),'small');
    pt_neg_big.(n_type{n}) = reach(pt_neg.(n_type{n}),'big');
    
    %COM p-vals
    [p_val_com_pos.(n_type{n}),~, stats_com_pos.(n_type{n})] = signrank(com_big_pos.(n_type{n}),com_small_pos.(n_type{n})); 
    [p_val_com_neg.(n_type{n}),~, stats_com_neg.(n_type{n})] = signrank(com_big_neg.(n_type{n}),com_small_neg.(n_type{n}));
    
    %PH p-vals
    [p_val_ph_pos.(n_type{n}),~, stats_ph_pos.(n_type{n})] = signrank(ph_pos_big.(n_type{n}), ph_pos_small.(n_type{n}));
    [p_val_ph_neg.(n_type{n}),~,stats_ph_neg.(n_type{n})] = signrank(ph_neg_big.(n_type{n}), ph_neg_small.(n_type{n}));
    
    %PT p-vals
    [p_val_pt_pos.(n_type{n}),~, stats_pt_pos.(n_type{n})] = signrank(pt_pos_big.(n_type{n}), pt_pos_small.(n_type{n}));
    [p_val_pt_neg.(n_type{n}),~, stats_pt_neg.(n_type{n})] = signrank(pt_neg_big.(n_type{n}), pt_neg_small.(n_type{n}));
    
    %IE p-val
    [p_val_ie.(n_type{n}),~, stats_ie.(n_type{n})] = signrank(ie_ratio_big.(n_type{n}), ie_ratio_small.(n_type{n}));
    
    %Total p-val
    [p_val_total_exc.(n_type{n}),~, stats_total_exc.(n_type{n})] = signrank(total_pos_big.(n_type{n}),total_pos_small.(n_type{n}));
    [p_val_total_inh.(n_type{n}),~, stats_total_inh.(n_type{n})] = signrank(total_neg_big.(n_type{n}), total_neg_small.(n_type{n}));
   
    
%     com_small_bf_neg.(n_type{n}) = reach(com_neg_bf.(n_type{n}),'small');
%     com_small_bf_pos.(n_type{n}) = reach(com_pos_bf.(n_type{n}),'small');
%     com_big_bf_neg.(n_type{n}) = reach(com_neg_bf.(n_type{n}),'big');
%     com_big_bf_pos.(n_type{n}) = reach(com_pos_bf.(n_type{n}),'big');
end
%COM +ve
meta.p_com_pos_exc = p_val_com_pos.exc; meta.p_com_pos_inh = p_val_com_pos.inh;
meta.z_com_pos_exc = stats_com_pos.exc.zval; meta.z_com_pos_inh = stats_com_pos.inh.zval;
meta.med_com_pos_exc = nanmedian(com_big_pos.exc - com_small_pos.exc); meta.med_com_pos_inh = nanmedian(com_big_pos.inh - com_small_pos.inh);

%COM -ve
meta.p_com_neg_exc = p_val_com_neg.exc; meta.p_com_neg_inh = p_val_com_neg.inh;
meta.z_com_neg_exc = stats_com_neg.exc.zval; meta.z_com_neg_inh = stats_com_neg.inh.zval;
meta.med_com_neg_exc = nanmedian(com_big_neg.exc - com_small_neg.exc); meta.med_com_neg_inh = nanmedian(com_big_neg.inh - com_small_neg.inh);

%PT +ve
meta.p_pt_pos_exc = p_val_pt_pos.exc; meta.p_pt_pos_inh = p_val_pt_pos.inh;
meta.z_pt_pos_exc = stats_pt_pos.exc.zval; meta.z_pt_pos_inh = stats_pt_pos.inh.zval;
meta.med_pt_pos_exc = nanmedian(pt_pos_big.exc - pt_pos_small.exc); meta.med_pt_pos_inh = nanmedian(pt_pos_big.inh - pt_pos_small.inh);

%PT -ve
meta.p_pt_neg_exc = p_val_pt_neg.exc; meta.p_pt_neg_inh = p_val_pt_neg.inh;
meta.z_pt_neg_exc = stats_pt_neg.exc.zval; meta.z_pt_neg_inh = stats_pt_neg.inh.zval;
meta.med_pt_neg_exc = nanmedian(pt_neg_big.exc - pt_neg_small.exc); meta.med_pt_neg_inh = nanmedian(pt_neg_big.inh - pt_neg_small.inh);

%PH +ve
meta.p_ph_pos_exc = p_val_ph_pos.exc; meta.p_ph_pos_inh = p_val_ph_pos.inh;
meta.z_ph_pos_exc = stats_ph_pos.exc.zval; meta.z_ph_pos_inh = stats_ph_pos.inh.zval;
meta.med_ph_pos_exc = nanmedian(log2(ph_pos_big.exc./ph_pos_small.exc)); meta.med_ph_pos_inh = nanmedian(log2(ph_pos_big.inh./ph_pos_small.inh));

%PH -ve
meta.p_ph_neg_exc = p_val_ph_neg.exc; meta.p_ph_neg_inh = p_val_ph_neg.inh;
meta.z_ph_neg_exc = stats_ph_neg.exc.zval; meta.z_ph_neg_inh = stats_ph_neg.inh.zval;
meta.med_ph_neg_exc = nanmedian(log2(ph_neg_big.exc./ph_neg_small.exc)); meta.med_ph_neg_inh = nanmedian(log2(ph_neg_big.inh./ph_neg_small.inh));

%% More stats

%COM
[meta.p_com_pos_small,~, stats_com_pos_small] = ranksum(com_small_pos.exc, com_small_pos.inh); meta.z_com_pos_small = stats_com_pos_small.zval;
[meta.p_com_pos_big,~, stats_com_pos_big] = ranksum(com_big_pos.exc, com_big_pos.inh); meta.z_com_pos_big = stats_com_pos_big.zval;
[meta.p_com_neg_small,~, stats_com_neg_small] = ranksum(com_small_neg.exc, com_small_neg.inh); meta.z_com_neg_small = stats_com_neg_small.zval;
[meta.p_com_neg_big,~, stats_com_neg_big] = ranksum(com_big_neg.exc, com_big_neg.inh); meta.z_com_neg_small = stats_com_neg_small.zval;

%PT
[meta.p_pt_pos_small,~, stats_pt_pos_small] = ranksum(pt_pos_small.exc, pt_pos_small.inh); meta.z_pt_pos_small = stats_pt_pos_small.zval;
[meta.p_pt_pos_big,~, stats_pt_pos_big] = ranksum(pt_pos_big.exc, pt_pos_big.inh); meta.z_pt_pos_big = stats_pt_pos_big.zval;
[meta.p_pt_neg_small,~, stats_pt_neg_small] = ranksum(pt_neg_small.exc, pt_neg_small.inh); meta.z_pt_neg_small = stats_pt_neg_small.zval;
[meta.p_pt_neg_big,~, stats_pt_neg_big] = ranksum(pt_neg_big.exc, pt_neg_big.inh); meta.z_pt_neg_big = stats_pt_neg_big.zval;

%PH
[meta.p_ph_pos_small,~, stats_ph_pos_small] = ranksum(ph_pos_small.exc, ph_pos_small.inh); meta.z_ph_pos_small = stats_ph_pos_small.zval;
[meta.p_ph_pos_big,~, stats_ph_pos_big] = ranksum(ph_pos_big.exc, ph_pos_big.inh); meta.z_ph_pos_big = stats_ph_pos_big.zval;
[meta.p_ph_neg_small,~, stats_ph_neg_small] = ranksum(ph_neg_small.exc, ph_neg_small.inh); meta.z_ph_neg_small = stats_ph_neg_small.zval;
[meta.p_ph_neg_big,~, stats_ph_neg_big] = ranksum(ph_neg_big.exc, ph_neg_big.inh); meta.z_ph_neg_big = stats_ph_neg_big.zval;

%% Write to excel file all the results

filename = ['EI_neurons_stats','.xlsx'];
save_name = fullfile(base_dir,filename);
writetable(struct2table(meta),save_name);