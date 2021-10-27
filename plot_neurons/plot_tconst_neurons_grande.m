%% Input variables
type = 'neuron'; % Which data to plot
%type -- 'neuron', 'sim'
%               'neuron' - Plot the neuronal data
%               'sim' - Plot the simulated model

switch type
    case 'neuron'
        kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
        suffix = 'Neuronal_data';
    case 'sim'
        kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_sim_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
        suffix = 'LNP';
end

get_neurons = 'all';
NPSP_th = 40;
plot_f_only = 1;
freq_up_bound = 17000;
freq_down_bound = 700;
save_dir_paper_fig4 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/fig_4/',suffix);
save_dir_paper_fig3 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/fig_3/',suffix);
save_dir_paper_fig2 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/fig_2/',suffix);
save_dir_paper_fig1_sup = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/fig_2_sup_1_2',suffix);

% save_dir_temp = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_sim_noneuron_norm/perfreq_noneuro/ridge/10ms/Plots';
% save_dir_paper_fig4 = save_dir_temp;
% save_dir_paper_fig3 = save_dir_temp;
% save_dir_paper_fig2 = save_dir_temp;
% save_dir_paper_fig1_sup = save_dir_temp;
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

bf_method = 'window_mean'; %The method use to extract the bf:
% bf_method -- 'max', 'window_max', 'window_mean'
%              'max' - Take the max across all history
%              'window_max' - Take the max in a given time window
%              'window_mean' - Take the mean in a given time window
bf_window_ms = 190; %The window for calcualting bf if this method is used

exc_small_color = [0.9804    0.4196    0.6431];
exc_big_color = [0.8 0 0];
inh_small_color = [0.0745 0.6235 1.0000];
inh_big_color = 'b';
font_type = 'Liberation Sans';
legend_font_sz = 32;
sz = 70;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

if ~exist('NPSP_th','var') || isempty(NPSP_th)
    NPSP_th = 40;
end

if ~exist('get_neurons','var') || isempty(get_neurons)
    get_neurons = 'all';
end

[save_dir,~] = fileparts(kernel_dir);
plot_dir = fullfile(save_dir,'Plots');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

save_dir_full = fullfile(plot_dir,get_neurons);
if ~exist(save_dir_full,'dir')
    mkdir(save_dir_full);
end

if ~strcmp(bf_method,'max')
    save_dir_bf = fullfile(save_dir_full,[bf_neurons,'_',bf_method,'_',num2str(bf_window_ms),'ms']);
else
    save_dir_bf = fullfile(save_dir_full,[bf_neurons,'_',bf_method]);
end

if ~exist(save_dir_bf,'dir')
    mkdir(save_dir_bf);
end

% r_type{1} = 'anech';
r_type{1} = 'small';
r_type{2} = 'big';
n_rooms = length(r_type);
%% Load the data
load(fullfile(kernel_dir,'info'),'info');
temp_files = dir([kernel_dir,'/*.mat']);
temp = load(fullfile(temp_files(1).folder,temp_files(1).name));
model = temp.kernel.model;

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

%% First load all the selected kernels
fprintf('== Loading the data ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'kernel');
    kernels{k,1} = kernel;
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% First compute the bf and com for every cluster
freqs = fliplr(kernels{1}.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_h = kernel.n_h;
dt_ms = round(kernel.dt_ms);
chop_ix = round(chop_ms/dt_ms);
bf_window_ix = round(bf_window_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;
fprintf('== Calcuating com and bf ==\n');tic;
for k = 1:n_clust
    sprintf('== Cluster %0.f/%0.f ==\n',k,n_clust);
    for r = 1:n_rooms
        room = r_type{r};
        switch model
            
            case {'sep','sep_kh'}
                [~,ix] = max(kernels{k}.(room).k_f);
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
                k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                
            case {'ridge','lasso','elastic'}
                k_fh = fliplr(kernels{k}.(room).main{end}.k_fh);
                k_fh = k_fh(:,1:chop_ix);
                k_fhn(:,:,k).(room) = k_fh;
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
                max_exc(k).(room) = max(k_h_pos_temp(:));
                max_inh(k).(room) = max(k_h_neg_temp(:));
                %Total Excitation and Inhibition
                total_exc(k).(room) = sum(k_h_pos_temp(:));
                total_inh(k).(room) = sum(k_h_neg_temp(:));
                %IE ratio
                ie_ratio(k).(room) = sum(k_h_neg_temp(:))/sum(k_h_pos_temp(:));
                %Peak Time (PT) values
                [~,max_ix_exc] = max(k_h_pos_temp(:));
                [~, max_ix_col_exc] = ind2sub(size(k_h_pos_temp),max_ix_exc);
                max_exc_time(k).(room) = h(max_ix_col_exc);
                
                [~,max_ix_inh] = max(k_h_neg_temp(:));
                [~, max_ix_col_inh] = ind2sub(size(k_h_neg_temp),max_ix_inh);
                max_inh_time(k).(room) = h(max_ix_col_inh);
  
                switch bf_method
                    case 'max'
                        k_f = max(k_fh_pos,[],2); %Take the max across history steps
                    case 'window_max'
                        k_f = max(k_fh_pos(:,1:bf_window_ix),[],2); %Take the max in a specified window
                    case 'window_mean'
                        k_f = mean(k_fh_pos(:,1:bf_window_ix),2); %Take the mean in a specified window
                end
                
                [~,ix] = max(k_f); %Find the index of the max frequency
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
        end
        
        %Get the center of mass (COM) values
        k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
        k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
        com_neg(k).(room) = (k_h_neg*h); %Compute a weighted sum of all values
        com_pos(k).(room) = (k_h_pos*h); %Compute a weighted sum of all values
    end
    bf_mean(k) = 2^(log2(bf(k).small*bf(k).big)/2);
    [~,ix_bf] = min(abs(bf_mean(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
    bf_closest(k) = freqs(ix_bf);
end
fprintf('== Done! This took %0.fs ==\n',toc);
%% Compute a mean STRF for each BF

for  r = 1:n_rooms
    room = r_type{r};
    k_fhn_temp = reach(k_fhn,room);
    for f = 1:length(freqs)
        curr_freq = freqs(f);
        ix = bf_closest == curr_freq;
        n_bf(f) = sum(ix);
        k_fh = mean(k_fhn_temp(:,:,ix),3);
        k_fh_bf.(room){f} = k_fh;
        k_fh_neg = abs(min(k_fh,0));
        k_fh_pos = abs(max(k_fh,0));
        %Get the IE ratio
        ie_ratio_bf(f).(room) = sum(k_fh_neg(:))/sum(k_fh_pos(:)); %Plot the max/sum of excitation inhibiton changes
        %Make sure don't divide by zero
        k_h_neg = mean(k_fh_neg);
        k_h_pos = mean(k_fh_pos);
        %Get the com
        k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
        k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
        com_neg_bf(f).(room) = (k_h_neg*h); %Compute a weighted sum of all values
        com_pos_bf(f).(room) = (k_h_pos*h); %Compute a weighted sum of all values
    end
end

params.freqs = freqs;
params.n_h = n_h-1;
params.dt_ms = dt_ms;
params.h_max_ms = chop_ms;
params.n_bf = n_bf;
params.font_type = font_type;
params.exc_small_color = exc_small_color; params.exc_big_color = exc_big_color;
params.inh_small_color = inh_small_color; params.inh_big_color = inh_big_color;
params.name = 'neurons';

%% Plot a mean STRF across all frequencies using averaging across all kernels
params.type = 'raw';
params.save_dir = save_dir_paper_fig2;
plot_mean_strf2(k_fhn, params);

%% Plot a mean STRF across all frequencies using pre-averaged kernels
params.type = 'bf';
plot_mean_strf(k_fh_bf,params);

%% Plot STRF and temporal profiles for each bf separately 
params.save_dir = save_dir_paper_fig1_sup;
plot_ind_bf(k_fh_bf,params);

%% Stats
%COM measures
com_small_neg = reach(com_neg,'small');
com_small_pos = reach(com_pos,'small');
com_big_neg = reach(com_neg,'big');
com_big_pos = reach(com_pos,'big');

%IE ratio
ie_ratio_small = reach(ie_ratio,'small');
ie_ratio_big = reach(ie_ratio,'big');

%Total Inhibiton and excitation
total_inh_small = reach(total_inh,'small');
total_inh_big = reach(total_inh,'big');
total_exc_small = reach(total_exc,'small');
total_exc_big = reach(total_exc,'big');

%Inhibition Max value
max_inh_small = reach(max_inh,'small');
max_inh_big = reach(max_inh,'big');
max_exc_small = reach(max_exc,'small');
max_exc_big = reach(max_exc,'big');

%Time of Max Inhibition
max_inh_time_small = reach(max_inh_time,'small');
max_inh_time_big = reach(max_inh_time,'big');
max_exc_time_small = reach(max_exc_time,'small');
max_exc_time_big = reach(max_exc_time,'big');

[p_val_neg,~,~] = signrank(com_big_neg,com_small_neg);
[p_val_pos,~,~] = signrank(com_big_pos,com_small_pos);
[p_val_ie,~,~] = signrank(ie_ratio_big,ie_ratio_small);
[p_val_total_inh,~,~] = signrank(total_inh_big,total_inh_small);
[p_val_max_inh,~,~] = signrank(max_inh_big,max_inh_small);
[p_val_max_inh_time,~,~] = signrank(max_inh_time_big,max_inh_time_small);
[p_val_total_exc,~,~] = signrank(total_exc_big,total_exc_small);
[p_val_max_exc,~,~] = signrank(max_exc_big,max_exc_small);
[p_val_max_exc_time,~,~] = signrank(max_exc_time_big,max_exc_time_small);

com_small_bf_neg = reach(com_neg_bf,'small');
com_small_bf_pos = reach(com_pos_bf,'small');
com_big_bf_neg = reach(com_neg_bf,'big');
com_big_bf_pos = reach(com_pos_bf,'big');

%% Compute some important difference measures
%COM
median_diff.com.neg = nanmedian(com_big_neg-com_small_neg);
median_diff.com.pos = nanmedian(com_big_pos-com_small_pos);
diff_vec.com.neg = com_big_neg-com_small_neg;
diff_vec.com.pos = com_big_pos-com_small_pos;

%PT
median_diff.pt.neg = nanmedian(max_inh_time_big-max_inh_time_small);
median_diff.pt.pos = nanmedian(max_exc_time_big-max_exc_time_small);
diff_vec.pt.neg = max_inh_time_big-max_inh_time_small;
diff_vec.pt.pos = max_exc_time_big-max_exc_time_small;

%PH
median_diff.ph.neg = nanmedian(log2(max_inh_big./max_inh_small));
median_diff.ph.pos = nanmedian(log2(max_exc_big./max_exc_small));
diff_vec.ph.neg = log2(max_inh_big./max_inh_small);
diff_vec.ph.pos = log2(max_exc_big./max_exc_small);

%% Plot histograms and scattergrams

if ~plot_f_only
    %% Define params for histogram for Excitation and Inhibition
    exc_color = 'r'; inh_color = 'b'; lw1 = 7; lw2 = 3;
    params_his_ei.fit = 'neurons'; params_his_ei.units = 'normal'; params_his_ei.calc_method = calc_method;
    %Plot general properties
    params_his_ei.hist_type = hist_type; params_his_ei.get_neurons = get_neurons; params_his_ei.NPSP_th = NPSP_th; params_his_ei.model = model; params_his_ei.font_type = font_type;
    %Plot colours
    params_his_ei.exc_color = exc_color; params_his_ei.inh_color = inh_color; 
    %Plot size
    params_his_ei.all_font_sz = all_font_sz; params_his_ei.lw1 = lw1; params_his_ei.lw2 = lw2;
    %Save dir
    params_his_ei.save_dir =  save_dir_paper_fig3;
    
    %% Define params for Single histogram
    %Plot general properties
    params_norm_his.fit = 'neurons'; params_norm_his.calc_method = calc_method;
    params_norm_his.hist_type = hist_type; params_norm_his.get_neurons = get_neurons; params_norm_his.NPSP_th = NPSP_th; params_norm_his.model = model; params_norm_his.font_type = font_type;
    %Plot colours
    params_norm_his.his_color = 'k';
    %Plot size
    params_norm_his.all_font_sz = all_font_sz; params_norm_his.lw1 = lw1; params_norm_his.lw2 = lw2;
    %Save dir
    params_norm_his.save_dir =  save_dir_paper_fig3;

    %% Define params for scatter plot for Excitation and Inhibition 
    %Plot general properties
    params_ei_scatter.fit = 'neurons'; params_ei_scatter.units = 'normal'; params_ei_scatter.calc_method = calc_method;
    params_ei_scatter.get_neurons = get_neurons; params_ei_scatter.NPSP_th = NPSP_th; params_ei_scatter.model = model; params_ei_scatter.font_type = font_type;
    %Plot colours
    params_ei_scatter.exc_color = exc_color; params_ei_scatter.inh_color = inh_color;
    %Plot size
    params_ei_scatter.lw = 3; params_ei_scatter.all_font_sz = all_font_sz; params_ei_scatter.sz = sz;
    %Save dir
    params_ei_scatter.save_dir = save_dir_paper_fig3;
    
    %% Define params for scatter plot with NPSP 
    %Plot general properties
    params_npsp_scatter.fit = 'neurons';
    params_npsp_scatter.get_neurons = get_neurons; params_npsp_scatter.NPSP_th = NPSP_th; params_npsp_scatter.model = model; params_npsp_scatter.font_type = font_type; 
    %Plot colours
    params_npsp_scatter.NPSPs = NPSPs;
    %Plot size
    params_npsp_scatter.lw = 3; params_npsp_scatter.all_font_sz = all_font_sz; params_npsp_scatter.sz = sz;
    %Save dir
    params_npsp_scatter.save_dir = save_dir_paper_fig3;
    
    %% COM Inhibition small vs big with NPSP Scatter
    params_npsp_scatter.specific_name = 'Inhibitory COM for ';
    params_npsp_scatter.npsp_on = 1;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com_big_neg; params_npsp_scatter.val_small = com_small_neg; 
    params_npsp_scatter.p_val = p_val_neg;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation with NPSP Scatter
    params_npsp_scatter.specific_name = 'Excitatory COM  ';
    params_npsp_scatter.val_big = com_big_pos; params_npsp_scatter.val_small = com_small_pos;
    params_npsp_scatter.p_val = p_val_pos;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = 'Inhibitory and excitatory com no npsp for ';
    params_ei_scatter.val_exc_big = com_big_pos; params_ei_scatter.val_exc_small = com_small_pos; params_ei_scatter.val_inh_big = com_big_neg; params_ei_scatter.val_inh_small = com_small_neg;
    params_ei_scatter.p_val_exc = p_val_pos; params_ei_scatter.p_val_inh = p_val_neg;
    params_ei_scatter.lim_val = lim_val_ms;
    plot_ei_scatter(params_ei_scatter);
    
    %% IE Ratio with NPSP Scatter
    params_npsp_scatter.specific_name = 'IE Ratio for ';
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie;
    params_npsp_scatter.lim_val = 2;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio small vs big Scatter
    params_npsp_scatter.specific_name = 'IE Ratio no NPSP for test ';
    params_npsp_scatter.npsp_on = 0;
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% Total Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = 'E and I total no NPSP for test ';
    params_ei_scatter.val_exc_big = total_exc_big; params_ei_scatter.val_exc_small = total_exc_small; params_ei_scatter.val_inh_big = total_inh_big; params_ei_scatter.val_inh_small = total_inh_small;
    params_ei_scatter.p_val_exc = p_val_total_exc; params_ei_scatter.p_val_inh = p_val_total_inh;
    params_ei_scatter.lim_val = 0.3;
    plot_ei_scatter(params_ei_scatter);
    
    %% PH Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = 'PH E and I no NPSP for ';
    params_ei_scatter.val_exc_big = max_exc_big; params_ei_scatter.val_exc_small = max_exc_small; params_ei_scatter.val_inh_big = max_inh_big; params_ei_scatter.val_inh_small = max_inh_small;
    params_ei_scatter.p_val_exc = p_val_max_exc; params_ei_scatter.p_val_inh = p_val_max_inh;
    params_ei_scatter.lim_val = 0.03;
    plot_ei_scatter(params_ei_scatter);
    
    %% PT Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = 'PT E and I no NPSP for test ';
    g_noise = 5*randn(1,length(max_inh_time_small));
    params_ei_scatter.val_exc_big = max_exc_time_big+g_noise; params_ei_scatter.val_exc_small = max_exc_time_small+g_noise; params_ei_scatter.val_inh_big = max_inh_time_big+g_noise; params_ei_scatter.val_inh_small = max_inh_time_small+g_noise;
    params_ei_scatter.p_val_exc = p_val_max_exc_time; params_ei_scatter.p_val_inh = p_val_max_inh_time;
    params_ei_scatter.lim_val = 200;
    plot_ei_scatter(params_ei_scatter);
    
    %% COM Excitation and Inhibiton small vs big Histogram
    params_his_ei.specific_name = ' COM ms change for ';
    params_his_ei.val_exc_big = com_big_pos; params_his_ei.val_exc_small = com_small_pos; params_his_ei.val_inh_big = com_big_neg; params_his_ei.val_inh_small = com_small_neg;  
    params_his_ei.p_val_exc = p_val_pos; params_his_ei.p_val_inh = p_val_neg;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);

    %% COM change ix Excitation and Inhibiton small vs big Histogram
    params_his_ei.units = 'change ix';
    params_his_ei.specific_name = ' COM change ix for test ';
    params_his_ei.val_exc_big = com_big_pos; params_his_ei.val_exc_small = com_small_pos; params_his_ei.val_inh_big = com_big_neg; params_his_ei.val_inh_small = com_small_neg;
    params_his_ei.p_val_exc = p_val_pos; params_his_ei.p_val_inh = p_val_neg;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% Total Excitation and Inhibiton small vs big room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' of Total amount change for ';
    params_his_ei.val_exc_big = total_exc_big; params_his_ei.val_exc_small = total_exc_small; params_his_ei.val_inh_big = total_inh_big; params_his_ei.val_inh_small = total_inh_small;
    params_his_ei.p_val_exc = p_val_total_exc; params_his_ei.p_val_inh = p_val_total_inh;
    params_his_ei.his_spacing = 0.1; params_his_ei.lim_val = 4;
    plot_ei_histogram(params_his_ei);
    
    %% PH Excitation and Inhibiton for small vs big room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' PH E and I change for ';
    params_his_ei.val_exc_big = max_exc_big; params_his_ei.val_exc_small = max_exc_small; params_his_ei.val_inh_big = max_inh_big; params_his_ei.val_inh_small = max_inh_small;
    params_his_ei.p_val_exc = p_val_max_exc; params_his_ei.p_val_inh = p_val_max_inh;
    params_his_ei.his_spacing = 0.125; params_his_ei.lim_val = 2.5;
    plot_ei_histogram(params_his_ei);
    
    %% PT Excitation and Inhibiton for small vs big room Histogram
    params_his_ei.units = 'normal';
    params_his_ei.specific_name = ' PT E and I change for ';
    params_his_ei.val_exc_big = max_exc_time_big; params_his_ei.val_exc_small = max_exc_time_small; params_his_ei.val_inh_big = max_inh_time_big; params_his_ei.val_inh_small = max_inh_time_small;
    params_his_ei.p_val_exc = p_val_max_exc_time; params_his_ei.p_val_inh = p_val_max_inh_time;
    params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
    plot_ei_histogram(params_his_ei);
      
    %% IE Ratio change ix small vs big room Histogram
    params_norm_his.specific_name = ' IE Ratio of com change ix for ';
    params_norm_his.units = 'change ix'; 
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie; 
    params_norm_his.his_spacing = 2.5; params_norm_his.lim_val = 70;
    plot_norm_histogram(params_norm_his);
    
    %% IE ratio for small vs big room Histogram
    params_norm_his.specific_name = ' of IE Ratio AU change for test ';
    params_norm_his.units = 'normal';
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie;
    params_norm_his.his_spacing = 0.05; params_norm_his.lim_val = 1;
    plot_norm_histogram(params_norm_his);
    
end

%% Plot no of neurons at BF as a function of freqeuncy
lw = 6;
figure('units','normalized','outerposition',[0 0 1 1]);
plot_freq = freqs;
plot(plot_freq,n_bf,'Color','k','LineWidth',lw);
xlabel('Best Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Number of neurons','FontSize',y_font_sz,'FontWeight','Normal');
title(['Number of neurons per BF  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir_bf,['Number of neurons per BF  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot inhibitory com vs frequency for all rooms
skip_f = 3;
lw = 5;
x_box_dis = 0.7;
%Get the BFs for all clusters for small and big rooms
switch bf_neurons
    case 'shared'
        bf_small = bf_closest; bf_small = bf_small(:);
        bf_big = bf_closest; bf_big = bf_big(:);
    case 'ind'
        bf_small = reach(bf,'small'); bf_small = bf_small(:);
        bf_big = reach(bf,'big'); bf_big = bf_big(:);
end


for f = 1:length(freqs)
    curr_freq = freqs(f);
    ix_small = bf_small == curr_freq;
    ix_big = bf_big == curr_freq;    
    com_small_mean(f) = nanmean(com_small_neg(ix_small));
    com_small_std(f) = nanstd(com_small_neg(ix_small))./sqrt(sum(ix_small));
    com_big_mean(f) = nanmean(com_big_neg(ix_big));
    com_big_std(f) = nanstd(com_big_neg(ix_big))./sqrt(sum(ix_big));
end

ix = freqs>freq_down_bound & freqs<freq_up_bound;
plot_freqs = fliplr(freqs(ix));
n_f = length(plot_freqs);

for f = 1:n_f
    f_labels{f} = num2str(plot_freqs(f)./1000,'%.1f');
end

plot_com_small_mean_inh = fliplr(com_small_mean(ix)); plot_com_small_std_inh = fliplr(com_small_std(ix));
plot_com_big_mean_inh = fliplr(com_big_mean(ix)); plot_com_big_std_inh = fliplr(com_big_std(ix));

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
shadedErrorBar(log10(plot_freqs),plot_com_small_mean_inh,plot_com_small_std_inh,{'LineWidth',lw,'Color',inh_small_color}); hold on;
shadedErrorBar(log10(plot_freqs),plot_com_big_mean_inh,plot_com_big_std_inh,{'LineWidth',lw,'Color',inh_big_color});
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('COM [ms]','FontSize',y_font_sz,'FontWeight','Normal');
% title(['Inhibitory COM vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room inh'),'LineStyle','none','Color',inh_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room inh'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
xticks([log10(plot_freqs(2:skip_f:n_f))]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
ylim([80 130]);
save_name = fullfile(save_dir_paper_fig4,['Inhibitory com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close;

%% Plot Agregated inhibitory com vs frequency for all rooms
f_step = 2;
f_start = [1:f_step+1:length(freqs)];
%Get the BFs for all clusters for small and big rooms
clear com_small_mean;
clear com_small_std;
clear com_big_mean;
clear com_big_std;
clear plot_freq;

for f = 1:length(f_start)
    curr_freqs = freqs(f_start(f):f_start(f)+f_step);
    ix_small = bf_small == curr_freqs(1) | bf_small == curr_freqs(2) | bf_small == curr_freqs(3);
    ix_big = bf_big == curr_freqs(1) | bf_big == curr_freqs(2) | bf_big == curr_freqs(3);   
    com_small_mean(f) = nanmean(com_small_neg(ix_small));
    com_small_std(f) = nanstd(com_small_neg(ix_small))./sqrt(sum(ix_small));
    com_big_mean(f) = nanmean(com_big_neg(ix_big));
    com_big_std(f) = nanstd(com_big_neg(ix_big))./sqrt(sum(ix_big));
    plot_freq(f) = 2^(log2(prod(curr_freqs))/length(curr_freqs));
end

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

shadedErrorBar(plot_freq,com_small_mean,com_small_std,{'LineWidth',lw,'Color',inh_small_color}); hold on;
shadedErrorBar(plot_freq,com_big_mean,com_big_std,{'LineWidth',lw,'Color','b'});
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Inhibitory \com vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room inh'),'LineStyle','none','Color',inh_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room inh'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
ylim([80 130]);
save_name = fullfile(save_dir_paper_fig4,['Agregated Inhibitory com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot excitatory com vs frequency for all rooms
for f = 1:length(freqs)    
    curr_freq = freqs(f);
    ix_small = bf_small == curr_freq;
    ix_big = bf_big == curr_freq;    
    com_small_mean(f) = nanmean(com_small_pos(ix_small));
    com_small_std(f) = nanstd(com_small_pos(ix_small))./sqrt(sum(ix_small));
    com_big_mean(f) = nanmean(com_big_pos(ix_big));
    com_big_std(f) = nanstd(com_big_pos(ix_big))./sqrt(sum(ix_big));
end

plot_com_small_mean_exc = fliplr(com_small_mean(ix)); plot_com_small_std_exc = fliplr(com_small_std(ix));
plot_com_big_mean_exc = fliplr(com_big_mean(ix)); plot_com_big_std_exc = fliplr(com_big_std(ix));

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); hold on;
shadedErrorBar(plot_freqs, plot_com_small_mean_exc, plot_com_small_std_exc,{'LineWidth',lw,'Color',exc_small_color});
shadedErrorBar(plot_freqs, plot_com_big_mean_exc, plot_com_big_std_exc,{'LineWidth',lw,'Color',exc_big_color});
hold off;

xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('COM [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Excitatory COM vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room exc'),'LineStyle','none','Color',exc_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room exc'),'LineStyle','none','Color',exc_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
ylim([35 85]);
save_name = fullfile(save_dir_paper_fig4,['Excitatory com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot excitatory and inhbiitory together com vs frequency for all rooms

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); hold on;
shadedErrorBar(plot_freqs, plot_com_small_mean_exc, plot_com_small_std_exc,{'LineWidth',lw,'Color',exc_small_color});
shadedErrorBar(plot_freqs, plot_com_big_mean_exc, plot_com_big_std_exc,{'LineWidth',lw,'Color',exc_big_color});
shadedErrorBar(plot_freqs, plot_com_small_mean_inh, plot_com_small_std_inh,{'LineWidth',lw,'Color',inh_small_color});
shadedErrorBar(plot_freqs, plot_com_big_mean_inh, plot_com_big_std_inh,{'LineWidth',lw,'Color', inh_big_color});
hold off;

% xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
% title(['Excitatory \com vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room inh'),'LineStyle','none','Color','b','FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room inh'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[x_box_dis 0.7 0.1 0.1],'String', sprintf('Large room exc'),'LineStyle','none','Color',exc_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[x_box_dis 0.65 0.1 0.1],'String', sprintf('Small room exc'),'LineStyle','none','Color',exc_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
ylim([35 130]);
save_name = fullfile(save_dir_paper_fig4,['Together com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;
%% Plot BF averaged inhibitory com vs frequency for all rooms
lw = 3;
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
plot(plot_freq,com_small_bf_neg,'Color',[0.31 0.75 0.405],'LineWidth',lw);
plot(plot_freq,com_big_bf_neg,'Color',[0.642 0.456 0.924],'LineWidth',lw);
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('COM [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['BF Averaged Inhibitory COM vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.75 0.8 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',inh_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
ylim([70 140]);
save_name = fullfile(save_dir_bf,['BF Averaged Inhibitory com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot BF averaged excitatory com vs frequency for all rooms
lw = 3;
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
plot(plot_freq,com_small_bf_pos,'Color',[0.31 0.75 0.405],'LineWidth',lw);
plot(plot_freq,com_big_bf_pos,'Color',[0.642 0.456 0.924],'LineWidth',lw);
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['BF Averaged Excitatory \com vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',exc_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',exc_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
ylim([20 100]);
save_name = fullfile(save_dir_bf,['BF Averaged Excitatory com vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot Inhibitory/Excitatory ratio vs frequency for all rooms
lw = 3;
for f = 1:length(freqs)    
    curr_freq = freqs(f);
    ix_small = bf_small == curr_freq;
    ix_big = bf_big == curr_freq;
    nn_f(f) = sum(ix_small);
    ie_ratio_small_mean(f) = nanmean(ie_ratio_small(ix_small));
    ie_ratio_small_std(f) = nanstd(ie_ratio_small(ix_small))./sqrt(sum(ix_small)); %Take sqrt here and ehck all other plots
    ie_ratio_big_mean(f) = nanmean(ie_ratio_big(ix_big));
    ie_ratio_big_std(f) = nanstd(ie_ratio_big(ix_big))./sqrt(sum(ix_big));
end
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
shadedErrorBar(plot_freq,ie_ratio_small_mean,ie_ratio_small_std,{'Color',[0.31 0.75 0.405]});
shadedErrorBar(plot_freq,ie_ratio_big_mean,ie_ratio_big_std,{'Color',[0.642 0.456 0.924]});
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('IE ratio [AU]','FontSize',y_font_sz,'FontWeight','Normal');
title(['IE Ratio vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
ylim([0 2]);
save_name = fullfile(save_dir_bf,['IE Ratio vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot bf consistency accross clusters with histograms
axis_sz = 15;
diff_bf = abs(log2(bf_big./bf_small));
figure('units','normalized','outerposition',[0 0 1 1]);
histogram(diff_bf,'Normalization','probability');
xlabel('\DeltaBF [Octaves]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Probability','FontSize',y_font_sz,'FontWeight','Normal');
title(['Histogram of BF consistency small<->big  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Histogram of BF consistency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot inhibitory com vs bf no binning small big separate
sz = 30;
com_small = reach(com_neg,'small');
com_big = reach(com_neg,'big');
figure('units','normalized','outerposition',[0 0 1 1]);
plot(bf_small,com_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
plot(bf_big,com_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);
xlabel('BF [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Inhibitory \com vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',inh_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gca, 'XScale', 'log');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Inhibitory com vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Compute and  plot regression lines for inhibitory coms

ix_small = bf_small>freq_down_bound & bf_small<freq_up_bound;
ix_big = bf_big>freq_down_bound & bf_big<freq_up_bound;
n_reg_small = sum(ix_small); 
n_reg_big = sum(ix_big); 

X_small = [ones(n_reg_small,1),log10(bf_small(ix_small))];
X_big = [ones(n_reg_big,1),log10(bf_big(ix_big))];
plot_bf_small = bf_small(ix_small);
plot_bf_large = bf_big(ix_big);

%Inhibitory regression
y_small_neg = com_small(ix_small); y_small_neg = y_small_neg(:);
y_big_neg = com_big(ix_big); y_big_neg = y_big_neg(:);

[b_small_neg,~,~,~,stats_small_neg] = regress(y_small_neg,X_small);
[b_big_neg,~,~,~,stats_big_neg] = regress(y_big_neg,X_big);
r_small_neg =  sign(b_small_neg(2))*sqrt(stats_small_neg(1)); p_small_neg = stats_small_neg(3);
r_big_neg = sign(b_big_neg(2))*sqrt(stats_big_neg(1)); p_big_neg = stats_big_neg(3);

figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(bf_big),com_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);hold on;
set(gca, 'XScale', 'log');
h1 = refline([b_big_neg(2),b_big_neg(1)]);
h1.Color = 'r';
h1.LineWidth = 3;
annotation('textbox', [0.75, 0.72, 0.3, 0.2], 'String', sprintf('\\beta_{1}=%.3f\n',b_big_neg(2)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.65, 0.3, 0.2], 'String', sprintf('r_{big}=%.2f\n',r_big_neg),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.58, 0.3, 0.2], 'String', sprintf('p=%.5f\n',p_big_neg),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
xlabel('BF [Hz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Inhibitory \com vs BF Large room  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Inhibitory big com vs BF regression for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;


figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(bf_small),com_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
set(gca, 'XScale', 'log');
h2 = refline([b_small_neg(2),b_small_neg(1)]);
h2.Color = 'r';
h2.LineWidth = 3;
annotation('textbox', [0.75, 0.72, 0.3, 0.2], 'String', sprintf('\\beta_{1}=%.3f\n',b_small_neg(2)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.65, 0.3, 0.2], 'String', sprintf('r_{small}=%.3f\n',r_small_neg),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.58, 0.3, 0.2], 'String', sprintf('p=%.5f\n',p_small_neg),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
xlabel('BF [Hz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Inhibitory \com vs BF Small room  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Inhibitory small com vs BF regression for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot excitatory com vs bf no binning small big separate
com_small = reach(com_pos,'small');
com_big = reach(com_pos,'big');
figure('units','normalized','outerposition',[0 0 1 1]);
plot(bf_small,com_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
plot(bf_big,com_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);
xlabel('BF [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Excitatory \com vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gca, 'XScale', 'log');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Excitatory com vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Compute and  plot regression lines for excitatory coms

%Excitatory regression
y_small_pos = com_small(ix_small); y_small_pos = y_small_pos(:);
y_big_pos = com_big(ix_big); y_big_pos = y_big_pos(:);

[b_small_pos,~,~,~,stats_small_pos] = regress(y_small_pos,X_small);
[b_big_pos,~,~,~,stats_big_pos] = regress(y_big_pos,X_big);
r_small_pos = sign(b_small_pos(2))*sqrt(stats_small_pos(1)); p_small_pos = stats_small_pos(3);
r_big_pos = sign(b_big_pos(2))*sqrt(stats_big_pos(1)); p_big_pos = stats_big_pos(3);


figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(bf_big),com_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);hold on;
set(gca, 'XScale', 'log');
h1 = refline([b_big_pos(2),b_big_pos(1)]);
h1.Color = 'r';
h1.LineWidth = 3;
annotation('textbox', [0.75, 0.72, 0.3, 0.2], 'String', sprintf('\\beta_{1}=%.3f\n',b_big_pos(2)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.65, 0.3, 0.2], 'String', sprintf('r_{big}=%.2f\n',r_big_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.58, 0.3, 0.2], 'String', sprintf('p=%.5f\n',p_big_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
xlabel('BF [Hz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Excitatory \com vs BF Large room  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Excitatory big com vs BF regression for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;


figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(bf_small),com_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
set(gca, 'XScale', 'log');
h2 = refline([b_small_pos(2),b_small_pos(1)]);
h2.Color = 'r';
h2.LineWidth = 3;
annotation('textbox', [0.75, 0.72, 0.3, 0.2], 'String', sprintf('\\beta_{1}=%.3f\n',b_small_pos(2)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.65, 0.3, 0.2], 'String', sprintf('r_{small}=%.3f\n',r_small_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox', [0.75, 0.58, 0.3, 0.2], 'String', sprintf('p=%.5f\n',p_small_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','Normal');
xlabel('BF [Hz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('\com [ms]','FontSize',y_font_sz,'FontWeight','Normal');
title(['Excitatory \com vs BF Small room  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_bf,['Excitatory small com vs BF regression for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot regression inhibitory and excitatory coms together
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); hold on;
plot(log10(plot_bf_small), y_small_neg,'.','Color',inh_small_color,'MarkerSize',sz);hold on;
plot(log10(plot_bf_large), y_big_neg,'.','Color',inh_big_color,'MarkerSize',sz);
plot(log10(plot_bf_small), y_small_pos,'.','Color',exc_small_color,'MarkerSize',sz);
plot(log10(plot_bf_large), y_big_pos,'.','Color',exc_big_color,'MarkerSize',sz);
set(gca, 'XScale', 'log');
h1 = refline([b_small_neg(2),b_small_neg(1)]); h1.Color = inh_small_color; h1.LineWidth = lw; h1.LineStyle = '--';
h2 = refline([b_big_neg(2),b_big_neg(1)]); h2.Color = inh_big_color; h2.LineWidth = lw; h2.LineStyle = '--';
h3 = refline([b_small_pos(2),b_small_pos(1)]); h3.Color = exc_small_color; h3.LineWidth = lw; h3.LineStyle = '--';
h4 = refline([b_big_pos(2),b_big_pos(1)]); h4.Color = exc_big_color; h4.LineWidth = lw; h4.LineStyle = '--';
hold off;

xticks([log10(plot_freqs(2:skip_f:n_f))]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
ylim([20 200]);
annotation('textbox',[0.54 0.88 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}\n',r_big_neg,p_star(p_big_neg)),'LineStyle','none','Color',inh_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.54 0.8 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}\n',r_small_neg,p_star(p_small_neg)),'LineStyle','none','Color',inh_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.75 0.88 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}\n',r_big_pos,p_star(p_big_pos)),'LineStyle','none','Color',exc_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.75 0.8 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}\n',r_small_pos,p_star(p_small_pos)),'LineStyle','none','Color',exc_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir_paper_fig4,['Regression for both excitatory and inhibitory com for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;

%% Plot regression inhibitory and excitatory coms together @ BF for both rooms
legend_font_sz = 50;
all_font_sz = 75;
n_reg = length(plot_freqs);
X = [ones(n_reg,1),log10(plot_freqs')];

%Inhibitory regression
y_small_neg = plot_com_small_mean_inh; y_small_neg = y_small_neg(:);
y_big_neg = plot_com_big_mean_inh; y_big_neg = y_big_neg(:);

[b_small_neg,~,~,~,stats_small_neg] = regress(y_small_neg,X);
[b_big_neg,~,~,~,stats_big_neg] = regress(y_big_neg,X);
r_small_neg =  sign(b_small_neg(2))*sqrt(stats_small_neg(1)); p_small_neg = stats_small_neg(3);
r_big_neg = sign(b_big_neg(2))*sqrt(stats_big_neg(1)); p_big_neg = stats_big_neg(3);

%Excitatory regression
y_small_pos = plot_com_small_mean_exc; y_small_pos = y_small_pos(:);
y_big_pos = plot_com_big_mean_exc; y_big_pos = y_big_pos(:);

[b_small_pos,~,~,~,stats_small_pos] = regress(y_small_pos,X);
[b_big_pos,~,~,~,stats_big_pos] = regress(y_big_pos,X);
r_small_pos =  sign(b_small_pos(2))*sqrt(stats_small_pos(1)); p_small_pos = stats_small_pos(3);
r_big_pos = sign(b_big_pos(2))*sqrt(stats_big_pos(1)); p_big_pos = stats_big_pos(3);

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); hold on;

shadedErrorBar(log10(plot_freqs), plot_com_small_mean_exc, plot_com_small_std_exc,{'LineWidth',lw,'Color',exc_small_color});
shadedErrorBar(log10(plot_freqs), plot_com_big_mean_exc, plot_com_big_std_exc,{'LineWidth',lw,'Color',exc_big_color});
shadedErrorBar(log10(plot_freqs), plot_com_small_mean_inh, plot_com_small_std_inh,{'LineWidth',lw,'Color',inh_small_color});
shadedErrorBar(log10(plot_freqs), plot_com_big_mean_inh, plot_com_big_std_inh,{'LineWidth',lw,'Color', inh_big_color});
set(gca, 'XScale', 'log');
h1 = refline([b_small_neg(2),b_small_neg(1)]); h1.Color = inh_small_color; h1.LineWidth = lw; h1.LineStyle = '--';
h2 = refline([b_big_neg(2),b_big_neg(1)]); h2.Color = inh_big_color; h2.LineWidth = lw; h2.LineStyle = '--';
h3 = refline([b_small_pos(2),b_small_pos(1)]); h3.Color = exc_small_color; h3.LineWidth = lw; h3.LineStyle = '--';
h4 = refline([b_big_pos(2),b_big_pos(1)]); h4.Color = exc_big_color; h4.LineWidth = lw; h4.LineStyle = '--';
hold off;

xticks([log10(plot_freqs(2:skip_f:n_f))]);
xticklabels(f_labels(2:skip_f:n_f));
ylim([30 135]);
yticks([30:20:130]);
set(gca,'XMinorTick','off');
annotation('textbox',[0.54 0.88 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}',r_big_neg,p_star(p_big_neg)),'LineStyle','none','Color',inh_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.54 0.8 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}',r_small_neg,p_star(p_small_neg)),'LineStyle','none','Color',inh_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.77 0.88 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}',r_big_pos,p_star(p_small_pos)),'LineStyle','none','Color',exc_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.77 0.8 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}',r_small_pos,p_star(p_big_pos)),'LineStyle','none','Color',exc_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir_paper_fig4,['Regression for both excitatory and inhibitory com at BF for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
% saveas(gcf,save_name);
% close all;
%% Plot IE ratio vs bf no binning small big separate
figure('units','normalized','outerposition',[0 0 1 1]);
plot(bf_small,ie_ratio_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
plot(bf_big,ie_ratio_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);
xlabel('BF [kHz]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('IE Ratio [AU]','FontSize',y_font_sz,'FontWeight','Normal');
title(['IE Ratio vs BF all neurons  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gca, 'XScale', 'log');
set(gcf,'color','w');
ylim([0 2]);
save_name = fullfile(save_dir_bf,['IE Ratio vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
saveas(gcf,save_name);
close all;