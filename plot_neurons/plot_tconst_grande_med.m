%% Input variables
type = 'neuron'; % Which data to plot
%type -- 'neuron', 'sim'
%               'neuron' - Plot the neuronal data
%               'sim' - Plot the simulated model

switch type
    
    case 'neuron'
        kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm_med/perfreq_noneuro/ridge/10ms/200ms';
        suffix = 'Neuronal_data';
        
    case 'sim'
        kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_sim/perfreq/ridge/10ms/200ms';
        suffix = 'LNP';
end

get_neurons = 'all';
NPSP_th = 40;
plot_f_only = 0;
freq_up_bound = 17000;
freq_down_bound = 700;
save_dir_paper_fig4 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/Med_test/fig_4/',suffix);
save_dir_paper_fig3 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/Med_test/fig_3/',suffix);
save_dir_paper_fig2 = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/Med_test/fig_2/',suffix);
save_dir_paper_fig1_sup = fullfile('/mnt/40086D4C086D41D0/Reverb_paper/Med_test/fig_1_sup/',suffix);

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

%Select the colors
exc_small_color = [0.9804 0.4196 0.6431];
exc_med_color = [0.8902 0.2098 0.3216];
exc_big_color = [0.8 0 0];

inh_small_color = [0.0745 0.6235 1.0000];
inh_med_color = [0.0372 0.3118 1];
inh_big_color = [0 0 1];

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
r_type{2} = 'med';
r_type{3} = 'big';
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
    bf_mean(k) = 2^(log2(bf(k).small*bf(k).med*bf(k).big)/3);
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

%% Plot a mean STRF across all frequencies using pre-averaged kernels
% params.type = 'bf';
% params.save_dir = save_dir_paper_fig2;
% plot_mean_strf(k_fh_bf,params);
% 
% %% Plot a mean STRF across all frequencies using averaging across all kernels
% params.type = 'raw';
% plot_mean_strf(k_fhn,params);
% 
% %% Plot STRF and temporal profiles for each bf separately 
% params.save_dir = save_dir_paper_fig1_sup;
% plot_ind_bf(k_fh_bf,params);

%% Stats
%COM measures
com.neg.small = reach(com_neg,'small');
com.pos.small = reach(com_pos,'small');
com.neg.med = reach(com_neg,'med');
com.pos.med = reach(com_pos,'med');
com.neg.large = reach(com_neg,'big');
com.pos.large = reach(com_pos,'big');

%IE ratio
ie_ratio_small = reach(ie_ratio,'small');
ie_ratio_med = reach(ie_ratio,'med');
ie_ratio_big = reach(ie_ratio,'big');

%Total Inhibiton and excitation
total_inh_small = reach(total_inh,'small');
total_inh_med = reach(total_inh,'med');
total_inh_big = reach(total_inh,'big');
total_exc_small = reach(total_exc,'small');
total_exc_med = reach(total_exc,'med');
total_exc_big = reach(total_exc,'big');

%Inhibition Max value
ph.neg.small = reach(max_inh,'small');
ph.neg.med = reach(max_inh,'med');
ph.neg.large = reach(max_inh,'big');
ph.pos.small = reach(max_exc,'small');
ph.pos.med = reach(max_exc,'med');
ph.pos.large = reach(max_exc,'big');

%Time of Max Inhibition
pt.neg.small = reach(max_inh_time,'small');
pt.neg.med = reach(max_inh_time,'med');
pt.neg.large = reach(max_inh_time,'big');
pt.pos.small = reach(max_exc_time,'small');
pt.pos.med = reach(max_exc_time,'med');
pt.pos.large = reach(max_exc_time,'big');

[p_val.com.neg.ls,~,~] = signrank(com.neg.large,com.neg.small);
[p_val.com.neg.lm,~,~] = signrank(com.neg.large,com.neg.med); 
[p_val.com.neg.ms,~,~] = signrank(com.neg.med,com.neg.small); 

[p_val.com.pos.ls,~,~] = signrank(com.pos.large,com.pos.small); 
[p_val.com.pos.lm,~,~] = signrank(com.pos.large,com.pos.med);
[p_val.com.pos.ms,~,~] = signrank(com.pos.med,com.pos.small);

[p_val_ie_ls,~,~] = signrank(ie_ratio_big,ie_ratio_small);
[p_val_ie_lm,~,~] = signrank(ie_ratio_big,ie_ratio_med);
[p_val_ie_ms,~,~] = signrank(ie_ratio_med,ie_ratio_small);

[p_val_total_inh_ls,~,~] = signrank(total_inh_big,total_inh_small);
[p_val_total_inh_lm,~,~] = signrank(total_inh_big,total_inh_med);
[p_val_total_inh_ms,~,~] = signrank(total_inh_med,total_inh_small);

[p_val.ph.neg.ls,~,~] = signrank(ph.neg.large,ph.neg.small);
[p_val.ph.neg.lm,~,~] = signrank(ph.neg.large,ph.neg.med);
[p_val.ph.neg.ms,~,~] = signrank(ph.neg.med,ph.neg.small);

[p_val.pt.neg.ls,~,~] = signrank(pt.neg.large,pt.neg.small);
[p_val.pt.neg.lm,~,~] = signrank(pt.neg.large,pt.neg.med);
[p_val.pt.neg.ms,~,~] = signrank(pt.neg.med,pt.neg.small);

[p_val_total_exc_ls,~,~] = signrank(total_exc_big,total_exc_small);
[p_val_total_exc_lm,~,~] = signrank(total_exc_big,total_exc_med);
[p_val_total_exc_ms,~,~] = signrank(total_exc_med,total_exc_small);

[p_val.ph.pos.ls,~,~] = signrank(ph.pos.large,ph.pos.small);
[p_val.ph.pos.lm,~,~] = signrank(ph.pos.large,ph.pos.med);
[p_val.ph.pos.ms,~,~] = signrank(ph.pos.med,ph.pos.small);

[p_val.pt.pos.ls,~,~] = signrank(pt.pos.large, pt.pos.small);
[p_val.pt.pos.lm,~,~] = signrank(pt.pos.large, pt.pos.med);
[p_val.pt.pos.ms,~,~] = signrank(pt.pos.med, pt.pos.small);

com_small_bf_neg = reach(com_neg_bf,'small');
com_small_bf_pos = reach(com_pos_bf,'small');
com_med_bf_neg = reach(com_neg_bf,'med');
com_med_bf_pos = reach(com_pos_bf,'med');
com_big_bf_neg = reach(com_neg_bf,'big');
com_big_bf_pos = reach(com_pos_bf,'big');

%% Compute some important difference measures
%COM
median_diff.com.neg.ls = nanmedian(com.neg.large-com.neg.small);
median_diff.com.pos.ls = nanmedian(com.pos.large-com.pos.small);
median_diff.com.neg.lm  = nanmedian(com.neg.large-com.neg.med);
median_diff.com.pos.lm = nanmedian(com.pos.large-com.pos.med);
median_diff.com.neg.ms  = nanmedian(com.neg.med-com.neg.small);
median_diff.com.pos.ms = nanmedian(com.pos.med-com.pos.small);

mean_diff.com.neg.ls = nanmean(com.neg.large-com.neg.small);
mean_diff.com.pos.ls = nanmean(com.pos.large-com.pos.small);
mean_diff.com.neg.lm  = nanmean(com.neg.large-com.neg.med);
mean_diff.com.pos.lm = nanmean(com.pos.large-com.pos.med);
mean_diff.com.neg.ms  = nanmean(com.neg.med-com.neg.small);
mean_diff.com.pos.ms = nanmean(com.pos.med-com.pos.small);

%PT
median_diff.pt.neg.ls = nanmedian(pt.neg.large-pt.neg.small);
median_diff.pt.pos.ls = nanmedian(pt.pos.large-pt.pos.small);
median_diff.pt.neg.lm = nanmedian(pt.neg.large-pt.neg.med);
median_diff.pt.pos.lm = nanmedian(pt.pos.large-pt.pos.med);
median_diff.pt.neg.ms = nanmedian(pt.neg.med-pt.neg.small);
median_diff.pt.pos.ms = nanmedian(pt.pos.med-pt.pos.small);

mean_diff.pt.neg.ls = nanmean(pt.neg.large-pt.neg.small);
mean_diff.pt.pos.ls = nanmean(pt.pos.large-pt.pos.small);
mean_diff.pt.neg.lm = nanmean(pt.neg.large-pt.neg.med);
mean_diff.pt.pos.lm = nanmean(pt.pos.large-pt.pos.med);
mean_diff.pt.neg.ms = nanmean(pt.neg.med-pt.neg.small);
mean_diff.pt.pos.ms = nanmean(pt.pos.med-pt.pos.small);

%PH
median_diff.ph.neg.ls = nanmedian(log2(ph.neg.large./ph.neg.small));
median_diff.ph.pos.ls = nanmedian(log2(ph.pos.large./ph.pos.small));
median_diff.ph.neg.lm = nanmedian(log2(ph.neg.large./ph.neg.med));
median_diff.ph.pos.lm = nanmedian(log2(ph.pos.large./ph.pos.med));
median_diff.ph.neg.ms = nanmedian(log2(ph.neg.med./ph.neg.small));
median_diff.ph.pos.ms = nanmedian(log2(ph.pos.med./ph.pos.small));

mean_diff.ph.neg.ls = nanmean(log2(ph.neg.large./ph.neg.small));
mean_diff.ph.pos.ls = nanmean(log2(ph.pos.large./ph.pos.small));
mean_diff.ph.neg.lm = nanmean(log2(ph.neg.large./ph.neg.med));
mean_diff.ph.pos.lm = nanmean(log2(ph.pos.large./ph.pos.med));
mean_diff.ph.neg.ms = nanmean(log2(ph.neg.med./ph.neg.small));
mean_diff.ph.pos.ms = nanmean(log2(ph.pos.med./ph.pos.small));

%% Perform Kruskal-Wallis and multcompare for COM, PT and PH
display_results = 'on';
mult_corr = 'lsd';
% 'tukey-kramer' or 'hsd'	
%     Tukey's honest significant difference criterion
% 'bonferroni'	
%     Bonferroni method
% 'dunn-sidak'	
%     Dunn and Sidák’s approach
% 'lsd'	
%     Fisher's least significant difference procedure
% 'scheffe'	
%     Scheffé's S procedure

%% COM

% - Inhibition
%Perform Kruskal-Wallis (KW) test to compare the distributions for COM
data_com_neg = [com.neg.small'; com.neg.med'; com.neg.large']; %Make a vector with the data for KW test
group_com_neg = [repmat({'Small'},length(com.neg.small),1); repmat({'Medium'},length(com.neg.med),1); repmat({'Large'},length(com.neg.large),1)]; %Make a vector with the group labels

[pval_kw.com.neg.single, pval_kw.com.neg.table, pval_kw.com.neg.stats] = kruskalwallis(data_com_neg, group_com_neg, display_results);
pval_kw.com.neg.multiple = multcompare(pval_kw.com.neg.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Excitation
%Perform Kruskal-Wallis (KW) test to compare the distributions for COM
data_com_pos = [com.pos.small'; com.pos.med'; com.pos.large']; %Make a vector with the data for KW test
group_com_pos = [repmat({'Small'},length(com.pos.small),1); repmat({'Medium'},length(com.pos.med),1); repmat({'Large'},length(com.pos.large),1)]; %Make a vector with the group labels

[pval_kw.com.pos.single, pval_kw.com.pos.table, pval_kw.com.pos.stats] = kruskalwallis(data_com_pos, group_com_pos, display_results);
pval_kw.com.pos.multiple = multcompare(pval_kw.com.pos.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

%% PT

% - Inhibition
%Perform Kruskal-Wallis (KW) test to compare the distributions for PT
data_pt_neg = [pt.neg.small'; pt.neg.med'; pt.neg.large']; %Make a vector with the data for KW test
group_pt_neg = [repmat({'Small'},length(pt.neg.small),1); repmat({'Medium'},length(pt.neg.med),1); repmat({'Large'},length(pt.neg.large),1)]; %Make a vector with the group labels

[pval_kw.pt.neg.single, pval_kw.pt.neg.table, pval_kw.pt.neg.stats] = kruskalwallis(data_pt_neg, group_pt_neg, display_results);
pval_kw.pt.neg.multiple = multcompare(pval_kw.pt.neg.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Excitation
%Perform Kruskal-Wallis (KW) test to compare the distributions for PT
data_pt_pos = [pt.pos.small'; pt.pos.med'; pt.pos.large']; %Make a vector with the data for KW test
group_pt_pos = [repmat({'Small'},length(pt.pos.small),1); repmat({'Medium'},length(pt.pos.med),1); repmat({'Large'},length(pt.pos.large),1)]; %Make a vector with the group labels

[pval_kw.pt.pos.single, pval_kw.pt.pos.table, pval_kw.pt.pos.stats] = kruskalwallis(data_pt_pos, group_pt_pos, display_results);
pval_kw.pt.pos.multiple = multcompare(pval_kw.pt.pos.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

%% PH

% - Inhibition
%Perform Kruskal-Wallis (KW) test to compare the distributions for PH
data_ph_neg = [ph.neg.small'; ph.neg.med'; ph.neg.large']; %Make a vector with the data for KW test
group_ph_neg = [repmat({'Small'},length(ph.neg.small),1); repmat({'Medium'},length(ph.neg.med),1); repmat({'Large'},length(ph.neg.large),1)]; %Make a vector with the group labels

[pval_kw.ph.neg.single, pval_kw.ph.neg.table, pval_kw.ph.neg.stats] = kruskalwallis(data_ph_neg, group_ph_neg, display_results);
pval_kw.ph.neg.multiple = multcompare(pval_kw.ph.neg.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Excitation
%Perform Kruskal-Wallis (KW) test to compare the distributions for PH
data_ph_pos = [ph.pos.small'; ph.pos.med'; ph.pos.large']; %Make a vector with the data for KW test
group_ph_pos = [repmat({'Small'},length(ph.pos.small),1); repmat({'Medium'},length(ph.pos.med),1); repmat({'Large'},length(ph.pos.large),1)]; %Make a vector with the group labels

[pval_kw.ph.pos.single, pval_kw.ph.pos.table, pval_kw.ph.pos.stats] = kruskalwallis(data_ph_pos, group_ph_pos, display_results);
pval_kw.ph.pos.multiple = multcompare(pval_kw.ph.pos.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);


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
    
    %% COM Inhibition large vs small with NPSP Scatter
    params_npsp_scatter.specific_name = 'Inhibitory LS COM for ';
    params_npsp_scatter.npsp_on = 1;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com.neg.large; params_npsp_scatter.val_small = com.neg.small; 
    params_npsp_scatter.p_val = p_val.com.neg.ls;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Inhibition large vs med with NPSP Scatter
    params_npsp_scatter.specific_name = 'Inhibitory LM COM for ';
    params_npsp_scatter.npsp_on = 1;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com.neg.large; params_npsp_scatter.val_small = com.neg.med; 
    params_npsp_scatter.p_val = p_val.com.neg.lm;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Inhibition med vs small with NPSP Scatter
    params_npsp_scatter.specific_name = 'Inhibitory MS COM for ';
    params_npsp_scatter.npsp_on = 1;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com.neg.med; params_npsp_scatter.val_small = com.neg.small; 
    params_npsp_scatter.p_val = p_val.com.neg.ms;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation large vs small with NPSP Scatter
    params_npsp_scatter.specific_name = 'Excitatory LS COM  ';
    params_npsp_scatter.val_big = com.pos.large; params_npsp_scatter.val_small = com.pos.small;
    params_npsp_scatter.p_val = p_val.com.pos.ls;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation large vs med with NPSP Scatter
    params_npsp_scatter.specific_name = 'Excitatory LS COM  ';
    params_npsp_scatter.val_big = com.pos.large; params_npsp_scatter.val_small = com.pos.med;
    params_npsp_scatter.p_val = p_val.com.pos.lm;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation med vs small with NPSP Scatter
    params_npsp_scatter.specific_name = 'Excitatory LS COM  ';
    params_npsp_scatter.val_big = com.pos.med; params_npsp_scatter.val_small = com.pos.small;
    params_npsp_scatter.p_val = p_val.com.pos.ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Excitation and Inhibiton large vs small Scatter No NPSP
    params_ei_scatter.specific_name = 'Inhibitory and excitatory LS com no npsp for ';
    params_ei_scatter.val_exc_big = com.pos.large; params_ei_scatter.val_exc_small = com.pos.small; params_ei_scatter.val_inh_big = com.neg.large; params_ei_scatter.val_inh_small = com.neg.small;
    params_ei_scatter.p_val_exc = p_val.com.pos.ls; params_ei_scatter.p_val_inh = p_val.com.neg.ls;
    params_ei_scatter.lim_val = lim_val_ms;
    plot_ei_scatter(params_ei_scatter);
    
    %% COM Excitation and Inhibiton large vs med Scatter No NPSP
    params_ei_scatter.specific_name = 'Inhibitory and excitatory LM com no npsp for ';
    params_ei_scatter.val_exc_big = com.pos.large; params_ei_scatter.val_exc_small = com.pos.med; params_ei_scatter.val_inh_big = com.neg.large; params_ei_scatter.val_inh_small = com.neg.med;
    params_ei_scatter.p_val_exc = p_val.com.pos.lm; params_ei_scatter.p_val_inh = p_val.com.neg.lm;
    params_ei_scatter.lim_val = lim_val_ms;
    plot_ei_scatter(params_ei_scatter);
    
    %% COM Excitation and Inhibiton med vs small Scatter No NPSP
    params_ei_scatter.specific_name = 'Inhibitory and excitatory MS com no npsp for ';
    params_ei_scatter.val_exc_big = com.pos.med; params_ei_scatter.val_exc_small = com.pos.small; params_ei_scatter.val_inh_big = com.neg.med; params_ei_scatter.val_inh_small = com.neg.small;
    params_ei_scatter.p_val_exc = p_val.com.pos.ms; params_ei_scatter.p_val_inh = p_val.com.neg.ms;
    params_ei_scatter.lim_val = lim_val_ms;
    plot_ei_scatter(params_ei_scatter);
    
    %% IE Ratio with NPSP large vs small Scatter
    params_npsp_scatter.specific_name = 'IE Ratio LS for ';
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie_ls;
    params_npsp_scatter.lim_val = 2;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio with NPSP large vs med Scatter
    params_npsp_scatter.specific_name = 'IE Ratio LM for ';
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_med;
    params_npsp_scatter.p_val = p_val_ie_lm;
    params_npsp_scatter.lim_val = 2;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio with NPSP med vs small Scatter
    params_npsp_scatter.specific_name = 'IE Ratio MS for ';
    params_npsp_scatter.val_big = ie_ratio_med; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie_ms;
    params_npsp_scatter.lim_val = 2;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio large vs small Scatter
    params_npsp_scatter.specific_name = 'IE Ratio LS no NPSP for test ';
    params_npsp_scatter.npsp_on = 0;
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie_ls;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio large vs med Scatter
    params_npsp_scatter.specific_name = 'IE Ratio LM no NPSP for test ';
    params_npsp_scatter.npsp_on = 0;
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_med;
    params_npsp_scatter.p_val = p_val_ie_lm;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% IE Ratio med vs small Scatter
    params_npsp_scatter.specific_name = 'IE Ratio MS no NPSP for test ';
    params_npsp_scatter.npsp_on = 0;
    params_npsp_scatter.val_big = ie_ratio_med; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% Total Excitation and Inhibiton large vs small Scatter
    params_ei_scatter.specific_name = 'E and I total LS no NPSP for test ';
    params_ei_scatter.val_exc_big = total_exc_big; params_ei_scatter.val_exc_small = total_exc_small; params_ei_scatter.val_inh_big = total_inh_big; params_ei_scatter.val_inh_small = total_inh_small;
    params_ei_scatter.p_val_exc = p_val_total_exc_ls; params_ei_scatter.p_val_inh = p_val_total_inh_ls;
    params_ei_scatter.lim_val = 0.3;
    plot_ei_scatter(params_ei_scatter);
    
    %% Total Excitation and Inhibiton large vs med Scatter
    params_ei_scatter.specific_name = 'E and I total LM no NPSP for test ';
    params_ei_scatter.val_exc_big = total_exc_big; params_ei_scatter.val_exc_small = total_exc_med; params_ei_scatter.val_inh_big = total_inh_big; params_ei_scatter.val_inh_small = total_inh_med;
    params_ei_scatter.p_val_exc = p_val_total_exc_lm; params_ei_scatter.p_val_inh = p_val_total_inh_lm;
    params_ei_scatter.lim_val = 0.3;
    plot_ei_scatter(params_ei_scatter);
    
    %% Total Excitation and Inhibiton med vs small Scatter
    params_ei_scatter.specific_name = 'E and I total MS no NPSP for test ';
    params_ei_scatter.val_exc_big = total_exc_med; params_ei_scatter.val_exc_small = total_exc_small; params_ei_scatter.val_inh_big = total_inh_med; params_ei_scatter.val_inh_small = total_inh_small;
    params_ei_scatter.p_val_exc = p_val_total_exc_ms; params_ei_scatter.p_val_inh = p_val_total_inh_ms;
    params_ei_scatter.lim_val = 0.3;
    plot_ei_scatter(params_ei_scatter);
    
    %% PH Excitation and Inhibiton large vs small Scatter
    params_ei_scatter.specific_name = 'PH E and I LS no NPSP for ';
    params_ei_scatter.val_exc_big = ph.pos.large; params_ei_scatter.val_exc_small = ph.pos.small; params_ei_scatter.val_inh_big = ph.neg.large; params_ei_scatter.val_inh_small = ph.neg.small;
    params_ei_scatter.p_val_exc = p_val.ph.pos.ls; params_ei_scatter.p_val_inh = p_val.ph.neg.ls;
    params_ei_scatter.lim_val = 0.03;
    plot_ei_scatter(params_ei_scatter);
    
    %% PH Excitation and Inhibiton large vs med Scatter
    params_ei_scatter.specific_name = 'PH E and I LM no NPSP for ';
    params_ei_scatter.val_exc_big = ph.pos.large; params_ei_scatter.val_exc_small = ph.pos.med; params_ei_scatter.val_inh_big = ph.neg.large; params_ei_scatter.val_inh_small = ph.neg.med;
    params_ei_scatter.p_val_exc = p_val.ph.pos.lm; params_ei_scatter.p_val_inh = p_val.ph.neg.lm;
    params_ei_scatter.lim_val = 0.03;
    plot_ei_scatter(params_ei_scatter);
    
    %% PH Excitation and Inhibiton med vs small Scatter
    params_ei_scatter.specific_name = 'PH E and I MS no NPSP for ';
    params_ei_scatter.val_exc_big = ph.pos.med; params_ei_scatter.val_exc_small = ph.pos.small; params_ei_scatter.val_inh_big = ph.neg.med; params_ei_scatter.val_inh_small = ph.neg.small;
    params_ei_scatter.p_val_exc = p_val.ph.pos.ms; params_ei_scatter.p_val_inh = p_val.ph.neg.ms;
    params_ei_scatter.lim_val = 0.03;
    plot_ei_scatter(params_ei_scatter);
    
    %% PT Excitation and Inhibiton large vs small Scatter
    params_ei_scatter.specific_name = 'PT E and I LS no NPSP for test ';
    g_noise = 5*randn(1,length(pt.neg.small));
    params_ei_scatter.val_exc_big = pt.pos.large+g_noise; params_ei_scatter.val_exc_small = pt.pos.small+g_noise; params_ei_scatter.val_inh_big = pt.neg.large+g_noise; params_ei_scatter.val_inh_small = pt.neg.small+g_noise;
    params_ei_scatter.p_val_exc = p_val.pt.pos.ls; params_ei_scatter.p_val_inh = p_val.pt.neg.ls;
    params_ei_scatter.lim_val = 200;
    plot_ei_scatter(params_ei_scatter);
    
    %% PT Excitation and Inhibiton large vs med Scatter
    params_ei_scatter.specific_name = 'PT E and I LM no NPSP for test ';
    g_noise = 5*randn(1,length(pt.neg.med));
    params_ei_scatter.val_exc_big = pt.pos.large+g_noise; params_ei_scatter.val_exc_small = pt.pos.med+g_noise; params_ei_scatter.val_inh_big = pt.neg.large+g_noise; params_ei_scatter.val_inh_small = pt.neg.med+g_noise;
    params_ei_scatter.p_val_exc = p_val.pt.pos.lm; params_ei_scatter.p_val_inh = p_val.pt.neg.lm;
    params_ei_scatter.lim_val = 200;
    plot_ei_scatter(params_ei_scatter);
    
    %% PT Excitation and Inhibiton large vs med Scatter
    params_ei_scatter.specific_name = 'PT E and I MS no NPSP for test ';
    g_noise = 5*randn(1,length(pt.neg.small));
    params_ei_scatter.val_exc_big = pt.pos.med+g_noise; params_ei_scatter.val_exc_small = pt.pos.small+g_noise; params_ei_scatter.val_inh_big = pt.neg.med+g_noise; params_ei_scatter.val_inh_small = pt.neg.small+g_noise;
    params_ei_scatter.p_val_exc = p_val.pt.pos.ms; params_ei_scatter.p_val_inh = p_val.pt.neg.ms;
    params_ei_scatter.lim_val = 200;
    plot_ei_scatter(params_ei_scatter);
    
    %% COM Excitation and Inhibiton large vs small Histogram
    params_his_ei.specific_name = ' COM LS ms change for ';
    params_his_ei.val_exc_big = com.pos.large; params_his_ei.val_exc_small = com.pos.small; params_his_ei.val_inh_big = com.neg.large; params_his_ei.val_inh_small = com.neg.small;  
    params_his_ei.p_val_exc = p_val.com.pos.ls; params_his_ei.p_val_inh = p_val.com.neg.ls;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% COM Excitation and Inhibiton large vs med Histogram
    params_his_ei.specific_name = ' COM LM ms change for ';
    params_his_ei.val_exc_big = com.pos.large; params_his_ei.val_exc_small = com.pos.med; params_his_ei.val_inh_big = com.neg.large; params_his_ei.val_inh_small = com.neg.med;  
    params_his_ei.p_val_exc = p_val.com.pos.lm; params_his_ei.p_val_inh = p_val.com.neg.lm;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);

    %% COM Excitation and Inhibiton med vs small Histogram
    params_his_ei.specific_name = ' COM MS ms change for ';
    params_his_ei.val_exc_big = com.pos.med; params_his_ei.val_exc_small = com.pos.small; params_his_ei.val_inh_big = com.neg.med; params_his_ei.val_inh_small = com.neg.small;  
    params_his_ei.p_val_exc = p_val.com.pos.ms; params_his_ei.p_val_inh = p_val.com.neg.ms;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% COM change ix Excitation and Inhibiton large vs small Histogram
    params_his_ei.units = 'change ix';
    params_his_ei.specific_name = ' COM LS change ix for test ';
    params_his_ei.val_exc_big = com.pos.large; params_his_ei.val_exc_small = com.pos.small; params_his_ei.val_inh_big = com.neg.large; params_his_ei.val_inh_small = com.neg.small;
    params_his_ei.p_val_exc = p_val.com.pos.ls; params_his_ei.p_val_inh = p_val.com.neg.ls;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% COM change ix Excitation and Inhibiton large vs med Histogram
    params_his_ei.units = 'change ix';
    params_his_ei.specific_name = ' COM LM change ix for test ';
    params_his_ei.val_exc_big = com.pos.large; params_his_ei.val_exc_small = com.pos.med; params_his_ei.val_inh_big = com.neg.large; params_his_ei.val_inh_small = com.neg.med;
    params_his_ei.p_val_exc = p_val.com.pos.lm; params_his_ei.p_val_inh = p_val.com.neg.lm;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% COM change ix Excitation and Inhibiton med vs small Histogram
    params_his_ei.units = 'change ix';
    params_his_ei.specific_name = ' COM MS change ix for test ';
    params_his_ei.val_exc_big = com.pos.med; params_his_ei.val_exc_small = com.pos.small; params_his_ei.val_inh_big = com.neg.med; params_his_ei.val_inh_small = com.neg.small;
    params_his_ei.p_val_exc = p_val.com.pos.ms; params_his_ei.p_val_inh = p_val.com.neg.ms;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% Total Excitation and Inhibiton large vs small room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' of Total amount change for LS ';
    params_his_ei.val_exc_big = total_exc_big; params_his_ei.val_exc_small = total_exc_small; params_his_ei.val_inh_big = total_inh_big; params_his_ei.val_inh_small = total_inh_small;
    params_his_ei.p_val_exc = p_val_total_exc_ls; params_his_ei.p_val_inh = p_val_total_inh_ls;
    params_his_ei.his_spacing = 0.1; params_his_ei.lim_val = 4;
    plot_ei_histogram(params_his_ei);
    
    %% Total Excitation and Inhibiton large vs med room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' of Total amount change for LM ';
    params_his_ei.val_exc_big = total_exc_big; params_his_ei.val_exc_small = total_exc_med; params_his_ei.val_inh_big = total_inh_big; params_his_ei.val_inh_small = total_inh_med;
    params_his_ei.p_val_exc = p_val_total_exc_lm; params_his_ei.p_val_inh = p_val_total_inh_lm;
    params_his_ei.his_spacing = 0.1; params_his_ei.lim_val = 4;
    plot_ei_histogram(params_his_ei);
    
    %% Total Excitation and Inhibiton med vs small room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' of Total amount change for MS ';
    params_his_ei.val_exc_big = total_exc_med; params_his_ei.val_exc_small = total_exc_small; params_his_ei.val_inh_big = total_inh_med; params_his_ei.val_inh_small = total_inh_small;
    params_his_ei.p_val_exc = p_val_total_exc_ms; params_his_ei.p_val_inh = p_val_total_inh_ms;
    params_his_ei.his_spacing = 0.1; params_his_ei.lim_val = 4;
    plot_ei_histogram(params_his_ei);
    
    %% PH Excitation and Inhibiton for large vs small room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' PH LS E and I change for ';
    params_his_ei.val_exc_big = ph.pos.large; params_his_ei.val_exc_small = ph.pos.small; params_his_ei.val_inh_big = ph.neg.large; params_his_ei.val_inh_small = ph.neg.small;
    params_his_ei.p_val_exc = p_val.ph.pos.ls; params_his_ei.p_val_inh = p_val.ph.neg.ls;
    params_his_ei.his_spacing = 0.125; params_his_ei.lim_val = 2.5;
    plot_ei_histogram(params_his_ei);
    
    %% PH Excitation and Inhibiton for large vs med room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' PH LM E and I change for ';
    params_his_ei.val_exc_big = ph.pos.large; params_his_ei.val_exc_small = ph.pos.med; params_his_ei.val_inh_big = ph.neg.large; params_his_ei.val_inh_small = ph.neg.med;
    params_his_ei.p_val_exc = p_val.ph.pos.lm; params_his_ei.p_val_inh = p_val.ph.neg.lm;
    params_his_ei.his_spacing = 0.125; params_his_ei.lim_val = 2.5;
    plot_ei_histogram(params_his_ei);
    
    %% PH Excitation and Inhibiton for med vs small room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' PH MS E and I change for ';
    params_his_ei.val_exc_big = ph.pos.med; params_his_ei.val_exc_small = ph.pos.small; params_his_ei.val_inh_big = ph.neg.med; params_his_ei.val_inh_small = ph.neg.small;
    params_his_ei.p_val_exc = p_val.ph.pos.ms; params_his_ei.p_val_inh = p_val.ph.neg.ms;
    params_his_ei.his_spacing = 0.125; params_his_ei.lim_val = 2.5;
    plot_ei_histogram(params_his_ei);
    
    %% PT Excitation and Inhibiton for large vs small room Histogram
    params_his_ei.units = 'normal';
    params_his_ei.specific_name = ' PT LS E and I change for ';
    params_his_ei.val_exc_big = pt.pos.large; params_his_ei.val_exc_small = pt.pos.small; params_his_ei.val_inh_big = pt.neg.large; params_his_ei.val_inh_small = pt.neg.small;
    params_his_ei.p_val_exc = p_val.pt.pos.ls; params_his_ei.p_val_inh = p_val.pt.neg.ls;
    params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
    plot_ei_histogram(params_his_ei);
    
    %% PT Excitation and Inhibiton for large vs med room Histogram
    params_his_ei.units = 'normal';
    params_his_ei.specific_name = ' PT LM E and I change for ';
    params_his_ei.val_exc_big = pt.pos.large; params_his_ei.val_exc_small = pt.pos.med; params_his_ei.val_inh_big = pt.neg.large; params_his_ei.val_inh_small = pt.neg.med;
    params_his_ei.p_val_exc = p_val.pt.pos.lm; params_his_ei.p_val_inh = p_val.pt.neg.lm;
    params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
    plot_ei_histogram(params_his_ei);
    
    %% PT Excitation and Inhibiton for med vs small room Histogram
    params_his_ei.units = 'normal';
    params_his_ei.specific_name = ' PT MS E and I change for ';
    params_his_ei.val_exc_big = pt.pos.med; params_his_ei.val_exc_small = pt.pos.small; params_his_ei.val_inh_big = pt.neg.med; params_his_ei.val_inh_small = pt.neg.small;
    params_his_ei.p_val_exc = p_val.pt.pos.ms; params_his_ei.p_val_inh = p_val.pt.neg.ms;
    params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
    plot_ei_histogram(params_his_ei);
      
    %% IE Ratio change ix large vs small room Histogram
    params_norm_his.specific_name = ' IE Ratio LS of com change ix for ';
    params_norm_his.units = 'change ix'; 
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie_ls; 
    params_norm_his.his_spacing = 2.5; params_norm_his.lim_val = 70;
    plot_norm_histogram(params_norm_his);
    
    %% IE Ratio change ix large vs med room Histogram
    params_norm_his.specific_name = ' IE Ratio LM of com change ix for ';
    params_norm_his.units = 'change ix';
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_med;
    params_norm_his.p_val = p_val_ie_lm;
    params_norm_his.his_spacing = 2.5; params_norm_his.lim_val = 70;
    plot_norm_histogram(params_norm_his);
    
    %% IE Ratio change ix med vs small room Histogram
    params_norm_his.specific_name = ' IE Ratio MS of com change ix for ';
    params_norm_his.units = 'change ix';
    params_norm_his.big_val = ie_ratio_med; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie_lm;
    params_norm_his.his_spacing = 2.5; params_norm_his.lim_val = 70;
    plot_norm_histogram(params_norm_his);
    
    %% IE ratio for large vs small room Histogram
    params_norm_his.specific_name = ' of IE LS Ratio AU change for test ';
    params_norm_his.units = 'normal';
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie_ls;
    params_norm_his.his_spacing = 0.05; params_norm_his.lim_val = 1;
    plot_norm_histogram(params_norm_his);
    
    %% IE ratio for large vs med room Histogram
    params_norm_his.specific_name = ' of IE LM Ratio AU change for test ';
    params_norm_his.units = 'normal';
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_med;
    params_norm_his.p_val = p_val_ie_lm;
    params_norm_his.his_spacing = 0.05; params_norm_his.lim_val = 1;
    plot_norm_histogram(params_norm_his);
    
    %% IE ratio for med vs small room Histogram
    params_norm_his.specific_name = ' of IE MS Ratio AU change for test ';
    params_norm_his.units = 'normal';
    params_norm_his.big_val = ie_ratio_med; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie_ms;
    params_norm_his.his_spacing = 0.05; params_norm_his.lim_val = 1;
    plot_norm_histogram(params_norm_his);
    
end

%% Violin plot of COM inhibition
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
bandwidth = 10;
int_names = {'Small','Medium','Large'};
violins = violinplot(com.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('COM inhibition');
ylabel(['COM (ms)']);
ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_COM_inhibition.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(com.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('COM inhibition');
ylabel(['COM (ms)']);
ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_COM_inhibition_ind.svg');
saveas(gcf, save_name);
close;

%% Violin plot of COM excitation
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Medium','Large'};
violins = violinplot(com.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('COM excitation');
ylabel(['COM (ms)']);
ylim([0 110]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_COM_excitation.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(com.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('COM excitation');
ylabel(['COM (ms)']);
ylim([0 110]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_COM_excitation_ind.svg');
saveas(gcf, save_name);
close;

%% Violin plot of PT inhibition
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Medium','Large'};
violins = violinplot(pt.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('PT inhibition');
ylabel(['PT (ms)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PT_inhibition.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(pt.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('PT inhibition');
ylabel(['PT (ms)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PT_inhibition_ind.svg');
saveas(gcf, save_name);
close;

%% Violin plot of PT excitation
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Medium','Large'};
violins = violinplot(pt.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width,'Bandwidth', bandwidth, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('PT excitation');
ylabel(['PT (ms)']);
% ylim([0 110]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PT_excitation.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(pt.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'Bandwidth', bandwidth, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('PT excitation');
ylabel(['PT (ms)']);
% ylim([0 110]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PT_excitation_ind.svg');
saveas(gcf, save_name);
close;

%% Violin plot of PH inhibition
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Medium','Large'};
violins = violinplot(ph.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('PH inhibition');
ylabel(['PH (AU)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PH_inhibition.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(ph.neg, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_med_color;
violins(3).ViolinColor = inh_big_color;
title('PH inhibition');
ylabel(['PH (AU)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PH_inhibition_ind.svg');
saveas(gcf, save_name);
close;

%% Violin plot of PH excitation
font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Medium','Large'};
violins = violinplot(ph.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('PH excitation');
ylabel(['PH (AU)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PH_excitation.svg');
saveas(gcf, save_name);
close;

%Make another violin plot with the individual points
figure('units','normalized','outerposition',[0 0 1 1]);
violins = violinplot(ph.pos, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_med_color;
violins(3).ViolinColor = exc_big_color;
title('PH excitation');
ylabel(['PH (AU)']);
% ylim([50 150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_bf, 'Violin_PH_excitation_ind.svg');
saveas(gcf, save_name);
close;

%% Box plot COM inhibition
axis_sz = 40;
n_points = length(com.neg.small);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('COM inhibition');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),com.neg.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),com.neg.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),com.neg.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [com.neg.small(n), com.neg.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [com.neg.med(n), com.neg.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([com.neg.small',com.neg.med',com.neg.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'COM_inhibition_box.svg');
saveas(gcf, full_name);
close all;

%% Box plot COM excitation
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('COM excitation');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),com.pos.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),com.pos.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),com.pos.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [com.pos.small(n), com.pos.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [com.pos.med(n), com.pos.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([com.pos.small',com.pos.med',com.pos.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'COM_excitation_box.svg');
saveas(gcf, full_name);
close all;

%% Box plot PT inhibition
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('PT inhibition');
ylabel('PT [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),pt.neg.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),pt.neg.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),pt.neg.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [pt.neg.small(n), pt.neg.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [pt.neg.med(n), pt.neg.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([pt.neg.small', pt.neg.med', pt.neg.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'PT_inhibition_box.svg');
saveas(gcf, full_name);
close all;

%% Box plot PT excitation
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('PT excitation');
ylabel('PT [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),pt.pos.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),pt.pos.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),pt.pos.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [pt.pos.small(n), pt.pos.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [pt.pos.med(n), pt.pos.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([pt.pos.small', pt.pos.med', pt.pos.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'PT_excitation_box.svg');
saveas(gcf, full_name);
close all;

%% Box plot PH inhibition
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('PH inhibition');
ylabel('PH [AU]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),ph.neg.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),ph.neg.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),ph.neg.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [ph.neg.small(n), ph.neg.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [ph.neg.med(n), ph.neg.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([ph.neg.small', ph.neg.med', ph.neg.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'PH_inhibition_box.svg');
saveas(gcf, full_name);
close all;

%% Box plot PH excitation
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('PH excitation');
ylabel('PH [AU]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),ph.pos.small,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
scatter(2*ones(n_points,1),ph.pos.med,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4]);
scatter(3*ones(n_points,1),ph.pos.large,'MarkerEdgeColor','k','MarkerFaceColor','k');


for n = 1:n_points
    plot([1,2], [ph.pos.small(n), ph.pos.med(n)],'Color',[0.6 0.6 0.6]);
    plot([2,3], [ph.pos.med(n), ph.pos.large(n)],'Color',[0.2 0.2 0.2]);
end

h = boxplot([ph.pos.small', ph.pos.med', ph.pos.large'], 'Notch', 'on', 'Labels', {'Small', 'Medium', 'Large'});
set(h,'LineWidth',4);
hold off;

full_name = fullfile(save_dir_bf,'PH_excitation_box.svg');
saveas(gcf, full_name);
close all;