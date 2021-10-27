%% Input variables

kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/Kernel_fits/kfolds/Sliding/perfreq_noneuro/ridge/10ms/200ms/window_4.0_s/stride_0.500_s';
suffix = 'Slide';

get_neurons = 'all';
NPSP_th = 40;
plot_f_only = 0;
freq_up_bound = 19000;
freq_down_bound = 400;
norm = true;

%% Params
chop_ms = 190; %The analysis window for the COM
hist_type = 'bar'; %The type of the histogram
calc_method = 'average'; %How to estimate important values
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
legend_font_sz = 32;
sz = 55;
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


r_type{1} = 'sl';
r_type{2} = 'ls';
n_rooms = length(r_type);
%% Load the data
load(fullfile(kernel_dir,'info'),'info');
temp_files = dir([kernel_dir,'/*.mat']);
temp = load(fullfile(temp_files(1).folder,temp_files(1).name));
model = temp.kernel{1}.model;

%Define the window, stride and center of the window at each tiem step
total_s = 8;
window_s = diff(temp.kernel{1}.params.sl(1).small);
stride_s = temp.kernel{2}.params.sl(1).small(1) - temp.kernel{1}.params.sl(1).small(1);
window_center_s = [window_s/2:stride_s:total_s-window_s/2];

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
n_steps = length(kernel); %Get the number of time steps
freqs = fliplr(kernel{1}.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_h = kernel{1}.n_h;
dt_ms = round(kernel{1}.dt_ms);
chop_ix = round(chop_ms/dt_ms);
bf_window_ix = round(bf_window_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;
fprintf('== Calcuating com and bf ==\n');tic;
for k = 1:n_clust
    sprintf('== Cluster %0.f/%0.f ==\n',k,n_clust);
    for r = 1:n_rooms
        room = r_type{r};
        for j = 1:n_steps
            switch model
                case {'sep','sep_kh'}
                    [~,ix] = max(kernels{k}.(room).k_f);
                    bf(k).(room) = freqs(ix); %Find the corresponding frequency
                    k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                    
                case {'ridge','lasso','elastic'}
                    k_fh = fliplr(kernels{k}{j}.(room).main{end}.k_fh);
                    k_fh = k_fh(:,1:chop_ix);
                    k_fhn(:,:,k).(room){j} = k_fh;
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
                    max_exc(k).(room)(j) = max(k_h_pos_temp(:));
                    max_inh(k).(room)(j) = max(k_h_neg_temp(:));
                    %Total Excitation and Inhibition
                    total_exc(k).(room)(j) = sum(k_h_pos_temp(:));
                    total_inh(k).(room)(j) = sum(k_h_neg_temp(:));
                    %IE ratio
                    ie_ratio(k).(room)(j) = sum(k_h_neg_temp(:))/sum(k_h_pos_temp(:));
                    %Peak Time (PT) values
                    [~,max_ix_exc] = max(k_h_pos_temp(:));
                    [~, max_ix_col_exc] = ind2sub(size(k_h_pos_temp),max_ix_exc);
                    max_exc_time(k).(room)(j) = h(max_ix_col_exc);
                    
                    [~,max_ix_inh] = max(k_h_neg_temp(:));
                    [~, max_ix_col_inh] = ind2sub(size(k_h_neg_temp),max_ix_inh);
                    max_inh_time(k).(room)(j) = h(max_ix_col_inh);
                    
                    switch bf_method
                        case 'max'
                            k_f = max(k_fh_pos,[],2); %Take the max across history steps
                        case 'window_max'
                            k_f = max(k_fh_pos(:,1:bf_window_ix),[],2); %Take the max in a specified window
                        case 'window_mean'
                            k_f = mean(k_fh_pos(:,1:bf_window_ix),2); %Take the mean in a specified window
                    end
                    
                    [~,ix] = max(k_f); %Find the index of the max frequency
                    bf(k).(room)(j) = freqs(ix); %Find the corresponding frequency
            end
            
            %Get the center of mass (COM) values
            k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
            k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
            com_neg(k).(room)(j) = (k_h_neg*h); %Compute a weighted sum of all values
            com_pos(k).(room)(j) = (k_h_pos*h); %Compute a weighted sum of all values
        end
    end
    bf_mean(k) = 2^(log2(prod(bf(k).sl)*prod(bf(k).ls)));
    [~,ix_bf] = min(abs(bf_mean(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
    bf_closest(k) = freqs(ix_bf);
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Compute a mean STRF for each BF

% for  r = 1:n_rooms
%     room = r_type{r};
%     k_fhn_temp = reach(k_fhn,room);
%     for f = 1:length(freqs)
%         curr_freq = freqs(f);
%         ix = bf_closest == curr_freq;
%         n_bf(f) = sum(ix);
%         k_fh = mean(k_fhn_temp(:,:,ix),3);
%         k_fh_bf.(room){f} = k_fh;
%         k_fh_neg = abs(min(k_fh,0));
%         k_fh_pos = abs(max(k_fh,0));
%         %Get the IE ratio
%         ie_ratio_bf(f).(room) = sum(k_fh_neg(:))/sum(k_fh_pos(:)); %Plot the max/sum of excitation inhibiton changes
%         %Make sure don't divide by zero
%         k_h_neg = mean(k_fh_neg);
%         k_h_pos = mean(k_fh_pos);
%         %Get the com
%         k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
%         k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
%         com_neg_bf(f).(room) = (k_h_neg*h); %Compute a weighted sum of all values
%         com_pos_bf(f).(room) = (k_h_pos*h); %Compute a weighted sum of all values
%     end
% end
% params.freqs = freqs;
% params.n_h = n_h-1;
% params.dt_ms = dt_ms;
% params.h_max_ms = chop_ms;
% params.n_bf = n_bf;
% params.font_type = font_type;
% params.exc_small_color = exc_small_color; params.exc_big_color = exc_big_color;
% params.inh_small_color = inh_small_color; params.inh_big_color = inh_big_color;
% params.name = 'neurons';
% 
% %% Plot a mean STRF across all frequencies using pre-averaged kernels
% params.type = 'bf';
% params.save_dir = save_dir;
% plot_mean_strf_switch(k_fh_bf,params);
% 
% %% Plot a mean STRF across all frequencies using averaging across all kernels
% params.type = 'raw';
% plot_mean_strf_switch(k_fhn,params);

%% Plot STRF and temporal profiles for each bf separately 
% params.save_dir = save_dir_paper_fig1_sup;
% plot_ind_bf(k_fh_bf,params);

%% Stats
%COM measures
com_neg_res.sl = reshape(reach(com_neg,'sl'),n_steps,n_clust)'; %Get all time steps for this room for every cluster and reshape to dimensions of n_clusters x n_time_steps
com_pos_res.sl = reshape(reach(com_pos,'sl'),n_steps,n_clust)';
com_neg_res.ls = reshape(reach(com_neg,'ls'),n_steps,n_clust)';
com_pos_res.ls = reshape(reach(com_pos,'ls'),n_steps,n_clust)';

%IE ratio
ie_ratio_res.sl = reshape(reach(ie_ratio,'sl'),n_steps,n_clust)';
ie_ratio_res.ls = reshape(reach(ie_ratio,'ls'),n_steps,n_clust)';

%Total Inhibiton and excitation
total_inh_res.sl = reshape(reach(total_inh,'sl'),n_steps,n_clust)';
total_exc_res.sl = reshape(reach(total_exc,'sl'),n_steps,n_clust)';
total_inh_res.ls = reshape(reach(total_inh,'ls'),n_steps,n_clust)';
total_exc_res.ls = reshape(reach(total_exc,'ls'),n_steps,n_clust)';

%Inhibition Max value
max_inh_res.sl = reshape(reach(max_inh,'sl'),n_steps,n_clust)';
max_exc_res.sl = reshape(reach(max_exc,'sl'),n_steps,n_clust)';
max_inh_res.ls = reshape(reach(max_inh,'ls'),n_steps,n_clust)';
max_exc_res.ls = reshape(reach(max_exc,'ls'),n_steps,n_clust)';

%Time of Max Inhibition
max_inh_time_res.sl = reshape(reach(max_inh_time,'sl'),n_steps,n_clust)';
max_exc_time_res.sl = reshape(reach(max_exc_time,'sl'),n_steps,n_clust)';
max_inh_time_res.ls = reshape(reach(max_inh_time,'ls'),n_steps,n_clust)';
max_exc_time_res.ls = reshape(reach(max_exc_time,'ls'),n_steps,n_clust)';

%% Optionally normalize
if norm
    %COM
    com_norm_neg = nanmean([com_neg_res.sl,com_neg_res.ls],2);
    com_norm_pos = nanmean([com_pos_res.sl,com_pos_res.ls],2);
    com_neg_res.sl = com_neg_res.sl - com_norm_neg;
    com_pos_res.sl = com_pos_res.sl - com_norm_pos;
    com_neg_res.ls = com_neg_res.ls - com_norm_neg;
    com_pos_res.ls = com_pos_res.ls - com_norm_pos;
    
    %PT
    pt_norm_neg = nanmean([max_inh_time_res.sl,max_inh_time_res.ls],2);
    pt_norm_pos = nanmean([max_exc_time_res.sl,max_exc_time_res.ls],2);
    max_inh_time_res.sl = max_inh_time_res.sl - pt_norm_neg;
    max_exc_time_res.sl = max_exc_time_res.sl - pt_norm_pos;
    max_inh_time_res.ls = max_inh_time_res.ls - pt_norm_neg;
    max_exc_time_res.ls = max_exc_time_res.ls - pt_norm_pos;
    
    %PH
    ph_norm_neg = nanmean([max_inh_res.sl,max_inh_res.ls],2);
    ph_norm_pos = nanmean([max_exc_res.sl,max_exc_res.ls],2);
    max_inh_res.sl = max_inh_res.sl - ph_norm_neg;
    max_exc_res.sl = max_exc_res.sl - ph_norm_pos;
    max_inh_res.ls = max_inh_res.ls - ph_norm_neg;
    max_exc_res.ls = max_exc_res.ls - ph_norm_pos;
end

% [p_val_neg_sl2_sl1,~,~] = signrank(com_sl2_neg,com_sl1_neg);
% [p_val_pos_sl2_sl1,~,~] = signrank(com_sl2_pos,com_sl1_pos);
% 
% [p_val_neg_ls2_ls1,~,~] = signrank(com_ls2_neg,com_ls1_neg);
% [p_val_pos_ls2_ls1,~,~] = signrank(com_ls2_pos,com_ls1_pos);
% 
% [p_val_neg_sl2_ls2,~,~] = signrank(com_sl2_neg,com_ls2_neg);
% [p_val_pos_sl2_ls2,~,~] = signrank(com_sl2_pos,com_ls2_pos);
% [p_val_ie,~,~] = signrank(ie_ratio_sl2,ie_ratio_sl1);
% [p_val_total_inh,~,~] = signrank(total_inh_sl2,total_inh_sl1);
% [p_val_max_inh,~,~] = signrank(max_inh_sl2,max_inh_sl1);
% [p_val_max_inh_time,~,~] = signrank(max_inh_time_sl2,max_inh_time_sl1);
% [p_val_total_exc,~,~] = signrank(total_exc_sl2,total_exc_sl1);
% [p_val_max_exc,~,~] = signrank(max_exc_sl2,max_exc_sl1);
% [p_val_max_exc_time,~,~] = signrank(max_exc_time_sl2,max_exc_time_sl1);

% com_sl_bf_neg = reach(com_neg_bf,'sl');
% com_sl_bf_pos = reach(com_pos_bf,'sl');
% com_ls_bf_neg = reach(com_neg_bf,'ls');
% com_ls_bf_pos = reach(com_pos_bf,'ls');
%% Plotting stuff
r_type_plot{1} = 'SL'; r_type_plot{2} = 'LS'; 
count = 0;
for r = 1:n_rooms
    for ii = 1:n_steps
        count = count+1;
        plot_lbl{count} = [r_type_plot{r},' ',num2str(window_center_s(ii),'%1.1f')];
    end
end
%% Plot box plot of COM inhibition
axis_sz = 25;
n_points = size(com_neg_res.sl,1);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory COM');
ylabel('COM [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;

count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),com_neg_res.(room)(:,ii),'k');
    end
end

boxplot([com_neg_res.sl, com_neg_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_COM_box.svg');
saveas(gcf, full_name);
close all;

%% Plot box plot of COM excitation
axis_sz = 25;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
boxplot([com_pos_res.sl, com_pos_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
title('Switching stimuli excitatory COM');
ylabel('COM [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');

hold on;
count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),com_pos_res.(room)(:,ii),'k');
    end
end

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_COM_box.svg');
saveas(gcf, full_name);
close all;

%% Plot box plot of PT inhibition
axis_sz = 25;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory PT');
ylabel('PT [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;

count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),max_inh_time_res.(room)(:,ii),'k');
    end
end

boxplot([max_inh_time_res.sl, max_inh_time_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_PT_box.svg');
saveas(gcf, full_name);
close all;

%% Plot box plot of PT excitation
axis_sz = 25;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
boxplot([max_exc_time_res.sl, max_exc_time_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
title('Switching stimuli excitatory PT');
ylabel('PT [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');

hold on;
count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),max_exc_time_res.(room)(:,ii),'k');
    end
end

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_PT_box.svg');
saveas(gcf, full_name);
close all;

%% Plot box plot of PH inhibition
axis_sz = 25;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory PH');
ylabel('PH [AU]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;

count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),max_inh_res.(room)(:,ii),'k');
    end
end

boxplot([max_inh_res.sl, max_inh_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_PH_box.svg');
saveas(gcf, full_name);
close all;

%% Plot box plot of PH excitation
axis_sz = 25;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
boxplot([max_exc_res.sl, max_exc_res.ls], 'Notch', 'on', 'Labels', plot_lbl);
title('Switching stimuli excitatory PH');
ylabel('PH [AU]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');

hold on;
count = 0;
for r = 1:n_rooms
    room = r_type{r};
    for ii = 1:n_steps
        count = count+1;
        scatter(count*ones(n_points,1),max_exc_res.(room)(:,ii),'k');
    end
end

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_PH_box.svg');
saveas(gcf, full_name);
close all;
%% Plot mean line and error bar COM of inhibition
axis_sz = 25;
lw = 3;
mean_com_neg = nanmean([com_neg_res.sl, com_neg_res.ls]);
std_com_neg = nanstd([com_neg_res.sl, com_neg_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory COM');
ylabel('COM [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1],mean_com_neg, std_com_neg, {'LineWidth',lw,'Color',inh_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_COM_line.svg');
saveas(gcf, full_name);
close all;

%% Plot mean line and error bar COM of excitation
axis_sz = 25;
lw = 3;
mean_com_pos = nanmean([com_pos_res.sl, com_pos_res.ls]);
std_com_pos = nanstd([com_pos_res.sl, com_pos_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli excitatory COM');
ylabel('COM [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1],mean_com_pos, std_com_pos, {'LineWidth',lw,'Color',exc_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_COM_line.svg');
saveas(gcf, full_name);
close all;

%% Plot mean line and error bar PT of inhibition
axis_sz = 25;
lw = 3;
mean_pt_neg = nanmean([max_inh_time_res.sl, max_inh_time_res.ls]);
std_pt_neg = nanstd([max_inh_time_res.sl, max_inh_time_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory PT');
ylabel('PT [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1],mean_pt_neg, std_pt_neg, {'LineWidth',lw,'Color',inh_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_PT_line.svg');
saveas(gcf, full_name);
close all;

%% Plot mean line and error bar PT of excitation
axis_sz = 25;
lw = 3;
mean_pt_pos = nanmean([max_exc_time_res.sl, max_exc_time_res.ls]);
std_pt_pos = nanstd([max_exc_time_res.sl, max_exc_time_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli excitatory PT');
ylabel('PT [ms]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1],mean_pt_pos, std_pt_pos, {'LineWidth',lw,'Color',exc_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_PT_line.svg');
saveas(gcf, full_name);
close all;

%% Plot mean line and error bar PH of inhibition
axis_sz = 25;
lw = 3;
mean_ph_neg = nanmean([max_inh_res.sl, max_inh_res.ls]);
std_ph_neg = nanstd([max_inh_res.sl, max_inh_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibitory PH');
ylabel('PH [AU]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1], mean_ph_neg, std_ph_neg, {'LineWidth',lw,'Color',inh_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_inhibition_PH_line.svg');
saveas(gcf, full_name);
close all;

%% Plot mean line and error bar PH of excitation
axis_sz = 25;
lw = 3;
mean_ph_pos = nanmean([max_exc_res.sl, max_exc_res.ls]);
std_ph_pos = nanstd([max_exc_res.sl, max_exc_res.ls])/sqrt(n_clust);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli excitatory PH');
ylabel('PH [AU]');
xlabel('STRF center [s]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
shadedErrorBar([0:1:2*n_steps-1], mean_ph_pos, std_ph_pos, {'LineWidth',lw,'Color',exc_big_color});
xticks([0:1:2*n_steps-1]);
xticklabels(plot_lbl);
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xline(n_steps-1+0.5,'k-',{'Transition point'},'LineWidth',3);
hold off;

full_name = fullfile(save_dir_bf,'Switch_excitation_PH_line.svg');
saveas(gcf, full_name);
close all;