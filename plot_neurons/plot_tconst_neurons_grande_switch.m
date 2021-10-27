%% Input variables

kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/Kernel_fits/kfolds/0.5_4_4.5_8s_version/perfreq_noneuro/ridge/10ms/200ms';
suffix = 'Neuronal_data';

get_neurons = 'all';
NPSP_th = 40;
CCnorm_on = 0;
CCnorm_th = 0.001;
plot_f_only = 0;
freq_up_bound = 17000;
freq_down_bound = 700;

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


r_type{1} = 'sl1';
r_type{2} = 'sl2';
r_type{3} = 'ls1';
r_type{4} = 'ls2';
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

ix_npsp = info.NPSP<NPSP_th; %Select all the neurons below certain NPSP


if CCnorm_on
    ix_CCnorm = info.CCnorm_mean>CCnorm_th & ~isnan(info.CCnorm_mean); %Select all the neurons that have at least a certain CCnorm
    ix = ix_qualia & ix_npsp & ix_CCnorm; %Find the intersection of the two
    CCnorm_mean = info.CCnorm_mean(ix);
    CCnorm_mean = CCnorm_mean(ix_select);
else
    ix = ix_qualia & ix_npsp; %Find the intersection of the two
end

NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);


%Sort in increasing NPSP
[NPSPs, ix_select] = sort(NPSPs,'ascend');
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
    bf_mean(k) = 2^(log2(bf(k).sl1*bf(k).sl2*bf(k).ls1*bf(k).ls2)/4);
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
params.save_dir = save_dir;
% plot_mean_strf_switch(k_fh_bf,params);

%% Plot a mean STRF across all frequencies using averaging across all kernels
params.type = 'raw';
plot_mean_strf_switch(k_fhn,params);

%% Plot STRF and temporal profiles for each bf separately 
% params.save_dir = save_dir_paper_fig1_sup;
% plot_ind_bf(k_fh_bf,params);

%% Stats
%COM measures
com_sl1_neg = reach(com_neg,'sl1');
com_sl1_pos = reach(com_pos,'sl1');
com_sl2_neg = reach(com_neg,'sl2');
com_sl2_pos = reach(com_pos,'sl2');

com_ls1_neg = reach(com_neg,'ls1');
com_ls1_pos = reach(com_pos,'ls1');
com_ls2_neg = reach(com_neg,'ls2');
com_ls2_pos = reach(com_pos,'ls2');

%IE ratio
ie_ratio_sl1 = reach(ie_ratio,'sl1');
ie_ratio_sl2 = reach(ie_ratio,'sl2');

ie_ratio_ls1 = reach(ie_ratio,'ls1');
ie_ratio_ls2 = reach(ie_ratio,'ls2');

%Total Inhibiton and excitation
total_inh_sl1 = reach(total_inh,'sl1');
total_inh_sl2 = reach(total_inh,'sl2');
total_exc_sl1 = reach(total_exc,'sl1');
total_exc_sl2 = reach(total_exc,'sl2');

total_inh_ls1 = reach(total_inh,'ls1');
total_inh_ls2 = reach(total_inh,'ls2');
total_exc_ls1 = reach(total_exc,'ls1');
total_exc_ls2 = reach(total_exc,'ls2');

%Inhibition Max value
max_inh_sl1 = reach(max_inh,'sl1');
max_inh_sl2 = reach(max_inh,'sl2');
max_exc_sl1 = reach(max_exc,'sl1');
max_exc_sl2 = reach(max_exc,'sl2');

max_inh_ls1 = reach(max_inh,'ls1');
max_inh_ls2 = reach(max_inh,'ls2');
max_exc_ls1 = reach(max_exc,'ls1');
max_exc_ls2 = reach(max_exc,'ls2');

%Time of Max Inhibition
max_inh_time_sl1 = reach(max_inh_time,'sl1');
max_inh_time_sl2 = reach(max_inh_time,'sl2');
max_exc_time_sl1 = reach(max_exc_time,'sl1');
max_exc_time_sl2 = reach(max_exc_time,'sl2');

max_inh_time_ls1 = reach(max_inh_time,'ls1');
max_inh_time_ls2 = reach(max_inh_time,'ls2');
max_exc_time_ls1 = reach(max_exc_time,'ls1');
max_exc_time_ls2 = reach(max_exc_time,'ls2');

%COM p-values
[p_val.com.neg.sl2_sl1,~,~] = signrank(com_sl2_neg, com_sl1_neg);
[p_val.com.pos.sl2_sl1,~,~] = signrank(com_sl2_pos, com_sl1_pos);

[p_val.com.neg.ls2_ls1,~,~] = signrank(com_ls2_neg, com_ls1_neg); 
[p_val.com.pos.ls2_ls1,~,~] = signrank(com_ls2_pos, com_ls1_pos); 

[p_val.com.neg.sl2_ls2,~,~] = signrank(com_sl2_neg, com_ls2_neg); 
[p_val.com.pos.sl2_ls2,~,~] = signrank(com_sl2_pos, com_ls2_pos); 

%PT p-values
[p_val.pt.neg.sl2_sl1,~,~] = signrank(max_inh_time_sl2, max_inh_time_sl1);
[p_val.pt.pos.sl2_sl1,~,~] = signrank(max_exc_time_sl2, max_exc_time_sl1);

[p_val.pt.neg.ls2_ls1,~,~] = signrank(max_inh_time_ls2, max_inh_time_ls1);
[p_val.pt.pos.ls2_ls1,~,~] = signrank(max_exc_time_ls2, max_exc_time_ls1);

[p_val.pt.neg.sl2_ls2,~,~] = signrank(max_inh_time_sl2, max_inh_time_ls2);
[p_val.pt.pos.sl2_ls2,~,~] = signrank(max_exc_time_sl2, max_exc_time_ls2);

%PH p-values
[p_val.ph.neg.sl2_sl1,~,~] = signrank(com_sl2_neg,com_sl1_neg);
[p_val.ph.pos.sl2_sl1,~,~] = signrank(com_sl2_pos,com_sl1_pos);

[p_val.ph.neg.ls2_ls1,~,~] = signrank(com_ls2_neg,com_ls1_neg);
[p_val.ph.pos.ls2_ls1,~,~] = signrank(com_ls2_pos,com_ls1_pos);

[p_val.ph.neg.sl2_ls2,~,~] = signrank(com_sl2_neg,com_ls2_neg);
[p_val.ph.pos.sl2_ls2,~,~] = signrank(com_sl2_pos,com_ls2_pos);


com_sl1_bf_neg = reach(com_neg_bf,'sl1');
com_sl1_bf_pos = reach(com_pos_bf,'sl1');
com_sl2_bf_neg = reach(com_neg_bf,'sl2');
com_sl2_bf_pos = reach(com_pos_bf,'sl2');

com_ls1_bf_neg = reach(com_neg_bf,'ls1');
com_ls1_bf_pos = reach(com_pos_bf,'ls1');
com_ls2_bf_neg = reach(com_neg_bf,'ls2');
com_ls2_bf_pos = reach(com_pos_bf,'ls2');

%% Plot box plot of inhibition no ind points
lw = 3;
axis_sz = 40;
n_points = length(com_sl1_neg);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibition');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),com_sl1_neg,'k');
scatter(2*ones(n_points,1),com_sl2_neg,'k');
scatter(3*ones(n_points,1),com_ls1_neg,'k');
scatter(4*ones(n_points,1),com_ls2_neg,'k');

boxplot([com_sl1_neg',com_sl2_neg',com_ls1_neg', com_ls2_neg'], 'Notch', 'on', 'Labels', {'SL1', 'SL2', 'LS1', 'LS2'});
ylim([50 150]);
hold off;

full_name = fullfile(save_dir_bf,['Switch_inhibition_noind_box','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;

%% Plot box plot of inhibition with ind points
lw = 3;
axis_sz = 40;
n_points = length(com_sl1_neg);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
title('Switching stimuli inhibition');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold on;
scatter(ones(n_points,1),com_sl1_neg,'k');
scatter(2*ones(n_points,1),com_sl2_neg,'k');
scatter(3*ones(n_points,1),com_ls1_neg,'k');
scatter(4*ones(n_points,1),com_ls2_neg,'k');

for n = 1:n_points
    plot([1,2,3,4], [com_sl1_neg(n), com_sl2_neg(n), com_ls1_neg(n), com_ls2_neg(n)],'k');
end

boxplot([com_sl1_neg',com_sl2_neg',com_ls1_neg', com_ls2_neg'], 'Notch', 'on', 'Labels', {'SL1', 'SL2', 'LS1', 'LS2'});
ylim([50 150]);
hold off;

full_name = fullfile(save_dir_bf,['Switch_inhibition_ind_box','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;

%% Plot Violin plot of inhibition
lw = 3;
axis_sz = 40;
n_points = length(com_sl1_neg);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
violinplot([com_sl1_neg',com_sl2_neg',com_ls1_neg', com_ls2_neg'], {'SL1', 'SL2', 'LS1', 'LS2'});
title('Switching stimuli inhibition');
ylabel('COM [ms]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
ylim([50 170]);

full_name = fullfile(save_dir_bf,['Switch_inhibition_violin_plot','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;
%% Plot box plot of excitation ind
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
boxplot([com_sl1_pos',com_sl2_pos',com_ls1_pos', com_ls2_pos'], 'Notch', 'on', 'Labels', {'SL1', 'SL2', 'LS1', 'LS2'});
title('Switching stimuli excitation');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');

hold on;
scatter(ones(n_points,1),com_sl1_pos,'k');
scatter(2*ones(n_points,1),com_sl2_pos,'k');
scatter(3*ones(n_points,1),com_ls1_pos,'k');
scatter(4*ones(n_points,1),com_ls2_pos,'k');

for n = 1:n_points
    plot([1,2,3,4], [com_sl1_pos(n), com_sl2_pos(n), com_ls1_pos(n), com_ls2_pos(n)],'k');
end

hold off;

full_name = fullfile(save_dir_bf,['Switch_excitation_ind_box','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;


%% Plot box plot of excitation no ind
axis_sz = 40;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
boxplot([com_sl1_pos',com_sl2_pos',com_ls1_pos', com_ls2_pos'], 'Notch', 'on', 'Labels', {'SL1', 'SL2', 'LS1', 'LS2'});
title('Switching stimuli excitation');
ylabel('COM [ms]');
xlabel('Condition');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');

hold on;
scatter(ones(n_points,1),com_sl1_pos,'k');
scatter(2*ones(n_points,1),com_sl2_pos,'k');
scatter(3*ones(n_points,1),com_ls1_pos,'k');
scatter(4*ones(n_points,1),com_ls2_pos,'k');
hold off;

full_name = fullfile(save_dir_bf,['Switch_excitation_noind_box','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;

%% Plot Violin plot of excitation
lw = 3;
axis_sz = 40;
n_points = length(com_sl1_neg);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
violinplot([com_sl1_pos',com_sl2_pos',com_ls1_pos', com_ls2_pos'], {'SL1', 'SL2', 'LS1', 'LS2'});
title('Switching stimuli inhibition');
ylabel('COM [ms]');
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
ylim([10 90]);

full_name = fullfile(save_dir_bf,['Switch_excitation_violin_plot','_NPSP_',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),'.svg']);
saveas(gcf, full_name);
close all;
%% Define params for histogram for Excitation and Inhibition
exc_color = 'r'; inh_color = 'b'; lw1 = 7; lw2 = 3;
params_his_ei.fit = 'neurons'; params_his_ei.units = 'normal'; params_his_ei.calc_method = calc_method;
%Plot general properties
params_his_ei.hist_type = hist_type; params_his_ei.get_neurons = get_neurons; params_his_ei.NPSP_th = NPSP_th; params_his_ei.CCnorm_th = CCnorm_th; params_his_ei.model = model; params_his_ei.font_type = font_type;
%Plot colours
params_his_ei.exc_color = exc_color; params_his_ei.inh_color = inh_color;
%Plot size
params_his_ei.all_font_sz = all_font_sz; params_his_ei.lw1 = lw1; params_his_ei.lw2 = lw2;
%Save dir
params_his_ei.save_dir =  save_dir_bf;

%% Define params for scatter plot for Excitation and Inhibition
%Plot general properties
lim_val_ms = 150;
params_ei_scatter.fit = 'neurons'; params_ei_scatter.units = 'normal'; params_ei_scatter.calc_method = calc_method;
params_ei_scatter.get_neurons = get_neurons; params_ei_scatter.NPSP_th = NPSP_th; params_ei_scatter.CCnorm_th = CCnorm_th; params_ei_scatter.model = model; params_ei_scatter.font_type = font_type;
%Plot colours
params_ei_scatter.exc_color = exc_color; params_ei_scatter.inh_color = inh_color;
%Plot size
params_ei_scatter.lw = 3; params_ei_scatter.all_font_sz = all_font_sz; params_ei_scatter.sz = sz;
%Save dir
params_ei_scatter.save_dir = save_dir_bf;

%% COM Excitation and Inhibiton small vs big Histogram SL2 - SL1
params_his_ei.specific_name = ' COM ms change SL2_SL1 for ';
params_his_ei.val_exc_big = com_sl2_pos; params_his_ei.val_exc_small = com_sl1_pos; params_his_ei.val_inh_big = com_sl2_neg; params_his_ei.val_inh_small = com_sl1_neg;
params_his_ei.p_val_exc = p_val.com.pos.sl2_sl1; params_his_ei.p_val_inh = p_val.com.neg.sl2_sl1;
params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
plot_ei_histogram(params_his_ei);

%% COM Excitation and Inhibiton small vs big Histogram LS2 - LS1
params_his_ei.specific_name = ' COM ms change LS2_LS1 for ';
params_his_ei.val_exc_big = com_ls2_pos; params_his_ei.val_exc_small = com_ls1_pos; params_his_ei.val_inh_big = com_ls2_neg; params_his_ei.val_inh_small = com_ls1_neg;
params_his_ei.p_val_exc = p_val.com.pos.ls2_ls1; params_his_ei.p_val_inh = p_val.com.neg.ls2_ls1;
params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
plot_ei_histogram(params_his_ei);

%% COM Excitation and Inhibiton small vs big Histogram SL2 - LS2
params_his_ei.specific_name = ' COM ms change SL2_LS2 for ';
params_his_ei.val_exc_big = com_sl2_pos; params_his_ei.val_exc_small = com_ls2_pos; params_his_ei.val_inh_big = com_sl2_neg; params_his_ei.val_inh_small = com_ls2_neg;
params_his_ei.p_val_exc = p_val.com.pos.sl2_ls2; params_his_ei.p_val_inh = p_val.com.neg.sl2_ls2;
params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
plot_ei_histogram(params_his_ei);

%% COM Excitation and Inhibiton small vs big Scatter SL2 - SL1
params_ei_scatter.specific_name = 'Inhibitory and excitatory SL2_SL1 com no npsp for ';
params_ei_scatter.val_exc_big = com_sl2_pos; params_ei_scatter.val_exc_small = com_sl1_pos; params_ei_scatter.val_inh_big = com_sl2_neg; params_ei_scatter.val_inh_small = com_sl1_neg;
params_ei_scatter.p_val_exc = p_val.com.pos.sl2_sl1; params_ei_scatter.p_val_inh = p_val.com.neg.sl2_sl1;
params_ei_scatter.lim_val = lim_val_ms;
plot_ei_scatter(params_ei_scatter);

%% COM Excitation and Inhibiton small vs big Scatter LS2 - LS1
params_ei_scatter.specific_name = 'Inhibitory and excitatory LS2_LS1 com no npsp for ';
params_ei_scatter.val_exc_big = com_ls2_pos; params_ei_scatter.val_exc_small = com_ls1_pos; params_ei_scatter.val_inh_big = com_ls2_neg; params_ei_scatter.val_inh_small = com_ls1_neg;
params_ei_scatter.p_val_exc = p_val.com.pos.ls2_ls1; params_ei_scatter.p_val_inh = p_val.com.neg.ls2_ls1;
params_ei_scatter.lim_val = lim_val_ms;
plot_ei_scatter(params_ei_scatter);

%% COM Excitation and Inhibiton small vs big Scatter SL2 - LS2
params_ei_scatter.specific_name = 'Inhibitory and excitatory SL2_LS2 com no npsp for ';
params_ei_scatter.val_exc_big = com_sl2_pos; params_ei_scatter.val_exc_small = com_ls2_pos; params_ei_scatter.val_inh_big = com_sl2_neg; params_ei_scatter.val_inh_small = com_ls2_neg;
params_ei_scatter.p_val_exc = p_val.com.pos.sl2_ls2; params_ei_scatter.p_val_inh = p_val.com.neg.sl2_ls2;
params_ei_scatter.lim_val = lim_val_ms;
plot_ei_scatter(params_ei_scatter);

%% Compute some important difference measures

%Median difference
%L2 - L1
%COM
median_diff.com.neg.l2_l1 = nanmedian(com_sl2_neg-com_sl1_neg);
median_diff.com.pos.l2_l1  = nanmedian(com_sl2_pos-com_sl1_pos);
%PT
median_diff.pt.neg.l2_l1  = nanmedian(max_inh_time_sl2-max_inh_time_sl1);
median_diff.pt.pos.l2_l1  = nanmedian(max_exc_time_sl2-max_exc_time_sl1);
%PH
median_diff.ph.neg.l2_l1  = nanmedian(log2(max_inh_sl2./max_inh_sl1));
median_diff.ph.pos.l2_l1  = nanmedian(log2(max_exc_sl2./max_exc_sl1));

%S2 - S1
%COM
median_diff.com.neg.s2_s1  = nanmedian(com_ls2_neg-com_ls1_neg);
median_diff.com.pos.s2_s1  = nanmedian(com_ls2_pos-com_ls1_pos);
%PT
median_diff.pt.neg.s2_s1  = nanmedian(max_inh_time_ls2-max_inh_time_ls1);
median_diff.pt.pos.s2_s1  = nanmedian(max_exc_time_ls2-max_exc_time_ls1);
%PH
median_diff.ph.neg.s2_s1  = nanmedian(log2(max_inh_ls2./max_inh_ls1));
median_diff.ph.pos.s2_s1  = nanmedian(log2(max_exc_ls2./max_exc_ls1));

%L2 - S2
%COM
median_diff.com.neg.l2_s2 = nanmedian(com_sl2_neg-com_ls2_neg);
median_diff.com.pos.l2_s2  = nanmedian(com_sl2_pos-com_ls2_pos);
%PT
median_diff.pt.neg.l2_s2  = nanmedian(max_inh_time_sl2-max_inh_time_ls2);
median_diff.pt.pos.l2_s2  = nanmedian(max_exc_time_sl2-max_exc_time_ls2);
%PH
median_diff.ph.neg.l2_s2  = nanmedian(log2(max_inh_sl2./max_inh_ls2));
median_diff.ph.pos.l2_s2  = nanmedian(log2(max_exc_sl2./max_exc_ls2));

%Mean difference
%L2 - L1
%COM
mean_diff.com.neg.l2_l1 = nanmean(com_sl2_neg-com_sl1_neg);
mean_diff.com.pos.l2_l1  = nanmean(com_sl2_pos-com_sl1_pos);
%PT
mean_diff.pt.neg.l2_l1  = nanmean(max_inh_time_sl2-max_inh_time_sl1);
mean_diff.pt.pos.l2_l1  = nanmean(max_exc_time_sl2-max_exc_time_sl1);
%PH
mean_diff.ph.neg.l2_l1  = nanmean(log2(max_inh_sl2./max_inh_sl1));
mean_diff.ph.pos.l2_l1  = nanmean(log2(max_exc_sl2./max_exc_sl1));

%S2 - S1
%COM
mean_diff.com.neg.s2_s1  = nanmean(com_ls2_neg-com_ls1_neg);
mean_diff.com.pos.s2_s1  = nanmean(com_ls2_pos-com_ls1_pos);
%PT
mean_diff.pt.neg.s2_s1  = nanmean(max_inh_time_ls2-max_inh_time_ls1);
mean_diff.pt.pos.s2_s1  = nanmean(max_exc_time_ls2-max_exc_time_ls1);
%PH
mean_diff.ph.neg.s2_s1  = nanmean(log2(max_inh_ls2./max_inh_ls1));
mean_diff.ph.pos.s2_s1  = nanmean(log2(max_exc_ls2./max_exc_ls1));

%L2 - S2
%COM
mean_diff.com.neg.l2_s2 = nanmean(com_sl2_neg-com_ls2_neg);
mean_diff.com.pos.l2_s2  = nanmean(com_sl2_pos-com_ls2_pos);
%PT
mean_diff.pt.neg.l2_s2  = nanmean(max_inh_time_sl2-max_inh_time_ls2);
mean_diff.pt.pos.l2_s2  = nanmean(max_exc_time_sl2-max_exc_time_ls2);
%PH
mean_diff.ph.neg.l2_s2  = nanmean(log2(max_inh_sl2./max_inh_ls2));
mean_diff.ph.pos.l2_s2  = nanmean(log2(max_exc_sl2./max_exc_ls2));