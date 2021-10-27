function plot_model_tconst_grande_cv_best(kernel_dir)
%% Params
save_dir = fullfile(kernel_dir,'Plots');
hist_type = 'bar';
calc_method = 'average'; %How to estimate impotant values
%calc_method -- 'raw', 'average'
%               'raw' - Use the whole receptive field 
%               'average' - Take an average across frequencies first

bf_method = 'window_mean'; %The method use to extract the bf:
% bf_method -- 'max', 'window_max', 'window_mean'
%              'max' - Take the max across all history
%              'window_max' - Take the max in a given time window
%              'window_mean' - Take the mean in a given time window
bf_window_ms = 190; %The window for calcualting bf if this method is used

plot_hyperparams = 0;
plot_histograms = 0;
chop_ms = 190; 
r_type{1} = 'small';
r_type{2} = 'med';
r_type{3} = 'big';
freq_up_bound = 17000;
freq_down_bound = 700;
font_type = 'Liberation Sans';
legend_font_sz = 32;
all_font_sz = 70;
sz = 40;
n_rooms = length(r_type);
exc_small_color = [0.9804    0.4196    0.6431];
exc_big_color = [0.8 0 0];
inh_small_color = [0.0745 0.6235 1.0000];
inh_big_color = 'b';
save_dir_paper_fig4 = '/mnt/40086D4C086D41D0/Reverb_paper/fig_4/Normative_model';
save_dir_paper_fig3 = '/mnt/40086D4C086D41D0/Reverb_paper/fig_3/Normative_model';
save_dir_paper_fig2 = '/mnt/40086D4C086D41D0/Reverb_paper/fig_2/Normative_model';
save_dir_paper_fig1_sup = '/mnt/40086D4C086D41D0/Reverb_paper/fig_2_sup_1_2/Normative_model';
%% Load the data
load(fullfile(kernel_dir,'kernel'));
best_fold = length(kernels.small{1}.main);
model = kernels.model;
%% First compute the bf and com for every cluster
freqs = fliplr(kernels.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_ker = length(freqs);
n_h = kernels.n_h-1;
dt_ms = round(kernels.dt_ms);
chop_ix = round(chop_ms/dt_ms);
bf_window_ix = round(bf_window_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;

fprintf('== Calcuating com and bf ==\n');tic;
for k = 1:n_ker
    for r = 1:n_rooms
        room = r_type{r};
        switch model
            case {'sep','sep_kh'}
                [~,ix] = max(kernels{k}.(room).k_f);
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
                k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                
            case {'ridge','lasso','manual_alpha','elastic'}
                k_fh = fliplr(kernels.(room){k}.main{best_fold}.k_fh);
                k_fh = k_fh(:,1:chop_ix);
                k_fh_bf.(room){k} = k_fh;
                k_fh_neg = abs(min(k_fh,0));
                k_fh_pos = abs(max(k_fh,0));        
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
        k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
        k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
        com_neg(k).(room) = (k_h_neg*h); %Compute a weighted sum of all values
        com_pos(k).(room) = (k_h_pos*h); %Compute a weighted sum of all values
        if plot_hyperparams
            alpha(k).(room) = kernels.(room){best_fold}.main.alpha;
            lambda(k).(room) = kernels.(room){best_fold}.main.lambda;
        end
    end
    bf_mean(k) = mean([bf(k).small,bf(k).big]);
    [~,ix_bf] = min(abs(bf_mean(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
    bf_closest(k) = freqs(ix_bf);
end

%% Make folders
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
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

%% Calculate some important differences
%COM
com_median_neg_diff = nanmedian(com_big_neg-com_small_neg);
com_median_pos_diff = nanmedian(com_big_pos-com_small_pos);
%PT
tmax_median_neg_diff = nanmedian(max_inh_time_big-max_inh_time_small);
tmax_median_pos_diff = nanmedian(max_exc_time_big-max_exc_time_small);
%PH
ph_median_neg_diff = nanmedian(log2(max_inh_big./max_inh_small));
ph_median_pos_diff = nanmedian(log2(max_exc_big./max_exc_small));

%% Plot mean STRF and k_h for the two rooms
params.freqs = freqs;
params.n_h = n_h;
params.dt_ms = dt_ms;
params.h_max_ms = chop_ms;
params.font_type = font_type;
params.exc_small_color = exc_small_color; params.exc_big_color = exc_big_color;
params.inh_small_color = inh_small_color; params.inh_big_color = inh_big_color;
params.save_dir = save_dir_paper_fig2;
params.name = 'normative model';
params.type = 'bf';
plot_mean_strf2(k_fh_bf,params);

%% Plot STRF and temporal profiles for each bf separately 
params.save_dir = save_dir_paper_fig1_sup;
plot_ind_bf(k_fh_bf,params);

if plot_histograms
    %% Define params for histogram for Excitation and Inhibition
    exc_color = 'r'; inh_color = 'b'; lw1 = 5; lw2 = 3;
    params_his_ei.units = 'normal'; params_his_ei.fit = 'normative'; params_his_ei.calc_method = calc_method;
    %Plot general properties
    params_his_ei.hist_type = hist_type; params_his_ei.model = model; params_his_ei.font_type = font_type;
    %Plot colours
    params_his_ei.exc_color = exc_color; params_his_ei.inh_color = inh_color; 
    %Plot size
    params_his_ei.all_font_sz = all_font_sz; params_his_ei.lw1 = lw1; params_his_ei.lw2 = lw2;
    %Save dir
    params_his_ei.save_dir =  save_dir_paper_fig3;
    %% Define params for Single histogram
    %Plot general properties
    params_norm_his.fit = 'normative'; params_norm_his.calc_method = calc_method;
    params_norm_his.hist_type = hist_type; params_norm_his.model = model; params_norm_his.font_type = font_type;
    %Plot colours
    params_norm_his.his_color = 'k';
    %Plot size
    params_norm_his.all_font_sz = all_font_sz; params_norm_his.lw1 = lw1; params_norm_his.lw2 = lw2;
    %Save dir
    params_norm_his.save_dir =  save_dir_paper_fig3;

    %% Define params for scatter plot for Excitation and Inhibition 
    %Plot general properties
    params_ei_scatter.units = 'normal'; params_ei_scatter.fit = 'normative'; params_ei_scatter.calc_method = calc_method;
    params_ei_scatter.model = model; params_ei_scatter.font_type = font_type;
    %Plot colours
    params_ei_scatter.exc_color = exc_color; params_ei_scatter.inh_color = inh_color;
    %Plot size
    params_ei_scatter.lw = 3; params_ei_scatter.all_font_sz = all_font_sz; params_ei_scatter.sz = sz;
    %Save dir
    params_ei_scatter.save_dir = save_dir_paper_fig3;
    
    %% Define params for scatter plot with NPSP 
    %Plot general properties
    params_npsp_scatter.fit = 'normative'; params_npsp_scatter.npsp_on = 0;
    params_npsp_scatter.model = model; params_npsp_scatter.font_type = font_type; 
    %Plot size
    params_npsp_scatter.lw = 3; params_npsp_scatter.all_font_sz = all_font_sz; params_npsp_scatter.sz = sz;
    %Save dir
    params_npsp_scatter.save_dir = save_dir_paper_fig3;
    
    %% COM Excitation and Inhibiton small vs big Scatter
    lim_val_ms = 150;
    params_ei_scatter.specific_name = ' Normative model Inhibitory and excitatory com no npsp for ';
    params_ei_scatter.val_exc_big = com_big_pos; params_ei_scatter.val_exc_small = com_small_pos; params_ei_scatter.val_inh_big = com_big_neg; params_ei_scatter.val_inh_small = com_small_neg;
    params_ei_scatter.p_val_exc = p_val_pos; params_ei_scatter.p_val_inh = p_val_neg;
    params_ei_scatter.lim_val = lim_val_ms;
    plot_ei_scatter(params_ei_scatter);
    
    %% IE Ratio small vs big Scatter
    params_npsp_scatter.specific_name = ' Normative model IE Ratio no NPSP for test ';
    params_npsp_scatter.val_big = ie_ratio_big; params_npsp_scatter.val_small = ie_ratio_small;
    params_npsp_scatter.p_val = p_val_ie;
    params_npsp_scatter.lim_val = 1;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% Total Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = ' Normative model E and I total no NPSP for test ';
    params_ei_scatter.val_exc_big = total_exc_big; params_ei_scatter.val_exc_small = total_exc_small; params_ei_scatter.val_inh_big = total_inh_big; params_ei_scatter.val_inh_small = total_inh_small;
    params_ei_scatter.p_val_exc = p_val_total_exc; params_ei_scatter.p_val_inh = p_val_total_inh;
    params_ei_scatter.lim_val = 0.3;
    plot_ei_scatter(params_ei_scatter);
    
    %% PH Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = ' Normative model PH E and I max no NPSP for ';
    params_ei_scatter.val_exc_big = max_exc_big; params_ei_scatter.val_exc_small = max_exc_small; params_ei_scatter.val_inh_big = max_inh_big; params_ei_scatter.val_inh_small = max_inh_small;
    params_ei_scatter.p_val_exc = p_val_max_exc; params_ei_scatter.p_val_inh = p_val_max_inh;
    params_ei_scatter.lim_val = 0.03;
    plot_ei_scatter(params_ei_scatter);
    
    %% PT Excitation and Inhibiton small vs big Scatter
    params_ei_scatter.specific_name = ' Normative model PT E and I no NPSP for test ';
    g_noise = 5*randn(1,length(max_inh_time_small));
    params_ei_scatter.val_exc_big = max_exc_time_big+g_noise; params_ei_scatter.val_exc_small = max_exc_time_small+g_noise; params_ei_scatter.val_inh_big = max_inh_time_big+g_noise; params_ei_scatter.val_inh_small = max_inh_time_small+g_noise;
    params_ei_scatter.p_val_exc = p_val_max_exc_time; params_ei_scatter.p_val_inh = p_val_max_inh_time;
    params_ei_scatter.lim_val = 100;
    plot_ei_scatter(params_ei_scatter);
    
    %% COM Excitation and Inhibiton small vs big Histogram
    params_his_ei.specific_name = ' Normative model of COM ms change for ';
    params_his_ei.val_exc_big = com_big_pos; params_his_ei.val_exc_small = com_small_pos; params_his_ei.val_inh_big = com_big_neg; params_his_ei.val_inh_small = com_small_neg;  
    params_his_ei.p_val_exc = p_val_pos; params_his_ei.p_val_inh = p_val_neg;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);

    %% COM change ix Excitation and Inhibiton small vs big Histogram
    params_his_ei.units = 'change ix';
    params_his_ei.specific_name = ' Normative model of COM change ix for test ';
    params_his_ei.val_exc_big = com_big_pos; params_his_ei.val_exc_small = com_small_pos; params_his_ei.val_inh_big = com_big_neg; params_his_ei.val_inh_small = com_small_neg;
    params_his_ei.p_val_exc = p_val_pos; params_his_ei.p_val_inh = p_val_neg;
    params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
    plot_ei_histogram(params_his_ei);
    
    %% Total Excitation and Inhibiton small vs big room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' Normative model of Total amount change for ';
    params_his_ei.val_exc_big = total_exc_big; params_his_ei.val_exc_small = total_exc_small; params_his_ei.val_inh_big = total_inh_big; params_his_ei.val_inh_small = total_inh_small;
    params_his_ei.p_val_exc = p_val_total_exc; params_his_ei.p_val_inh = p_val_total_inh;
    params_his_ei.his_spacing = 0.1; params_his_ei.lim_val = 4;
    plot_ei_histogram(params_his_ei);
    
    %% PH Excitation and Inhibiton for small vs big room Histogram
    params_his_ei.units = 'ratio';
    params_his_ei.specific_name = ' Normative model PH E and I change for ';
    params_his_ei.val_exc_big = max_exc_big; params_his_ei.val_exc_small = max_exc_small; params_his_ei.val_inh_big = max_inh_big; params_his_ei.val_inh_small = max_inh_small;
    params_his_ei.p_val_exc = p_val_max_exc; params_his_ei.p_val_inh = p_val_max_inh;
    params_his_ei.his_spacing = 0.125; params_his_ei.lim_val = 2.5;
    plot_ei_histogram(params_his_ei);
    
    %% PT Excitation and Inhibiton for small vs big room Histogram
    params_his_ei.units = 'normal';
    params_his_ei.specific_name = ' Normative model PT E and I change for ';
    params_his_ei.val_exc_big = max_exc_time_big; params_his_ei.val_exc_small = max_exc_time_small; params_his_ei.val_inh_big = max_inh_time_big; params_his_ei.val_inh_small = max_inh_time_small;
    params_his_ei.p_val_exc = p_val_max_exc_time; params_his_ei.p_val_inh = p_val_max_inh_time;
    params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
    plot_ei_histogram(params_his_ei);
      
    %% IE Ratio change ix small vs big room Histogram
    params_norm_his.specific_name = ' Normative model IE Ratio of com change ix for ';
    params_norm_his.units = 'change ix'; 
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie; 
    params_norm_his.his_spacing = 2.5; params_norm_his.lim_val = 20;
    plot_norm_histogram(params_norm_his);
    
    %% IE ratio for small vs big room Histogram
    params_norm_his.specific_name = ' Normative model of IE Ratio AU change for test ';
    params_norm_his.units = 'normal';
    params_norm_his.big_val = ie_ratio_big; params_norm_his.small_val = ie_ratio_small;
    params_norm_his.p_val = p_val_ie;
    params_norm_his.his_spacing = 0.05; params_norm_his.lim_val = 1;
    plot_norm_histogram(params_norm_his);
end

%% Plot inhibitory COM vs frequency for all rooms

skip_f = 3;
lw = 6;
x_box_dis = 0.75;
ix = freqs>freq_down_bound & freqs<freq_up_bound;
plot_freqs = fliplr(freqs(ix));
n_f = length(plot_freqs);

for f = 1:n_f
    f_labels{f} = num2str(plot_freqs(f)./1000,'%.1f');
end

plot_com_small_neg = fliplr(com_small_neg(ix));
plot_com_large_neg = fliplr(com_big_neg(ix));
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

% semilogx(plot_freqs,plot_com_small_neg,'Color',[0.31 0.75 0.405],'LineWidth',lw);hold on;
semilogx(plot_freqs,plot_com_small_neg,'Color',inh_small_color,'LineWidth',lw);hold on;
% semilogx(plot_freq,com_med_neg,'Color',[0.872 0.496 0.192],'LineWidth',lw);
% semilogx(plot_freqs,plot_com_large_neg,'Color',[0.642 0.456 0.924],'LineWidth',lw); hold off;
semilogx(plot_freqs,plot_com_large_neg,'Color','b','LineWidth',lw); hold off;
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
% xlabel('Frequency [kHz]','FontSize',all_font_sz,'FontWeight','Normal');
ylabel('COM [ms]','FontSize',all_font_sz,'FontWeight','Normal');
title(['Inhibitory COM vs frequency for normative ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room inh'),'LineStyle','none','Color',inh_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room inh'),'LineStyle','none','Color',inh_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
ylim([50 100]);
save_name = fullfile(save_dir_paper_fig4,['Inhibitory COM vs freq for normative ',model,' kernel.svg']);
saveas(gcf,save_name);
close;

%% Plot excitatory com vs frequency for all rooms
plot_com_small_pos = fliplr(com_small_pos(ix));
plot_com_large_pos = fliplr(com_big_pos(ix));

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
% semilogx(plot_freqs,plot_com_small_pos,'Color',[0.31 0.75 0.405],'LineWidth',lw);hold on;
semilogx(plot_freqs,plot_com_small_pos,'Color',exc_small_color,'LineWidth',lw);hold on;
% semilogx(plot_freqs,com_med_pos,'Color',[0.872 0.496 0.192],'LineWidth',lw);
% semilogx(plot_freqs,plot_com_large_pos,'Color',[0.642 0.456 0.924],'LineWidth',lw);hold off;
semilogx(plot_freqs,plot_com_large_pos,'Color',exc_big_color,'LineWidth',lw);hold off;
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
% xlabel('Frequency [kHz]','FontSize',all_font_sz,'FontWeight','Normal');
ylabel('COM [ms]','FontSize',all_font_sz,'FontWeight','Normal');
title(['Excitatory COM vs frequency for normative ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room exc'),'LineStyle','none','Color',exc_big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room exc'),'LineStyle','none','Color',exc_small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',axis_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',axis_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',axis_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
ylim([25 75]);
save_name = fullfile(save_dir_paper_fig4,['Excitatory COM vs freq for normative ',model,' kernel.svg']);
saveas(gcf,save_name);
close;

%% Plot excitatory and inhibitory com vs frequency together for all rooms
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
semilogx(plot_freqs,plot_com_small_neg,'Color',inh_small_color,'LineWidth',lw);hold on;
semilogx(plot_freqs,plot_com_large_neg,'Color',inh_big_color,'LineWidth',lw); 
semilogx(plot_freqs,plot_com_small_pos,'Color',exc_small_color,'LineWidth',lw);
semilogx(plot_freqs,plot_com_large_pos,'Color',exc_big_color,'LineWidth',lw);hold off;
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
% xlabel('Frequency [kHz]','FontSize',all_font_sz,'FontWeight','Normal');
% ylabel('\com [ms]','FontSize',all_font_sz,'FontWeight','Normal');
% title(['Excitatory \com vs frequency for normative ',model,' kernel']);
annotation('textbox',[x_box_dis 0.8 0.1 0.1],'String', sprintf('Large room inh'),'LineStyle','none','Color',inh_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.75 0.1 0.1],'String', sprintf('Small room inh'),'LineStyle','none','Color',inh_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.7 0.1 0.1],'String', sprintf('Large room exc'),'LineStyle','none','Color',exc_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[x_box_dis 0.65 0.1 0.1],'String', sprintf('Small room exc'),'LineStyle','none','Color',exc_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
ylim([40 90]);
save_name = fullfile(save_dir_paper_fig4,['Both com vs freq for normative ',model,' kernel.svg']);
saveas(gcf,save_name);
close;

%% Fit regression for the two Rooms
legend_font_sz = 50;
all_font_sz = 75;
plot_freqs = plot_freqs(:);
n_reg = length(plot_freqs);
X_small = [ones(n_reg,1),log10(plot_freqs)];
X_big = [ones(n_reg,1),log10(plot_freqs)];

%Inhibitory regression
y_small_neg = plot_com_small_neg; y_small_neg = y_small_neg(:);
y_big_neg = plot_com_large_neg; y_big_neg = y_big_neg(:);

[b_small_neg,~,~,~,stats_small_neg] = regress(y_small_neg,X_small);
[b_big_neg,~,~,~,stats_big_neg] = regress(y_big_neg,X_big);
r_small_neg = sign(b_small_neg(2))*sqrt(stats_small_neg(1)); p_small_neg = stats_small_neg(3);
r_big_neg = sign(b_big_neg(2))*sqrt(stats_big_neg(1)); p_big_neg = stats_big_neg(3);

%Excitatory regression
y_small_pos = plot_com_small_pos; y_small_pos = y_small_pos(:);
y_big_pos = plot_com_large_pos; y_big_pos = y_big_pos(:);

[b_small_pos,~,~,~,stats_small_pos] = regress(y_small_pos,X_small);
[b_big_pos,~,~,~,stats_big_pos] = regress(y_big_pos,X_big);
r_small_pos = sign(b_small_pos(2))*sqrt(stats_small_pos(1)); p_small_pos = stats_small_pos(3);
r_big_pos = sign(b_big_pos(2))*sqrt(stats_big_pos(1)); p_big_pos = stats_big_pos(3);

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

plot(log10(plot_freqs), plot_com_small_neg,'-','Color',inh_small_color,'LineWidth',lw);hold on;
plot(log10(plot_freqs), plot_com_large_neg,'-','Color',inh_big_color,'LineWidth',lw);
plot(log10(plot_freqs), plot_com_small_pos,'-','Color',exc_small_color,'LineWidth',lw);
plot(log10(plot_freqs), plot_com_large_pos,'-','Color',exc_big_color,'LineWidth',lw);
set(gca, 'XScale', 'log');
h1 = refline([b_small_neg(2),b_small_neg(1)]); h1.Color = inh_small_color; h1.LineWidth = lw; h1.LineStyle = '--';
h2 = refline([b_big_neg(2),b_big_neg(1)]); h2.Color = inh_big_color; h2.LineWidth = lw; h2.LineStyle = '--';
h3 = refline([b_small_pos(2),b_small_pos(1)]); h3.Color = exc_small_color; h3.LineWidth = lw; h3.LineStyle = '--';
h4 = refline([b_big_pos(2),b_big_pos(1)]); h4.Color = exc_big_color; h4.LineWidth = lw; h4.LineStyle = '--';
hold off;
xticks([log10(plot_freqs(2:skip_f:n_f))]);
xticklabels(f_labels(2:skip_f:n_f));
ylim([40 90]);
yticks([40:10:90]);

set(gca,'XMinorTick','off');
annotation('textbox',[0.38 0.91 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}',r_big_neg,p_star(p_big_neg)),'LineStyle','none','Color',inh_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.38 0.83 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}\n',r_small_neg,p_star(p_small_neg)),'LineStyle','none','Color',inh_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.66 0.91 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}\n',r_big_pos,p_star(p_big_pos)),'LineStyle','none','Color',exc_big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.66 0.83 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}\n',r_small_pos,p_star(p_small_pos)),'LineStyle','none','Color',exc_small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir_paper_fig4,['Regression  vs BF individual rooms.svg']);
saveas(gcf,save_name);
close;
%% Plot bf consistency accross clusters with histograms
bf_small = reach(bf,'small')/1000;
bf_big = reach(bf,'big')/1000;
diff_bf = abs(log2(bf_big./bf_small));
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
histogram(diff_bf,'Normalization','probability');
xlabel('\DeltaBF [Octaves]','FontSize',16,'FontWeight','Normal');
ylabel('Probability','FontSize',16,'FontWeight','Normal');
title(['Histogram of BF consistency small<->big  for normative ',model,' kernel']);
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir,['Histogram of BF consistency for normative ',model,' kernel.png']);
export_fig(save_name);
close all;

if plot_hyperparams
    %% Plot alpha for the 3 different conditions to make sure there is no bias
    lw = 3;
    alpha_small = reach(alpha,'small');
    alpha_med = reach(alpha,'med');
    alpha_big = reach(alpha,'big');
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(freqs,alpha_big,'r','LineWidth', lw);hold on;
    plot(freqs,alpha_med,'Color',[0.85, 0.325, 0.098],'LineWidth', lw)
    plot(freqs,alpha_small,'b','LineWidth', lw);
    xlabel('BF [kHz]','FontSize',16,'FontWeight','Normal');
    ylabel('\alpha [AU]','FontSize',16,'FontWeight','Normal');
    title(['\alpha parameter vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Normal');
    set(gca, 'XScale', 'log');
    legend('Big Room','Medium Room','Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['alpha vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
    %% Plot the difference between big and small, big and medium, medium and small
    alpha_dif_big_small = alpha_big - alpha_small;
    alpha_dif_big_med = alpha_big - alpha_med;
    alpha_dif_med_small = alpha_med - alpha_small;
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(freqs,alpha_dif_big_small,'r','LineWidth', lw);hold on;
    plot(freqs,alpha_dif_big_med,'Color',[0.85, 0.325, 0.098],'LineWidth', lw)
    plot(freqs,alpha_dif_med_small,'b','LineWidth', lw);
    xlabel('BF [kHz]','FontSize',16,'FontWeight','Normal');
    ylabel('\alpha difference [AU]','FontSize',16,'FontWeight','Normal');
    title(['\alpha difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Normal');
    set(gca, 'XScale', 'log');
    legend('Big - Small Room','Big - Medium Room','Medium - Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['alpha difference vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
    
    %% Plot lambda for the 3 different conditions to make sure there is no bias
    lw = 3;
    lambda_small = reach(lambda,'small');
    lambda_med = reach(lambda,'med');
    lambda_big = reach(lambda,'big');
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(freqs,lambda_big,'r','LineWidth', lw);hold on;
    plot(freqs,lambda_med,'Color',[0.85, 0.325, 0.098],'LineWidth', lw)
    plot(freqs,lambda_small,'b','LineWidth', lw);
    xlabel('BF [kHz]','FontSize',16,'FontWeight','Normal');
    ylabel('\lambda [AU]','FontSize',16,'FontWeight','Normal');
    title(['\lambda parameter vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    set(gca, 'XScale', 'log');
    legend('Big Room','Medium Room','Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['lambda vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
    %% Plot the difference between big and small, big and medium, medium and small
    lambda_dif_big_small = lambda_big - lambda_small;
    lambda_dif_big_med = lambda_big - lambda_med;
    lambda_dif_med_small = lambda_med - lambda_small;
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(freqs,lambda_dif_big_small,'r','LineWidth', lw);hold on;
    plot(freqs,lambda_dif_big_med,'Color',[0.85, 0.325, 0.098],'LineWidth', lw)
    plot(freqs,lambda_dif_med_small,'b','LineWidth', lw);
    xlabel('BF [kHz]','FontSize',16,'FontWeight','Normal');
    ylabel('\lambda difference [AU]','FontSize',16,'FontWeight','Normal');
    title(['\lambda difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Normal');
    set(gca, 'XScale', 'log');
    legend('Big - Small Room','Big - Medium Room','Medium - Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['lambda difference vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
    
end