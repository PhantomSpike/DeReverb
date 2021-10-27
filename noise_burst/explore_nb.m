%% Define params
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Noise_burst_analysis';

normalize = 0; %Whether to normalize the nb COM calculation by removing the minimum value for each neuron  

%Which plots to plot
plot_nb = 0;
plot_nb_all_neurons = 0;
plot_histograms = 1;
plot_scatter = 0;
plot_violins = 0;

%Plotting params
xlim_start_ms = -10; 
xlim_end_ms = 110;
ylim_start_ms = 2; 
ylim_end_ms = 20;

%NPSP threshold for natural stimuli and the noise
NPSP_th = 40;
NPSP_nb_th = 40;

%Time window for COM calculation in ms
com_start_ms = 0;
com_end_ms = 100;

%Time window for calculating the baseline
baseline_start_ms = -50;
baseline_end_ms = 0;

get_neurons = 'all';
dt_ms = 10; %The bin size for the psth
psth_start_ms = -100; %The starting point of the psth relative to the noiseburst onset
psth_end_ms = 600; %The ending point of the psth relative to the noiseburst onset
nb_len_ms = 500; %The actual lenght of the noiseburst in ms
room_type{1} = 'small'; room_type{2} = 'large';
n_rooms = length(room_type);
n_stim = 2; %The number of stimuli for each type of room
n_nb = 7; %The number of noise bursts for each room and stimulus
font_type = 'Liberation Sans';
small_col = [0.592, 0.737, 0.384];
large_col = [0.173, 0.373, 0.176];
t_edges_ms = [psth_start_ms:dt_ms:psth_end_ms]; %Edges for the psth in s
t_edges_s = t_edges_ms/1000; %Convert to s
%% Find the indices and timing of the nosie bursts
meta_dir1 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Experimental_stimuli_old/Ronnie_23_01_2018_stimuli/reverb_with_noise_stimuli/wav_files/reverb_with_noise_same/alexreverbmetadata';
meta_dir2 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Experimental_stimuli_old/Derry_Kilkenny_Cork_stims/reverb_with_noise_stimuli/wav_files/alexreverbmetadata.mat';
meta_dir3  = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Experimental_stimuli_new/Noah_Derekah_stims/reverb_with_noise_stimuli/wav_files/alexreverbmetadata.mat';
n_animal_cases = 3;

%Define which stimuli are small and large for the original data
select{1}.small(1) = 2; select{1}.small(2) = 5;
select{1}.large(1) = 3; select{1}.large(2) = 6;

select{2}.small(1) = 2; select{2}.small(2) = 5;
select{2}.large(1) = 3; select{2}.large(2) = 6;

select{3}.small(1) = 1; select{3}.small(2) = 4;
select{3}.large(1) = 3; select{3}.large(2) = 6;

%Define which stimuli are small and large for the way the neural data is
%organized

select_y{1}.small(1) = 3; select_y{1}.small(2) = 4;
select_y{1}.large(1) = 5; select_y{1}.large(2) = 6;

select_y{2}.small(1) = 3; select_y{2}.small(2) = 4;
select_y{2}.large(1) = 5; select_y{2}.large(2) = 6;

select_y{3}.small(1) = 1; select_y{3}.small(2) = 2;
select_y{3}.large(1) = 5; select_y{3}.large(2) = 6;

%Load the metadata
meta{1} = load(meta_dir1,'alexreverbmetadata');
meta{2} = load(meta_dir2,'alexreverbmetadata');
meta{3} = load(meta_dir3,'alexreverbmetadata');

%Find out which presentations have a nosie burst
for a = 1:n_animal_cases
    hasnoise = [];
    hasnoise = reach(meta{a}.alexreverbmetadata,'hasnoise');
    hasnoise = hasnoise([1,10,2:9],:); %Rearrange the rows to match the neural data
    for r = 1:n_rooms
        room = room_type{r};
        for s = 1:n_stim
            nb_ix{a}.(room){s} = find(hasnoise(:,select{a}.(room)(s))==1);
        end
    end
end

%Get the exact timing of the noise burst
for a = 1:n_animal_cases
    noisepos = [];
    noisepos = reach(meta{a}.alexreverbmetadata,'noisepos');
    noisepos = noisepos([1,10,2:9],:); %Rearrange the rows to match the neural data
    for r = 1:n_rooms
        room = room_type{r};
        for s = 1:n_stim
            nb_pos_s{a}.(room){s} = noisepos(nb_ix{a}.(room){s},select{a}.(room)(s));
        end
    end
end

%% Select the data based on NPSP
load(fullfile(cluster_dir,'info'),'info');
switch get_neurons
    case 'all'
        ix_qualia = ones(length(info.cluster_id),1); %Get all the neurons
    case 'good'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'good'),info.quality,'UniformOutput',false)); %Find all the good units
    case 'mua'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'mua'),info.quality,'UniformOutput',false)); %Find all the mua units
end

ix_npsp = info.NPSP<NPSP_th; %Select all the neurons below certain NPSP
ix = ix_qualia & ix_npsp; %Find the intersection of the two

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

%% First load all the selected clusters
fprintf('== Processing clusters ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(cluster_dir, strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'data');
    
    if ismember(data.params.animal_name,["Ronnie","PLP"])
        a = 1;
    elseif ismember(data.params.animal_name,["Cork","Kilkenny","Derry"])
        a = 2;
    elseif ismember(data.params.animal_name,["Noah","Derekah"])
        a = 3;
    else
        error('Unrecognized animal');
    end
    
    %Make the psths for the noisebusrts for each room and stimulus
    for r = 1:n_rooms
        room = room_type{r}; %Seelct the room
        
        for s = 1:n_stim
            ix_room = select_y{a}.(room)(s); %Select the right index for the room and stim
            room_name = [room,'_',num2str(s)];
            for n = 1:n_nb
                ix_rep = nb_ix{a}.(room){s}(n); %Get the rep that had a noiseburst
                nb_start_time_s = nb_pos_s{a}.(room){s}(n); %Get the onset time of the noiseburst in s
                psth_edges_s = nb_start_time_s + t_edges_s; %Make the edges of the histogram
                y = data.stim(ix_room).repeat(ix_rep).spiketimes; %Get the spiketimes corresponding to the right cluster, room and repeat
                psth_all(k,1).(room_name)(n,:) = histc(y,psth_edges_s);%Generate the psth
            end
        end
    end
    
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Combine the results for all neurons

%Get the 4 different conditions across all neurons
small_1 = reach(psth_all,'small_1'); small_2 = reach(psth_all,'small_2');
large_1 = reach(psth_all,'large_1'); large_2 = reach(psth_all,'large_2');

%Combine them together
small_psth = [small_1; small_2];
large_psth = [large_1; large_2];

%Convert firing rate to Hz
small_psth = small_psth*(1000/dt_ms);
large_psth = large_psth*(1000/dt_ms);

n_measure = length(small_psth); %Find the total number of points used for the psths

%Compute meand and std
small_psth_mean = nanmean(small_psth);
large_psth_mean = nanmean(large_psth);

small_psth_std = nanstd(small_psth);
large_psth_std = nanstd(large_psth);

%% Compute NPSP based on the noise burst
npsp_start_ms = 0;
npsp_end_ms = 500;
[~, start_ix] = min(abs(npsp_start_ms-t_edges_ms));
[~, end_ix] = min(abs(npsp_end_ms-t_edges_ms));
fprintf('== Calculating NPSP for clusters ==\n');tic;
for k = 1:n_clust
    psth_sahani = [psth_all(k).small_1(:,start_ix:end_ix), psth_all(k).small_2(:,start_ix:end_ix), psth_all(k).large_1(:,start_ix:end_ix), psth_all(k).large_2(:,start_ix:end_ix)];
    [SP, NP, TP, SP_std_error] = sahani_quick(psth_sahani);
    npsp_nb(k) = NP/SP;
end
fprintf('== Done! This took %0.fs ==\n',toc);

[~, ix_select_nb] = sort(npsp_nb,'ascend');
psth_all_plot = psth_all(ix_select_nb);

%% Plot the mean firing rate awith error bars for the two different conditions
if plot_nb
    
    %Plot the mean firing rate +/- SEM
    lw = 3;
    axis_sz = 55;
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    shadedErrorBar(t_edges_ms, small_psth_mean, small_psth_std/sqrt(n_measure), {'LineWidth',lw,'Color', small_col}); hold on;
    shadedErrorBar(t_edges_ms, large_psth_mean, large_psth_std/sqrt(n_measure), {'LineWidth',lw,'Color', large_col});
    annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Large'),'LineStyle','none','Color',large_col,'FontSize',axis_sz,'FontWeight','bold');
    annotation('textbox',[0.8 0.64 0.1 0.1],'String', sprintf('Small'),'LineStyle','none','Color',small_col,'FontSize',axis_sz,'FontWeight','bold');
    xline(0,'k-',{'NB ON'},'LineWidth',lw);
    xline(500,'k-',{'NB OFF'},'LineWidth',lw);
    hold off;
    title('Noise burst all neurons all reps');
    ylabel('Firing rate [Hz]');
    xlabel('Time [ms]');
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    save_name = fullfile(save_dir,['Noise_burst_NPSP<',num2str(NPSP_th),'.svg']);
    saveas(gcf, save_name)
    close;
    
    %% Plot the same but with smaller time window
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    shadedErrorBar(t_edges_ms, small_psth_mean, small_psth_std/sqrt(n_measure), {'LineWidth',lw,'Color', small_col}); hold on;
    shadedErrorBar(t_edges_ms, large_psth_mean, large_psth_std/sqrt(n_measure), {'LineWidth',lw,'Color', large_col});
    xlim([xlim_start_ms xlim_end_ms]);
    ylim([ylim_start_ms ylim_end_ms]);
    annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large'),'LineStyle','none','Color',large_col,'FontSize',axis_sz,'FontWeight','bold');
    annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small'),'LineStyle','none','Color',small_col,'FontSize',axis_sz,'FontWeight','bold');
    xline(0,'k-',{'NB ON'},'LineWidth',lw);
    hold off;
    title('Noise burst all neurons all reps');
    ylabel('Firing rate [Hz]');
    xlabel('Time [ms]');
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    save_name = fullfile(save_dir,['Noise_burst_short_NPSP<',num2str(NPSP_th),'.svg']);
    saveas(gcf, save_name)
    close;
    
end

%% Plot the analysis for all the neurons
if plot_nb_all_neurons
    
    save_dir_ind = fullfile(save_dir,'Ind_neurons');
    
    if ~exist(save_dir_ind, 'dir')
        mkdir(save_dir_ind);
    end
    
    plot_nb_all_neurons(psth_all_plot, t_edges_ms, save_dir_ind)
    
end


%% Compute COM for every cluster for small and large
save_dir_stats = fullfile(save_dir,'Stats');

if ~exist(save_dir_stats, 'dir')
    mkdir(save_dir_stats);
end


[~, start_ix] = min(abs(com_start_ms-t_edges_ms));
[~, end_ix] = min(abs(com_end_ms-t_edges_ms));

[~, baseline_start_ix] = min(abs(baseline_start_ms-t_edges_ms));
[~, baseline_end_ix] = min(abs(baseline_end_ms-t_edges_ms));

ix_select_com = npsp_nb<NPSP_nb_th; %Select the clusters with good NPSP only
npsp_com = npsp_nb(ix_select_com);
psth_all_com = psth_all(ix_select_com);

t_ms = t_edges_ms(start_ix:end_ix); %Get the time to use for the COM calculation
t_ms = t_ms(:); %Force a column vector

n_clust = length(psth_all_com);

com_small = [];
com_large = [];
for k = 1:n_clust
    
    psth_small = [psth_all_com(k).small_1(:,start_ix:end_ix); psth_all_com(k).small_2(:,start_ix:end_ix)];
    psth_large = [psth_all_com(k).large_1(:,start_ix:end_ix); psth_all_com(k).large_2(:,start_ix:end_ix)];
    
    psth_small_mean = mean(psth_small);
    psth_large_mean = mean(psth_large);
    
    if normalize
        psth_baseline_small = mean(mean([psth_all_com(k).small_1(:,baseline_start_ix:baseline_end_ix); psth_all_com(k).small_2(:,baseline_start_ix:baseline_end_ix)]));
        psth_baseline_large = mean(mean([psth_all_com(k).large_1(:,baseline_start_ix:baseline_end_ix); psth_all_com(k).large_2(:,baseline_start_ix:baseline_end_ix)],2));
        
        psth_small_mean = psth_small_mean - psth_baseline_small;
        psth_large_mean = psth_large_mean - psth_baseline_large;
    end
    
    y_small_t = psth_small_mean./sum(psth_small_mean(:)); %Scale the values to sum to 1 for the small room
    y_small_t = y_small_t(:)'; %Force a row vector
    y_large_t = psth_large_mean./sum(psth_large_mean(:)); %Scale the values to sum to 1 for the small room
    y_large_t = y_large_t(:)'; %Force a row vector
    
    com_small(k) = (y_small_t*t_ms); %Compute a weighted sum of all values for small room
    com_large(k) = (y_large_t*t_ms); %Compute a weighted sum of all values for large room
end

[p_val_com,~,~] = signrank(com_large,com_small);

%Compute the median and mean differences
com_diff.median = nanmedian(com_large - com_small);
com_diff.mean = nanmean(com_large - com_small);
%% Define params for histogram plot
hist_type = 'bar'; %The type of the histogram
model = 'noise_burst';
all_font_sz = 55;
lw1 = 7; lw2 = 3;
%Plot general properties
params_norm_his.fit = 'neurons'; params_norm_his.calc_method = 'bf';
params_norm_his.hist_type = hist_type; params_norm_his.get_neurons = get_neurons; params_norm_his.NPSP_th = NPSP_nb_th; params_norm_his.model = model; params_norm_his.font_type = font_type;
%Plot colours
params_norm_his.his_color = 'k';
%Plot size
params_norm_his.all_font_sz = all_font_sz; params_norm_his.lw1 = lw1; params_norm_his.lw2 = lw2;
%Save dir
params_norm_his.save_dir =  save_dir_stats;

%% Define params for scatter plot with NPSP
%Plot general properties
params_npsp_scatter.fit = 'neurons';
params_npsp_scatter.get_neurons = get_neurons; params_npsp_scatter.NPSP_th = NPSP_nb_th; params_npsp_scatter.model = model; params_npsp_scatter.font_type = font_type;
%Plot colours
params_npsp_scatter.NPSPs = npsp_com;
%Plot size
params_npsp_scatter.lw = 3; params_npsp_scatter.all_font_sz = all_font_sz; params_npsp_scatter.sz = all_font_sz;
%Save dir
params_npsp_scatter.save_dir = save_dir_stats;

%% COM for small vs large room Histogram
if plot_histograms
    
    params_norm_his.specific_name = [' Noise burst COM small vs large ','_time_',num2str(com_start_ms),'_',num2str(com_end_ms),'_'];
    params_norm_his.units = 'normal';
    params_norm_his.big_val = com_large; params_norm_his.small_val = com_small;
    params_norm_his.p_val = p_val_com;
    params_norm_his.his_spacing = 8; params_norm_his.lim_val = 60;
    plot_norm_histogram(params_norm_his);
    
end

if plot_scatter
    %% COM Inhibition small vs big with NPSP Scatter
    params_npsp_scatter.specific_name = ['Noise burst COM scatter ','_time_',num2str(com_start_ms),'_',num2str(com_end_ms),'_'];
    params_npsp_scatter.npsp_on = 1;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com_large; params_npsp_scatter.val_small = com_small;
    params_npsp_scatter.p_val = p_val_com;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
    %% COM Inhibition small vs big without NPSP Scatter
    params_npsp_scatter.specific_name = ['Noise burst COM scatter no npsp  ','_time_',num2str(com_start_ms),'_',num2str(com_end_ms),'_'];
    params_npsp_scatter.npsp_on = 0;
    lim_val_ms = 150;
    params_npsp_scatter.val_big = com_large; params_npsp_scatter.val_small = com_small;
    params_npsp_scatter.p_val = p_val_com;
    params_npsp_scatter.lim_val = lim_val_ms;
    plot_npsp_scatter(params_npsp_scatter);
    
end

%% Plot Violin plot of COM
if plot_violins
    
    lw = 3;
    axis_sz = 40;
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    violinplot([com_small', com_large'], {'Small Room', 'Large Room'});
    title('COM noise burst small vs large');
    ylabel('COM [ms]');
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    % ylim([10 200]);
    
    full_name = fullfile(save_dir_stats,['Violin_plot','_NPSP<',num2str(NPSP_nb_th),'_time_',num2str(com_start_ms),'_',num2str(com_end_ms),'.svg']);
    saveas(gcf, full_name);
    close all;
    
end