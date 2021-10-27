%% Define params
data_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/For_analysis/All_data';
coch_file = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch/Switch_Noah_Derekah/specpower/10ms/coch_all_conditions_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/Check_psth';
get_neurons = 'all';
NPSP_th = 100;
dt_ms = 100;

t_edges_s = [0:dt_ms/1000:40];
fr_convert = 1000/dt_ms;
n_stim = 8;
n_reps = 10;
load(coch_file,'coch');
% t_edges_s = coch(1).t;
%% Select the data
load(fullfile(data_dir,'info'),'info');

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

%% Load the desired data and compute PSTH
fprintf('== Loading the data ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(data_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'data');
    for s = 1:n_stim
        for r = 1:n_reps
            temp_psth(r,:) = histc(data.stim(s).repeat(r).spiketimes, t_edges_s);
        end
       mean_temp(s,:)  = mean(temp_psth);
    end
    psth_clusters{k,1}.small_start = mean_temp(1:4,:);
    psth_clusters{k,1}.large_start = mean_temp(5:8,:);
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Compute global psths and plot them
small_start_all_psth = reach(cell2mat(psth_clusters),'small_start');
large_start_all_psth = reach(cell2mat(psth_clusters),'large_start');

lw = 3;
lw_ref = 2;
all_font_sz = 30;
figure('units','normalized','outerposition',[0 0 1 1]);
stairs(t_edges_s, mean(small_start_all_psth)*fr_convert,'LineWidth',lw);
hold on;
xline(8,'r-',{'SL'},'LineWidth',lw_ref);
xline(16,'r-',{'LS'},'LineWidth',lw_ref);
xline(24,'r-',{'SL'},'LineWidth',lw_ref);
xline(32,'r-',{'LS'},'LineWidth',lw_ref);
xlabel('Time [s]');
ylabel('Firing rate [Hz]');
ylim([2.5 9]);
hold off;
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir, 'Small_start.svg');
saveas(gcf,save_name);
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
stairs(t_edges_s, mean(large_start_all_psth)*fr_convert,'LineWidth',lw);
hold on;
xline(8,'r-',{'LS'},'LineWidth',lw_ref);
xline(16,'r-',{'SL'},'LineWidth',lw_ref);
xline(24,'r-',{'LS'},'LineWidth',lw_ref);
xline(32,'r-',{'SL'},'LineWidth',lw_ref);
xlabel('Time [s]');
ylabel('Firing rate [Hz]');
ylim([2.5 9]);
hold off;
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
save_name = fullfile(save_dir, 'Large_start.svg');
saveas(gcf,save_name);
close all;