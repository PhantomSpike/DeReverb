%% Define params

dt_ms = 10; %The binning for the psths
NPSP_th = 40;
t_start_s = 0;
t_end_s = 36;

save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis';
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Directory with all the clusters

% if ~exist('n_cores','var') || isempty(n_cores)
%     n_cores = 20;
% end

dt_s = dt_ms/1000;
t_edges_s = [t_start_s:dt_s:t_end_s]; %Make the time edges for the histogram

col.small = [0.592, 0.737, 0.384];
col.large = [0.173, 0.373, 0.176];
all_font_sz = 35;
axis_sz = all_font_sz;

%% Select the clusters to plot
load(fullfile(cluster_dir,'info.mat'),'info'); %Load the info file with the meta data
ix = info.NPSP<NPSP_th; %Find only the neurons that are below the NPSP th
%Get the necessary fields
cluster_ids = info.cluster_id(ix);
pen_names = info.pen_name(ix);
animal_names = info.animal_name(ix);
qualities = info.quality(ix);
NPSPs = info.NPSP(ix);
n_clusters = sum(ix);

%% Get the data
% delete(gcp('nocreate'));
% parpool('local',n_cores);

psth_all.small.stim1 = zeros(n_clusters, length(t_edges_s));
psth_all.small.stim2 = zeros(n_clusters, length(t_edges_s));
psth_all.large.stim1 = zeros(n_clusters, length(t_edges_s));
psth_all.large.stim2 = zeros(n_clusters, length(t_edges_s));

fprintf('== Making the PSTHs ==\n'); tic;

for c = 1:n_clusters
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
    %Load the data
    temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    animal_name = temp.data.params.animal_name;
    y = temp.data;
    n_stim = length(y.stim);
    n_rep = length(y.stim(1).repeat);
    
    for s = 1:n_stim
        psth_temp = [];
        for r = 1:n_rep
            psth_temp(r,:) = histc(y.stim(s).repeat(r).spiketimes, t_edges_s);
        end
        y_temp(s,:) = mean(psth_temp);
    end
    
    if ismember(animal_name, ["Ronnie","PLP","Cork","Kilkenny","Derry"])
        
        y_t_small_1 = y_temp(3,:);
        y_t_small_2 = y_temp(4,:);
        y_t_large_1 = y_temp(5,:);
        y_t_large_2 = y_temp(6,:);
        %         y_t_anech = [y_temp(1,:),y_temp(2,:)];
        
    elseif ismember(animal_name, ["Noah","Derekah"])
        
        y_t_small_1 = y_temp(1,:);
        y_t_small_2 = y_temp(2,:);
        y_t_large_1 = y_temp(5,:);
        y_t_large_2 = y_temp(6,:);
        %         y_t_med = [y_temp(3,:), y_temp(4,:)];
        
    else
        error('Unrecognized animal');
    end
    
    psth_all.small.stim1(c,:) = y_t_small_1;
    psth_all.small.stim2(c,:) = y_t_small_2;
    psth_all.large.stim1(c,:) = y_t_large_1;
    psth_all.large.stim2(c,:) = y_t_large_2;

end

fprintf('== Done! Fitting all kernels took %0.fs ==\n',toc);

%% Average firing rates

mean_fr.small.stim1 = mean(psth_all.small.stim1,2)*(1000/dt_ms);
mean_fr.small.stim2 = mean(psth_all.small.stim2,2)*(1000/dt_ms);
mean_fr.large.stim1 = mean(psth_all.large.stim1,2)*(1000/dt_ms);
mean_fr.large.stim2 = mean(psth_all.large.stim2,2)*(1000/dt_ms);

%% Histogram of firing rates

%Stim 1
bin_width = 0.5;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(mean_fr.small.stim1,'BinWidth',bin_width,'FaceColor',col.small,'EdgeColor',col.small); hold on;
histogram(mean_fr.large.stim1,'BinWidth',bin_width,'FaceColor',col.large,'EdgeColor',col.large);
xlabel('Average firing rate (Hz)','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of average firing rates Stim 1');

annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram average firing rates stim 1','.svg']);
saveas(gcf, file_name);
close;

%Stim 2
bin_width = 0.5;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(mean_fr.small.stim2,'BinWidth',bin_width,'FaceColor',col.small,'EdgeColor',col.small); hold on;
histogram(mean_fr.large.stim2,'BinWidth',bin_width,'FaceColor',col.large,'EdgeColor',col.large);
xlabel('Average firing rate (Hz)','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of average firing rates Stim 2');

annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram average firing rates stim 2','.svg']);
saveas(gcf, file_name);
close;

%% Violin plots of firing rates

%Stim 1
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Large'};
violins = violinplot([mean_fr.small.stim1, mean_fr.large.stim1], int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = col.small;
violins(2).ViolinColor = col.large;
ylim([0 20]);
title('Average firing rates Stim 1');
ylabel(['Average firing rate (Hz)']);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',all_font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin plots firing rates stim1.svg');
saveas(gcf, save_name);
close;

%Stim 2
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Large'};
violins = violinplot([mean_fr.small.stim2, mean_fr.large.stim2], int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = col.small;
violins(2).ViolinColor = col.large;
ylim([0 20]);
title('Average firing rates Stim 2');
ylabel(['Average firing rate (Hz)']);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',all_font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin plots firing rates stim2.svg');
saveas(gcf, save_name);
close;

%% Compute stats

[pval_stim1,~,~] = signrank(mean_fr.large.stim1, mean_fr.small.stim1);
[pval_stim2,~,~] = signrank(mean_fr.large.stim2, mean_fr.small.stim2);
median_diff.stim1 = nanmedian(mean_fr.large.stim1 - mean_fr.small.stim1);
median_diff.stim2 = nanmedian(mean_fr.large.stim2 - mean_fr.small.stim2);
mean_diff.stim1 = nanmean(mean_fr.large.stim1 - mean_fr.small.stim1);
mean_diff.stim2 = nanmean(mean_fr.large.stim2 - mean_fr.small.stim2);

%% Histogram of differences in firing rates

%Stim 1
bin_width = 0.25;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(mean_fr.large.stim1 - mean_fr.small.stim1,'BinWidth',bin_width,'FaceColor','k','EdgeColor','k'); hold on;
xlabel('Difference in average firing rate Large - Small (Hz)','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of the difference in average firing rates Stim 1');
xline(0,'r','LineWidth',3);
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('%s',p_star(pval_stim1)),'LineStyle','none','FontSize',all_font_sz,'FontWeight','bold');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram of differences in firing rates stim 1','.svg']);
saveas(gcf, file_name);
close;

%Stim 2
bin_width = 0.25;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(mean_fr.large.stim2 - mean_fr.small.stim2,'BinWidth',bin_width,'FaceColor','k','EdgeColor','k'); hold on;
xlabel('Difference in average firing rate Large - Small (Hz)','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of the difference in average firing rates Stim 2');
xline(0,'r','LineWidth',3);
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('%s',p_star(pval_stim1)),'LineStyle','none','FontSize',all_font_sz,'FontWeight','bold');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram of differences in firing rates stim 2','.svg']);
saveas(gcf, file_name);
close;