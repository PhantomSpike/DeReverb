%% Params
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis';
col.small = [0.592, 0.737, 0.384];
col.large = [0.173, 0.373, 0.176];
all_font_sz = 35;
axis_sz = all_font_sz;

%% Load the data

file_name = fullfile(kernel_dir, 'info');
load(file_name);
CCnorm_small = info.CCnorm_mean.small;
CCnorm_large = info.CCnorm_mean.big;
hasNaN_small = info.CCnorm_hasNaN.small;
hasNaN_large = info.CCnorm_hasNaN.big;

%Select the indices of neurons which do not have any NaNs in either the
%Small and Large room

ix = ~(hasNaN_small | hasNaN_large);

%% Histogram of CCnorm

%No selection
bin_width = 0.05;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(CCnorm_small,'BinWidth',bin_width,'FaceColor',col.small,'EdgeColor',col.small); hold on;
histogram(CCnorm_large,'BinWidth',bin_width,'FaceColor',col.large,'EdgeColor',col.large);
xlabel('Average CCnorm','FontWeight','normal');
ylabel('Number of neurons','FontWeight','normal');
title('Histogram of average CCnorm');

annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram average CCnorm','.svg']);
saveas(gcf, file_name);
close;

%Selected without NaNs
bin_width = 0.025;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(CCnorm_small(ix),'BinWidth',bin_width,'FaceColor',col.small,'EdgeColor',col.small); hold on;
histogram(CCnorm_large(ix),'BinWidth',bin_width,'FaceColor',col.large,'EdgeColor',col.large);
xlabel('Average CCnorm','FontWeight','normal');
ylabel('Number of neurons','FontWeight','normal');
title('Histogram of average CCnorm no NaNs');

annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram average CCnorm no NaNs','.svg']);
saveas(gcf, file_name);
close;

%% Violin plots of CCnorm

%Stim 1
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Large'};
violins = violinplot([CCnorm_small, CCnorm_large], int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = col.small;
violins(2).ViolinColor = col.large;
% ylim([0 20]);
title('Average CCnorm violin plots');
ylabel(['Average CCnorm']);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',all_font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin plots CCnorm.svg');
saveas(gcf, save_name);
close;

%Stim 2
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Large'};
violins = violinplot([CCnorm_small(ix), CCnorm_large(ix)], int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = col.small;
violins(2).ViolinColor = col.large;
% ylim([0 20]);
title('Average CCnorm violin plots no NaNs');
ylabel(['Average CCnorm']);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',all_font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin plots CCnorm no NaNs.svg');
saveas(gcf, save_name);
close;

%% Compute stats

[pval_all,~,~] = signrank(CCnorm_large, CCnorm_small);
[pval_no_nans,~,~] = signrank(CCnorm_large(ix), CCnorm_small(ix));
median_diff.all = nanmedian(CCnorm_large - CCnorm_small);
median_diff.no_nans = nanmedian(CCnorm_large(ix) - CCnorm_small(ix));
mean_diff.all = nanmean(CCnorm_large - CCnorm_small);
mean_diff.no_nans  = nanmean(CCnorm_large(ix) - CCnorm_small(ix));

%% Histogram of differences in CCnorm

%Stim 1
bin_width = 0.05;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(CCnorm_large - CCnorm_small,'BinWidth',bin_width,'FaceColor','k','EdgeColor','k'); hold on;
xlabel('Difference in CCnorm Large - Small','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of the difference in CCnorm');
xline(0,'r','LineWidth',3);
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('%s',p_star(pval_all)),'LineStyle','none','FontSize',all_font_sz,'FontWeight','bold');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram of differences in CCnorm','.svg']);
saveas(gcf, file_name);
close;

%Stim 2
bin_width = 0.04;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(CCnorm_large(ix) - CCnorm_small(ix),'BinWidth',bin_width,'FaceColor','k','EdgeColor','k'); hold on;
xlabel('Difference in CCnorm Large - Small','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of the difference in CCnorm no NaNs');
xline(0,'r','LineWidth',3);
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('%s',p_star(pval_no_nans)),'LineStyle','none','FontSize',all_font_sz,'FontWeight','bold');

set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(save_dir,['Histogram of differences in CCnorm no NaNs','.svg']);
saveas(gcf, file_name);
close;