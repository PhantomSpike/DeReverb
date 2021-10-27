%% Params
data_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis'; %Absolute path to the kernel used for the prediction

col.small = [0.592, 0.737, 0.384];
col.large = [0.173, 0.373, 0.176];
all_font_sz = 35;
axis_sz = all_font_sz;
lw = 3;

%% Load the data
file_name = fullfile(data_dir, 'gain_data.mat');
load(file_name);

%% Histogram of Slopes

%No selection
bin_width = 0.1;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(model_ln.small.slope,'BinWidth',bin_width,'FaceColor',col.small,'EdgeColor',col.small); hold on;
histogram(model_ln.big.slope,'BinWidth',bin_width,'FaceColor',col.large,'EdgeColor',col.large);
xlabel('Slope at midpoint','FontWeight','normal');
ylabel('Number of neurons','FontWeight','normal');
title('Histogram of Slope at midpoint');

annotation('textbox',[0.8 0.80 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.72 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(data_dir,['Histogram of Slope at midpoint','.svg']);
saveas(gcf, file_name);
close;

%% Violin plots of CCnorm

figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'Small','Large'};
violins = violinplot([model_ln.small.slope(:), model_ln.big.slope(:)], int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', true, 'ShowNotches', false);
violins(1).ViolinColor = col.small;
violins(2).ViolinColor = col.large;
ylim([0 8]);
title('Slope at midpoint violin plots');
ylabel(['Slope at midpoint']);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',all_font_sz,'FontWeight','Normal');
save_name = fullfile(data_dir, 'Violin plots Slope at midpoint.svg');
saveas(gcf, save_name);
close;

%% Compute stats

[pval_all,~,~] = signrank(model_ln.big.slope, model_ln.small.slope);

median_diff.all = nanmedian(model_ln.big.slope - model_ln.small.slope);

mean_diff.all = nanmean(model_ln.big.slope - model_ln.small.slope);


%% Histogram of differences in CCnorm

bin_width = 0.2;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(model_ln.big.slope - model_ln.small.slope,'BinWidth',bin_width,'FaceColor','k','EdgeColor','k'); hold on;
xlabel('Difference in Slope at midpoint Large - Small','FontWeight','normal');
ylabel('Number of neurons','FontWeight','bold');
title('Histogram of the difference in Slope at midpoint');
xline(0,'r','LineWidth',3);
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('%s',p_star(pval_all)),'LineStyle','none','FontSize',all_font_sz,'FontWeight','bold');
xlim([-6 6])
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
hold off;
file_name = fullfile(data_dir,['Histogram of differences in Slope at midpoint','.svg']);
saveas(gcf, file_name);
close;

%% Plot mean +/- SEM Slopes
lw = 5;
n_tp = length(y_hat_lin_t_small{1});
n_clust = length(y_hat_lin_t_small);

%Get the data all together
y_hat_lin_t_small_all = reshape(cell2mat(y_hat_lin_t_small),n_tp,n_clust)';
y_hat_lin_t_big_all = reshape(cell2mat(y_hat_lin_t_big),n_tp,n_clust)';

y_hat_ln_t_small_all = reshape(cell2mat(y_hat_ln_t_small),n_tp,n_clust)';
y_hat_ln_t_big_all = reshape(cell2mat(y_hat_ln_t_big),n_tp,n_clust)';

y_t_small_all = reshape(cell2mat(y_t_small),n_tp,n_clust)';
y_t_big_all = reshape(cell2mat(y_t_big),n_tp,n_clust)';

%Compute the mean
y_hat_lin_t_small_mean = mean(y_hat_lin_t_small_all);
y_hat_lin_t_big_mean = mean(y_hat_lin_t_big_all);

y_hat_ln_t_small_mean = mean(y_hat_ln_t_small_all);
y_hat_ln_t_big_mean = mean(y_hat_ln_t_big_all);

y_t_small_all_mean = mean(y_t_small_all);
y_t_big_all_mean = mean(y_t_big_all);

%Compute the SEM
y_hat_lin_t_small_sem = std(y_hat_lin_t_small_all)/sqrt(n_clust);
y_hat_lin_t_big_sem = std(y_hat_lin_t_big_all)/sqrt(n_clust);

y_hat_ln_t_small_sem = std(y_hat_ln_t_small_all)/sqrt(n_clust);
y_hat_ln_t_big_sem = std(y_hat_ln_t_big_all)/sqrt(n_clust);

y_t_small_all_sem = std(y_t_small_all)/sqrt(n_clust);
y_t_big_all_sem = std(y_t_big_all)/sqrt(n_clust);

figure('units','normalized','outerposition',[0 0 1 1]);

[x, y] = binplot(y_hat_lin_t_small_mean, y_hat_ln_t_small_mean, 100); plot(x(1:end-1),y(1:end-1),'Color',col.small, 'LineWidth',lw); hold on;
[x, y] = binplot(y_hat_lin_t_big_mean, y_hat_ln_t_big_mean, 100); plot(x(1:end-1),y(1:end-1),'Color',col.large, 'LineWidth',lw);
% [x, y] = binplot(y_hat_lin_t_small_mean+y_hat_lin_t_small_sem, y_hat_ln_t_small_mean+y_hat_ln_t_small_sem, 100); plot(x(1:end-1),y(1:end-1),'--','Color',col.small, 'LineWidth',lw); 
% [x, y] = binplot(y_hat_lin_t_small_mean-y_hat_lin_t_small_sem, y_hat_ln_t_small_mean-y_hat_ln_t_small_sem, 100); plot(x(1:end-1),y(1:end-1),'--','Color',col.small, 'LineWidth',lw); 
% [x, y] = binplot(y_hat_lin_t_big_mean+y_hat_lin_t_big_sem, y_hat_ln_t_big_mean+y_hat_ln_t_big_sem, 100); plot(x(1:end-1),y(1:end-1),'--','Color',col.large, 'LineWidth',lw); 
% [x, y] = binplot(y_hat_lin_t_big_mean-y_hat_lin_t_big_sem, y_hat_ln_t_big_mean-y_hat_ln_t_big_sem, 100); plot(x(1:end-1),y(1:end-1),'--','Color',col.large, 'LineWidth',lw); 
annotation('textbox',[0.7 0.2 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',col.large,'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.7 0.12 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',col.small,'FontSize',axis_sz,'FontWeight','bold');
xlabel('yhat');
ylabel('y');
title('Output non-linearity averaged across all neurons')
set(gcf,'color','w');
hold off;
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
file_name = fullfile(data_dir,['Output non-lineairty across all neurons','.svg']);
saveas(gcf, file_name);
close;
