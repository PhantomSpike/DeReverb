%% Params
slope_file = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis/gain_data.mat';
npsp_file = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms/info.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis/Slope';
save_dir_ind = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis/Plots_LN_ordered';
axis_sz = 30;
col_npsp = 'k';
lw = 4;
col.small = [0.592, 0.737, 0.384];
col.large = [0.173, 0.373, 0.176];

%% Load the data
load(slope_file);
load(npsp_file);
npsp = info.NPSP(:);
slope.small = model_ln.small.slope(:);
slope.big = model_ln.big.slope(:);

%% Calculate the average firing rates
mean_fr.small = cell2mat(cellfun(@mean, y_t_small, 'UniformOutput', false));
mean_fr.big = cell2mat(cellfun(@mean, y_t_big, 'UniformOutput', false));

%% Calculate corr and best line
gain_ratio = slope.big./slope.small;

%NPSP
%Find corr coef 
[r.npsp, p.npsp] = corrcoef(npsp,gain_ratio);

%Fit the best line
X = [ones(size(npsp)), npsp];
y = gain_ratio;
b.npsp = regress(y,X);

%Firing rate
%Find corr coef 
[r.fr, p.fr] = corrcoef(mean_fr.big, gain_ratio);

%Fit the best line
X = [ones(size(mean_fr.big)), mean_fr.big];
y = gain_ratio;
b.fr = regress(y,X);

%% Plot NPSP
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(npsp, gain_ratio,'MarkerEdgeColor',col_npsp,'MarkerFaceColor',col_npsp); hold on;
h1 = refline(b.npsp(2), b.npsp(1));
h1.Color = 'r'; h1.LineWidth = lw;
yline(1,'--b','LineWidth',lw);
xlabel('NPSP');
ylabel('Gain Ratio');
title('Gain Ratio vs NPSP');
set(gcf,'color','w');

annotation('textbox',[0.75 0.8 0.1 0.1],'String', sprintf('r=%0.3f',r.npsp(1,2)),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('p=%0.3f',p.npsp(1,2)),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');

set(gca,'FontSize',axis_sz,'FontWeight','Normal');
file_name = fullfile(save_dir,'Gain ratio vs NPSP.svg');
saveas(gcf,file_name);
close;


%% Plot Firing rate
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(mean_fr.big, gain_ratio,'MarkerEdgeColor',col_npsp,'MarkerFaceColor',col_npsp); hold on;
h1 = refline(b.fr(2), b.fr(1));
h1.Color = 'r'; h1.LineWidth = lw;
yline(1,'--b','LineWidth',lw);
xlabel('Firing rate Large room');
ylabel('Gain Ratio');
title('Gain Ratio vs Firing rate Large room');
set(gcf,'color','w');

annotation('textbox',[0.75 0.8 0.1 0.1],'String', sprintf('r=%0.3f',r.fr(1,2)),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('p=%0.3f',p.fr(1,2)),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');

set(gca,'FontSize',axis_sz,'FontWeight','Normal');
file_name = fullfile(save_dir,'Gain ratio vs Firing rate in Large Room.svg');
saveas(gcf,file_name);
close;

%% Sort by firing rate of the large room

[~,ix] = sort(mean_fr.big,'descend');
y_hat_lin_t_small_ord = y_hat_lin_t_small(ix);
y_hat_lin_t_big_ord = y_hat_lin_t_big(ix);
y_hat_ln_t_small_ord = y_hat_ln_t_small(ix);
y_hat_ln_t_big_ord = y_hat_ln_t_big(ix);
y_t_small_ord = y_t_small(ix);
y_t_big_ord = y_t_big(ix);

%% Plot ind neurons 
lim_val = 0.4;

rows = 5;
cols = 6;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = 0.05; space_h = 0.06; space_v = 0.09;
[pos]=subplot_pos(rows,cols,edgel,edger,edgeh,edgeb,space_h,space_v);
clust_per_batch = rows*cols;
axis_sz = 20;
n_clust = length(y_t_small);
num_batch = ceil(n_clust/clust_per_batch);
curr_n_clust = clust_per_batch;

count = 0;


for b = 1:num_batch
    figure('units','normalized','outerposition',[0 0 1 1]);
    save_name = fullfile(save_dir_ind,['batch_',num2str(b),'_LN_model.svg']);
    
    if b == num_batch
        curr_n_clust = n_clust - (num_batch-1)*clust_per_batch;
    end
    
    for j = 1:curr_n_clust
        count = count+1;
        subplot('position',pos{j});
        [x, y] = binplot(y_hat_lin_t_small_ord{count}, y_t_small_ord{count}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.small,'MarkerFaceColor',col.small); hold on;
        [x, y] = binplot(y_hat_lin_t_small_ord{count}, y_hat_ln_t_small_ord{count}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.small, 'LineWidth',lw);
        [x, y] = binplot(y_hat_lin_t_big_ord{count}, y_t_big_ord{count}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.large,'MarkerFaceColor',col.large);
        [x, y] = binplot(y_hat_lin_t_big_ord{count}, y_hat_ln_t_big_ord{count}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.large, 'LineWidth',lw);
        xlim([0 lim_val]);
        ylim([0 lim_val]);
        xlabel('yhat');
        ylabel('y');
        set(gcf,'color','w');
        hold off;
        set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    end
    saveas(gcf,save_name);
    close;
 
end