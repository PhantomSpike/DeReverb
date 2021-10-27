function plot_model_tconst_grande_cv_all(kernel_dir)
%% Params
temp_dir = fullfile(kernel_dir,'Plots');
plot_hyperparams = 0;
chop_ms = 190;
hist_type = 'bar';
r_type{1} = 'small';
r_type{2} = 'med';
r_type{3} = 'big';
n_rooms = length(r_type);
%% Load the data
load(fullfile(kernel_dir,'kernel'));
model = kernels.model;
k_folds = length(kernels.small{1}.main)-1;
if ~exist(temp_dir,'dir')
    mkdir(temp_dir);
end
save_dir_full = fullfile(temp_dir,'all_parts');
if ~exist(save_dir_full,'dir')
    mkdir(save_dir_full);
end

fprintf('== Calcuating tau and bf ==\n');tic;
for kf = 1:k_folds
    %% First compute the bf and tau for every cluster
    freqs = fliplr(kernels.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
    n_ker = length(freqs);
    n_h = kernels.n_h-1;
    dt_ms = round(kernels.dt_ms);
    chop_ix = round(chop_ms/dt_ms);
    h = (1:1:chop_ix)';
    h = dt_ms*h;
    for k = 1:n_ker
        for r = 1:n_rooms
            room = r_type{r};
            switch model
                case {'sep','sep_kh'}
                    [~,ix] = max(kernels{k}.(room).k_f);
                    bf(k).(room) = freqs(ix); %Find the corresponding frequency
                    k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                    
                case {'ridge','lasso','manual_alpha','elastic'}
                    k_fh = fliplr(kernels.(room){k}.main{kf}.k_fh);
                    k_fh = k_fh(:,1:chop_ix);
                    k_fh_neg = abs(min(k_fh,0));
                    k_fh_pos = abs(max(k_fh,0));
                    k_h_neg = mean(k_fh_neg);
                    k_h_pos = mean(k_fh_pos);
                    k_f = mean(k_fh_pos,2); %Take the mean across history steps
                    [~,ix] = max(k_f); %Find the index of the max frequency
                    bf(k,kf).(room) = freqs(ix); %Find the corresponding frequency
            end
            k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
            k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
            tau_neg(k,kf).(room) = (k_h_neg*h); %Compute a weighted sum of all values
            tau_pos(k,kf).(room) = (k_h_pos*h); %Compute a weighted sum of all values
            if plot_hyperparams
                alpha(k,kf).(room) = kernels.(room){k}.main{kf}.alpha;
                lambda(k,kf).(room) = kernels.(room){k}.main{kf}.lambda;
            end
        end
        bf_mean(k,kf) = mean([bf(k,kf).small,bf(k,kf).big]);
        [~,ix_bf] = min(abs(bf_mean(k,kf) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
        bf_closest(k,kf) = freqs(ix_bf);
    end
    
end

%% Stats
tau_small_neg = reach(tau_neg,'small');
tau_small_pos = reach(tau_pos,'small');
tau_med_neg = reach(tau_neg,'med');
tau_med_pos = reach(tau_pos,'med');
tau_big_neg = reach(tau_neg,'big');
tau_big_pos = reach(tau_pos,'big');

%Agregate everything together for each room
tau_small_neg_all = tau_small_neg(:); tau_big_neg_all = tau_big_neg(:);
tau_med_pos_all = tau_med_pos(:); tau_med_pos_all = tau_med_pos(:);
tau_small_pos_all = tau_small_pos(:); tau_big_pos_all = tau_big_pos(:);

%Mean and std for all rooms
tau_small_neg_mean = nanmean(tau_small_neg,2); tau_med_neg_mean = nanmean(tau_med_neg,2); tau_big_neg_mean = nanmean(tau_big_neg,2);
tau_small_neg_std = nanstd(tau_small_neg,[],2); tau_med_neg_std = nanstd(tau_med_neg,[],2); tau_big_neg_std = nanstd(tau_big_neg,[],2);
tau_small_pos_mean = nanmean(tau_small_pos,2); tau_med_pos_mean = nanmean(tau_med_pos,2); tau_big_pos_mean = nanmean(tau_big_pos,2);
tau_small_pos_std = nanstd(tau_small_pos,[],2); tau_med_pos_std = nanstd(tau_med_pos,[],2); tau_big_pos_std = nanstd(tau_big_pos,[],2);

%Compute stats
[p_val_neg,~,~] = signrank(tau_big_neg_all,tau_small_neg_all);
[p_val_pos,~,~] = signrank(tau_big_pos_all,tau_small_pos_all);
n_total = numel(tau_small_neg_all);
%% Plot the small vs big inhibitory tau w/o NPSP
lim_val_ms = 150;
font_type = 'Liberation Sans';
sz = 35;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_neg_all,tau_big_neg_all,sz,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
scatter(tau_small_pos_all,tau_big_pos_all,sz,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
% title(['Inhibitory \tau_{big} vs \tau_{small}  for normative ',model,' kernel']);
% annotation('textbox',[0.65 0.2 0.1 0.1],'String', sprintf('n=%0.f',n_total),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.65 0.15 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.65 0.1 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
hold off;
save_name = fullfile(save_dir_full,['Scatter tau no npsp for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot a histogram of inhibitory and excitatory tau_big vs tau_small as COM change ix
d_val = 2;
max_val = 60;
lw = 3;
axis_sz = 20;
%Neagtive
tau_diff_neg = tau_big_neg_all-tau_small_neg_all;
tau_sum_neg = tau_big_neg_all + tau_small_neg_all;
com_change_ix_neg = 100*(tau_diff_neg./tau_sum_neg);
edges = [-max_val:d_val:max_val];
counts_neg = histcounts(com_change_ix_neg,edges);
tau_diff_pos = tau_big_pos_all-tau_small_pos_all;
tau_sum_pos = tau_big_pos_all + tau_small_pos_all;
com_change_ix = 100*(tau_diff_pos./tau_sum_pos);
counts_pos = histcounts(com_change_ix,edges);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

switch hist_type
    case 'stairs'
        fc_neg = 'none';
        fc_pos = 'none';
    case 'bar'
        fc_neg = 'b';
        fc_pos = 'r';
end
        
hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor',fc_neg,'EdgeColor','b','DisplayStyle',hist_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor',fc_pos,'EdgeColor','r','DisplayStyle',hist_type,'LineWidth', lw);
xlabel('100*(\tau_{large}- \tau_{small})/(\tau_{big} + \tau_{small}) [Cetner of mass change index]','FontSize',16,'FontWeight','bold');
ylabel('Number of kernels','FontSize',16,'FontWeight','bold');
title(['Percentage change of Inhibitory \tau  for normative',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_total),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
hold off;
save_name = fullfile(save_dir_full,['Histogram ',hist_type,' of per change for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot a histogram of inhibitory tau_big vs tau_small in ms
lw = 3;
axis_sz = 20;
ms_spacing = 2;
lim_val_ms = 60;
edges = [-lim_val_ms:ms_spacing:lim_val_ms];
counts_neg = histcounts(tau_diff_neg,edges);
counts_pos = histcounts(tau_diff_pos,edges);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor',fc_neg,'EdgeColor','b','DisplayStyle',hist_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor',fc_pos,'EdgeColor','r','DisplayStyle',hist_type,'LineWidth', lw);
xlabel('\tau_{large}- \tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('Number of model kernels','FontSize',y_font_sz,'FontWeight','bold');
% title(['Change of inhibitory and excitatory \tau in ms for normative ',model,' kernel']);
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_total),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
xlim([-lim_val_ms, lim_val_ms]);
hold off;
save_name = fullfile(save_dir_full,['Histogram ',hist_type,' of tau ms change for normative ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot inhibitory tau vs frequency for all rooms
lw = 3;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);hold on;
plot_freq = freqs;
shadedErrorBar(plot_freq,tau_small_neg_mean,tau_small_neg_std./sqrt(k_folds),{'Color',[0.31 0.75 0.405]});
shadedErrorBar(plot_freq,tau_med_neg_mean,tau_med_neg_std./sqrt(k_folds),{'Color',[0.872 0.496 0.192]});
shadedErrorBar(plot_freq,tau_big_neg_mean,tau_big_neg_std./sqrt(k_folds),{'Color',[0.642 0.456 0.924]});
hold off;
xlabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
title(['Inhibitory \tau vs frequency for normative ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_total),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir_full,['Inhibitory tau vs freq for normative ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot excitatory tau vs frequency for all rooms
lw = 3;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); hold on;
plot_freq = freqs;
shadedErrorBar(plot_freq,tau_small_pos_mean,tau_small_pos_std./sqrt(k_folds),{'Color',[0.31 0.75 0.405]});
shadedErrorBar(plot_freq,tau_med_pos_mean,tau_med_pos_std./sqrt(k_folds),{'Color',[0.872 0.496 0.192]});
shadedErrorBar(plot_freq,tau_big_pos_mean,tau_big_pos_std./sqrt(k_folds),{'Color',[0.642 0.456 0.924]});
hold off;
xlabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
title(['Excitatory \tau vs frequency for normative ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_total),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir_full,['Excitatory tau vs freq for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot bf consistency accross clusters with histograms
bf_small = reach(bf,'small')/1000;
bf_big = reach(bf,'big')/1000;
diff_bf = abs(bf_big - bf_small);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
histogram(diff_bf,'Normalization','probability');
xlabel('\DeltaBF [kHz]','FontSize',16,'FontWeight','bold');
ylabel('Probability','FontSize',16,'FontWeight','bold');
title(['Histogram of BF consistency small<->big  for normative ',model,' kernel']);
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Histogram of BF consistency for normative ',model,' kernel.png']);
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\alpha [AU]','FontSize',16,'FontWeight','bold');
    title(['\alpha parameter vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gca, 'XScale', 'log');
    legend('Big Room','Medium Room','Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir_full,['alpha vs BF individual rooms  for normative ',model,' kernel.png']);
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\alpha difference [AU]','FontSize',16,'FontWeight','bold');
    title(['\alpha difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gca, 'XScale', 'log');
    legend('Big - Small Room','Big - Medium Room','Medium - Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir_full,['alpha difference vs BF individual rooms  for normative ',model,' kernel.png']);
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\lambda [AU]','FontSize',16,'FontWeight','bold');
    title(['\lambda parameter vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gca, 'XScale', 'log');
    legend('Big Room','Medium Room','Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir_full,['lambda vs BF individual rooms  for normative ',model,' kernel.png']);
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\lambda difference [AU]','FontSize',16,'FontWeight','bold');
    title(['\lambda difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gca, 'XScale', 'log');
    legend('Big - Small Room','Big - Medium Room','Medium - Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir_full,['lambda difference vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
end
