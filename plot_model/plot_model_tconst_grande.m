function plot_model_tconst_grande(kernel_dir)
%% Params
save_dir = fullfile(kernel_dir,'Plots');
plot_hyperparams = 0;
chop_ms = 190;
r_type{1} = 'small';
r_type{2} = 'med';
r_type{3} = 'big';
n_rooms = length(r_type);
%% Load the data
load(fullfile(kernel_dir,'kernel'));
model = kernels.model;
%% First compute the bf and tau for every cluster
freqs = fliplr(kernels.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_ker = length(freqs);
n_h = kernels.n_h-1;
dt_ms = round(kernels.dt_ms);
chop_ix = round(chop_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;
fprintf('== Calcuating tau and bf ==\n');tic;
for k = 1:n_ker
    for r = 1:n_rooms
        room = r_type{r};
        switch model
            case {'sep','sep_kh'}      
                [~,ix] = max(kernels{k}.(room).k_f);
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
                k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                
            case {'ridge','lasso','manual_alpha','elastic'}
                k_fh = fliplr(kernels.(room){k}.main.k_fh);
                k_fh = k_fh(:,1:chop_ix);
                k_fh_neg = abs(min(k_fh,0));
                k_fh_pos = abs(max(k_fh,0));
                k_h_neg = mean(k_fh_neg);
                k_h_pos = mean(k_fh_pos);
                k_f = mean(k_fh_pos,2); %Take the mean across history steps
                [~,ix] = max(k_f); %Find the index of the max frequency
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
        end
        k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
        k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
        tau_neg(k).(room) = (k_h_neg*h); %Compute a weighted sum of all values
        tau_pos(k).(room) = (k_h_pos*h); %Compute a weighted sum of all values
        if plot_hyperparams
            alpha(k).(room) = kernels.(room){k}.main.alpha;
            lambda(k).(room) = kernels.(room){k}.main.lambda;
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
tau_small_neg = reach(tau_neg,'small');
tau_small_pos = reach(tau_pos,'small');
tau_med_neg = reach(tau_neg,'med');
tau_med_pos = reach(tau_pos,'med');
tau_big_neg = reach(tau_neg,'big');
tau_big_pos = reach(tau_pos,'big');
[p_val_neg,~,~] = signrank(tau_big_neg,tau_small_neg);
[p_val_pos,~,~] = signrank(tau_big_pos,tau_small_pos);

%% Plot the small vs big inhibitory tau w/o NPSP
lim_val = 100;
axis_sz = 20;
sz = 30;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
scatter(tau_small_neg,tau_big_neg,sz,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(tau_small_pos,tau_big_pos,sz,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');
xlabel('\tau_{small} [ms]','FontSize',16,'FontWeight','bold');
ylabel('\tau_{big} [ms]','FontSize',16,'FontWeight','bold');
title(['Inhibitory \tau_{big} vs \tau_{small}  for normative ',model,' kernel']);
annotation('textbox',[0.65 0.2 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.65 0.15 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.65 0.1 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
axis equal;
l = [0 lim_val];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
hold off;
save_name = fullfile(save_dir,['Scatter tau no npsp for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot a histogram of inhibitory and excitatory tau_big vs tau_small as COM change ix
plot_type = 'bar';
max_val = 50;
d_val = 2.5;
lw = 3;
axis_sz = 20;
%Neagtive
tau_diff_neg = tau_big_neg-tau_small_neg;
tau_sum_neg = tau_big_neg + tau_small_neg;
com_change_ix_neg = 100*(tau_diff_neg./tau_sum_neg);
edges = [-max_val:d_val:max_val];
counts_neg = histcounts(com_change_ix_neg,edges);
tau_diff_pos = tau_big_pos-tau_small_pos;
tau_sum_pos = tau_big_pos + tau_small_pos;
com_change_ix = 100*(tau_diff_pos./tau_sum_pos);
counts_pos = histcounts(com_change_ix,edges);
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor','b','EdgeColor','b','DisplayStyle',plot_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor','r','EdgeColor','r','DisplayStyle',plot_type,'LineWidth', lw);
xlabel('100*(\tau_{big}- \tau_{small})/(\tau_{big} + \tau_{small}) [Cetner of mass change index]','FontSize',16,'FontWeight','bold');
ylabel('Number of kernels','FontSize',16,'FontWeight','bold');
title(['Percentage change of Inhibitory \tau  for normative',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
hold off;
save_name = fullfile(save_dir,['Histogram of per change for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot a histogram of inhibitory tau_big vs tau_small in ms
lw = 3;
axis_sz = 20;
extra_ms = 20;
ms_spacing = 2;
lim_val = max(abs([tau_diff_neg(:);tau_diff_pos(:)])) + extra_ms;
edges = [-lim_val:ms_spacing:lim_val];
counts_neg = histcounts(tau_diff_neg,edges);
counts_pos = histcounts(tau_diff_pos,edges);
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor','b','EdgeColor','b','DisplayStyle',plot_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor','r','EdgeColor','r','DisplayStyle',plot_type,'LineWidth', lw);
xlabel('\tau_{big}- \tau_{small} [ms]','FontSize',16,'FontWeight','bold');
ylabel('Number of neurons','FontSize',16,'FontWeight','bold');
title(['Change of inhibitory \tau in ms for normative ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
hold off;
save_name = fullfile(save_dir,['Histogram of tau ms change for normative ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot inhibitory tau vs frequency for all rooms
lw = 3;
% for f = 1:length(freqs)
%     curr_freq = freqs(f);
%     ix = bf_closest==curr_freq;
%     tau_small_mean(f) = nanmean(tau_small_neg(ix));
%     tau_medium_mean(f) = nanmean(tau_med_neg(ix));
%     tau_big_mean(f) = nanmean(tau_big_neg(ix));
%     tau_small_std(f) = nanstd(tau_small_neg(ix));  
%     tau_medium_std(f) = nanstd(tau_med_neg(ix));
%     tau_big_std(f) = nanstd(tau_big_neg(ix));
% end
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
plot(plot_freq,tau_small_neg,'Color',[0.31 0.75 0.405],'LineWidth',lw);hold on;
plot(plot_freq,tau_med_neg,'Color',[0.872 0.496 0.192],'LineWidth',lw);
plot(plot_freq,tau_big_neg,'Color',[0.642 0.456 0.924],'LineWidth',lw);
xlabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
title(['Inhibitory \tau vs frequency for normative ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir,['Inhibitory tau vs freq for normative ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot excitatory tau vs frequency for all rooms
lw = 3;
% for f = 1:length(freqs)
%     curr_freq = freqs(f);
%     ix = bf_closest==curr_freq;
%     tau_small_mean(f) = nanmean(tau_small_pos(ix));
%     tau_small_std(f) = nanstd(tau_small_pos(ix));
%     tau_medium_mean(f) = nanmean(tau_med_pos(ix));
%     tau_medium_std(f) = nanstd(tau_med_pos(ix));
%     tau_big_mean(f) = nanmean(tau_big_pos(ix));
%     tau_big_std(f) = nanstd(tau_big_pos(ix));
% end
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
plot(plot_freq,tau_small_pos,'Color',[0.31 0.75 0.405],'LineWidth',lw);hold on;
plot(plot_freq,tau_med_pos,'Color',[0.872 0.496 0.192],'LineWidth',lw);
plot(plot_freq,tau_big_pos,'Color',[0.642 0.456 0.924],'LineWidth',lw);
xlabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
title(['Excitatory \tau vs frequency for normative ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('med room'),'LineStyle','none','Color',[0.872 0.496 0.192],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',axis_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.65 0.1 0.1],'String', sprintf('n=%0.f',n_ker),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir,['Excitatory tau vs freq for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot bf consistency accross clusters with histograms
bf_small = reach(bf,'small')/1000;
bf_big = reach(bf,'big')/1000;
diff_bf = abs(bf_big - bf_small);
figure('units','normalized','outerposition',[0 0 1 1]);
histogram(diff_bf,'Normalization','probability');
xlabel('\DeltaBF [kHz]','FontSize',16,'FontWeight','bold');
ylabel('Probability','FontSize',16,'FontWeight','bold');
title(['Histogram of BF consistency small<->big  for normative ',model,' kernel']);
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir,['Histogram of BF consistency for normative ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot inhibitory tau vs bf no binning small big separate
% tau_small = reach(tau_neg,'small');
% tau_big = reach(tau_neg,'big');
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(bf_mean,tau_small,'b');hold on;
% plot(bf_mean,tau_big,'r');
% xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
% ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
% title(['Inhibitory \tau vs BF individual rooms  for normative ',model,' kernel']);
% set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
% set(gca, 'XScale', 'log');
% set(gcf,'color','w');
% save_name = fullfile(save_dir,['Inhibitory tau vs BF individual rooms  for normative ',model,' kernel.png']);
% export_fig(save_name);
% close all;
% %% Plot excitatory tau vs bf no binning small big separate
% tau_small = reach(tau_pos,'small');
% tau_big = reach(tau_pos,'big');
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(bf_mean,tau_small,'b');hold on;
% plot(bf_mean,tau_big,'r');
% xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
% ylabel('\tau [ms]','FontSize',16,'FontWeight','bold');
% title(['Excitatory \tau vs BF individual rooms  for normative ',model,' kernel']);
% set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
% set(gca, 'XScale', 'log');
% set(gcf,'color','w');
% save_name = fullfile(save_dir,['Excitatory tau vs BF individual rooms  for normative ',model,' kernel.png']);
% export_fig(save_name);
% close all;

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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\alpha difference [AU]','FontSize',16,'FontWeight','bold');
    title(['\alpha difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\lambda [AU]','FontSize',16,'FontWeight','bold');
    title(['\lambda parameter vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
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
    xlabel('BF [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('\lambda difference [AU]','FontSize',16,'FontWeight','bold');
    title(['\lambda difference vs BF individual rooms  for normative ',model,' kernel']);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gca, 'XScale', 'log');
    legend('Big - Small Room','Big - Medium Room','Medium - Small Room');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['lambda difference vs BF individual rooms  for normative ',model,' kernel.png']);
    export_fig(save_name);
    close all;
    
end
