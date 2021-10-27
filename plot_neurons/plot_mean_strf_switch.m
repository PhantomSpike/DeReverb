function plot_mean_strf_switch(k_fhn,params)
%% Plotting params
freqs = params.freqs;
n_f = length(freqs);
n_h = params.n_h+1;
dt = params.dt_ms;
h_max_ms = params.h_max_ms+dt;
save_dir = params.save_dir;
h_steps = [0:dt:h_max_ms];
h_plot = [0:dt:h_max_ms];
font_type = params.font_type;
exc_small_color = params.exc_small_color; exc_big_color = params.exc_big_color;
inh_small_color = params.inh_small_color; inh_big_color = params.inh_big_color;
lw = 5;
axis_sz = 91.5;
legend_font_sz = 32;
title_sz = 15;
h_lim_ms = h_max_ms;
h_lim = round(h_lim_ms/dt);
skip_f = 5;
skip_h = 5;
skip_h_ms = 50;

for f = 1:n_f
    f_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    y_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for h = 1:skip_h:n_h+1
    count = count+1;
    h_labels{count} = num2str(h_steps(h),'%.0f');
end

rooms{1} = 'sl1';
rooms{2} = 'sl2';
rooms{3} = 'ls1';
rooms{4} = 'ls2';

switch params.type
    case 'bf'
        k_fhn_sl1 =  reshape(cell2mat(k_fhn.sl1),n_f,n_h-1,n_f);
        k_fhn_sl2 =  reshape(cell2mat(k_fhn.sl2),n_f,n_h-1,n_f);
        k_fhn_ls1 =  reshape(cell2mat(k_fhn.ls1),n_f,n_h-1,n_f);
        k_fhn_ls2 =  reshape(cell2mat(k_fhn.ls2),n_f,n_h-1,n_f);
        n = n_f;
    case 'raw'
        k_fhn_sl1 = reach(k_fhn,'sl1');
        k_fhn_sl2 = reach(k_fhn,'sl2');
        k_fhn_ls1 = reach(k_fhn,'ls1');
        k_fhn_ls2 = reach(k_fhn,'ls2');
        n = size(k_fhn_sl1,3);
end

%Get the mean strf
k_fh_mean.sl1 = nanmean(k_fhn_sl1,3);
k_fh_mean.sl2 = nanmean(k_fhn_sl2,3);
k_fh_mean.ls1 = nanmean(k_fhn_ls1,3);
k_fh_mean.ls2 = nanmean(k_fhn_ls2,3);
%Get the -ve and +ve part separately
k_fhn_sl1_neg = min(k_fhn_sl1,0);
k_fhn_sl1_pos = max(k_fhn_sl1,0);
k_fhn_sl2_neg = min(k_fhn_sl2,0);
k_fhn_sl2_pos = max(k_fhn_sl2,0);

k_fhn_ls1_neg = min(k_fhn_ls1,0);
k_fhn_ls1_pos = max(k_fhn_ls1,0);
k_fhn_ls2_neg = min(k_fhn_ls2,0);
k_fhn_ls2_pos = max(k_fhn_ls2,0);
%Compute the mean across frequency
k_hn_neg_mean_sl1 = squeeze(nanmean(k_fhn_sl1_neg,1));
k_hn_pos_mean_sl1 = squeeze(nanmean(k_fhn_sl1_pos,1));
k_hn_neg_mean_sl2 = squeeze(nanmean(k_fhn_sl2_neg,1));
k_hn_pos_mean_sl2 = squeeze(nanmean(k_fhn_sl2_pos,1));

k_hn_neg_mean_ls1 = squeeze(nanmean(k_fhn_ls1_neg,1));
k_hn_pos_mean_ls1 = squeeze(nanmean(k_fhn_ls1_pos,1));
k_hn_neg_mean_ls2 = squeeze(nanmean(k_fhn_ls2_neg,1));
k_hn_pos_mean_ls2 = squeeze(nanmean(k_fhn_ls2_pos,1));

max_val = max(abs([k_hn_neg_mean_sl1(:);k_hn_pos_mean_sl1(:);k_hn_neg_mean_sl2(:);k_hn_pos_mean_sl2(:);k_hn_neg_mean_ls1(:);k_hn_pos_mean_ls1(:);k_hn_neg_mean_ls2(:);k_hn_pos_mean_ls2(:)]));
k_hn_neg_mean_sl1 = k_hn_neg_mean_sl1/max_val;
k_hn_pos_mean_sl1 = k_hn_pos_mean_sl1/max_val;
k_hn_neg_mean_sl2 = k_hn_neg_mean_sl2/max_val;
k_hn_pos_mean_sl2 = k_hn_pos_mean_sl2/max_val;

k_hn_neg_mean_ls1 = k_hn_neg_mean_ls1/max_val;
k_hn_pos_mean_ls1 = k_hn_pos_mean_ls1/max_val;
k_hn_neg_mean_ls2 = k_hn_neg_mean_ls2/max_val;
k_hn_pos_mean_ls2 = k_hn_pos_mean_ls2/max_val;
%Compute the mean for the k_h part for +ve and -ve part
k_h_neg_mean_sl1 = nanmean(k_hn_neg_mean_sl1,2);
k_h_neg_mean_sl2 = nanmean(k_hn_neg_mean_sl2,2);
k_h_pos_mean_sl1 = nanmean(k_hn_pos_mean_sl1,2);
k_h_pos_mean_sl2 = nanmean(k_hn_pos_mean_sl2,2);

k_h_neg_mean_ls1 = nanmean(k_hn_neg_mean_ls1,2);
k_h_neg_mean_ls2 = nanmean(k_hn_neg_mean_ls2,2);
k_h_pos_mean_ls1 = nanmean(k_hn_pos_mean_ls1,2);
k_h_pos_mean_ls2 = nanmean(k_hn_pos_mean_ls2,2);

%Compute the SEM
k_h_neg_sem_sl1 = nanstd(k_hn_neg_mean_sl1,[],2)/sqrt(n);
k_h_neg_sem_sl2 = nanstd(k_hn_neg_mean_sl2,[],2)/sqrt(n);
k_h_pos_sem_sl1 = nanstd(k_hn_pos_mean_sl1,[],2)/sqrt(n);
k_h_pos_sem_sl2 = nanstd(k_hn_pos_mean_sl2,[],2)/sqrt(n);

k_h_neg_sem_ls1 = nanstd(k_hn_neg_mean_ls1,[],2)/sqrt(n);
k_h_neg_sem_ls2 = nanstd(k_hn_neg_mean_ls2,[],2)/sqrt(n);
k_h_pos_sem_ls1 = nanstd(k_hn_pos_mean_ls1,[],2)/sqrt(n);
k_h_pos_sem_ls2 = nanstd(k_hn_pos_mean_ls2,[],2)/sqrt(n);

max_val = max(abs([k_fh_mean.sl1(:);k_fh_mean.sl2(:);k_fh_mean.ls1(:);k_fh_mean.ls2(:)]));

%% STRF - k_fh

for k = 1:4
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    room = rooms{k};
    k_fh = [k_fh_mean.(room),zeros(n_f,1)];
    k_fh = k_fh./max_val;
    
    switch params.name
        case 'normative model'
            imagesc(k_fh,[-0.7 0.7]);
        case 'neurons'
            imagesc(k_fh,[-0.7 0.7]);
    end
    colormap('redblue');
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([0:skip_h:h_lim]);
    xticklabels(h_labels);
    xlim([0 h_lim]);
    set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    save_name = fullfile(save_dir,strjoin({'Mean','STRF',room,'room',params.name,params.type,'.svg'},'_'));
    saveas(gcf,save_name);
    close;
end

%% Temporal part of kernel - k_h

%Excitation
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
shadedErrorBar(h_plot,[0;k_h_pos_mean_sl1;0],[0;k_h_pos_sem_sl1;0],{'LineWidth',lw,'Color',exc_big_color}); hold on;
shadedErrorBar(h_plot,[0;k_h_pos_mean_sl2;0],[0;k_h_pos_sem_sl2;0],{'LineWidth',lw,'Color',exc_big_color});
shadedErrorBar(h_plot,[0;k_h_pos_mean_ls1;0],[0;k_h_pos_sem_ls1;0],{'LineWidth',lw,'Color',exc_small_color});
shadedErrorBar(h_plot,[0;k_h_pos_mean_ls2;0],[0;k_h_pos_sem_ls2;0],{'LineWidth',lw,'Color',exc_small_color});
xticks([0:skip_h_ms:h_max_ms]);
set(gca,'Xticklabel',[]);
xlim([0 h_lim_ms]);
switch params.name
    case 'normative model'
        yticks([0:0.5:1]);
    case 'neurons'
        yticks([0:0.1:0.2]);
end
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold off;
save_name = fullfile(save_dir,strjoin({'Temporal','excitatory','profile',params.name,params.type,'kh.svg'},'_'));
saveas(gcf,save_name);
close;

%Inhibition
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
shadedErrorBar(h_plot,[0;k_h_neg_mean_sl1;0],[0;k_h_neg_sem_sl1;0],{'LineWidth',lw,'Color','c'}); hold on;
shadedErrorBar(h_plot,[0;k_h_neg_mean_sl2;0],[0;k_h_neg_sem_sl2;0],{'LineWidth',lw,'Color',inh_big_color});
shadedErrorBar(h_plot,[0;k_h_neg_mean_ls1;0],[0;k_h_neg_sem_ls1;0],{'LineWidth',lw,'Color','m'});
shadedErrorBar(h_plot,[0;k_h_neg_mean_ls2;0],[0;k_h_neg_sem_ls2;0],{'LineWidth',lw,'Color',inh_small_color});
legend('L1','L2','S1','S2');
xticks([0:skip_h_ms:h_max_ms]);
xticklabels(h_labels);
xlim([0 h_lim_ms]);
ylim([-0.125 0]);
switch params.name
    case 'normative model'
        yticks([-0.3:0.15:0]);
    case 'neurons'
        yticks([-0.15:0.05:0]);
end
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold off;
save_name = fullfile(save_dir,strjoin({'Temporal','inhibitory','profile',params.name,params.type,'kh.svg'},'_'));
saveas(gcf,save_name);
close;


% annotation('textbox',[0.54 0.62 0.1 0.1],'String', sprintf('Large room inhibition'),'LineStyle','none','Color',inh_big_color,'FontSize',legend_font_sz,'FontWeight','bold');
% annotation('textbox',[0.54 0.8 0.1 0.1],'String', sprintf('Small room excitation'),'LineStyle','none','Color',exc_small_color,'FontSize',legend_font_sz,'FontWeight','bold');
% annotation('textbox',[0.54 0.74 0.1 0.1],'String', sprintf('Large room excitation'),'LineStyle','none','Color',exc_big_color,'FontSize',legend_font_sz,'FontWeight','bold');
% annotation('textbox',[0.54 0.68 0.1 0.1],'String', sprintf('Small room inhibition'),'LineStyle','none','Color',inh_small_color,'FontSize',legend_font_sz,'FontWeight','bold');