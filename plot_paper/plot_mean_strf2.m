function plot_mean_strf2(k_fhn,params)
%% Plotting params
freqs = params.freqs;
n_f = length(freqs);
n_h = params.n_h+1;
dt = params.dt_ms;
h_max_ms = params.h_max_ms+dt;
save_dir = params.save_dir;
h_steps = [0:dt:h_max_ms];
h_com = [10:dt:h_max_ms-dt];
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

plot_strf = 0;

%Arrow stuff
arrow_sz = 20; 
tip_angle = 30;
base_angle = 40;
% tip_length = 30;
tip_length = 42;

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

rooms{1} = 'large';
rooms{2} = 'small';

switch params.type
    case 'bf'
        k_fhn_small =  reshape(cell2mat(k_fhn.small),n_f,n_h-1,n_f);
        k_fhn_large =  reshape(cell2mat(k_fhn.big),n_f,n_h-1,n_f);
        n = n_f;
    case 'raw'
        k_fhn_small = reach(k_fhn,'small');
        k_fhn_large = reach(k_fhn,'big');
        n = size(k_fhn_small,3);
end

%Get the mean strf
k_fh_mean.small = nanmean(k_fhn_small,3);
k_fh_mean.large = nanmean(k_fhn_large,3);

%Get the -ve and +ve part separately
k_fhn_small_neg = min(k_fhn_small,0);
k_fhn_small_pos = max(k_fhn_small,0);
k_fhn_large_neg = min(k_fhn_large,0);
k_fhn_large_pos = max(k_fhn_large,0);

%Compute the mean across frequency 
k_hn_neg_mean_small = squeeze(nanmean(k_fhn_small_neg,1));
k_hn_pos_mean_small = squeeze(nanmean(k_fhn_small_pos,1));
k_hn_neg_mean_large = squeeze(nanmean(k_fhn_large_neg,1));
k_hn_pos_mean_large = squeeze(nanmean(k_fhn_large_pos,1));

%Normalise all the values to scale to -1/+1 for the inh/exc
max_val_neg_small = max(abs(k_hn_neg_mean_small));
max_val_pos_small = max(abs(k_hn_pos_mean_small));
max_val_neg_large = max(abs(k_hn_neg_mean_large));
max_val_pos_large = max(abs(k_hn_pos_mean_large));

k_hn_neg_mean_small = k_hn_neg_mean_small./max_val_neg_small;
k_hn_pos_mean_small = k_hn_pos_mean_small./max_val_pos_small;
k_hn_neg_mean_large = k_hn_neg_mean_large./max_val_neg_large;
k_hn_pos_mean_large = k_hn_pos_mean_large./max_val_pos_large;

%Compute the mean for the k_h part for +ve and -ve part
k_h_neg_mean_small = nanmean(k_hn_neg_mean_small,2);
k_h_neg_mean_large = nanmean(k_hn_neg_mean_large,2);
k_h_pos_mean_small = nanmean(k_hn_pos_mean_small,2);
k_h_pos_mean_large = nanmean(k_hn_pos_mean_large,2);

%Compute the SEM
k_h_neg_sem_small = nanstd(k_hn_neg_mean_small,[],2)/sqrt(n);
k_h_neg_sem_large = nanstd(k_hn_neg_mean_large,[],2)/sqrt(n);
k_h_pos_sem_small = nanstd(k_hn_pos_mean_small,[],2)/sqrt(n);
k_h_pos_sem_large = nanstd(k_hn_pos_mean_large,[],2)/sqrt(n);
max_val = max(abs([k_fh_mean.small(:);k_fh_mean.large(:)]));

%% Compute the COM

%Inhibition
k_h_neg_mean_norm_small = k_h_neg_mean_small./sum(k_h_neg_mean_small(:)); %Scale the values to sum to 1 small room
k_h_neg_mean_norm_large = k_h_neg_mean_large./sum(k_h_neg_mean_large(:)); %large room

%Compute COM using a weighted sum of all values
com_neg_small = (h_com*k_h_neg_mean_norm_small); 
com_neg_large = (h_com*k_h_neg_mean_norm_large); 

%Excitation
k_h_pos_mean_norm_small = k_h_pos_mean_small./sum(k_h_pos_mean_small(:)); %Scale the values to sum to 1 small room
k_h_pos_mean_norm_large = k_h_pos_mean_large./sum(k_h_pos_mean_large(:)); %large room

%Compute COM using a weighted sum of all values
com_pos_small = (h_com*k_h_pos_mean_norm_small); 
com_pos_large = (h_com*k_h_pos_mean_norm_large);

%% STRF - k_fh

if plot_strf
    
    for k = 1:2
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
    
end
%% Temporal part of kernel - k_h

%Excitation
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
shadedErrorBar(h_plot,[0;k_h_pos_mean_small;0],[0;k_h_pos_sem_small;0],{'LineWidth',lw,'Color',exc_small_color}); hold on;
shadedErrorBar(h_plot,[0;k_h_pos_mean_large;0],[0;k_h_pos_sem_large;0],{'LineWidth',lw,'Color',exc_big_color});

xticks([0:skip_h_ms:h_max_ms]);
set(gca,'Xticklabel',[]);
xlim([0 h_lim_ms]);
switch params.name
    case 'normative model'
        yticks([0:0.5:1]);
        min_pos = 0.001;
        max_pos = 0.333;
    case 'neurons'
        ylim([0 0.8]);
        yticks([0:0.4:0.8]);
        min_pos = 0.0001;
        max_pos = 0.2683;
end

%Plot the COM as arrows
arrow([com_pos_large max_pos],[com_pos_large min_pos],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_big_color);
arrow([com_pos_small max_pos],[com_pos_small min_pos],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_small_color);

set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',axis_sz,'FontWeight','Normal');
hold off;
save_name = fullfile(save_dir,strjoin({'Temporal','excitatory','profile',params.name,params.type,'kh.svg'},'_'));
saveas(gcf,save_name);
close;

%Inhibition
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
shadedErrorBar(h_plot,[0;k_h_neg_mean_small;0],[0;k_h_neg_sem_small;0],{'LineWidth',lw,'Color',inh_small_color}); hold on;
shadedErrorBar(h_plot,[0;k_h_neg_mean_large;0],[0;k_h_neg_sem_large;0],{'LineWidth',lw,'Color',inh_big_color});
    
xticks([0:skip_h_ms:h_max_ms]);
xticklabels(h_labels);
xlim([0 h_lim_ms]);

switch params.name
    case 'normative model'
        yticks([-1:0.5:0]);
        min_neg = -0.999;
        max_neg = -0.666;
    case 'neurons'
        yticks([-0.6:0.3:0]);
        ylim([-0.6 0]);
        min_neg = -0.599;
        max_neg = -0.399;
end

%Plot the COM as arrows
arrow([com_neg_large max_neg],[com_neg_large min_neg],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_big_color);
arrow([com_neg_small max_neg],[com_neg_small min_neg],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_small_color);

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

