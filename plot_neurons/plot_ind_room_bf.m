function plot_ind_room_bf(k_fh_bf,params)
%This function will plot all the f-band kernels for a given room size
start_f = 10000;
stop_f = 700;
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
lw = 3;
axis_sz = 17;
legend_font_sz = 32;
title_sz = 25;
h_lim_ms = h_max_ms;
h_lim = round(h_lim_ms/dt);
skip_f = 5;
skip_h = 5;
skip_h_ms = 50;

rooms{1} = 'big';
rooms{2} = 'small';


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

row = 2;
col = 10;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = 0.05; space_h = per; space_v = 0.08;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

%Select frequencies for plotting
ix_start = find(freqs<start_f); ix_start = ix_start(1);
ix_stop = find(freqs>stop_f); ix_stop = ix_stop(end);
ix_freq = [ix_start:ix_stop];
n_fplot = length(ix_freq);
%% STRF - k_fh
for r = 1:2
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    room = rooms{r};
    for k = 1:n_fplot
        subplot('position',pos{k});
        ix = ix_freq(k);
        k_fh_small_temp = k_fh_bf.big{ix};
        k_fh_large_temp = k_fh_bf.small{ix};
        max_val = max(abs([k_fh_small_temp(:);k_fh_large_temp(:)]));
        k_fh = k_fh_bf.(room){ix};
        k_fh = k_fh./max_val;
        k_fh = [k_fh,zeros(n_f,1)];
        k_fh_neg_temp = min(k_fh,0);
        k_fh_pos_temp = max(k_fh,0);
        k_h_neg.(room){k} = mean(k_fh_neg_temp);
        k_h_pos.(room){k} = mean(k_fh_pos_temp);
        switch params.name
            case 'normative model'
                imagesc(k_fh,[-0.7 0.7]);
                title([f_labels{ix},' kHz']);
            case 'neurons'
                imagesc(k_fh,[-1 1]);
%                 title([f_labels{ix},' kHz',' n=',num2str(params.n_bf(ix))]);
                title([f_labels{ix},' kHz']);
        end
        colormap('redblue');
       
        yticks([1:skip_f:n_f]);
        yticklabels(y_labels);
        xticks([0:skip_h:h_lim]);
        xticklabels(h_labels);
        xlim([0 h_lim]);
        set(gca,'FontSize',axis_sz,'FontWeight','Bold');
        set(gcf,'color','w');
    end
%     sgtitle([params.name,' ',room,' room'],'FontSize',title_sz,'Color','k','FontWeight','Bold');
    save_name = fullfile(save_dir,strjoin({'STRFs',params.name,'kfh.svg'},'_'));
    saveas(gcf,save_name);
%     print('-r300',gcf,save_name,'-dpng');
    close;
end

%% Temporal part of kernel - k_h
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type); 

for k = 1:n_fplot
    subplot('position',pos{k});
    ix = ix_freq(k);
    plot(h_plot,[0 k_h_neg.small{k}],'Color',inh_small_color,'LineWidth',lw); hold on;
    plot(h_plot,[0 k_h_neg.big{k}],'Color',inh_big_color,'LineWidth',lw);
    plot(h_plot,[0 k_h_pos.small{k}],'Color',exc_small_color,'LineWidth',lw);
    plot(h_plot,[0 k_h_pos.big{k}],'Color',exc_big_color,'LineWidth',lw);
    yline(0,'--','LineWidth',lw,'Color','k');
    xticks([0:skip_h_ms:h_max_ms]);
    xticklabels(h_labels);
    switch params.name
        case 'normative model'
            title([f_labels{ix},' kHz']);
        case 'neurons'
%             title([f_labels{ix},' kHz',' n=',num2str(params.n_bf(ix))]);
            title([f_labels{ix},' kHz']);
    end
    set(gcf,'color','w');
    set(gca,'FontSize',axis_sz,'FontWeight','Bold');
    hold off;
    xlim([0 h_lim_ms]);
end
% 
% sgtitle(['Temporal kernel profiles ',params.name],'FontSize',title_sz,'Color','k','FontWeight','Bold');
save_name = fullfile(save_dir,strjoin({'Temporal','profile','all_freq',params.name,'kh.svg'},'_'));
saveas(gcf,save_name);
% print('-r300',gcf,save_name,'-dpng');
close;
