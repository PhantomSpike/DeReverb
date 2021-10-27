function plot_ind_bf(k_fh_bf,params)
%This function will plot all the f-band kernels for a given room size
start_f = 700;
stop_f = 17000;
%% Plotting params
freqs = params.freqs;
freqs_plot = freqs;
%Flip the data for plotting
freqs = fliplr(freqs);
k_fh_bf.small = fliplr(k_fh_bf.small);
k_fh_bf.big = fliplr(k_fh_bf.big);

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
lw = 4;
axis_sz = 60;
title_sz = 60;
h_lim_ms = h_max_ms;
h_lim = round(h_lim_ms/dt);
skip_f = 8;
skip_h = 10;
skip_h_ms = 100;

%Arrow params
arrow_sz = 10;
tip_angle = 30;
base_angle = 40;
tip_length = 30;
epsilon = 0.01;
per_arrow = 0.4;

for f = 1:n_f
    f_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    y_labels{count} = num2str(freqs_plot(f)./1000,'%.1f');
end

count = 0;
for h = 1:skip_h:n_h+1
    count = count+1;
    h_labels{count} = num2str(h_steps(h),'%.0f');
end

row = 2;
col = 3;
per = 0.03;
edgel = 0.095; edger = 0.01; edgeh = 0.096; edgeb = 0.098; space_h = 0.02; space_v = 0.065;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

%Select frequencies for plotting
ix_start = find(freqs>start_f); ix_start = ix_start(1);
ix_stop = find(freqs<stop_f); ix_stop = ix_stop(end);
ix_freq = [ix_start:ix_stop];
n_f_select = length(ix_freq);


%% STRF - k_fh
count_big = 0;
count_small = 0;
n_batch = 8;

for b = 1:n_batch
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    for k = 1:row*col
        subplot('position',pos{k});
        if ismember(k,[1:col])
            count_big = count_big+1;
            ix = ix_freq(count_big);
            k_fh = k_fh_bf.big{ix};
        else
            count_small = count_small+1;
            ix = ix_freq(count_small);
            k_fh = k_fh_bf.small{ix};
        end
%         k_fh_small_temp = k_fh_bf.big{ix}; %Before this meant that the
%         magnitude of the two STRFs (small and large) are shown relative
%         to each other. Now they are normalized only within themselves
%         k_fh_large_temp = k_fh_bf.small{ix};
%         max_val = max(abs([k_fh_small_temp(:);k_fh_large_temp(:)]));

        max_val = max(abs(k_fh(:))); %
        k_fh = k_fh./max_val;
        k_fh = [k_fh,zeros(n_f,1)];
        k_fh_neg_temp = min(k_fh,0);
        k_fh_pos_temp = max(k_fh,0);
        
        if ismember(k,[1:col])
            k_h_neg.big{count_big} = mean(k_fh_neg_temp);
            k_h_pos.big{count_big} = mean(k_fh_pos_temp);
        else
            k_h_neg.small{count_small} = mean(k_fh_neg_temp);
            k_h_pos.small{count_small} = mean(k_fh_pos_temp);
        end
        
        switch params.name
            case 'normative model'
                imagesc(k_fh,[-0.8 0.8]);
            case 'neurons'
                imagesc(k_fh,[-1.2 1.2]);
        end
        colormap('redblue');
        
        if ismember(k,[1:col])
            title([f_labels{ix},' kHz'],'fontsize',title_sz,'FontWeight','Normal');
        end
        
        yticks([1:skip_f:n_f]);
        xticks([0:skip_h:h_lim]);
        
        set(gca,'Xticklabel',[]);
        set(gca,'Yticklabel',[]);
        if ismember(b,[1:2:n_batch]) && k==(row-1)*col + 1
            xticklabels(h_labels);
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',axis_sz,'FontWeight','Normal')
            yticklabels(y_labels);
            c = get(gca,'YTickLabel');
            set(gca,'YTickLabel',c,'fontsize',axis_sz,'FontWeight','Normal')
        end
        
        xlim([0 h_lim]);
        set(gcf,'color','w');
        set(gca,'TickDir','out');
    end
    save_name = fullfile(save_dir,strjoin({'STRFs','batch',num2str(b),params.name,'kfh.svg'},'_'));
    saveas(gcf,save_name);
    close;
end

%% Temporal part of kernel - k_h

switch params.name
    case 'normative model'
        y_down_lim = -0.5;
        y_tick_inh = [y_down_lim:abs(y_down_lim/2):0];
        y_tick_exc = [0:0.5:1];
    case 'neurons'
        y_down_lim = -0.8;
        y_tick_inh = [y_down_lim:abs(y_down_lim/2):0];
        y_tick_exc = [0:0.5:1];
end
count_exc = 0;
count_inh = 0;

for b = 1:n_batch
    figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
    for k = 1:row*col
        subplot('position',pos{k});
        if ismember(k,[1:col])
            count_exc = count_exc+1;
            ix = ix_freq(count_exc);
            k_h_big = k_h_pos.big{count_exc};
            k_h_small = k_h_pos.small{count_exc};
            big_color = exc_big_color;
            small_color = exc_small_color;
            max_val_small = max(k_h_pos.small{count_exc}(:));
            max_val_big = max(k_h_pos.big{count_exc}(:));
            yticks(y_tick_exc);
        else
            count_inh = count_inh+1;
            ix = ix_freq(count_inh);
            k_h_big = k_h_neg.big{count_inh};
            k_h_small = k_h_neg.small{count_inh};
            big_color = inh_big_color;
            small_color = inh_small_color;
            max_val_small = max(k_h_pos.small{count_inh}(:));
            max_val_big = max(k_h_pos.big{count_inh}(:));
            yticks(y_tick_inh);
        end
        
        %Compute COM
        k_h_big_com = k_h_big./sum(k_h_big(:)); %Scale the values to sum to 1 for inhibition
        k_h_small_com = k_h_small./sum(k_h_small(:)); %Scale the values to sum to 1 for excitation
        com_big = k_h_big_com(:)'*h_steps(1:end-1)';
        com_small = k_h_small_com(:)'*h_steps(1:end-1)';
        
        %Plot the temporal profiles
        hold on;
%         k_h_small_plot = k_h_small./max(abs(k_h_small(:)));
%         k_h_big_plot = k_h_big./max(abs(k_h_big(:)));
        k_h_small_plot = k_h_small./max_val_small;
        k_h_big_plot = k_h_big./max_val_big;
        plot(h_plot,[0 k_h_big_plot],'Color',big_color,'LineWidth',lw);
        plot(h_plot,[0 k_h_small_plot],'Color',small_color,'LineWidth',lw);
        
        xticks([0:skip_h_ms:h_lim_ms]);
        
        set(gca,'Xticklabel',[]);
        set(gca,'Yticklabel',[]);
        
        if ismember(b,[1:2:n_batch]) && k==(row-2)*col+1
            set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto','fontsize',axis_sz,'FontWeight','Normal')
        end
        
        if ismember(b,[1:2:n_batch]) && k==(row-1)*col+1
            xticklabels(h_labels);
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',axis_sz,'FontWeight','Normal')
            set(gca,'YTickLabelMode', 'auto','fontsize',axis_sz,'FontWeight','Normal')
        end

        
        set(gcf,'color','w');
        set(gca,'TickDir','out');
        xlim([0 h_lim_ms]);
        
        if ismember(k,[1:col])
            y_up_lim = 1;
            ylim([0 y_up_lim]);
            arrow([com_big per_arrow*y_up_lim],[com_big epsilon],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_big_color);
            arrow([com_small per_arrow*y_up_lim],[com_small epsilon],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_small_color);
        else
            ylim([y_down_lim 0]);
            arrow([com_big y_down_lim + per_arrow*abs(y_down_lim)],[com_big y_down_lim + epsilon],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_big_color);
            arrow([com_small y_down_lim + per_arrow*abs(y_down_lim)],[com_small y_down_lim + epsilon],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_small_color);
        end
        
        if ismember(k,[1:col])
            title([f_labels{ix},' kHz'],'fontsize',title_sz,'FontWeight','Normal');
        end
        
        hold off;

        
    end
    
    save_name = fullfile(save_dir,strjoin({'Temporal','profile','batch',num2str(b),'all_freq',params.name,'kh.svg'},'_'));
    saveas(gcf,save_name);
    close;
end