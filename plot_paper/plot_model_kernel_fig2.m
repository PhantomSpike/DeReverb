clear all; close all;
%% Global params
model = 'ridge';
kernel_path = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms/kernel.mat';
bfs(1) = 16; %Define the normative kernels to be used
bfs(2) = 6.5;
bfs(3) = 2;
rooms{1} = 'big';
rooms{2} = 'small';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_2';
%% Detailed params
font_type = 'Liberation Sans';
sz = 70;
arrow_sz = 15;
tip_angle = 30;
base_angle = 40;
tip_length = 55;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

h_lim_ms = 210;
bfs = bfs*1e+3;
n_ker = length(bfs);
row = 2;
col = 3;
lw = 5;
per = 0.01;
edgel = 0.107; edger = 0.01; edgeh = 0.035; edgeb = 0.12; space_h = 0.0375; space_v = 0.077;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

exc_small_color = [0.9804    0.4196    0.6431];
exc_big_color = [0.8 0 0];
inh_small_color = [0.0745 0.6235 1.0000];
inh_big_color = 'b';

load(kernel_path,'kernels');
freqs = fliplr(kernels.freqs);
n_f = length(freqs);
n_h = kernels.n_h;
dt = round(kernels.dt_ms);
h_max_ms = kernels.h_max_ms;
h_steps = [0:dt:h_max_ms-dt];
h_plot = [0:dt:h_max_ms];
h_lim_ms = h_max_ms;
h_lim = round(h_lim_ms/dt);
skip_f = 8;
skip_h = 10;
skip_h_ms = 100;

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    f_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for t = 1:skip_h:n_h+1
    count = count+1;
    h_labels{count} = num2str(h_plot(t),'%.0f');
end
%% Plot the STRF model kernels
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
for k = 1:n_ker
    bf = bfs(k);
    [~,ix] = min(abs(bf-freqs));
    for r = 1:2
        room = rooms{r};
        switch model
            case {'sep','sep_kh'}
                k_fh = kernels.kernel.(room).k_f*fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = k_h.(room)/max(abs(k_h.(room)));
            case {'lasso','ridge'}
                k_fh = fliplr(kernels.(room){ix}.main{end}.k_fh);
        end
        subplot('position',pos{(r-1)*col + k});
        k_fh = k_fh./max(abs(k_fh(:)));
        imagesc(k_fh);
        caxis([-0.7 0.7]);
        colormap('redblue');
        set(gca,'Xticklabel',[],'Yticklabel',[]);
        yticks([1:skip_f:n_f]);
        xticks([0:skip_h:h_lim]);
        if k == 1 && r == 2
            yticklabels(f_labels);
            xticklabels(h_labels);
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',x_font_sz,'FontWeight','Normal')
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'fontsize',y_font_sz,'FontWeight','Normal')
        end
        set(gcf,'color','w');
        xlim([0 h_lim]); 
        set(gca,'TickDir','out');
    end
end

save_name = fullfile(save_dir,['Model_kernels_STRF_examples.svg']);
saveas(gcf,save_name);
close all;

%% Plot the temporal kernels
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

for k = 1:n_ker
    bf = bfs(k);
    [~,ix] = min(abs(bf-freqs));
    for r = 1:2
        room = rooms{r};
        switch model
            case {'sep','sep_kh'}
                k_fh = kernels.kernel.(room).k_f*fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = k_h.(room)/max(abs(k_h.(room)));
            case {'lasso','ridge'}
                k_fh = fliplr(kernels.(room){ix}.main{end}.k_fh);
                k_fh_neg = abs(min(k_fh,0));
                k_fh_pos = abs(max(k_fh,0));
                k_h_neg.(room) = mean(k_fh_neg);
                k_h_pos.(room) = mean(k_fh_pos);
                k_h_neg_f.(room) = k_h_neg.(room)./sum(k_h_neg.(room)(:)); %Scale the values to sum to 1 for inhibition
                k_h_pos_f.(room) = k_h_pos.(room)./sum(k_h_pos.(room)(:)); %Scale the values to sum to 1 for excitation
        end
        tau_neg.(room) = (k_h_neg_f.(room)*h_steps'); %Compute a weighted sum of all values
        tau_pos.(room) = (k_h_pos_f.(room)*h_steps'); %Compute a weighted sum of all values
    end
    %Excitation
    subplot('position',pos{k});
    set(gca,'Xticklabel',[],'Yticklabel',[]);
    if k == 1
        set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto','fontsize',y_font_sz,'FontWeight','Normal')
    end
    xticks([0:skip_h_ms:h_max_ms]);
    yticks([0:0.5:1]);
    hold on;
%     max_val = max([k_h_pos.big(:);k_h_pos.small(:);k_h_neg.big(:);k_h_neg.small(:)]);
    max_val_small = max([k_h_pos.small(:);abs(k_h_neg.small(:))]);
    max_val_big = max([k_h_pos.big(:);abs(k_h_neg.big(:))]);
    plot(h_plot,[0 k_h_pos.big./max_val_big],'Color',exc_big_color,'LineWidth',lw);
    plot(h_plot,[0 k_h_pos.small./max_val_small],'Color',exc_small_color,'LineWidth',lw);
    arrow([tau_pos.big 0.4],[tau_pos.big 0.01],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_big_color);
    arrow([tau_pos.small 0.4],[tau_pos.small 0.01],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',exc_small_color);
    hold off;
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    xlim([0 h_lim_ms]);

    
    %Inhibition
    subplot('position',pos{col+k});
    set(gca,'Xticklabel',[],'Yticklabel',[]);
    if k == 1
        xticks([0:skip_h_ms:h_max_ms]);
        xticklabels(h_labels);
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',x_font_sz,'FontWeight','Normal')
        set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto','fontsize',y_font_sz,'FontWeight','Normal')
    end
    xticks([0:skip_h_ms:h_max_ms]);
    yticks([-0.5:0.25:0]);
    hold on;
    plot(h_plot,[0 -k_h_neg.big./max_val_big],'Color',inh_big_color,'LineWidth',lw);
    plot(h_plot,[0 -k_h_neg.small./max_val_small],'Color',inh_small_color,'LineWidth',lw);
    arrow([tau_neg.big -0.3],[tau_neg.big -0.499],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_big_color);
    arrow([tau_neg.small -0.3],[tau_neg.small -0.499],'Width',arrow_sz,'Length',tip_length,'BaseAngle',base_angle,'TipAngle',tip_angle,'EdgeColor','none','FaceColor',inh_small_color);
    hold off;
    set(gcf,'color','w');
    set(gca,'TickDir','out');
    xlim([0 h_lim_ms]);
    ylim([-0.501 0]);

end

save_name = fullfile(save_dir,['Model_kernels_EI_examples.svg']);
saveas(gcf,save_name);
close all;