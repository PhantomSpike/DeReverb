%% Input params
ir_small_path = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Full_IRs_calc/ferret_small_full.mat';
ir_large_path = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Full_IRs_calc/ferret_large_full.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_4';
ear = 'right';
font_type = 'Liberation Sans';
freq_up_bound = 17000;
freq_down_bound = 700;
all_font_sz = 75;
legend_font_sz = 50;
skip_f = 3;
plot_coch = 0;
plot_decay_profiles = 0;
small_color = [0.592, 0.737, 0.384];
big_color = [0.173, 0.373, 0.176];

%% Load the IRs
data_small = load(ir_small_path);
data_large = load(ir_large_path);
fs = data_small.Fs;
switch ear
    case 'left'
        ir.small = data_small.data(:,1); ir.large = data_large.data(:,1);
    case 'right'
        ir.small = data_small.data(:,2); ir.large = data_large.data(:,2);
end

dt_s = 1/fs;
l_small_s = length(ir.small)/fs;
l_large_s = length(ir.large)/fs;
t_small = [dt_s:dt_s:l_small_s];
t_large = [dt_s:dt_s:l_large_s];

%% Make cochleagrams
dt_ms = 10; %Time bin in ms
coch_type = 'specpower'; %The type of cochleagram
freq_spacing = 'log'; %The spacing of the frequencies
f_min = 400; %The minimal frequency
f_max = 19000; %The maximal frequency
n_f = 30; %The total number of freqeuncies
color_map = 'inferno'; %Options are 'inferno', 'magma', 'plasma', 'viridis'

coch(1).reverb_cond = 'small';
coch(2).reverb_cond = 'large';
coch(1).type = coch_type;
n_cond = length(coch);

for s = 1:n_cond
    room = coch(s).reverb_cond;
    switch coch_type
        case 'specpower'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_power(ir.(room), fs, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'speclog'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_log(ir.(room), fs, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'spechill'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_Hill(ir.(room), fs, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'bencoch'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram(ir.(room), fs, dt_ms, freq_spacing, f_min, f_max, n_f);
    end
    coch(s).X_ft = flipud(coch(s).X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
end

freqs = fliplr(coch(1).params.freqs);

%% Plot the cochleagrams
if plot_coch
    num_freq = length(freqs);
    skip_f = 2;
    freqs_plot = freqs(1:skip_f:end);
    n_tlab = 5;
    
    for f = 1:length(freqs_plot)
        f_plot_labels{f} = num2str(freqs_plot(f)./1000,'%.1f');
    end
    
    for s = 1:n_cond
        
        skip_t = round(length(coch(s).t)/n_tlab);
        for tm = 1:n_tlab
            x_labels{tm} = num2str(coch(s).t((tm-1)*skip_t +1),'%.2f');
        end
        
        max_val = max(coch(s).X_ft(:));
        coch(s).X_ft = coch(s).X_ft - max_val;
        figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
        imagesc(coch(s).X_ft);
        colorbar;
        colormap(color_map);
        caxis([-80 0]);
        yticks([1:skip_f:length(freqs)]);
        yticklabels(f_plot_labels);
        xticks([1:skip_t:length(coch(s).t)]);
        xticklabels(x_labels);
        title(coch(s).reverb_cond);
        set(gca,'FontName','Arial','FontSize',15,'FontWeight','Normal');
        xlabel('Time [s]','FontSize',16,'FontWeight','Normal');
        ylabel('Freqeuncy [kHz]','FontSize',16,'FontWeight','Normal');
        set(gcf,'color','w');
        set(gca,'TickDir','out');
        save_name = fullfile(save_dir,['Cochleagram_',coch(s).reverb_cond,'_',coch_type,'.svg']);
        saveas(gcf,save_name);
        close;
    end
end

%% Compute the decay rates from the differen cochlear subbands
end_s.small = 0.3;
end_s.large = 1.5;
for s = 1:n_cond
    room = coch(s).reverb_cond;
    t_s = coch(s).t;
    ix = t_s<end_s.(room);
    t_fit_s = t_s(ix);
    
    for f = 1:n_f
        x_f = coch(s).X_ft(f,ix); %Get the necessary frequency channel
        X = [ones(length(t_fit_s),1),t_fit_s(:)]; %Make design matrix
        y = x_f(:);
        
        [b.(room)(f,:),~,~,~,stats.(room)(f,:)] = regress(y,X);
        r.(room)(f) = stats.(room)(f,1); p.(room) = stats.(room)(f,3);
    end
    avg_power.(room) = mean(coch(s).X_ft(:,ix),2);
end

rt_60.small = -60./b.small(:,2);
rt_60.large = -60./b.large(:,2);
rt_60.freqs = freqs;

%% Plot the decay profiles of each freqeuncy band
if plot_decay_profiles
    %Params for subplots title textbox
    tpos_x = 0.4; tpos_y = 0.9; tsz_x = 0.3; tsz_y = 0.3;
    %Params for RT60, R^2 textbox
    Apos_x = 0.7; Apos_y = 0.7; Asz_x = 0.28; Asz_y = 0.35;
    
    axis_sz = 10;
    
    lw = 3;
    row = 6;
    col = 5;
    per = 0.03;
    edgel = 0.05; edger = 0.02; edgeh = per; edgeb = 0.06; space_h = per; space_v = 0.05;
    [pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
    xdiff = (pos{2}(1) - pos{1}(1)) - space_h;
    ydiff = (pos{1}(2) - pos{col+1}(2)) - space_v;
    
    for f = 1:n_f
        f_labels{f} = num2str(freqs(f)./1000,'%.1f');
    end
    
    
    
    for s = 1:size(coch,2)
        figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
        room = coch(s).reverb_cond;
        
        for f = 1:n_f
            subplot('position',pos{f});
            plot(coch(s).t,coch(s).X_ft(f,:),'LineWidth',lw);
            hold on;
            h1 = refline([b.(room)(f,2),b.(room)(f,1)]);
            h1.Color = 'r';
            h1.LineWidth = 3;
            h1.LineStyle = '--';
            hold off;
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Normal');
            
            if ismember(f,[1:col:row*col])
                ylabel('Power [dB]','FontSize',16,'FontWeight','Normal');
            end
            
            if ismember(f,[(row-1)*col+1:row*col])
                xlabel('Time [ms]','FontSize',16,'FontWeight','Normal');
            end
            
            %Title
            annotation(gcf,'textbox',pos{f} + [xdiff*tpos_x, ydiff*tpos_y -pos{f}(3)*(1-tsz_x) -pos{f}(4)*(1-tsz_y)],...
                'String', sprintf([f_labels{f},'kHz']),'LineStyle','none','FontSize',axis_sz,'FontWeight','Normal');
            %RT60, R^2 textbox
            annotation(gcf,'textbox',pos{f} + [xdiff*Apos_x, ydiff*Apos_y -pos{f}(3)*(1-Asz_x) -pos{f}(4)*(1-Asz_y)],...
                'String', sprintf('RT_{60}=%.1fs\nR^{2}=%.1f',rt_60.(room)(f),r.(room)(f)),'LineStyle','none','FontSize',axis_sz,'FontWeight','Normal','Color','r');
            
            set(gcf,'color','w');
            set(gca,'TickDir','out');
            ylim([-90 0])
        end
        
        save_name = fullfile(save_dir,['Cochleagram_decay_rate_',coch(s).reverb_cond,'_',coch_type,'.svg']);
        saveas(gcf,save_name);
        close;
    end
end

%% Plot RT60 vs freqeuncy
lw = 6;
ix = rt_60.freqs>freq_down_bound & rt_60.freqs<freq_up_bound;
plot_freqs = fliplr(rt_60.freqs(ix));
n_f = length(plot_freqs);

for f = 1:n_f
    f_labels{f} = num2str(plot_freqs(f)./1000,'%.1f');
end

plot_rt60_small = flipud(rt_60.small(ix));
plot_rt60_large = flipud(rt_60.large(ix));

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
semilogx(plot_freqs,plot_rt60_small,'Color',small_color,'LineWidth',lw);hold on;
semilogx(plot_freqs,plot_rt60_large,'Color',big_color,'LineWidth',lw);
hold off;
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
yticks([0:1:4]);
set(gca,'XMinorTick','off');
% xlabel('Frequency [kHz]','FontSize',all_font_sz,'FontWeight','Normal');
% ylabel('RT_{60} [s]','FontSize',all_font_sz,'FontWeight','Normal');
% title(['RT_{60} vs Frequency for individual rooms']);
annotation('textbox',[0.72 0.8 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.72 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,['RT60 vs BF individual rooms.svg']);
saveas(gcf,save_name);
close;

%% Fit regression for the two Rooms
sz = 30;
plot_freqs = plot_freqs(:);
n_reg = length(plot_freqs);
X_small = [ones(n_reg,1),log10(plot_freqs)];
X_big = [ones(n_reg,1),log10(plot_freqs)];

y_small_neg = plot_rt60_small; y_small_neg = y_small_neg(:);
y_big_neg = plot_rt60_large; y_big_neg = y_big_neg(:);

[b_small,~,~,~,stats_small] = regress(y_small_neg,X_small);
[b_big,~,~,~,stats_big] = regress(y_big_neg,X_big);
r_small = sign(b_small(2))*sqrt(stats_small(1)); p_small = stats_small(3);
r_big = sign(b_big(2))*sqrt(stats_big(1)); p_big = stats_big(3);

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
plot(log10(plot_freqs),plot_rt60_small,'-','Color',small_color,'LineWidth',lw);hold on;
plot(log10(plot_freqs),plot_rt60_large,'-','Color',big_color,'LineWidth',lw);
set(gca, 'XScale', 'log');
h1 = refline([b_small(2),b_small(1)]); h1.Color = small_color; h1.LineWidth = lw; h1.LineStyle = '--';
h2 = refline([b_big(2),b_big(1)]); h2.Color = big_color; h2.LineWidth = lw; h2.LineStyle = '--';
hold off;
xticks([log10(plot_freqs(2:skip_f:n_f))]);
xticklabels(f_labels(2:skip_f:n_f));
ylim([0 4]);
yticks([0:1:4]);
set(gca,'XMinorTick','off');
annotation('textbox',[0.6 0.82 0.1 0.1],'String', sprintf('r_{large}= %.2f^{%s}',r_big,p_star(p_small)),'LineStyle','none','Color',big_color,'FontSize',legend_font_sz,'FontWeight','Normal');
annotation('textbox',[0.6 0.72 0.1 0.1],'String', sprintf('r_{small}= %.2f^{%s}',r_small,p_star(p_big)),'LineStyle','none','Color',small_color,'FontSize',legend_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,['Regression RT60 vs BF individual rooms.svg']);
saveas(gcf,save_name);
close;

%% Plot Avg Power vs freqeuncy
plot_avg_power_small = flipud(avg_power.small(ix));
plot_avg_power_large = flipud(avg_power.large(ix));
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
semilogx(plot_freqs,plot_avg_power_small,'Color',small_color,'LineWidth',lw);hold on;
semilogx(plot_freqs,plot_avg_power_large,'Color',big_color,'LineWidth',lw);
xticks([plot_freqs(2:skip_f:n_f)]);
xticklabels(f_labels(2:skip_f:n_f));
set(gca,'XMinorTick','off');
% xlabel('Frequency [kHz]','FontSize',all_font_sz,'FontWeight','Normal');
ylabel('Average Power [dB]','FontSize',all_font_sz,'FontWeight','Normal');
title(['Average Power vs BF for individual rooms']);
annotation('textbox',[0.75 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',big_color,'FontSize',all_font_sz,'FontWeight','Normal');
annotation('textbox',[0.75 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',small_color,'FontSize',all_font_sz,'FontWeight','Normal');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,['Average Power vs BF individual rooms.svg']);
saveas(gcf,save_name);
close;


