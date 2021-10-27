%% Define vars
kernel_path = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms/kernel.mat'; %Absolute path to the kernel used for the prediction
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_1';
bf = 3.8; %Define the normative kernel to be used
bf = bf*1e+3;
load(kernel_path);

%% Plotting params
h_lim_ms = 210;
freqs = fliplr(kernels.freqs);
n_f = length(freqs);
n_h = kernels.n_h;
dt = round(kernels.dt_ms);
h_max_ms = kernels.h_max_ms;
h_steps = [0:dt:h_max_ms];
h_lim = round(h_lim_ms/dt);
skip_f = 5;
skip_h = 5;

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    f_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for t = 1:skip_h:n_h+1
    count = count+1;
    h_labels{count} = num2str(h_steps(t),'%.0f');
end

[~,ix] = min(abs(bf-freqs)); %Find the kernel with the closest frequency

%% Plot the kernel
font_type = 'Liberation Sans';
sz = 100;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

k_fh = fliplr(kernels.small{ix}.main{end}.k_fh);
k_fh = k_fh./max(abs(k_fh(:)));
k_fh = [k_fh,zeros(n_f,1)];
imagesc(k_fh);
caxis([-1 1]);
colormap('redblue');
% colorbar('Ticks',[-1 0 1]);
% ht = text(22.5,9,'Weight [AU]');
% set(ht,'Rotation',-90,'FontSize',sz,'FontWeight','Normal','FontName',font_type);
yticks([1:skip_f:n_f]);
yticklabels(f_labels);
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
xticks([1:skip_h:h_max_ms+dt]);
xticklabels(h_labels);
% xlabel('History [ms]','FontSize',x_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Example_model_kernel.svg');
saveas(gcf,save_name);
close all;