%% Define and load vars
kernel_path = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms/kernel.mat'; %Absolute path to the kernel used for the prediction
coch_path = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/specpower/10ms/coch_all_conditions_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_1'; %Directory where to save the plots
fprintf('== Loading the data ==\n');tic;
load(kernel_path);
load(coch_path);
room_names{1} = 'anech_orig'; room_names{2} = 'small_orig'; room_names{3} = 'big_orig';
room_names{4} = 'anech'; room_names{5} = 'small';  room_names{6} = 'big';

%Get the cochleagrams
anech1_X_ft = coch(1).X_ft;
anech2_X_ft = coch(2).X_ft;
small1_X_ft = coch(3).X_ft;
small2_X_ft = coch(4).X_ft;
big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;

%Get the kernels
n_ker = length(kernels.small);
for k = 1:n_ker
    kernel_small{k} = kernels.small{k}.main{end};
    kernel_med{k} = kernels.med{k}.main{end};
    kernel_big{k} = kernels.big{k}.main{end};
end
n_h = kernels.n_h;
fprintf('== Done! This took %0.fs ==\n',toc);

%% Normalize the cochleagrams
mean_small = mean([small1_X_ft, small2_X_ft],2); std_small = std([small1_X_ft, small2_X_ft],[],2);
mean_big = mean([big1_X_ft, big2_X_ft],2); std_big = std([big1_X_ft, big2_X_ft],[],2);
        
small1_X_ft_norm = (small1_X_ft - mean_small)./std_small;
small2_X_ft_norm = (small2_X_ft - mean_small)./std_small;
big1_X_ft_norm = (big1_X_ft - mean_big)./std_big;
big2_X_ft_norm = (big2_X_ft - mean_big)./std_big;


%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small1_X_fht = tensorize(small1_X_ft_norm,n_h);
small2_X_fht = tensorize(small2_X_ft_norm,n_h);
big1_X_fht = tensorize(big1_X_ft_norm,n_h);
big2_X_fht = tensorize(big2_X_ft_norm,n_h);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Concatenate together
small_X_fht = cat(3,small1_X_fht,small2_X_fht);
big_X_fht = cat(3,big1_X_fht,big2_X_fht);
n_t = size(small_X_fht,3);

%% Make predictions
fprintf('== Making predicted cochleagrams ==\n');tic;
small_anech_X_ft_hat = zeros(n_ker,n_t);
big_anech_X_ft_hat = zeros(n_ker,n_t);

for k = 1:n_ker
    small_anech_X_ft_hat(k,:) = kernelconv(small_X_fht, kernel_small{k});
    big_anech_X_ft_hat(k,:) = kernelconv(big_X_fht, kernel_big{k});
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Undo normalization to convert back to original values
anech_X_ft_hat_final.anech = [anech1_X_ft, anech2_X_ft];
anech_X_ft_hat_final.small = (small_anech_X_ft_hat.*std_small) + mean_small;
anech_X_ft_hat_final.big = (big_anech_X_ft_hat.*std_big) + mean_big;

anech_X_ft_hat_final.anech_orig = [anech1_X_ft, anech2_X_ft];
anech_X_ft_hat_final.small_orig = [small1_X_ft, small2_X_ft];
anech_X_ft_hat_final.big_orig = [big1_X_ft, big2_X_ft];

%% Plot original cochleagrams and predictions
max_val = 40; %How much to subtract from the cochleagram for normalization
up_lim_dB = 0; %The upper limit in dB
low_lim_dB = -80; %The lower limit in dB
start_s = 0;
end_s = 4.01;
dt_s = 1;
skip_f = 5;

t = coch(1).t;
dt = mean(diff(t));
dt_samples = round(dt_s/dt);
t(end:2*length(t)) = [t(end):dt:t(end)*2];
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);
n_steps = round(length(t)/dt_samples)+1;
color_map = 'inferno';
freqs = fliplr(coch(1).params.freqs); 
n_f = length(freqs);
freqs = freqs(1:skip_f:end);
row = 2;
col = 3;
per = 0.03;
edgel = 0.05; edger = 0.02; edgeh = 0.05; edgeb = 0.07; space_h = 0.06; space_v = 0.09;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_steps
    x_labels{tm} = num2str(t((tm-1)*dt_samples +1),'%.0f');
end

figure('units','normalized','outerposition',[0 0 1 1]);

for r = 1:6
    room = room_names{r};
    subplot('position',pos{r});
    imagesc(anech_X_ft_hat_final.(room)(:,ix_start:ix_end)-max_val);
    colorbar;
    colormap(color_map);
    caxis([low_lim_dB up_lim_dB]);
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:dt_samples:length(t)]);
    xticklabels(x_labels);
    title(room);
    set(gca,'FontName','Arial','FontSize',15,'FontWeight','Normal');
    xlabel('Time [s]','FontSize',14,'FontWeight','Normal');
    ylabel('Freqeuncy [kHz]','FontSize',14,'FontWeight','Normal');
end
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Predicted_coch_all.png');
export_fig(save_name);
close all;

%% Plot individual cochleagrams
sz = 100;
all_font_sz = sz;
x_font_sz = sz;
y_font_sz = sz;
font_type = 'Liberation Sans';

%1 - Original anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.anech_orig(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%2 - Small reverb
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.small_orig(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Small_coch_example.svg');
saveas(gcf,save_name);
close all;

%3 - Big reverb
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.big_orig(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Big_coch_example.svg');
saveas(gcf,save_name);
close all;

%4 - Small reverb predicted anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.small(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Small_pred_anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%5 - Big reverb predicted anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.big(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Big_pred_anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%% 6 - Colorbar
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(anech_X_ft_hat_final.big(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
c = colorbar;
ht = text(455,10,'Power (dB)');
set(ht,'Rotation',-90,'FontSize',sz,'FontWeight','Normal','FontName',font_type);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Colorbar.svg');
saveas(gcf,save_name);
close all;