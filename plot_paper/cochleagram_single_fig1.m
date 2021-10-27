clear all; close all;
coch_file = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Ronnie_PLP/spechill/10ms/coch_all_conditions_spechill.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_1';
cond = 1;
%% Load the data
load(coch_file);
%% Params
start_s = 0;
end_s = 4;
color_map = 'inferno';
font_sz = 30;
%% Make legend
freqs = fliplr(coch(1).params.freqs);
t = coch(1).t;
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);
X_ft = coch(cond).X_ft(:,ix_start:ix_end);
n_f = length(freqs);
skip_f = 2;
freqs = freqs(1:skip_f:end);
n_tlab = 10;
skip_t = round(length(t)/n_tlab);
figure('units','normalized','outerposition',[0 0 1 1]);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(t((tm-1)*skip_t +1),'%.1f');
end
%% Plot the cochleagram
imagesc(X_ft);
colorbar;
colormap(color_map);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:skip_t:length(coch(1).t)]);
xticklabels(x_labels);
% title(coch(cond).reverb_cond);
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Bold');
xlabel('Time [s]','FontSize',font_sz,'FontWeight','bold');
ylabel('Freqeuncy [kHz]','FontSize',font_sz,'FontWeight','bold');

set(gcf,'color','w');
save_name = fullfile(save_dir,[coch(cond).reverb_cond,' Cochleagram','.png']);
export_fig(save_name);