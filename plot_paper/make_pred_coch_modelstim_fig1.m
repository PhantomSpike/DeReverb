%% Define and load vars
kernel_path = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms/kernel.mat'; %Absolute path to the kernel used for the prediction
coch_path = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat'; %Absolute path to the cochleagram used for the prediction
fprintf('== Loading the data ==\n');tic;
load(kernel_path);
load(coch_path);
room_names{1} = 'anech_orig'; room_names{2} = 'small_orig'; room_names{3} = 'med_orig'; room_names{4} = 'big_orig';
room_names{5} = 'anech'; room_names{6} = 'small'; room_names{7} = 'med'; room_names{8} = 'big';

%Get the cochleagrams
anech_X_ft = coch(1).X_ft;
small_X_ft = coch(2).X_ft;
med_X_ft = coch(3).X_ft;
big_X_ft = coch(4).X_ft;
%Get the kernels
n_ker = length(kernels.small);
for k = 1:n_ker
    kernel_small{k} = kernels.small{k}.main{end};
    kernel_med{k} = kernels.med{k}.main{end};
    kernel_big{k} = kernels.big{k}.main{end};
end
n_h = kernels.n_h;
n_t = size(small_X_ft,2);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Normalize the cochleagrams
mean_small = mean(small_X_ft,2); std_small = std(small_X_ft,[],2);
mean_med = mean(med_X_ft,2); std_med = std(med_X_ft,[],2);
mean_big = mean(big_X_ft,2); std_big = std(big_X_ft,[],2);

small_X_ft_norm = (small_X_ft - mean_small)./std_small;
med_X_ft_norm = (med_X_ft - mean_med)./std_med;
big_X_ft_norm = (big_X_ft - mean_big)./std_big;

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small_X_fht = tensorize(small_X_ft_norm,n_h);
med_X_fht = tensorize(med_X_ft_norm,n_h);
big_X_fht = tensorize(big_X_ft_norm,n_h);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Make predictions
fprintf('== Making predicted cochleagrams ==\n');tic;
small_anech_X_ft_hat = zeros(n_ker,n_t);
med_anech_X_ft_hat = zeros(n_ker,n_t);
big_anech_X_ft_hat = zeros(n_ker,n_t);

for k = 1:n_ker
    small_anech_X_ft_hat(k,:) = kernelconv(small_X_fht, kernel_small{k});
    med_anech_X_ft_hat(k,:) = kernelconv(med_X_fht, kernel_med{k});
    big_anech_X_ft_hat(k,:) = kernelconv(big_X_fht, kernel_big{k});
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Undo normalization to convert back to original values
anech_X_ft_hat_final.anech = coch(1).X_ft;
anech_X_ft_hat_final.small = (small_anech_X_ft_hat.*std_small) + mean_small;
anech_X_ft_hat_final.med = (med_anech_X_ft_hat.*std_med) + mean_med;
anech_X_ft_hat_final.big = (big_anech_X_ft_hat.*std_big) + mean_big;

anech_X_ft_hat_final.anech_orig = anech_X_ft;
anech_X_ft_hat_final.small_orig = small_X_ft;
anech_X_ft_hat_final.med_orig = med_X_ft;
anech_X_ft_hat_final.big_orig = big_X_ft;

%% Plot original and predictions
%Find the max value
max_val = 40;
start_s = 500;
end_s = 504;
t = coch(1).t;
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);
color_map = 'inferno';
freqs = fliplr(coch(1).params.freqs); 
n_f = length(freqs);
n_tlab = 8;
skip_f = 2;
skip_t = round(length(t)/n_tlab);
row = 2;
col = 4;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = per; space_h = 0.05; space_v = 0.06;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(t((tm-1)*skip_t +1),'%.1f');
end

figure('units','normalized','outerposition',[0 0 1 1]);

for r = 1:8
    room = room_names{r};
    subplot('position',pos{r});
    imagesc(anech_X_ft_hat_final.(room)(:,ix_start:ix_end)-max_val);
    colorbar;
    colormap(color_map);
    caxis([-80 0]);
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_t:length(coch(1).t)]);
    xticklabels(x_labels);
    title(room);
    set(gca,'FontName','Arial','FontSize',15,'FontWeight','Normal');
%     xlabel('Time [s]','FontSize',14,'FontWeight','Normal');
%     ylabel('Freqeuncy [kHz]','FontSize',14,'FontWeight','Normal');
end
set(gcf,'color','w');
set(gca,'TickDir','out');