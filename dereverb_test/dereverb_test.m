%% Define and load vars
kernel_path = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms/kernel.mat'; %Absolute path to the kernel used for the prediction
coch_path_train = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat'; %Absolute path to the cochleagram used for the prediction
coch_path_test = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/100s_dereverb_new_exp_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/Dereverb_test';

norm = 1; %Optionally normalize the cohcleagrams by removing the mean and dividing by the std

fprintf('== Loading the data ==\n');tic;
load(kernel_path);
coch_train = load(coch_path_train);
% room_names{1} = 'anech_orig'; room_names{2} = 'small_orig'; room_names{3} = 'med_orig'; room_names{4} = 'big_orig';
% room_names{5} = 'anech'; room_names{6} = 'small'; room_names{7} = 'med'; room_names{8} = 'big';

%Get the cochleagrams of the training data
anech_X_ft_train = coch_train.coch(1).X_ft;
small_X_ft_train = coch_train.coch(2).X_ft;
med_X_ft_train = coch_train.coch(3).X_ft;
large_X_ft_train = coch_train.coch(4).X_ft;

%Get the cochleagrams of the test data
coch_test = load(coch_path_test);

anech_X_ft_test = coch_test.coch(1).X_ft;
small_X_ft_test = coch_test.coch(2).X_ft;
med_X_ft_test = coch_test.coch(3).X_ft;
large_X_ft_test = coch_test.coch(4).X_ft;

%Get the kernels
n_ker = length(kernels.small);
for k = 1:n_ker
    kernel_small{k} = kernels.small{k}.main{end};
    kernel_med{k} = kernels.med{k}.main{end};
    kernel_large{k} = kernels.big{k}.main{end};
end
n_h = kernels.n_h;
n_t = size(small_X_ft_test,2);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Normalize the cochleagrams
mean_anech = mean(anech_X_ft_train,2); std_anech = std(anech_X_ft_train,[],2);
mean_small = mean(small_X_ft_train,2); std_small = std(small_X_ft_train,[],2);
mean_med = mean(med_X_ft_train,2); std_med = std(med_X_ft_train,[],2);
mean_large = mean(large_X_ft_train,2); std_large = std(large_X_ft_train,[],2);

anech_X_ft_norm = (anech_X_ft_test - mean_anech)./std_anech;
small_X_ft_norm = (small_X_ft_test - mean_small)./std_small;
med_X_ft_norm = (med_X_ft_test - mean_med)./std_med;
large_X_ft_norm = (large_X_ft_test - mean_large)./std_large;

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small_X_fht = tensorize(small_X_ft_norm,n_h);
med_X_fht = tensorize(med_X_ft_norm,n_h);
large_X_fht = tensorize(large_X_ft_norm,n_h);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Make predictions
fprintf('== Making predicted cochleagrams ==\n');tic;
small_anech_X_ft_hat = zeros(n_ker,n_t);
med_anech_X_ft_hat = zeros(n_ker,n_t);
large_anech_X_ft_hat = zeros(n_ker,n_t);

for k = 1:n_ker
    small_anech_X_ft_hat(k,:) = kernelconv(small_X_fht, kernel_small{k});
    med_anech_X_ft_hat(k,:) = kernelconv(med_X_fht, kernel_med{k});
    large_anech_X_ft_hat(k,:) = kernelconv(large_X_fht, kernel_large{k});
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Undo normalization to convert back to original values
X_ft_hat.small = (small_anech_X_ft_hat.*std_small) + mean_small;
X_ft_hat.med = (med_anech_X_ft_hat.*std_med) + mean_med;
X_ft_hat.large = (large_anech_X_ft_hat.*std_large) + mean_large;

X_ft_orig.anech = anech_X_ft_test;
X_ft_orig.small = small_X_ft_test;
X_ft_orig.med = med_X_ft_test;
X_ft_orig.large = large_X_ft_test;

pred_rooms = fieldnames(X_ft_hat);
orig_rooms = fieldnames(X_ft_orig);

if norm
    
    for j = 1:length(pred_rooms)
        room = pred_rooms{j};
        X_ft_hat.(room) = (X_ft_hat.(room) - mean(X_ft_hat.(room)(:)))./std(X_ft_hat.(room)(:));
    end
    
    for i = 1:length(orig_rooms)
        room = orig_rooms{i};
        X_ft_orig.(room) = (X_ft_orig.(room) - mean(X_ft_orig.(room)(:)))./std( X_ft_orig.(room)(:));
    end
    
end
%% Compute the CC and MSE between the original cochleagram and the anechoic
kfold = 10; %How many folds to use for the testing set

%Find out how many folds need to be ran
n_t = size(X_ft_orig.anech,2);
chunk_sz = floor(n_t/kfold);

%Run the model k times
for k = 1:kfold
    
    start_ix = (k-1)*chunk_sz + 1;
    end_ix = start_ix + chunk_sz -1;
    val_ix = [start_ix:end_ix];
    
    c = 0;
    for i = 1:3
        %Compare original anechoic and reverberant cochleagrams
        c = c+1;
        pred_room = pred_rooms{i};
        A = []; B = []; C = [];
        A = X_ft_orig.anech(:,val_ix); B = X_ft_orig.(pred_room)(:,val_ix);
        r = corrcoef(A(:), B(:));
        pred.CC(k,c) = r(1,2);
        mse = mean((A(:) - B(:)).^2);
        pred.MSE(k,c) = mse;
        
        %Compare original anechoic and predicted anechoic cochleagrams
        c = c+1;
        C = X_ft_hat.(pred_room)(:,val_ix);
        r = corrcoef(A(:), C(:));
        pred.CC(k,c) = r(1,2);
        mse = mean((A(:) - C(:)).^2);
        pred.MSE(k,c) = mse;
    end
    
end
pred.comparisons = {'Anech-Small Orig','Anech-Small Pred', 'Anech-Med Orig','Anech-Med Pred','Anech-Large Orig','Anech-Large Pred'};

%% Compute pvals for the comparisons

%Small room
pval.CC.small = signrank(pred.CC(:,1), pred.CC(:,2));
pval.MSE.small = signrank(pred.MSE(:,1), pred.MSE(:,2));

%Medium room
pval.CC.med = signrank(pred.CC(:,3), pred.CC(:,4));
pval.MSE.med = signrank(pred.MSE(:,3), pred.MSE(:,4));

%Large room
pval.CC.large = signrank(pred.CC(:,5), pred.CC(:,6));
pval.MSE.large = signrank(pred.MSE(:,5), pred.MSE(:,6));

%% Compute the median and mean differences

%Median
    %Small room
    differ.med.CC.small = nanmedian(pred.CC(:,2) - pred.CC(:,1));
    differ.med.MSE.small = nanmedian(pred.MSE(:,2) - pred.MSE(:,1));

    %Medium room
    differ.med.CC.med = nanmedian(pred.CC(:,4) - pred.CC(:,3));
    differ.med.MSE.med = nanmedian(pred.MSE(:,4) - pred.MSE(:,3));

    %Large room
    differ.med.CC.large = nanmedian(pred.CC(:,6) - pred.CC(:,5));
    differ.med.MSE.large = nanmedian(pred.MSE(:,6) - pred.MSE(:,5));
    
 %Mean
    %Small room
    differ.mean.CC.small = nanmean(pred.CC(:,2) - pred.CC(:,1));
    differ.mean.MSE.small = nanmean(pred.MSE(:,2) - pred.MSE(:,1));

    %Medium room
    differ.mean.CC.med = nanmean(pred.CC(:,4) - pred.CC(:,3));
    differ.mean.MSE.med = nanmean(pred.MSE(:,4) - pred.MSE(:,3));

    %Large room
    differ.mean.CC.large = nanmean(pred.CC(:,6) - pred.CC(:,5));
    differ.mean.MSE.large = nanmean(pred.MSE(:,6) - pred.MSE(:,5));
    
%% Plot the CC

data.mean = mean(pred.CC);
data.sem = std(pred.CC)/sqrt(n_ker);
errhigh = data.sem;
errlow  = data.sem;


X{1} = categorical({'Anech-Small Original'});
X{2} = categorical({'Anech-Small Predicted'});
X{3} = categorical({'Anech-Med Original'});
X{4} = categorical({'Anech-Med Predicted'});
X{5} = categorical({'Anech-Large Original'});
X{6} = categorical({'Anech-Large Predicted'});

X_final = categorical({'Anech-Small Original','Anech-Small Predicted', 'Anech-Med Original','Anech-Med Predicted','Anech-Large Original','Anech-Large Predicted'});
X_final = reordercats(X_final,{'Anech-Small Original','Anech-Small Predicted', 'Anech-Med Original','Anech-Med Predicted','Anech-Large Original','Anech-Large Predicted'});
font_sz = 30;

small_color = [0.592, 0.737, 0.384];
med_color = [0.3825, 0.5550, 0.2800];
large_color = [0.173, 0.373, 0.176];

colors = {small_color, small_color, med_color, med_color, large_color, large_color};


figure('units','normalized','outerposition',[0 0 1 1]);


hold on;

for i = 1:6
    hb(i) = bar(X{i}, data.mean(i));
    hb(i).FaceColor = colors{i};
end

ylim([0.7 1]);
er = errorbar(X_final, data.mean, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

annotation('textbox',[0.8 0.85 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',large_color,'FontSize',font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.79 0.1 0.1],'String', sprintf('Medium room'),'LineStyle','none','Color',med_color,'FontSize',font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.73 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',small_color,'FontSize',font_sz,'FontWeight','bold');

% title('Corr coeff of predicted with actual cochleagrams');
hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir,['Pred_test_CC','.svg']);
saveas(gcf, save_name);
close;

%% Plot the MSE
data.mean = mean(pred.MSE);
data.sem = std(pred.MSE)/sqrt(n_ker);
errhigh = data.sem;
errlow  = data.sem;

figure('units','normalized','outerposition',[0 0 1 1]);


hold on;

for i = 1:6
    hb(i) = bar(X{i}, data.mean(i));
    hb(i).FaceColor = colors{i};
end

er = errorbar(X_final, data.mean, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

% title('MSE of predicted with actual cochleagrams');

annotation('textbox',[0.8 0.9 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',large_color,'FontSize',font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.84 0.1 0.1],'String', sprintf('Medium room'),'LineStyle','none','Color',med_color,'FontSize',font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.78 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',small_color,'FontSize',font_sz,'FontWeight','bold');

hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir,['Pred_test_MSE','.svg']);
saveas(gcf, save_name);
close;

%% Plot original and predictions
%Find the max value

if norm
    max_val = 3; %How much to subtract from the cochleagram for normalization
    up_lim_dB = 0; %The upper limit in dB
    low_lim_dB = -5; %The lower limit in dB
else
    max_val = 40; %How much to subtract from the cochleagram for normalization
    up_lim_dB = 0; %The upper limit in dB
    low_lim_dB = -80; %The lower limit in dB
end
start_s = 50;
end_s = 54;
t = coch_test.coch(1).t;
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);

dt_s = 1;
dt = mean(diff(t));
dt_samples = round(dt_s/dt);

color_map = 'inferno';
freqs = fliplr(coch_test.coch(1).params.freqs); 
n_f = length(freqs);
n_tlab = 8;
skip_f = 2;
skip_t = round(length(t)/n_tlab);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(t((tm-1)*skip_t +1),'%.1f');
end

%% Plot individual cochleagrams
sz = 40;
all_font_sz = sz;
x_font_sz = sz;
y_font_sz = sz;
font_type = 'Liberation Sans';

%1 - Original anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_orig.anech(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time(s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%2 - Small reverb
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_orig.small(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Small_coch_example.svg');
saveas(gcf,save_name);
close all;

%3 - Medium reverb
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_orig.med(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Med_coch_example.svg');
saveas(gcf,save_name);
close all;

%4 - Large reverb
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_orig.large(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Large_coch_example.svg');
saveas(gcf,save_name);
close all;

%5 - Small reverb predicted anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_hat.small(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Small_pred_anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%6 - Medium reverb predicted anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_hat.med(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Med_pred_anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;

%7 - Big reverb predicted anechoic
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
imagesc(X_ft_hat.large(:,ix_start:ix_end)-max_val);
colormap(color_map);
caxis([low_lim_dB up_lim_dB]);
yticks([1:skip_f:n_f]);
yticklabels(y_labels);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
xlabel('Time (s)','FontSize',x_font_sz,'FontWeight','Normal');
ylabel('Freqeuncy (kHz)','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'Large_pred_anechoic_coch_example.svg');
saveas(gcf,save_name);
close all;