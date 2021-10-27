%In this scipt I will:
% - Load the STRFs of all neurons with certain NPSP
% - Unroll into a vector
% - Combine in a matrix - H
% - Load the PSTHs of the same neurons
% - Combine into a matrix - R
% - Compute the pseudo-inverse of H - F
% - S = FR to make a predicted cochleagram

clear all; close all;
%% Params
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms/'; %Directory with the STRFs
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Directory with all the clusters
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/Decoding/Flat_prior';
normz = 'global'; % 'perfreq', 'global', 'none'
method = 'STA'; %  'NRC', 'reg', 'STA'
coch_type = 'specpower';
dt_ms = 10;
% normz = 'perfreq';
coch_file = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
get_neurons = 'all';
NPSP_th = 40;
n_f = 30;
chop_ms = 190; 
h_max_ms = 200;


r_type{1} = 'small';
r_type{2} = 'large';
n_rooms = length(r_type);

%% Load the STRF data
load(fullfile(kernel_dir,'info'),'info');
temp_files = dir([kernel_dir,'/*.mat']);
temp = load(fullfile(temp_files(1).folder,temp_files(1).name));
model = temp.kernel.model;

%% Select the data
switch get_neurons
    case 'all'
        ix_qualia = ones(length(info.cluster_id),1); %Get all the neurons
    case 'good'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'good'),info.quality,'UniformOutput',false)); %Find all the good units
    case 'mua'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'mua'),info.quality,'UniformOutput',false)); %Find all the mua units
end
ix_npsp = info.NPSP<NPSP_th; %Find all the neurons below certain NPSP
ix = ix_qualia & ix_npsp; %Find the intersection of the two

NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);

%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
qualities = qualities(ix_select);
n_clust = length(cluster_ids);

%% Load all the selected kernels
% fprintf('== Loading the data ==\n');tic;
% 
% for k = 1:n_clust
%     c_name = fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
%     load(c_name,'kernel');
%     
%     %Get the small and large STRFs
%     k_fh_small = kernel.small.main{end}.k_fh(:,2:end); 
%     k_fh_large = kernel.big.main{end}.k_fh(:,2:end); 
%     
%     %Get the bias
%     bias_small = kernel.small.main{end}.c;  
%     bias_large = kernel.big.main{end}.c;  
%     
%     %Make the H matrix for small and large - (n_f*n_h) x N
%     % n_f = no of frequencies; n_h = no of history steps; N = no of neurons
% %     k_fh_small = [bias_small;k_fh_small(:)];
% %     k_fh_large = [bias_large;k_fh_large(:)];
%     
%     k_fh_small = (k_fh_small-nanmean(k_fh_small(:)))./nanstd(k_fh_small(:));
%     k_fh_large = (k_fh_large-nanmean(k_fh_large(:)))./nanstd(k_fh_large(:));
%     
%     H.small(:,k) = k_fh_small(:);
%     H.large(:,k) = k_fh_large(:);
% %     H.small(:,k) = [bias_small; k_fh_small(:)];
% %     H.large(:,k) = [bias_large; k_fh_large(:)];
%     
% end
% 
% freqs = kernel.freqs;
% n_f = length(freqs);
% n_h = kernel.n_h;
% 
% fprintf('== Done! This took %0.fs ==\n',toc);


%% Load the cochleagrams

coch_in = load(coch_file,'coch');
coch = coch_in.coch;

anech1_X_ft = coch(1).X_ft;
anech2_X_ft = coch(2).X_ft;
small1_X_ft = coch(3).X_ft;
small2_X_ft = coch(4).X_ft;
large1_X_ft = coch(5).X_ft;
large2_X_ft = coch(6).X_ft;

switch normz
    
    case 'perfreq'
        mean_anech = mean([anech1_X_ft, anech2_X_ft],2); std_anech = std([anech1_X_ft, anech2_X_ft],[],2);
        mean_small = mean([small1_X_ft, small2_X_ft],2); std_small = std([small1_X_ft, small2_X_ft],[],2);
        mean_large = mean([large1_X_ft, large2_X_ft],2); std_large = std([large1_X_ft, large2_X_ft],[],2);
        anech1_X_ft = (anech1_X_ft - mean_anech)./std_anech; anech2_X_ft = (anech2_X_ft - mean_anech)./std_anech;
        small1_X_ft = (small1_X_ft - mean_small)./std_small; small2_X_ft = (small2_X_ft - mean_small)./std_small;
        large1_X_ft = (large1_X_ft - mean_large)./std_large; large2_X_ft = (large2_X_ft - mean_large)./std_large;
        
    case 'global'
        mean_anech = mean([anech1_X_ft(:); anech2_X_ft(:)]); std_anech = std([anech1_X_ft(:); anech2_X_ft(:)]);
        mean_small = mean([small1_X_ft(:); small2_X_ft(:)]); std_small = std([small1_X_ft(:); small2_X_ft(:)]);
        mean_large = mean([large1_X_ft(:); large2_X_ft(:)]); std_large = std([large1_X_ft(:); large2_X_ft(:)]);
        anech1_X_ft = (anech1_X_ft - mean_anech)./std_anech; anech2_X_ft = (anech2_X_ft - mean_anech)./std_anech;
        small1_X_ft = (small1_X_ft - mean_small)./std_small; small2_X_ft = (small2_X_ft - mean_small)./std_small;
        large1_X_ft = (large1_X_ft - mean_large)./std_large; large2_X_ft = (large2_X_ft - mean_large)./std_large;
        
    case 'none'
        
end

X_ft_orig.anech = [anech1_X_ft, anech2_X_ft];
X_ft_orig.small = [small1_X_ft, small2_X_ft];
X_ft_orig.large = [large1_X_ft, large2_X_ft];
        
dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used  
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram

save_dir_full = fullfile(save_dir,normz);
if ~exist(save_dir_full,'dir')
    mkdir(save_dir_full)
end

%% Tensorize and make design matrix

small1_X_fht = tensorize(small1_X_ft,n_h);
small2_X_fht = tensorize(small2_X_ft,n_h);
large1_X_fht = tensorize(large1_X_ft,n_h);
large2_X_fht = tensorize(large1_X_ft,n_h);

%Combine the two stimuli together
small_X_fht = cat(3, small1_X_fht, small2_X_fht);
large_X_fht = cat(3, large1_X_fht, large2_X_fht);

T = length(small_X_fht);

%Reshape into a design matrix S - (n_f*n_h) x T
small_X_fh_t = reshape(small_X_fht,n_f*n_h,T);
large_X_fh_t = reshape(large_X_fht,n_f*n_h,T);

%% Load the spike times, make the PSTHs and the R matrix

sprintf('== Making the PSTHs ==\n');tic;

for c = 1:n_clust

    temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    animal_name = temp.data.params.animal_name;
    
    %Get the PSTH for the two conditions
    psth = make_psth(temp.data, t_edges_s, animal_name);
    
    psth_small = [psth.small_1, psth.small_2];
    psth_large = [psth.large_1, psth.large_2];
    
    %Make the R matrix for the two conditions - N x T
    % N = no of neurons; T = number of time points
    R.small(c,:) = psth_small; 
    R.large(c,:) = psth_large; 
    
end

fprintf('== Done! This took %0.fs ==\n',toc);

%% Fit the STRFs usign normalized reverse correlation

switch method
    case 'reg'
        lambda = 10000;
        H.small =  inv(small_X_fh_t*small_X_fh_t' + lambda.*eye(n_f*n_h))*small_X_fh_t*R.small';
        H.large =  inv(large_X_fh_t*large_X_fh_t' + lambda.*eye(n_f*n_h))*large_X_fh_t*R.large';
        
    case 'NRC'
        lambda = 0;
        H.small =  inv(small_X_fh_t*small_X_fh_t')*small_X_fh_t*R.small';
        H.large =  inv(large_X_fh_t*large_X_fh_t')*large_X_fh_t*R.large';
        
    case 'STA'
        lambda = 0;
        H.small =  small_X_fh_t*R.small';
        H.large =  large_X_fh_t*R.large';
        
end

save_dir_full2 = fullfile(save_dir_full,method);
if ~exist(save_dir_full2,'dir')
    mkdir(save_dir_full2)
end

%% Make the predicted cochleagrams

%Make the flat prior reconstruction matrix F
F.small = inv(H.small*H.small')*H.small;
F.large = inv(H.large*H.large')*H.large;

% F.small = pinv(H.small');
% F.large = pinv(H.large');

%Make the predicted cochleagram for small and large
X_fh_t.small = F.small*R.small;
X_fh_t.large = F.large*R.large;

%Convert to the size of the original cochleagrams

for t = 1:2*length(t_edges_s)
    
    %Get the Cochleagram w/ history X_fh for each time point t but ignore
    %the bias
    X_fh_small = X_fh_t.small(1:end,t); 
    X_fh_large = X_fh_t.large(1:end,t);
    
    %Reshape to the same form as the STRF
    X_fh_small = reshape(X_fh_small,n_f,n_h);
    X_fh_large = reshape(X_fh_large,n_f,n_h);
    
    %Take all the frequency channels and the current time bin from the
    %history
    X_ft_hat.small(:,t) = X_fh_small(:,end);
    X_ft_hat.large(:,t) = X_fh_large(:,end);
    
end

%% Get the CC
pred_rooms = {'small', 'large'};
orig_rooms = {'anech','small','large'};
c = 0;
for p = 1:length(pred_rooms)
    pred_room = pred_rooms{p};
    for o = 1:length(orig_rooms)
        c = c + 1;
        orig_room = orig_rooms{o};
        r = corrcoef(X_ft_hat.(pred_room)(:), X_ft_orig.(orig_room)(:));
        CC(c).corrcoef = r(1,2);
        CC(c).predicted = pred_room;
        CC(c).original = orig_room;
    end
end

%% Get the MSE

c = 0;
for p = 1:length(pred_rooms)
    pred_room = pred_rooms{p};
    for o = 1:length(orig_rooms)
        c = c + 1;
        orig_room = orig_rooms{o};
        mse = mean((X_ft_hat.(pred_room)(:) - X_ft_orig.(orig_room)(:)).^2);
        MSE(c).mse = mse;
        MSE(c).predicted = pred_room;
        MSE(c).original = orig_room;
    end
end

%% Plot the CC


%Plot the data comparing EXC and all GAD only
X{1} = categorical({'Small - Anech'});
X{2} = categorical({'Small - Small'});
X{3} = categorical({'Small - Large'});
X{4} = categorical({'Large - Anech'});
X{5} = categorical({'Large - Small'});
X{6} = categorical({'Large - Large'});
X_final = categorical({'Small - Anech', 'Small - Small','Small - Large','Large - Anech','Large - Small','Large - Large'});
X_final = reordercats(X_final,{'Small - Anech', 'Small - Small','Small - Large','Large - Anech','Large - Small','Large - Large'});
font_sz = 30;
anechoic_color = [197, 225, 165]/255;
small_color = [0.592, 0.737, 0.384];
big_color = [0.173, 0.373, 0.176];

colors = {anechoic_color, small_color, big_color,anechoic_color, small_color, big_color};


figure('units','normalized','outerposition',[0 0 1 1]);


hold on;

for i = 1:6
    hb(i) = bar(X{i}, CC(i).corrcoef);
    hb(i).FaceColor = colors{i};
end


title('Corr coeff of predicted with actual cochleagrams');
hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_full2,['CC','.svg']);
saveas(gcf, save_name);
close;

%% Plot the MSE
font_sz = 30;
anechoic_color = [197, 225, 165]/255;
small_color = [0.592, 0.737, 0.384];
big_color = [0.173, 0.373, 0.176];

colors = {anechoic_color, small_color, big_color,anechoic_color, small_color, big_color};


figure('units','normalized','outerposition',[0 0 1 1]);


hold on;

for i = 1:6
    hb(i) = bar(X{i}, MSE(i).mse);
    hb(i).FaceColor = colors{i};
end

title('MSE of predicted with actual cochleagrams');
hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir_full2,['MSE','.svg']);
saveas(gcf, save_name);
close;

%% Plot the predicted cochleagrams
title_names = {'Original Anechoic','Original Small','Original Large',...
    'Predicted anechoic - Small','Predicted anechoic - Large'};

start_s = 50;
end_s = 54;

t = [dt_ms/1000:dt_ms/1000:T*(dt_ms/1000)];
n_tlab = 10;
skip_f = 4;
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);
skip_t = round(length(t)/n_tlab);

row = 2;
col = 3;
per = 0.03;
edgel = 0.07; edger = 0.02; edgeh = 0.05; edgeb = 0.09; space_h = 0.05; space_v = 0.1;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

c = 0;
for f = 1:skip_f:numel(freqs)
    c = c+1;
    y_labels{c} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(t((tm-1)*skip_t +1),'%.0f');
end

figure('units','normalized','outerposition',[0 0 1 1]);

for j = 1:5
    subplot('position',pos{j});
    title_name = title_names{j};
    
    if ismember(j,[1:3])
        orig_room = orig_rooms{j};
        imagesc(X_ft_orig.(orig_room)(:,ix_start:ix_end));
    else
        pred_room = pred_rooms{j-3};
        imagesc(X_ft_hat.(pred_room)(:,ix_start:ix_end));
    end

    colormap('inferno');
%     caxis([-40 20]);
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_t:length(t)]);
    xticklabels(x_labels);
    title(title_name);
    set(gca,'FontName','Arial','FontSize',20,'FontWeight','Normal');
    %     xlabel('Time [s]','FontSize',14,'FontWeight','Normal');
    %     ylabel('Freqeuncy [kHz]','FontSize',14,'FontWeight','Normal');
    
    if j==4
        xlabel('Time (s)');
        ylabel('Frequency (kHz)');
    end
    
    if j==5
        colorbar;
    end
    
end
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir_full2,'BF_predicted_cochleagrams.svg');
saveas(gcf,save_name);
close;


%% Plot the STRFs

ker_per_plot = 25;
% n_groups = floor(n_clust/ker_per_plot);
n_groups = 2;
row = 5;
col = 5;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = per; space_h = 0.05; space_v = 0.06;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
count = 0;

for i  = 1:n_groups
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    for j = 1:ker_per_plot
        count = count + 1;
        subplot('position',pos{j});
        k_fh = reshape(H.small(:,count),n_f,n_h);
        k_fh = k_fh./max(abs(k_fh(:)));
        imagesc(k_fh);caxis([-1 1]); colormap('redblue');
        set(gcf,'color','w');
    end
    save_name = fullfile(save_dir_full2,['Small_STRF_group_',num2str(i),'_lambda_',num2str(lambda),'.svg']);
    saveas(gcf,save_name);
    close;
end

count = 0;
for i  = 1:n_groups
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    for j = 1:ker_per_plot
        count = count + 1;
        subplot('position',pos{j});
        k_fh = reshape(H.large(:,count),n_f,n_h);
        k_fh = k_fh./max(abs(k_fh(:)));
        imagesc(k_fh);caxis([-1 1]); colormap('redblue');
        set(gcf,'color','w');
    end
    save_name = fullfile(save_dir_full2,['Large_STRF_group_',num2str(i),'_lambda_',num2str(lambda),'.svg']);
    saveas(gcf,save_name);
    close;
end
