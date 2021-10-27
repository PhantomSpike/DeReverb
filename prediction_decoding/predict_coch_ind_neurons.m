%In this scipt I will:
% - Load the STRFs of all neurons with certain NPSP
% - Compute the BF for each STRF
% - Combine all STRFs with the same BF
% - Convovle the STRF with the respctive cochleagram, small or large
% - Compute MSE and CC for the predicted clean anech cochleagram vs the
% reverberant ones

%% Params
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms/';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/Decoding';
coch_type = 'specpower';
dt_ms = 10;
normz = 'perfreq';
coch_file = fullfile(fullfile('/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/',coch_type),['/',num2str(dt_ms),'ms','/coch_all_conditions_',coch_type]);
get_neurons = 'all';
NPSP_th = 40;
chop_ms = 190; 
h_max_ms = 200;
freq_up_bound = 17000;
freq_down_bound = 700;

calc_method = 'average'; %How to estimate impotant values
%calc_method -- 'raw', 'average'
%               'raw' - Use the whole receptive field 
%               'average' - Take an average across frequencies first

bf_neurons = 'shared'; %How to get the BF
% bf_method -- 'ind', 'shared'
%              'ind' - Treat the BFs separately 
%              'shared' - Take an log weighted mean BF between the two

bf_method = 'max'; %The method use to extract the bf:
% bf_method -- 'max', 'window_max', 'window_mean'
%              'max' - Take the max across all history
%              'window_max' - Take the max in a given time window
%              'window_mean' - Take the mean in a given time window
bf_window_ms = 190; %The window for calcualting bf if this method is used

r_type{1} = 'small';
r_type{2} = 'big';
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

%% First load all the selected kernels
fprintf('== Loading the data ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'kernel');
    kernels{k,1} = kernel;
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Load the cochleagrams
coch_in = load(coch_file,'coch');
coch = coch_in.coch;

%% Setup and params

anech1_X_ft = coch(1).X_ft;
anech2_X_ft = coch(2).X_ft;
small1_X_ft = coch(3).X_ft;
small2_X_ft = coch(4).X_ft;
big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;

dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used  
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram

%% Get the original cochleagrams with selected frequencies
X_ft_orig.anech = [anech1_X_ft, anech2_X_ft];
X_ft_orig.small = [small1_X_ft, small2_X_ft];
X_ft_orig.big = [big1_X_ft, big2_X_ft];

%% Optionally normalize the data
switch normz
    case {'perfreq','perfreq_noneuro'}

        mean_small = mean([small1_X_ft, small2_X_ft],2); std_small = std([small1_X_ft, small2_X_ft],[],2);
        mean_big = mean([big1_X_ft, big2_X_ft],2); std_big = std([big1_X_ft, big2_X_ft],[],2);

        small1_X_ft_norm = (small1_X_ft - mean_small)./std_small; small2_X_ft_norm = (small2_X_ft - mean_small)./std_small;
        big1_X_ft_norm = (big1_X_ft - mean_big)./std_big; big2_X_ft_norm = (big2_X_ft - mean_big)./std_big;
    case 'none'
end

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;

small1_X_fht = tensorize(small1_X_ft_norm,n_h);
small2_X_fht = tensorize(small2_X_ft_norm,n_h);
big1_X_fht = tensorize(big1_X_ft_norm,n_h);
big2_X_fht = tensorize(big2_X_ft_norm,n_h);

%Combine the two stimuli together
X_fht.small = cat(3,small1_X_fht,small2_X_fht);
X_fht.big = cat(3,big1_X_fht,big2_X_fht);
n_t = size(X_fht.small,3);
fprintf('== Done! This took %0.fs ==\n',toc);

pred_rooms = fieldnames(X_fht);
orig_rooms = fieldnames(X_ft_orig);

%% First compute the BF for every cluster
freqs = fliplr(kernels{1}.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_h = kernel.n_h;
dt_ms = round(kernel.dt_ms);
% chop_ix = round(chop_ms/dt_ms);
bf_window_ix = round(bf_window_ms/dt_ms);
% h = (1:1:chop_ix)';
% h = dt_ms*h;

fprintf('== Calcuating BF ==\n');tic;
for k = 1:n_clust
    sprintf('== Cluster %0.f/%0.f ==\n',k,n_clust);
    for r = 1:n_rooms
        room = r_type{r};
        
        switch model    
            case {'sep','sep_kh'}
                [~,ix] = max(kernels{k}.(room).k_f);
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
                k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                
            case {'ridge','lasso','elastic'}
                
                k_fh = fliplr(kernels{k}.(room).main{end}.k_fh);
%                 k_fh = k_fh(:,1:chop_ix);
                kernel = kernels{k}.(room).main{end};
                
                %Get all +ve and -ve values separately
                k_fh_pos = abs(max(k_fh,0));
                
                %Get the mean across freqeuncies
                k_h_pos = mean(k_fh_pos);
  
                switch bf_method
                    case 'max'
                        k_f = max(k_fh_pos,[],2); %Take the max across history steps
                    case 'window_max'
                        k_f = max(k_fh_pos(:,1:bf_window_ix),[],2); %Take the max in a specified window
                    case 'window_mean'
                        k_f = mean(k_fh_pos(:,1:bf_window_ix),2); %Take the mean in a specified window
                end
                
                [~,ix] = max(k_f); %Find the index of the max frequency
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
        end
        
        %Make a predicted cochleagram
        X_ft_hat_temp.(room) = kernelconv(X_fht.(room), kernel);
        
    end
    bf_mean(k) = 2^(log2(bf(k).small*bf(k).big)/2);
    [~,ix_bf] = min(abs(bf_mean(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
    bf_closest(k) = freqs(ix_bf);
    
    %Undo transformations
    X_ft_hat.small(k,:) = (X_ft_hat_temp.small*std_small(ix_bf)) + mean_small(ix_bf);
    X_ft_hat.big(k,:) = (X_ft_hat_temp.big*std_big(ix_bf)) + mean_big(ix_bf);
    
    
    %Compute the CC and MSE between the predicted anechoic and original
    %cochleagrams
    c = 0;
    for p = 1:length(pred_rooms)
        pred_room = pred_rooms{p};
        for o = 1:length(orig_rooms)
            c = c+1;
            orig_room = orig_rooms{o};
            r = corrcoef(X_ft_hat.(pred_room)(k,:), X_ft_orig.(orig_room)(ix_bf,:));
            pred.CC(k,c) = r(1,2);
            mse = mean((X_ft_hat.(pred_room)(k,:) - X_ft_orig.(orig_room)(ix_bf,:)).^2);
            pred.MSE(k,c) = mse;
        end
    end
    pred.BF(k) = bf_closest(k);
end
pred.comparisons = {'Small-Anech','Small-Small','Small-Large','Large-Anech','Large-Small','Large-Large'};
fprintf('== Done! This took %0.fs ==\n',toc);



%% Plot the CC

data.mean = mean(pred.CC);
data.sem = std(pred.CC)/sqrt(n_clust);
errhigh = data.sem;
errlow  = data.sem;


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
    hb(i) = bar(X{i}, data.mean(i));
    hb(i).FaceColor = colors{i};
end

er = errorbar(X_final, data.mean, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

title('Corr coeff of predicted with actual cochleagrams');
hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir,['Ind_neurons_CC','.svg']);
saveas(gcf, save_name);
close;

%% Plot the MSE

data.mean = mean(pred.MSE);
data.sem = std(pred.MSE)/sqrt(n_clust);
errhigh = data.sem;
errlow  = data.sem;

%Plot the data comparing EXC and all GAD only
font_sz = 30;
anechoic_color = [197, 225, 165]/255;
small_color = [0.592, 0.737, 0.384];
big_color = [0.173, 0.373, 0.176];

colors = {anechoic_color, small_color, big_color,anechoic_color, small_color, big_color};


figure('units','normalized','outerposition',[0 0 1 1]);


hold on;

for i = 1:6
    hb(i) = bar(X{i}, data.mean(i));
    hb(i).FaceColor = colors{i};
end

er = errorbar(X_final, data.mean, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

title('MSE of predicted with actual cochleagrams');
hold off;

set(gcf,'color','w');
set(gca,'TickDir','out'); 
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir,['Ind_neurons_MSE','.svg']);
saveas(gcf, save_name);
close;