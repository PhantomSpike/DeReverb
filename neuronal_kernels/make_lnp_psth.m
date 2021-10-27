%% Define and load vars
%Real data
cluster_dir_real = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Absolute path to the real neuronal data containing the psths
% LNP kernels
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms'; %Absolute path to the kernel used for the prediction
%Get the cochleagrams
coch_file1 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Ronnie_PLP/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch_file2 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch_file3 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Noah_Derekah/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch{1} = load(coch_file1,'coch');
coch{2} = load(coch_file2,'coch');
coch{3} = load(coch_file3,'coch');

%Where to save the data
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_noneuron_norm/perfreq_noneuro/ridge/10ms/Predicted_PSTHs_noneuro_norm';


rnd_seed = 7; %Set the random seed for reproducibility 
neuron_norm = 0; %Logical flag whether to normalize the neuronal psth before fitting the NL
NPSP_th = 40;
n_h = 20; %Number of history steps
n_trials = 10; %How many trials to simulate
t_edges_s = coch{1}.coch(1).t; %Get the edges in sec as used in the cochleagrams
plot_check = 1; %Flag whether to plot the fits as they appear

%% Select the clusters to fit
load(fullfile(cluster_dir,'info.mat'),'info'); %Load the info file with the meta data
ix = info.NPSP<NPSP_th; %Find only the neurons that are below the NPSP th
%Get the necessary fields
cluster_ids = info.cluster_id(ix);
pen_names = info.pen_name(ix);
animal_names = info.animal_name(ix);
qualities = info.quality(ix);
NPSPs = info.NPSP(ix);
n_clust = sum(ix);

for j = 1:3
    
    %Initialise
    small1_X_ft = [];
    small2_X_ft = [];
    big1_X_ft = [];
    big2_X_ft = [];
    
    %% Normalize the cochleagrams
    switch j
        
        case{1,2}
            
            small1_X_ft = coch{j}.coch(3).X_ft;
            small2_X_ft = coch{j}.coch(4).X_ft;
            
        case {3}
            
            small1_X_ft = coch{j}.coch(1).X_ft;
            small2_X_ft = coch{j}.coch(2).X_ft;
    end
    
    big1_X_ft = coch{j}.coch(5).X_ft;
    big2_X_ft = coch{j}.coch(6).X_ft;
    
    mean_reverb = mean([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],2); std_reverb = std([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],[],2);
    
    small1_X_ft_norm = (small1_X_ft - mean_reverb)./std_reverb;
    small2_X_ft_norm = (small2_X_ft - mean_reverb)./std_reverb;
    big1_X_ft_norm = (big1_X_ft - mean_reverb)./std_reverb;
    big2_X_ft_norm = (big2_X_ft - mean_reverb)./std_reverb;
    
    
    %% Tensorize
    fprintf('== Tensorizing the cochleagrams ==\n');tic;
    small1_X_fht = tensorize(small1_X_ft_norm,n_h);
    small2_X_fht = tensorize(small2_X_ft_norm,n_h);
    big1_X_fht = tensorize(big1_X_ft_norm,n_h);
    big2_X_fht = tensorize(big2_X_ft_norm,n_h);
    fprintf('== Done! This took %0.fs ==\n',toc);
    
    %% Concatenate together
    small_X_fht{j} = cat(3,small1_X_fht,small2_X_fht);
    big_X_fht{j} = cat(3,big1_X_fht,big2_X_fht);
    n_t = size(small_X_fht{j},3);
    
end

%% Make simulated LNP neuron
temp_psth = cell(n_clust,1);
count = 0;
tic;
rng(rnd_seed);

for c = 1:n_clust
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clust);
    %Initialise vars
    y_t_small_pois = [];
    y_t_big_pois = [];
    temp_real = [];
    y_temp = [];
    temp = [];

    %Get the psth (y_real) for every cluster
    temp_real = load(fullfile(cluster_dir_real,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
    y = temp_real.data;
    
    %Get the animal name
    animal_name = temp_real.data.params.animal_name;
    
    % Make histograms for the 2 reverb conditions: stim
    %For "Ronnie", "PLP", "Cork", "Kilkenny", "Derry" this is:
    % [1,2]->anech, [3,4]->small reverb, [5,6]->large reverb
    %For "Derekah" and "Noah" this is:
    % [1,2]->small reverb, [3,4] ->medium reverb,[5,6]->large reverb
    
    n_stim = length(y.stim);
    n_rep = length(y.stim(1).repeat);
    
    for s = 1:n_stim
        psth_temp = [];
        for r = 1:n_rep
            psth_temp(r,:) = histc(y.stim(s).repeat(r).spiketimes,t_edges_s);
        end
        y_temp(s,:) = mean(psth_temp);
    end
    
    switch animal_name
        
        case {"Ronnie","PLP"}
            y_t_small = [y_temp(3,:),y_temp(4,:)];
            k = 1;
            
        case {"Cork","Kilkenny","Derry"}
            y_t_small = [y_temp(3,:),y_temp(4,:)];
            k = 2;
            
        case {"Noah","Derekah"}
            y_t_small = [y_temp(1,:),y_temp(2,:)];
            k = 3;
    end
    
    y_t_big = [y_temp(5,:),y_temp(6,:)];
    y_t_total = [y_t_small,y_t_big];
    
    %Optionally normalize the neuronal data before fitting NL
    if neuron_norm
        mean_y_t_total = mean(y_t_total);
        std_y_t_total = std(y_t_total);
        y_t_total = (y_t_total-mean_y_t_total)./std_y_t_total;
    end
    
    %Get the single lin kernel for every cluster
    temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')));
    ker = temp.kernel.reverb.main{end};
    
    %Make the linear prediction
    y_hat_lin_t_small = kernelconv(small_X_fht{k}, ker);
    y_hat_lin_t_big = kernelconv(big_X_fht{k}, ker);
    
    %Combine the predicted data together
    y_hat_lin_t_total = [y_hat_lin_t_small,y_hat_lin_t_big];
    
    %Fit a single NL for both conditions
    model = getlnmodel3(y_hat_lin_t_total, y_t_total);
    
    %Get the output of the LN models
    y_hat_ln_t_small = lnmodel(model.params, y_hat_lin_t_small);
    y_hat_ln_t_big = lnmodel(model.params, y_hat_lin_t_big);
    
    %Save the predicted values for later plotting
    y_hat_lin_t_total_save{c} = y_hat_lin_t_total;
    y_t_total_save{c} = y_t_total;
    y_hat_ln_t_total_save{c} = [y_hat_ln_t_small,y_hat_ln_t_big];
    
    %Remove -ve values (ReLU)
    y_hat_ln_t_small = max(y_hat_ln_t_small,0);
    y_hat_ln_t_big = max(y_hat_ln_t_big,0);
    
    %Add Poisson noise in 'virtual trials'
    for tr = 1:n_trials
        y_t_small_pois(tr,:) = poissrnd(y_hat_ln_t_small);
        y_t_big_pois(tr,:) = poissrnd(y_hat_ln_t_big);
    end
    
    %Average the virtual trials
    temp_psth{c}.y_t_small_mean = mean(y_t_small_pois);
    temp_psth{c}.y_t_big_mean = mean(y_t_big_pois);
    
    %Save the poisson ones too
    y_hat_lnp_t_total_save{c} = [temp_psth{c}.y_t_small_mean,temp_psth{c}.y_t_big_mean];
    
    %Record if there are any NaN values
    if sum(isnan(temp_psth{c}.y_t_small_mean))>0 || sum(isnan(temp_psth{c}.y_t_big_mean))>0
        count = count + 1;
        save_nan.small_count(count) = sum(isnan(temp_psth{c}.y_t_small_mean));
        save_nan.big_count(count) = sum(isnan(temp_psth{c}.y_t_big_mean));
        save_nan.cluster{count} = strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_');
    end
    
    %Optionally, plot the fits to check if they are sensible
    if plot_check
        subplot(1,2,1);
        [x, y] = binplot(y_hat_lin_t_total, y_t_total, 100); scatter(x,y); hold on;
        [x, y] = binplot(y_hat_lin_t_total, [y_hat_ln_t_small,y_hat_ln_t_big], 100); plot(x,y,'r');
        title('LN model');
        xlabel('yhat');
        ylabel('y');
        hold off;
        subplot(1,2,2);
        [x, y] = binplot(y_hat_lin_t_total, y_t_total, 100); scatter(x,y); hold on;
        [x, y] = binplot(y_hat_lin_t_total, [temp_psth{c}.y_t_small_mean,temp_psth{c}.y_t_big_mean], 100); plot(x,y,'k');
        title('LNP model');
        xlabel('yhat');
        ylabel('y');
        hold off;
        pause(0.25);
    end
end

%Raise a warning if NaN values exist
if count>0
    warning('There were NaN values when generating the data');
end

%% Plot the results as subplots
[temp_dir,curr_dir] = fileparts(save_dir);
save_dir_plots = fullfile(temp_dir,'Plots');
if ~exist(save_dir_plots, 'dir')
    mkdir(save_dir_plots);
end
save_dir_specific = fullfile(save_dir_plots,curr_dir);
if ~exist(save_dir_specific, 'dir')
    mkdir(save_dir_specific);
end
row = 5;
col = 6;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = 0.05; space_h = 0.06; space_v = 0.09;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
clust_per_batch = row*col;
axis_sz = 20;

num_batch = ceil(n_clust/clust_per_batch);
curr_n_clust = clust_per_batch;

count1 = 0;
count2 = 0;

for b = 1:num_batch
    figure('units','normalized','outerposition',[0 0 1 1]);
    save_name = fullfile(save_dir_specific,['batch_',num2str(b),'_LN_model.svg']);
    
    if b == num_batch
        curr_n_clust = n_clust - (num_batch-1)*clust_per_batch;
    end
    
    for j = 1:curr_n_clust
        count1 = count1+1;
        subplot('position',pos{j});
        [x, y] = binplot(y_hat_lin_t_total_save{count1}, y_t_total_save{count1}, 100); scatter(x,y); hold on;
        [x, y] = binplot(y_hat_lin_t_total_save{count1},y_hat_ln_t_total_save{count1} , 100); plot(x,y,'r');
        xlabel('yhat');
        ylabel('y');
        set(gcf,'color','w');
        hold off;
        set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    end
    saveas(gcf,save_name);
    close;
end

curr_n_clust = clust_per_batch;
for b = 1:num_batch
    figure('units','normalized','outerposition',[0 0 1 1]);
    save_name = fullfile(save_dir_specific,['batch_',num2str(b),'_LNP_model.svg']);
    
    if b == num_batch
        curr_n_clust = n_clust - (num_batch-1)*clust_per_batch;
    end
    
    for j = 1:curr_n_clust
        count2 = count2+1;
        subplot('position',pos{j});
        [x, y] = binplot(y_hat_lin_t_total_save{count2}, y_t_total_save{count2}, 100); scatter(x,y); hold on;
        [x, y] = binplot(y_hat_lin_t_total_save{count2},y_hat_lnp_t_total_save{count2} , 100); plot(x,y,'k');
        xlabel('yhat');
        ylabel('y');
        set(gcf,'color','w');
        hold off;
        set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    end
    saveas(gcf,save_name);
    close;
end

%% Save the results
fprintf('== Saving the results ==\n');
for cl = 1:n_clust
    psth.y_t_small_mean = temp_psth{cl}.y_t_small_mean;
    psth.y_t_big_mean = temp_psth{cl}.y_t_small_mean;
    save(fullfile(save_dir,strjoin({animal_names{cl},pen_names{cl},num2str(cluster_ids(cl))},'_')),'psth');
end

clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.quality = qualities;
info.pen_name = pen_names;
save(fullfile(save_dir,'info'),'info');

fprintf('== Done! This took %0.fs ==\n',toc);
