%% Params
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
data_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/For_analysis/All_data'; %Absolute path to the kernel used for the prediction
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Gain_effects_analysis';

%Get the cochleagrams
coch_file1 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Ronnie_PLP/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch_file2 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch_file3 = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Noah_Derekah/specpower/10ms/coch_all_conditions_specpower.mat'; %Absolute path to the cochleagrams
coch{1} = load(coch_file1,'coch');
coch{2} = load(coch_file2,'coch');
coch{3} = load(coch_file3,'coch');


neuron_norm = 0; %Logical flag whether to normalize the neuronal psth before fitting the NL
NPSP_th = 40;
n_h = 20; %Number of history steps

t_edges_s = coch{1}.coch(1).t; %Get the edges in sec as used in the cochleagrams

plot_check = 0;
col.small = [0.592, 0.737, 0.384];
col.large = [0.173, 0.373, 0.176];
lw = 3;

%% Select the clusters to fit
load(fullfile(kernel_dir,'info.mat'),'info'); %Load the info file with the meta data
ix = info.NPSP<NPSP_th; %Find only the neurons that are below the NPSP th
%Get the necessary fields
cluster_ids = info.cluster_id(ix);
pen_names = info.pen_name(ix);
animal_names = info.animal_name(ix);
qualities = info.quality(ix);
NPSPs = info.NPSP(ix);
n_clust = sum(ix);

%% Make the cochleagrams

for j = 1:3
    
    %Initialise
    small1_X_ft = [];
    small2_X_ft = [];
    big1_X_ft = [];
    big2_X_ft = [];
    
    % Normalize the cochleagrams
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
    
    
    % Tensorize
    fprintf('== Tensorizing the cochleagrams ==\n');tic;
    small1_X_fht = tensorize(small1_X_ft_norm,n_h);
    small2_X_fht = tensorize(small2_X_ft_norm,n_h);
    big1_X_fht = tensorize(big1_X_ft_norm,n_h);
    big2_X_fht = tensorize(big2_X_ft_norm,n_h);
    fprintf('== Done! This took %0.fs ==\n',toc);
    
    % Concatenate together
    small_X_fht{j} = cat(3,small1_X_fht,small2_X_fht);
    big_X_fht{j} = cat(3,big1_X_fht,big2_X_fht);
    n_t = size(small_X_fht{j},3);
    
end

%% Make linear prediction and fit NL
temp_psth = cell(n_clust,1);
y_t_small = cell(n_clust,1);
y_t_big = cell(n_clust,1);
tic;


parfor c = 1:n_clust
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clust);
    %Initialise vars

    temp_real = [];
    y_temp = [];
    temp = [];

    %Get the psth (y_real) for every cluster
    temp_real = load(fullfile(data_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')),'data');
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
            y_t_small{c} = [y_temp(3,:),y_temp(4,:)];
            X_fht_small = small_X_fht{1};
            X_fht_big = big_X_fht{1};
            
        case {"Cork","Kilkenny","Derry"}
            y_t_small{c} = [y_temp(3,:),y_temp(4,:)];
            X_fht_small = small_X_fht{2};
            X_fht_big = big_X_fht{2};
            
        case {"Noah","Derekah"}
            y_t_small{c} = [y_temp(1,:),y_temp(2,:)];
            X_fht_small = small_X_fht{3};
            X_fht_big = big_X_fht{3};
    end
    
    y_t_big{c} = [y_temp(5,:),y_temp(6,:)];
    
    %Optionally normalize the neuronal data before fitting NL
    if neuron_norm
        mean_y_t_total = mean([y_t_small{c}(:);y_t_big{c}(:)]);
        std_y_t_total = std([y_t_small{c}(:);y_t_big{c}(:)]);
        y_t_small{c} = (y_t_small{c}-mean_y_t_total)./std_y_t_total;
        y_t_big{c} = (y_t_big{c}-mean_y_t_total)./std_y_t_total;
    end
    
    %Get the single lin kernel for every cluster
    temp = load(fullfile(kernel_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')));
    ker = temp.kernel.reverb.main{end};
    
    %Make the linear prediction
    y_hat_lin_t_small{c} = kernelconv(X_fht_small, ker);
    y_hat_lin_t_big{c} = kernelconv(X_fht_big, ker);
    
    
    %Fit a NL for the Small and Large rooms
    lnmodel_small{c} = getlnmodel3(y_hat_lin_t_small{c}, y_t_small{c});
    lnmodel_big{c} = getlnmodel3(y_hat_lin_t_big{c}, y_t_big{c});
    
    %Get the output of the LN models
    y_hat_ln_t_small{c} = lnmodel(lnmodel_small{c}.params, y_hat_lin_t_small{c});
    y_hat_ln_t_big{c} = lnmodel(lnmodel_big{c}.params, y_hat_lin_t_big{c});
    
    
    
    %Optionally, plot the fits to check if they are sensible
    if plot_check
        [x, y] = binplot(y_hat_lin_t_small{c}, y_t_small{c}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.small,'MarkerFaceColor',col.small); hold on;
        [x, y] = binplot(y_hat_lin_t_small{c}, y_hat_ln_t_small{c}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.small, 'LineWidth',lw);
        [x, y] = binplot(y_hat_lin_t_big{c}, y_t_big{c}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.large,'MarkerFaceColor',col.large);
        [x, y] = binplot(y_hat_lin_t_big{c}, y_hat_ln_t_big{c}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.large, 'LineWidth',lw);
        title('LN model');
        xlabel('yhat');
        ylabel('y');
        hold off;
        pause;
    end
end

%% Save the Results
model_ln.small.slope = reach(cell2mat(lnmodel_small),'slope_at_midpoint');
model_ln.big.slope = reach(cell2mat(lnmodel_big),'slope_at_midpoint');
model_ln.small.params = reach(cell2mat(lnmodel_small),'params');
model_ln.small.params = reshape(model_ln.small.params',4,n_clust)';
model_ln.big.params = reach(cell2mat(lnmodel_big),'params');
model_ln.big.params = reshape(model_ln.small.params',4,n_clust)';
var_name = fullfile(save_dir, 'gain_data.mat');
save(var_name, 'y_hat_lin_t_small', 'y_hat_lin_t_big', 'y_hat_ln_t_small', 'y_hat_ln_t_big', 'y_t_small', 'y_t_big', 'model_ln');


%% Plot the individual fits

save_dir_plots = fullfile(save_dir,'Plots_LN');
if ~exist(save_dir_plots, 'dir')
    mkdir(save_dir_plots);
end

rows = 5;
cols = 6;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.05; edgeb = 0.05; space_h = 0.06; space_v = 0.09;
[pos]=subplot_pos(rows,cols,edgel,edger,edgeh,edgeb,space_h,space_v);
clust_per_batch = rows*cols;
axis_sz = 20;

num_batch = ceil(n_clust/clust_per_batch);
curr_n_clust = clust_per_batch;

count = 0;


for b = 1:num_batch
    figure('units','normalized','outerposition',[0 0 1 1]);
    save_name = fullfile(save_dir_plots,['batch_',num2str(b),'_LN_model.svg']);
    
    if b == num_batch
        curr_n_clust = n_clust - (num_batch-1)*clust_per_batch;
    end
    
    for j = 1:curr_n_clust
        count = count+1;
        subplot('position',pos{j});
        [x, y] = binplot(y_hat_lin_t_small{count}, y_t_small{count}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.small,'MarkerFaceColor',col.small); hold on;
        [x, y] = binplot(y_hat_lin_t_small{count}, y_hat_ln_t_small{count}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.small, 'LineWidth',lw);
        [x, y] = binplot(y_hat_lin_t_big{count}, y_t_big{count}, 100); scatter(x(1:end-1),y(1:end-1),'MarkerEdgeColor',col.large,'MarkerFaceColor',col.large);
        [x, y] = binplot(y_hat_lin_t_big{count}, y_hat_ln_t_big{count}, 100); plot(x(1:end-1),y(1:end-1),'Color',col.large, 'LineWidth',lw);
        xlabel('yhat');
        ylabel('y');
        set(gcf,'color','w');
        hold off;
        set(gca,'FontSize',axis_sz,'FontWeight','Normal');
    end
    saveas(gcf,save_name);
    close;
end