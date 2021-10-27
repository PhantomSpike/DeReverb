function CCnorm = get_CCnorm(coch_in, y_in, kernel, normz)

%% Setup and params
coch = coch_in.coch;
X_ft(1).small = coch(1).X_ft;
X_ft(2).small = coch(2).X_ft;
X_ft(3).small = coch(3).X_ft;
X_ft(4).small = coch(4).X_ft;
X_ft(1).large = coch(5).X_ft;
X_ft(2).large = coch(6).X_ft;
X_ft(3).large = coch(7).X_ft;
X_ft(4).large = coch(8).X_ft;

n_reps = length(y_in.stim(1).repeat);

y(1).small = y_in.stim(1).repeat;
y(2).small = y_in.stim(2).repeat;
y(3).small = y_in.stim(3).repeat;
y(4).small = y_in.stim(4).repeat;
y(1).large = y_in.stim(5).repeat;
y(2).large = y_in.stim(6).repeat;
y(3).large = y_in.stim(7).repeat;
y(4).large = y_in.stim(8).repeat;

n_h = kernel.n_h; %Find the number of history steps necessary given the bin size and max history 
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram
n_stim = length(coch)/2; %Get the number of different stimuli
type_start{1} = 'small'; type_start{2} = 'large';
n_types = length(type_start);
n_pres = 2; %How many times each condition (sl1, sl2, ls1, ls2) appears in each presentation
n_total = n_types*n_stim*n_pres; %Number of total occurences for each condition (sl1, sl2, ls1, ls2)

%% Define the periods of the 4 stimuli

%sl1 - small-to-large transition first 4s
%sl2 - small-to-large transition last 4s
%ls1 - large-to-small transition first 4s
%ls2 - large-to-small transition last 4s

% This is a 4s version
%Small start (R1)
% sl1(1).small = [8, 12]; sl1(2).small = [24, 28];
% sl2(1).small = [12, 16]; sl2(2).small = [28, 32];
% ls1(1).small = [16, 20]; ls1(2).small = [32, 36];
% ls2(1).small = [20, 24]; ls2(2).small = [36, 40];
% 
% %Large start (R2)
% sl1(1).large = [16, 20]; sl1(2).large = [32, 36];
% sl2(1).large = [20, 24]; sl2(2).large = [36, 40];
% ls1(1).large = [8, 12]; ls1(2).large = [24, 28];
% ls2(1).large = [12, 16]; ls2(2).large = [28, 32];

%This is a 3s version
%Small start (R1)
% sl1(1).small = [8, 11]; sl1(2).small = [24, 27];
% sl2(1).small = [12, 15]; sl2(2).small = [28, 31];
% ls1(1).small = [16, 19]; ls1(2).small = [32, 35];
% ls2(1).small = [20, 23]; ls2(2).small = [36, 39];
% 
% %Large start (R2)
% sl1(1).large = [16, 19]; sl1(2).large = [32, 35];
% sl2(1).large = [20, 23]; sl2(2).large = [36, 39];
% ls1(1).large = [8, 11]; ls1(2).large = [24, 27];
% ls2(1).large = [12, 15]; ls2(2).large = [28, 31];

%This is a 2s version
%Small start (R1)
% sl1(1).small = [8, 10]; sl1(2).small = [24, 26];
% sl2(1).small = [12, 14]; sl2(2).small = [28, 30];
% ls1(1).small = [16, 18]; ls1(2).small = [32, 34];
% ls2(1).small = [20, 22]; ls2(2).small = [36, 38];
% 
% %Large start (R2)
% sl1(1).large = [16, 18]; sl1(2).large = [32, 34];
% sl2(1).large = [20, 22]; sl2(2).large = [36, 38];
% ls1(1).large = [8, 10]; ls1(2).large = [24, 26];
% ls2(1).large = [12, 14]; ls2(2).large = [28, 30];

%This is 1-4s and 5-8s version
%Small start (R1)
% sl1(1).small = [9, 12]; sl1(2).small = [25, 28];
% sl2(1).small = [13, 16]; sl2(2).small = [29, 32];
% ls1(1).small = [17, 20]; ls1(2).small = [33, 36];
% ls2(1).small = [21, 24]; ls2(2).small = [37, 40];
% 
% %Large start (R2)
% sl1(1).large = [17, 20]; sl1(2).large = [33, 36];
% sl2(1).large = [21, 24]; sl2(2).large = [37, 40];
% ls1(1).large = [9, 12]; ls1(2).large = [25, 28];
% ls2(1).large = [13, 16]; ls2(2).large = [29, 32];

%This is 0.5-4s and 4.5-8s version
%Small start (R1)
sl1(1).small = [8.5, 12]; sl1(2).small = [24.5, 28];
sl2(1).small = [12.5, 16]; sl2(2).small = [28.5, 32];
ls1(1).small = [16.5, 20]; ls1(2).small = [32.5, 36];
ls2(1).small = [20.5, 24]; ls2(2).small = [36.5, 40];

%Large start (R2)
sl1(1).large = [16.5, 20]; sl1(2).large = [32.5, 36];
sl2(1).large = [20.5, 24]; sl2(2).large = [36.5, 40];
ls1(1).large = [8.5, 12]; ls1(2).large = [24.5, 28];
ls2(1).large = [12.5, 16]; ls2(2).large = [28.5, 32];

%% Colect the four conditions where they belong
count = 0; %Counter variable

for t = 1:n_types
    type = type_start{t};
    
    for p = 1:n_pres
        %Define the indices
        ix_sl1 = t_edges_s>=sl1(p).(type)(1) & t_edges_s<sl1(p).(type)(2);
        ix_sl2 = t_edges_s>=sl2(p).(type)(1) & t_edges_s<sl2(p).(type)(2);
        ix_ls1 = t_edges_s>=ls1(p).(type)(1) & t_edges_s<ls1(p).(type)(2);
        ix_ls2 = t_edges_s>=ls2(p).(type)(1) & t_edges_s<ls2(p).(type)(2);
        
        %Get the corresponding time edges for the psths
        t_edges_s_sl1 = t_edges_s(ix_sl1);
        t_edges_s_sl2 = t_edges_s(ix_sl2);
        t_edges_s_ls1 = t_edges_s(ix_ls1);
        t_edges_s_ls2 = t_edges_s(ix_ls2);
        
        for s = 1:n_stim
            count = count+1;
            %Get the right data from the cochleagrams
            X_ft_sl1{1,count} = X_ft(s).(type)(:,ix_sl1);
            X_ft_sl2{1,count} = X_ft(s).(type)(:,ix_sl2);
            X_ft_ls1{1,count} = X_ft(s).(type)(:,ix_ls1);
            X_ft_ls2{1,count} = X_ft(s).(type)(:,ix_ls2);
            
            psth_temp_sl1 = [];
            psth_temp_sl2 = [];
            psth_temp_ls1 = [];
            psth_temp_ls2 = [];
            
            %Get the psths of the neurons
            for  r = 1:n_reps
                psth_temp_sl1(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_sl1);
                psth_temp_sl2(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_sl2);
                psth_temp_ls1(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_ls1);
                psth_temp_ls2(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_ls2);
            end
            
            %Average across reps
            y_sl1{1,count} = psth_temp_sl1;
            y_sl2{1,count} = psth_temp_sl2;
            y_ls1{1,count} = psth_temp_ls1;
            y_ls2{1,count} = psth_temp_ls2;
            
        end
    end
end
%% Optionally normalize the data
switch normz
    
    case {'perfreq','perfreq_noneuro'}
        %Calculate the mean and std for each frequency across all stimuli of the
        %same kind
        X_ft_sl1_all = cell2mat(X_ft_sl1); X_ft_sl2_all = cell2mat(X_ft_sl2); 
        X_ft_ls1_all = cell2mat(X_ft_ls1); X_ft_ls2_all = cell2mat(X_ft_ls2);
        
        mean_sl1 = mean(X_ft_sl1_all, 2); mean_sl2 = mean(X_ft_sl2_all, 2);
        mean_ls1 = mean(X_ft_ls1_all, 2); mean_ls2 = mean(X_ft_ls2_all, 2);
        
        std_sl1 = std(X_ft_sl1_all,[],2); std_sl2 = std(X_ft_sl2_all,[],2);
        std_ls1 = std(X_ft_ls1_all,[],2); std_ls2 = std(X_ft_ls2_all,[],2);

        for j = 1:n_total
            
            X_ft_sl1{j} = (X_ft_sl1{j} - mean_sl1)./std_sl1;
            X_ft_sl2{j} = (X_ft_sl2{j} - mean_sl2)./std_sl2;
            X_ft_ls1{j} = (X_ft_ls1{j} - mean_ls1)./std_ls1;
            X_ft_ls2{j} = (X_ft_ls2{j} - mean_ls2)./std_ls2;
        end
        
    case 'none'
end

%% Tensorize and remove history that includes the previous condition
fprintf('== Tensorizing the cochleagrams ==\n');tic;

for j = 1:n_total
    
    %Tensorize the cochleagram
    X_fht_sl1{j,1} = tensorize(X_ft_sl1{j}, n_h);
    X_fht_sl2{j,1} = tensorize(X_ft_sl2{j}, n_h);
    X_fht_ls1{j,1} = tensorize(X_ft_ls1{j}, n_h);
    X_fht_ls2{j,1} = tensorize(X_ft_ls2{j}, n_h);
    
    %Remove the parts that include history from previous condition
    X_fht_sl1{j}(:,:,1:n_h-1) = [];
    X_fht_sl2{j}(:,:,1:n_h-1) = [];
    X_fht_ls1{j}(:,:,1:n_h-1) = [];
    X_fht_ls2{j}(:,:,1:n_h-1) = [];
    
    %Remove the same parts from the psths
    y_sl1{j}(:,1:n_h-1) = [];
    y_sl2{j}(:,1:n_h-1) = [];
    y_ls1{j}(:,1:n_h-1) = [];
    y_ls2{j}(:,1:n_h-1) = [];
    
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Combine
%Combine all stimuli of the same kind together for GLM fitting along the
%time dimension (actual order doesn't matter as long as neuronal data is the same)

X_fht_sl1_all = cat(3, X_fht_sl1{:});
X_fht_sl2_all = cat(3, X_fht_sl2{:});
X_fht_ls1_all = cat(3, X_fht_ls1{:});
X_fht_ls2_all = cat(3, X_fht_ls2{:});

%Combine the neuronal data too
y_nt_sl1 = cat(2, y_sl1{:});
y_nt_sl2 = cat(2, y_sl2{:});
y_nt_ls1 = cat(2, y_ls1{:});
y_nt_ls2 = cat(2, y_ls2{:});

%Optionally Normalize neuronal firing too
switch normz
    case 'perfreq'
        y_nt_sl1 = (y_nt_sl1 - mean(y_nt_sl1))./std(y_nt_sl1);
        y_nt_sl2 = (y_nt_sl2 - mean(y_nt_sl2))./std(y_nt_sl2);
        y_nt_ls1 = (y_nt_ls1 - mean(y_nt_ls1))./std(y_nt_ls1);
        y_nt_ls2 = (y_nt_ls2 - mean(y_nt_ls2))./std(y_nt_ls2);
        
    case 'perfreq_noneuro'
        
    case 'none'
        
end


%% Compute CCnorm for the different folds and conditions
kernel_sl1 = kernel.sl1.allKernels;
kernel_sl2 = kernel.sl2.allKernels;
kernel_ls1 = kernel.ls1.allKernels;
kernel_ls2 = kernel.ls2.allKernels;

%SL1
CCnorm.sl1 = kfold_CCnorm(X_fht_sl1_all, y_nt_sl1, kernel_sl1);
%SL2
CCnorm.sl2 = kfold_CCnorm(X_fht_sl2_all, y_nt_sl2, kernel_sl2);
%LS1
CCnorm.ls1 = kfold_CCnorm(X_fht_ls1_all, y_nt_ls1, kernel_ls1);
%LS2
CCnorm.ls2 = kfold_CCnorm(X_fht_ls2_all, y_nt_ls2, kernel_ls2);

CCnorm.mean = mean([CCnorm.sl1.CCnorm_mean, CCnorm.sl2.CCnorm_mean, CCnorm.ls1.CCnorm_mean, CCnorm.ls2.CCnorm_mean]);

end