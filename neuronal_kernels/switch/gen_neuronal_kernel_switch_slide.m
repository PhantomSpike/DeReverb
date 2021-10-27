function kernel =  gen_neuronal_kernel_switch_slide(coch_in,y_in,params)

%% Unroll the params
h_max_ms = params.h_max_ms;
model = params.model;
normz = params.normz;
kfold = params.kfold;
sl = params.sl;
ls = params.ls;

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


dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used  
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram
n_stim = length(coch)/2; %Get the number of different stimuli
type_start{1} = 'small'; type_start{2} = 'large';
n_types = length(type_start);
n_pres = 2; %How many times each condition (sl, ls) appears in each presentation
n_total = n_types*n_stim*n_pres; %Number of total occurences for each condition (sl1, sl2, ls1, ls2)



%% Colect the four conditions where they belong
count = 0; %Counter variable

for t = 1:n_types
    type = type_start{t};
    
    for p = 1:n_pres
        %Define the indices
        ix_sl = t_edges_s>=sl(p).(type)(1) & t_edges_s<sl(p).(type)(2);
        ix_ls = t_edges_s>=ls(p).(type)(1) & t_edges_s<ls(p).(type)(2);

        
        %Get the corresponding time edges for the psths
        t_edges_s_sl = t_edges_s(ix_sl);
        t_edges_s_ls = t_edges_s(ix_ls);
        
        for s = 1:n_stim
            count = count+1;
            %Get the right data from the cochleagrams
            X_ft_sl{1,count} = X_ft(s).(type)(:,ix_sl);
            X_ft_ls{1,count} = X_ft(s).(type)(:,ix_ls);

            
            psth_temp_sl = [];
            psth_temp_ls = [];
  
            
            %Get the psths of the neurons
            for  r = 1:n_reps
                psth_temp_sl(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_sl);
                psth_temp_ls(r,:) = histc(y(s).(type)(r).spiketimes, t_edges_s_ls);
            end
            
            %Average across reps
            y_sl{1,count} = mean(psth_temp_sl);
            y_ls{1,count} = mean(psth_temp_ls); 
            
        end
    end
end

%% Optionally normalize the data
switch normz
    
    case {'perfreq','perfreq_noneuro'}
        %Calculate the mean and std for each frequency across all stimuli of the
        %same kind
        X_ft_sl_all = cell2mat(X_ft_sl); 
        X_ft_ls_all = cell2mat(X_ft_ls); 
        
        mean_sl = mean(X_ft_sl_all, 2);
        mean_ls = mean(X_ft_ls_all, 2); 
        
        std_sl = std(X_ft_sl_all,[],2); 
        std_ls = std(X_ft_ls_all,[],2); 

        for j = 1:n_total
            
            X_ft_sl{j} = (X_ft_sl{j} - mean_sl)./std_sl;
            X_ft_ls{j} = (X_ft_ls{j} - mean_ls)./std_ls;

        end
        
    case 'none'
end

%% Tensorize and remove history that includes the previous condition
fprintf('== Tensorizing the cochleagrams ==\n');tic;

for j = 1:n_total
    
    %Tensorize the cochleagram
    X_fht_sl{j,1} = tensorize(X_ft_sl{j}, n_h);
    X_fht_ls{j,1} = tensorize(X_ft_ls{j}, n_h);
    
    %Remove the parts that include history from previous condition
    X_fht_sl{j}(:,:,1:n_h-1) = [];
    X_fht_ls{j}(:,:,1:n_h-1) = [];
    
    %Remove the same parts from the psths
    y_sl{j}(1:n_h-1) = [];
    y_ls{j}(1:n_h-1) = [];
    
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Combine
%Combine all stimuli of the same kind together for GLM fitting along the
%time dimension (actual order doesn't matter as long as neuronal data is the same)

X_fht_sl_all = cat(3, X_fht_sl{:});
X_fht_ls_all = cat(3, X_fht_ls{:});

%Combine the neuronal data too
y_t_sl = cat(2, y_sl{:});
y_t_ls = cat(2, y_ls{:});

%% Fit the model for the different rooms

%Normalize neuronal firing too
switch normz
    case 'perfreq'
        y_t_sl = (y_t_sl - mean(y_t_sl))./std(y_t_sl);
        y_t_ls = (y_t_ls - mean(y_t_ls))./std(y_t_ls);
        
    case 'perfreq_noneuro'
        
    case 'none'
        
end

if kfold
    %SL
    kernel.sl = run_model_kfold(X_fht_sl_all, y_t_sl, model);
    %LS
    kernel.ls = run_model_kfold(X_fht_ls_all, y_t_ls, model);
    
else
    %SL1
    kernel.sl = run_model(X_fht_sl_all, y_t_sl, model, k_f);
    %LS1
    kernel.ls = run_model(X_fht_ls_all, y_t_ls, model, k_f);
    
end

kernel(1).coch_type = coch(1).type;
kernel(1).dt_ms = dt_ms;
kernel(1).h_max_ms = h_max_ms;
kernel(1).n_h = n_h;
kernel(1).model = model;
kernel(1).freqs = freqs;
end