function kernel =  gen_neuronal_kernel2(coch_in,y,h_max_ms,model,normz,kfold)

%% Setup and params
coch = coch_in.coch;
small1_X_ft = coch(1).X_ft;
small2_X_ft = coch(2).X_ft;
med1_X_ft = coch(3).X_ft;
med2_X_ft = coch(4).X_ft;
big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;

dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram

%% Optionally normalize the data
switch normz
    case {'perfreq','perfreq_noneuro'}
        mean_small = mean([small1_X_ft, small2_X_ft],2); std_small = std([small1_X_ft, small2_X_ft],[],2);
        mean_med = mean([med1_X_ft, med2_X_ft],2); std_med = std([med1_X_ft, med2_X_ft],[],2);
        mean_big = mean([big1_X_ft, big2_X_ft],2); std_big = std([big1_X_ft, big2_X_ft],[],2);
        small1_X_ft = (small1_X_ft - mean_small)./std_small; small2_X_ft = (small2_X_ft - mean_small)./std_small;
        med1_X_ft = (med1_X_ft - mean_med)./std_med; med2_X_ft = (med2_X_ft - mean_med)./std_med;
        big1_X_ft = (big1_X_ft - mean_big)./std_big; big2_X_ft = (big2_X_ft - mean_big)./std_big;
    case 'none'
end

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small1_X_fht = tensorize(small1_X_ft,n_h);
small2_X_fht = tensorize(small2_X_ft,n_h);
med1_X_fht = tensorize(med1_X_ft,n_h);
med2_X_fht = tensorize(med2_X_ft,n_h);
big1_X_fht = tensorize(big1_X_ft,n_h);
big2_X_fht = tensorize(big2_X_ft,n_h);

%Combine the two stimuli together
small_X_fht = cat(3, small1_X_fht, small2_X_fht);
med_X_fht = cat(3, med1_X_fht, med2_X_fht);
big_X_fht = cat(3, big1_X_fht, big2_X_fht);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Make histograms for this cluster for all 3 conditions
n_stim = length(y.stim);
n_rep = length(y.stim(1).repeat);

for s = 1:n_stim
    psth_temp = [];
    for r = 1:n_rep
        psth_temp(r,:) = histc(y.stim(s).repeat(r).spiketimes,t_edges_s);
    end
    y_temp(s,:) = mean(psth_temp);
end

y_t_small = [y_temp(1,:), y_temp(2,:)];
y_t_med = [y_temp(3,:), y_temp(4,:)];
y_t_big = [y_temp(5,:), y_temp(6,:)];

%% Fit the model for the different rooms

%Get a single k_f for all conditions
if strcmp(model,'sep_kh_all')
    all_X_fht = cat(3, small_X_fht, med_X_fht, big_X_fht); %Put all cochleagrams together
    y_t_all = [y_t_small, y_t_med, y_t_big]; %Put all neuronal psths together
    kernel_temp = run_model(all_X_fht, y_t_all, 'sep');
    k_f = kernel_temp.k_f;
elseif strcmp(model,'sep_kh')
    all_X_fht = cat(3,small_X_fht,big_X_fht); %Put all cochleagrams together
    y_t_all = [y_t_small,y_t_big]; %Put all neuronal psths together
    kernel_temp = run_model(all_X_fht,y_t_all,'sep');
    k_f = kernel_temp.k_f;
else
    k_f = [];
end

%Normalize neuronal firing too
switch normz
    case 'perfreq' 
        y_t_small = (y_t_small - mean(y_t_small))./std(y_t_small);
        y_t_med = (y_t_med - mean(y_t_med))./std(y_t_med);
        y_t_big = (y_t_big - mean(y_t_big))./std(y_t_big);
        
    case 'perfreq_noneuro'
        
    case 'none'
        
end

if kfold
    %Small Room
    kernel.small = run_model_kfold(small_X_fht, y_t_small, model);
    %Medium Room
    kernel.med = run_model_kfold(med_X_fht, y_t_med, model);
    %Big Room
    kernel.big = run_model_kfold(big_X_fht, y_t_big, model);
else
    %Small Room
    kernel.small = run_model(small_X_fht, y_t_small, model, k_f);
    %Medium Room
    kernel.med = run_model(med_X_fht, y_t_med, model, k_f);
    %Big Room
    kernel.big = run_model(big_X_fht, y_t_big, model, k_f);
end

kernel(1).coch_type = coch(1).type;
kernel(1).dt_ms = dt_ms;
kernel(1).h_max_ms = h_max_ms;
kernel(1).n_h = n_h;
kernel(1).model = model;
kernel(1).freqs = freqs;

end