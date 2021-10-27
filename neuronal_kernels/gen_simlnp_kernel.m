function kernel =  gen_simlnp_kernel(coch_in , psth, h_max_ms, model, normz, kfold, animal_name)

%% Setup and params
coch = coch_in.coch;

switch animal_name
    
    case {"Ronnie","PLP","Cork","Kilkenny","Derry"}
        small1_X_ft = coch(3).X_ft;
        small2_X_ft = coch(4).X_ft;
        
    case {"Noah","Derekah"}
        small1_X_ft = coch(1).X_ft;
        small2_X_ft = coch(2).X_ft;
end

big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;

dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used  
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram

%% Optionally normalize the data
switch normz
    
    case {'perfreq','perfreq_noneuro'}
        mean_small =  mean([small1_X_ft, small2_X_ft],2); std_small = std([small1_X_ft, small2_X_ft],[],2);
        mean_big = mean([big1_X_ft, big2_X_ft],2); std_big = std([big1_X_ft, big2_X_ft],[],2);
        
        small1_X_ft = (small1_X_ft - mean_small)./std_small; small2_X_ft = (small2_X_ft - mean_small)./std_small;
        big1_X_ft = (big1_X_ft - mean_big)./std_big; big2_X_ft = (big2_X_ft - mean_big)./std_big;
        
    case 'none'
        
end

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small1_X_fht = tensorize(small1_X_ft,n_h);
small2_X_fht = tensorize(small2_X_ft,n_h);
big1_X_fht = tensorize(big1_X_ft,n_h);
big2_X_fht = tensorize(big2_X_ft,n_h);

%Combine the stimuli together
small_X_fht = cat(3,small1_X_fht,small2_X_fht);
big_X_fht = cat(3,big1_X_fht,big2_X_fht);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Make histograms for this cluster for all 3 conditions

y_t_small = psth.y_t_small_mean;
y_t_big = psth.y_t_big_mean;

%% Fit the model for the different rooms

%Get a single k_f for all conditions
if strcmp(model,'sep_kh_all')
    all_X_fht = cat(3,anech_X_fht,small_X_fht,big_X_fht); %Put all cochleagrams together
    y_t_all = [y_t_anech,y_t_small,y_t_big]; %Put all neuronal psths together
    kernel_temp = run_model(all_X_fht,y_t_all,'sep');
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
        y_t_big = (y_t_big - mean(y_t_big))./std(y_t_big);
        
    case 'perfreq_noneuro'
        
    case 'none'
    
end

if kfold
    %Small Room
    kernel.small = run_model_kfold(small_X_fht, y_t_small, model);
    %Big Room
    kernel.big = run_model_kfold(big_X_fht, y_t_big, model);

else
    %Small Room
    kernel.small = run_model(small_X_fht, y_t_small, model);
    %Big Room
    kernel.big = run_model(big_X_fht, y_t_big, model);
end

kernel(1).coch_type = coch(1).type;
kernel(1).dt_ms = dt_ms;
kernel(1).h_max_ms = h_max_ms;
kernel(1).n_h = n_h;
kernel(1).model = model;
kernel(1).freqs = freqs;

end