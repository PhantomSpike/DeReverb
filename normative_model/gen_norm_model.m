function kernel =  gen_norm_model(coch,h_max_ms,model,normz,kfolds)
norm2 = false; %Whether to do the the normalization with fitting and then throwing away the last bin
%% Load the cochleagrams
switch coch(1).pos
    case 'single'
        small_X_ft = coch(2).X_ft;
        med_X_ft = coch(3).X_ft;
        big_X_ft = coch(4).X_ft;
        small_anech_X_ft = coch(1).X_ft;
        med_anech_X_ft = coch(1).X_ft;
        big_anech_X_ft = coch(1).X_ft;
        
    case 'many'
        small_X_ft = coch(1).reverb.X_ft;
        med_X_ft = coch(2).reverb.X_ft;
        big_X_ft = coch(3).reverb.X_ft;
        small_anech_X_ft = coch(1).anechoic.X_ft;
        med_anech_X_ft = coch(2).anechoic.X_ft;
        big_anech_X_ft = coch(3).anechoic.X_ft;
end

dt_s = coch(1).params.dt_sec; %Get the cochleagram time bin in s
dt_ms = dt_s*1000; %Convert the dt used for the cochlea into ms
n_h = round(h_max_ms/dt_ms); %Find the number of history steps necessary given the bin size and max history
freqs = coch(1).params.freqs; %Get the freqeuncies that were used  
%% Normalize the cochleagrams before fitting the model
%This could be no normalization, global z-scoring or per freqeuncy
%z-scoring

fprintf('== Normalizing the cochleagrams ==\n');tic;
if strcmp(normz,'none')
    small_norm_X_ft = small_X_ft; small_anech_norm_X_ft = small_anech_X_ft;
    med_norm_X_ft = med_X_ft; med_anech_norm_X_ft = med_anech_X_ft;
    big_norm_X_ft = big_X_ft; big_anech_norm_X_ft = big_anech_X_ft;
else
    [small_norm_X_ft,small_anech_norm_X_ft] = norm_data(small_X_ft,small_anech_X_ft,normz);
    [med_norm_X_ft,med_anech_norm_X_ft] = norm_data(med_X_ft,med_anech_X_ft,normz);
    [big_norm_X_ft,big_anech_norm_X_ft] = norm_data(big_X_ft,big_anech_X_ft,normz);
end

fprintf('== Done! This took %0.fs ==\n',toc);
%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small_X_fht = tensorize(small_norm_X_ft,n_h);
med_X_fht = tensorize(med_norm_X_ft,n_h);
big_X_fht = tensorize(big_norm_X_ft,n_h);

if norm2
    mn_small = mean(small_X_fht(:,h_chop:end,:),2);
    mn_med = mean(med_X_fht(:,h_chop:end,:),2);
    mn_big = mean(big_X_fht(:,h_chop:end,:),2);
    small_X_fht = cat(2,small_X_fht(:,1:h_chop-1,:),mn_small);
    med_X_fht = cat(2,med_X_fht(:,1:h_chop-1,:),mn_med);
    big_X_fht = cat(2,big_X_fht(:,1:h_chop-1,:),mn_big);
end
    
fprintf('== Done! This took %0.fs ==\n',toc);
%% Fit the model for the different rooms

if kfolds
    %Small Room
    kernel.small = run_model_normcv(small_X_fht,small_anech_norm_X_ft,model,'small');
    %Medium Room
    kernel.med = run_model_normcv(med_X_fht,med_anech_norm_X_ft,model,'medium');
    %Big Room
    kernel.big = run_model_normcv(big_X_fht,big_anech_norm_X_ft,model,'big');
else
    %Small Room
    kernel.small = run_model_norm(small_X_fht,small_anech_norm_X_ft,model,'small');
    %Medium Room
    kernel.med = run_model_norm(med_X_fht,med_anech_norm_X_ft,model,'medium');
    %Big Room
    kernel.big = run_model_norm(big_X_fht,big_anech_norm_X_ft,model,'big');
end

kernel(1).coch_type = coch(1).type;
kernel(1).dt_ms = dt_ms;
kernel(1).h_max_ms = h_max_ms;
kernel(1).n_h = n_h;
if isnumeric(model)
    model_name = 'manual_alpha';
else
    model_name = model;
end
kernel(1).model = model_name;
kernel(1).normz = normz;
kernel(1).freqs = freqs;

end