function kernel =  gen_norm_model_test(coch,h_max_ms,model,normz,per)

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

dt_ms = coch(1).params.dt_sec*1000; %Convert the dt used for the cochlea into ms
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
fprintf('== Done! This took %0.fs ==\n',toc);

%% Select Y% of the data sequentially
per = per/100;
sz_t = length(small_norm_X_ft);
ix_select = [1:round(sz_t*per)];
%The tensoirzed reverb cochleagrams
chop_small_X_fht = small_X_fht(:,:,ix_select);
chop_med_X_fht = med_X_fht(:,:,ix_select);
chop_big_X_fht = big_X_fht(:,:,ix_select);
%The anechoic cochleagram that is predicted
chop_small_anech_norm_X_ft = small_anech_norm_X_ft(:,ix_select);
chop_med_anech_norm_X_ft = med_anech_norm_X_ft(:,ix_select);
chop_big_anech_norm_X_ft = big_anech_norm_X_ft(:,ix_select);
%% Fit the model for the different rooms

%Small Room
kernel.small = run_model_norm_test(chop_small_X_fht,chop_small_anech_norm_X_ft,model,'small');
% %Medium Room
% kernel.med = run_model_norm_test(chop_med_X_fht,chop_med_anech_norm_X_ft,model,'medium');
% %Big Room
% kernel.big = run_model_norm_test(chop_big_X_fht,chop_big_anech_norm_X_ft,model,'big');

kernel(1).coch_type = coch(1).type;
kernel(1).dt_ms = dt_ms;
kernel(1).h_max_ms = h_max_ms;
kernel(1).n_h = n_h;
kernel(1).model = model;
kernel(1).normz = normz;
kernel(1).freqs = freqs;

end