count = 0;
%% Model 1
coch_file = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_new_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_new_400_19k_specpower_kfolds';
model = 'lasso';
normz = 'none';
kfolds = 1;
h_max_ms = 200;
n_cores = 15;

try
    make_model(coch_file,h_max_ms,model,normz,n_cores,save_dir,kfolds);
catch
    count = count + 1;
    [~,fname] = fileparts(coch_file);
    model_fail{count} = [num2str(model),'_',num2str(h_max_ms),'ms_',fname];
end

%% Model 2
coch_file = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_new_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_new_400_19k_specpower_kfolds';
model = 'ridge';
normz = 'none';
kfolds = 1;
h_max_ms = 200;
n_cores = 15;

try
    make_model(coch_file,h_max_ms,model,normz,n_cores,save_dir,kfolds);
catch
    count = count + 1;
    [~,fname] = fileparts(coch_file);
    model_fail{count} = [num2str(model),'_',num2str(h_max_ms),'ms_',fname];
end

%% Model 3
coch_file = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_new_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_new_400_19k_specpower_kfolds';
model = 'lasso';
normz = 'perfreq';
kfolds = 1;
h_max_ms = 200;
n_cores = 15;

try
    make_model(coch_file,h_max_ms,model,normz,n_cores,save_dir,kfolds);
catch
    count = count + 1;
    [~,fname] = fileparts(coch_file);
    model_fail{count} = [num2str(model),'_',num2str(h_max_ms),'ms_',fname];
end

%% Model 4
coch_file = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s_exp_nojap_new_400_19k_specpower/closed/single/chopped/10ms/specpower/coch_all_conditions_closed_singlepos_specpower.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_new_400_19k_specpower_kfolds';
model = 'ridge';
normz = 'perfreq';
kfolds = 1;
h_max_ms = 200;
n_cores = 15;

try
    make_model(coch_file,h_max_ms,model,normz,n_cores,save_dir,kfolds);
catch
    count = count + 1;
    [~,fname] = fileparts(coch_file);
    model_fail{count} = [num2str(model),'_',num2str(h_max_ms),'ms_',fname];
end