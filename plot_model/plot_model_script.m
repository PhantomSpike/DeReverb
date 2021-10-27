%% Define models
dir_name = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms'; %Name of the parent directory 
cv = 1;
files = dir([dir_name,'/**/*kernel.mat']);
model_dir = {files.folder};
%% Run the models
n_model = length(model_dir);

for m = 1:n_model
   plot_norm_model(model_dir{m},cv);
end

