dir_name = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/600s_exp_nojap_400_19k_specpower_kfolds/closed/single/chopped/perfreq/ridge/10ms/200ms'; %Name of the parent directory 
cv = 'best';
files = dir([dir_name,'/**/*kernel.mat']);
kernel_dir = {files.folder};
n_dir = length(kernel_dir);

for d = 1:n_dir
    switch cv
        case 'all'
            plot_model_tconst_grande_cv_all(kernel_dir{d});
        case 'ind'
            plot_model_tconst_grande_cv(kernel_dir{d});
        case 'best'
            plot_model_tconst_grande_cv_best(kernel_dir{d});
        case 'none'
            plot_model_tconst_grande(kernel_dir{d});
    end
end