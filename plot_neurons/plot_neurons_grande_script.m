dir_name{1} = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/perfreq/ridge/10ms/200ms';
NPSP_th = 40;
get_neurons = 'all';
n_dir = length(dir_name);

for d = 1:n_dir
    plot_tconst_neurons_grande(dir_name{d},NPSP_th,get_neurons);
end