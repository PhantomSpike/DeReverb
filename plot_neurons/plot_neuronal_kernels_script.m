kernel_dir{1} = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model_sim_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
type = 'lnp';
NPSP_th = 40;
ker_per_plot = 10;

for d = 1:length(kernel_dir)
    plot_all_neuronal_kernels(kernel_dir{d},NPSP_th,type,ker_per_plot);
end