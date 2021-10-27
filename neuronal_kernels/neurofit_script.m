%% Neuronal data
time_list_ms = {25, 50, 100, 150, 300, 350};
time_list_string = cellfun(@num2str, time_list_ms, 'UniformOutput', false);
n_cond = length(time_list_ms);
model = 'ridge';
normz = 'perfreq_noneuro';
kfold = 1;
dt_ms = 10;
h_max_ms = 200;
n_cores = 20;

for j = 1:n_cond
    skip_ms = time_list_ms{j};
    save_dir = ['/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/Kernel_fits/kfolds/',num2str(skip_ms/1000),'_4_',num2str(4 + skip_ms/1000),'_8s_version'];
    
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    
    make_neuronal_models_switch(model,h_max_ms,dt_ms,normz,skip_ms,save_dir,kfold,n_cores);
end

