function plot_norm_model(model_dir,cv)
%Here I should make a bigger function and an inner function that plots just
%one set of kernels
%% Define params
save_dir = fullfile(model_dir,'/Plots');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
%% Load the data
kernel_list = dir(fullfile(model_dir,'*kernel.mat'));
kernel_list = kernel_list(~ismember({kernel_list.name},{'.','..'}));
fname = fullfile(kernel_list(1).folder,kernel_list(1).name);
load(fname,'kernels');
params.freqs = kernels.freqs;
params.dt_ms = kernels.dt_ms;
params.n_h = kernels.n_h;
params.model = kernels.model; 
params.h_max_ms = kernels.h_max_ms;
params.coch_type = kernels.coch_type;
params.normz = kernels.normz;
params.save_dir = save_dir;
%% Plot the kernels for a given room

if cv
    plot_ind_room_cv(kernels.small,params,'small');
    plot_ind_room_cv(kernels.med,params,'med');
    plot_ind_room_cv(kernels.big,params,'big');
else
    plot_ind_room(kernels.small,params,'small');
    plot_ind_room(kernels.med,params,'med');
    plot_ind_room(kernels.big,params,'big');
end

end
