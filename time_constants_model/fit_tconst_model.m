function fit_tconst_model(kernel_name)
%% Params
r_type{1} = 'small';
r_type{2} = 'med';
r_type{3} = 'big';
alpha_init_ms = 20;
ub_A_B_all = [5,25];
beta_init_ms = [];
interpolate = 0;
upsamp_factor = 100;
%% Load the kernel fits
save_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/test';
[a,model_name] = fileparts(kernel_name);
[~,room_gen_name] = fileparts(a);
model_name = strrep(model_name,'_',' ');
room_gen_name = strrep(room_gen_name,'_',' ');
load(kernel_name);
n_ker = length(kernels.small); %The number of kernels
n_rooms = length(r_type);
%% Loop through the all the f-bands (kernels) for every room
for u = 1:length(ub_A_B_all)
    ub_A_B = ub_A_B_all(u);
    input_params.time_bin_ms = kernels.dt_ms;
    input_params.alpha_init = alpha_init_ms;
    input_params.beta_init = beta_init_ms;
    input_params.ub_A_B = ub_A_B;
    input_params.interpolate = interpolate;
    if interpolate
        input_params.upsamp_factor = upsamp_factor;
        int_name = ['_inter_',num2str(upsamp_factor),'fold'];
    else
        int_name = [];
    end
    
    fprintf('== Generating dexp fits for the normative models ==\n');tic;
    for r = 1:n_rooms
        fprintf('== Working on %s room ==\n',r_type{r});
        for k = 1:n_ker
            fprintf('== kernel %0.f/%0.f ==\n',k,n_ker);tic
            k_h = flipud(kernels.(r_type{r}){k}.k_h); %Get the k_h from the norm model
            k_h = k_h./max(abs(k_h(:))); %Scale the k_h between [-1;1] so it doesn't have super small values
            input_params.A_init = max(k_h(:)); %Set the A_init to the maximum of the kernel
            input_params.B_init = abs(min(k_h(:))); %Set the B_init to the minimum of the kernel
            dexp_fit = fitdexp(k_h,input_params);
            kh_fits.(r_type{r}){k} = dexp_fit;
            fprintf('== Done! This took %.0fs ==\n',toc);
        end
    end
    
    save_dir_full = fullfile(save_dir,[room_gen_name,'_',model_name,'_ub_A_B_',num2str(ub_A_B),int_name]);
    if ~exist(save_dir_full, 'dir')
        mkdir(save_dir_full)
    end
    var_name = fullfile(save_dir_full,'kh_fits.mat');
    save(var_name,'kh_fits');
    %Plot the results
    params.r_type = r_type;
    params.save_dir_full = save_dir_full;
    params.room_gen_name = room_gen_name;
    params.model_name = model_name;
    params.int_name = int_name;
    params.freqs = kernels.freqs;
    params.n_rooms = n_rooms;
    params.n_ker = n_ker;
    params.ub_A_B = ub_A_B;
    plot_tconst(kh_fits,params)   
end