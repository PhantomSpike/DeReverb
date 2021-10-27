function kh_fits = fit_tconst_neuron(kernel_name)
%% Params
r_type{1} = 'anech';
r_type{2} = 'small';
r_type{3} = 'big';
alpha_init_ms = 20;
ub_A_B_all = [25];
beta_init_ms = [];
interpolate = 0;
upsamp_factor = 100;
%% Load the kernel fits
load(kernel_name,'kernel');
n_rooms = length(r_type);
%% Fit a dexp model (beta) for every room
for u = 1:length(ub_A_B_all)
    ub_A_B = ub_A_B_all(u);
    dt_ms = round(kernel.dt_ms);
    input_params.time_bin_ms = dt_ms;
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
    
    fprintf('== Generating dexp fits ==\n');tic;
    for r = 1:n_rooms
        fprintf('== Working on %s room ==\n',r_type{r});
        k_h = flipud(kernel.(r_type{r}).k_h); %Get the k_h from the neuronal kernel fits
        k_f = kernel.(r_type{r}).k_f;
        k_h_full  = k_h;
        [~,ix_max] = max(k_h); %Find the index of the maximum value
        k_h = k_h(ix_max:end); %Take only the values from the max onwards
        k_h = k_h./max(abs(k_h(:))); %Scale the k_h between [-1;1] so it doesn't have super small values
        input_params.A_init = max(k_h(:)); %Set the A_init to the maximum of the kernel
        input_params.B_init = abs(min(k_h(:))); %Set the B_init to the minimum of the kernel
        dexp_fit = fitdexp(k_h,input_params);
        dexp_fit.alpha = dexp_fit.alpha + (ix_max-1)*dt_ms; %Add the number of time points from the start to the max to compensate for the removal
        dexp_fit.beta = dexp_fit.beta + (ix_max-1)*dt_ms; %Add the number of time points from the start to the max to compensate for the removal
        kh_fits.(r_type{r}) = dexp_fit;
        kh_fits.(r_type{r}).k_h_full = k_h_full; %Save k_h before truncation
        kh_fits.(r_type{r}).k_f = k_f; %Save the k_f for later
    end
    fprintf('== Done! This took %0.fs ==\n',toc); 
end
kh_fits.params = kernel.params;
kh_fits.freqs = kernel.freqs;
kh_fits.dt_ms = dt_ms;
%% Find the center of mass (tau) for every room
n_h = kernel.n_h;
h = (0:1:n_h-1)';
h = dt_ms*h;
for r = 1:n_rooms
    k_h = flipud(kernel.(r_type{r}).k_h); %Get the k_h from the norm model
    k_h = min(k_h,0); %Remove all positive values
    k_h = abs(k_h); %Make all inhibitory values positive
    k_h = k_h./sum(k_h(:)); %Scale the values to sum to 1
    kh_fits.(r_type{r}).tau = (k_h'*h); %Compute a weighted sum of all values
end