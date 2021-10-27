%With this function we get the center of mass of the inhibition 
dir_name = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits/Current/900s_exp_400_19k_final'; %Name of the parent directory 
model_kernels = dir([dir_name,'/**/*kernel.mat']);


r_type{1} = 'big';
r_type{2} = 'med';
r_type{3} = 'small';
n_rooms = length(r_type);

for m = 1:length(model_kernels)
    load(fullfile(model_kernels(m).folder,model_kernels(m).name));
    model = kernels.model;
    dt_ms = kernels.dt_ms;
    freqs = kernels.freqs;
    n_ker = length(freqs);
    n_h = kernels.n_h;
    h = (0:1:n_h-1)'; 
    h = dt_ms*h;
    for r = 1:n_rooms
        for k = 1:n_ker
            switch model
                case {'lasso','ridge','elastic'}
                    k_fh = fliplr(kernels.(r_type{r}){k}.main.k_fh);
                    k_fh = abs(min(k_fh,0));
                    k_h = mean(k_fh);
                    k_h = k_h./sum(k_h(:));
                    tau.(r_type{r})(k) = k_h*h;
                case 'sep'
                    k_h = flipud(kernels.(r_type{r}){k}.k_h); %Get the k_h from the norm model
                    k_h = min(k_h,0); %Remove all positive values
                    k_h = abs(k_h); %Make all inhibitory values positive
                    k_h = k_h./sum(k_h(:)); %Scale the values to sum to 1
                    tau.(r_type{r})(k) = (k_h'*h); %Compute a weighted sum of all values
            end
        end
        tau.(r_type{r}) = flipud(tau.(r_type{r})(:));
    end
    model_dir = model_kernels(m).folder;
    save_dir = fullfile(model_dir,'Plots');
    %Plot the betas vs frequency for all rooms
    line_sz = 4;
    axis_sz = 25;
    freqs = fliplr(freqs);
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    col{1} = 'r';
    col{2} = [0.91, 0.41, 0.17];
    col{3} = 'b';
    for j = 1:n_rooms
        plot(freqs,flipud(tau.(r_type{j})(:)),'Linewidth',line_sz,'Color',col{j});
    end
    xlabel('Frequency [kHz]');
    ylabel('\tau values [ms]');
    title('\tau (center of mass) for inhibition vs freqeuncy for all rooms');
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gcf,'color','w');
    set(gca, 'XScale', 'log');
    legend('big','med','small');
    save_name = fullfile(save_dir,['All_rooms_tau_vs_freq.png']);
    export_fig(save_name);
    close all;
end
