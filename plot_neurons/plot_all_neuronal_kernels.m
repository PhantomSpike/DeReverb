function plot_all_neuronal_kernels(kernel_dir,NPSP_th,type,ker_per_plot)

%% Default params
if ~exist('ker_per_plot','var') || isempty(ker_per_plot)
    ker_per_plot = 10;
end

if ~exist('NPSP_th','var') || isempty(NPSP_th)
    NPSP_th = 10;
end

%% Load the data
load(fullfile(kernel_dir,'info'),'info');

%% Select the data
ix = info.NPSP<NPSP_th;
NPSPs = info.NPSP(ix);
clusters = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);

%Sort in increasing NPSP
[~,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
clusters = clusters(ix_select);

%% Load the kernels
n_ker = length(ix_select);
kernels = cell(n_ker,1);
if strcmp(type,'slide')
    temp_ker = load(fullfile(kernel_dir,strjoin({animal_names{1},pen_names{1},num2str(clusters(1))},'_')));
    n_steps = length(temp_ker.kernel);
    middle = ceil(n_steps/2);
    for k = 1:n_ker
        temp = load(fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(clusters(k))},'_')));
        kernels{k,1}.kernel = temp.kernel(1,[1,middle,n_steps]);
    end
else
    for k = 1:n_ker
        kernels{k,1} = load(fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(clusters(k))},'_')));
    end
end


%% Plot them
[root_dir,type_ker] = fileparts(kernel_dir);
plot_dir = fullfile(root_dir,'Plots');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
plot_dir = fullfile(plot_dir,type_ker);
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

switch type
    case 'neurons'
        plot_neuronal_kernels(kernels,plot_dir,ker_per_plot, animal_names);
    case 'switch'
        plot_neuronal_kernels_switch(kernels,plot_dir,ker_per_plot);
    case 'slide'
        plot_neuronal_kernels_slide(kernels,plot_dir,ker_per_plot);
    case 'lnp'
        plot_lnp_kernels(kernels,plot_dir,ker_per_plot)
end

end

