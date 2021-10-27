function make_model(coch_file,h_max_ms,model,normz,n_cores,save_dir,kfolds)
%% Define the params for all possible models
% save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Model_fits_600s';
if ~exist('model','var') || isempty(model)
    model = 'sep';
end

if ~exist('normz','var') || isempty(normz)
    normz = 'none';
end

if ~exist('h_max_ms','var') || isempty(h_max_ms)
    h_max_ms = 150;
end

if ~exist('n_cores','var') || isempty(n_cores)
    n_cores = 6;
end

if ~exist('kfolds','var') || isempty(n_cores)
    kfolds = true;
end
%% Fit the model
delete(gcp('nocreate'));
parpool('local',n_cores);

load(coch_file,'coch');
r_type = coch(1).r_type;
pos = coch(1).pos;
dt_ms = coch(1).params.dt_sec*1000;
chop = coch(1).chop;

if isnumeric(model)
    model_name = ['alpha_',num2str(model)];
else
    model_name = model;
end

fprintf('== Generating %s ==\n',['Room: ',r_type,pos,' positions ','Model:',h_max_ms,model_name,' Normz:',normz]);
kernels =  gen_norm_model(coch,h_max_ms,model,normz,kfolds);


room_dir = fullfile(save_dir,r_type);
if ~exist(room_dir, 'dir')
    mkdir(room_dir);
end

pos_dir = fullfile(room_dir,pos);
if ~exist(pos_dir, 'dir')
    mkdir(pos_dir);
end

chop_dir = fullfile(pos_dir,chop);
if ~exist(chop_dir, 'dir')
    mkdir(chop_dir);
end

norm_dir = fullfile(chop_dir,normz);
if ~exist(norm_dir, 'dir')
    mkdir(norm_dir);
end

model_dir = fullfile(norm_dir,model_name);
if ~exist(model_dir, 'dir')
    mkdir(model_dir);
end

bin_dir = fullfile(model_dir,[num2str(dt_ms,'%1.0f'),'ms']);
if ~exist(bin_dir, 'dir')
    mkdir(bin_dir);
end

hist_dir = fullfile(bin_dir,[num2str(h_max_ms,'%1.0f'),'ms']);
if ~exist(hist_dir, 'dir')
    mkdir(hist_dir);
end


save(fullfile(hist_dir,'kernel.mat'),'kernels','-v7.3');

