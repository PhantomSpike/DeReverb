neural_data_path = '/mnt/40086D4C086D41D0/Reverb_paper/fig_3_sup_1/Neuronal_data/neuronal_data.mat';
model_data_path = '/mnt/40086D4C086D41D0/Reverb_paper/fig_3_sup_1/LNP/lnp_data.mat';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_3_sup_1';

neurons = load(neural_data_path);
model = load(model_data_path);

%% Define params for scatter plot for Excitation and Inhibition
hist_type = 'bar'; get_neurons = 'all'; NPSP_th = 40; model_name = 'ridge'; font_type = 'Liberation Sans';

exc_color = 'r'; inh_color = 'b'; lw1 = 7; lw2 = 3;
all_font_sz = 55;
params_his_ei.fit = 'neurons'; params_his_ei.units = 'normal'; params_his_ei.calc_method = 'average';
%Plot general properties
params_his_ei.hist_type = hist_type; params_his_ei.get_neurons = get_neurons; params_his_ei.NPSP_th = NPSP_th; params_his_ei.model = model_name; params_his_ei.font_type = font_type;
%Plot colours
params_his_ei.exc_color = exc_color; params_his_ei.inh_color = inh_color;
%Plot size
params_his_ei.all_font_sz = all_font_sz; params_his_ei.lw1 = lw1; params_his_ei.lw2 = lw2;
%Save dir
params_his_ei.save_dir =  save_dir;


%% COM Excitation and Inhibiton stats
[p_val.com.neg,~,~] = signrank(neurons.diff_vec.com.neg, model.diff_vec.com.neg);
[p_val.com.pos,~,~] = signrank(neurons.diff_vec.com.pos, model.diff_vec.com.pos);
median_diff.com.neg = nanmedian(neurons.diff_vec.com.neg - model.diff_vec.com.neg);
median_diff.com.pos = nanmedian(neurons.diff_vec.com.pos - model.diff_vec.com.pos);
lim_val_ms = 60;

%% COM Excitation and Inhibiton small vs big Histogram
params_his_ei.specific_name = ' COM ms change for ';
params_his_ei.val_exc_big = neurons.diff_vec.com.pos; params_his_ei.val_exc_small = model.diff_vec.com.pos; params_his_ei.val_inh_big = neurons.diff_vec.com.neg; params_his_ei.val_inh_small = model.diff_vec.com.neg;
params_his_ei.p_val_exc = p_val.com.pos; params_his_ei.p_val_inh = p_val.com.neg;
params_his_ei.his_spacing = 2.5; params_his_ei.lim_val = 60;
plot_ei_histogram(params_his_ei);

%% PT Excitation and Inhibiton stats
[p_val.pt.neg,~,~] = signrank(neurons.diff_vec.pt.neg, model.diff_vec.pt.neg);
[p_val.pt.pos,~,~] = signrank(neurons.diff_vec.pt.pos, model.diff_vec.pt.pos);
median_diff.pt.neg = nanmedian(neurons.diff_vec.pt.neg - model.diff_vec.pt.neg);
median_diff.pt.pos = nanmedian(neurons.diff_vec.pt.pos - model.diff_vec.pt.pos);

%% PT Excitation and Inhibiton for small vs big room Histogram
params_his_ei.units = 'normal';
params_his_ei.specific_name = ' PT E and I change for ';
params_his_ei.val_exc_big = neurons.diff_vec.com.pos; params_his_ei.val_exc_small = model.diff_vec.com.pos; params_his_ei.val_inh_big = neurons.diff_vec.com.neg; params_his_ei.val_inh_small = model.diff_vec.com.neg;
params_his_ei.p_val_exc = p_val.pt.pos; params_his_ei.p_val_inh = p_val.pt.neg;
params_his_ei.his_spacing = 10; params_his_ei.lim_val = 120;
plot_ei_histogram(params_his_ei);
