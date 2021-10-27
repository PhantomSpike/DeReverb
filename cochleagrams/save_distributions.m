function save_distributions(th)
%% Cochleagram params
stim_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Conditions/900s_exp_mixed/Closed_single_pos_chopped_900s_exp_mixed';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/Troubleshoot';
r_type = 'closed';
pos = 'single';
chop = 'chopped';
fs_desired = 44100; %The sampling rate
dt_ms = 10; %Time bin in ms
coch_type = 'spechill'; %The type of cochleagram 
freq_spacing = 'log'; %The spacing of the frequencies
f_min = 400; %The minimal frequency
f_max = 19000; %The maximal frequency
n_f = 30; %The total number of freqeuncies
color_map = 'inferno'; %Options are 'inferno', 'magma', 'plasma', 'viridis'
%% Load the data
fname{1} = 'anech_room';
fname{2} = 'small_room';
fname{3} = 'med_room';
fname{4} = 'big_room';
n_stim = numel(fname);

fprintf('== Loading the data ==\n');
store = cell(4,1);
for s = 1:n_stim
    store{s} = load(fullfile(stim_dir,[fname{s},'.mat']));
end
%% Make the cochleagrams
coch(1).reverb_cond = 'anech';
coch(2).reverb_cond = 'small';
coch(3).reverb_cond = 'med';
coch(4).reverb_cond = 'big';
coch(1).type = coch_type;
coch(1).r_type = r_type; 
coch(1).pos = pos;
coch(1).chop = chop;
num_cond = length(coch);

fprintf('== Making the cochleagrams ==\n');
for s = 1:n_stim
    fprintf('== Cochleagram %.0f/%0.f ==\n',s,n_stim);tic;
    switch coch_type
        case 'spechill'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_Hill_th(store{s}.data, fs_desired, dt_ms, freq_spacing,th, f_min, f_max, n_f);
        case 'specpower'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_power_th(store{s}.data, fs_desired, dt_ms, freq_spacing,th, f_min, f_max, n_f);   
    end
    coch(s).X_ft = flipud(coch(s).X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
    fprintf('== Done! This took %0.fs ==\n',toc);
end
save(fullfile(save_dir,'coch'),'coch');