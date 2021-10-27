function data = prepare_trial_data(clust_info,sorted_dir,cluster_id_list)
%% Define parameters

min_trig_length_s = 39.9; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1; %The minmum inter trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
actual_stimlength_s = 39.9; %The length of each presentation without a noise burst in s
num_stim = 6; 
num_repeats = 10; %Number of different stimuli. They are 6 of them (3 reverb conditions and 2 sound stim sets) with 10 repeat each. Therefore the indices are like this
%For all other ferrets:
% 1 – anechoic stim 1 -> [1:10]
% 2 – anechoic stim 2 -> [11:20]
% 3 – small reverb stim 1 -> [21:30]
% 4 – small reverb stim 2 -> [31:40]
% 5 – big reverb stim 1 -> [41:50]
% 6 – big reverb stim 2 -> [51:60]

%For Noah and new one
% 1 – small stim 1 -> [1:10]
% 2 – small stim 2 -> [11:20]
% 3 – med reverb stim 1 -> [21:30]
% 4 – med reverb stim 2 -> [31:40]
% 5 – big reverb stim 1 -> [41:50]
% 6 – big reverb stim 2 -> [51:60]
%% Load synch_ch and grid_info
load(fullfile(sorted_dir,'synch_ch.mat'),'synch_ch');

grid_list = dir(fullfile(sorted_dir,'Meta'));
grid_list = grid_list(~ismember({grid_list.name},{'.','..'}));
grid_root = fullfile(grid_list(1).folder,grid_list(1).name);
dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
load(grid_filename,'grid');
%% Get the triggers
fs = clust_info.fs;
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
start_time_s = start_time_ms/1000; % Convert to s
%% Loop through every cluster and save the spiketimes for the particular stim and repeat
num_clust = length(cluster_id_list); %Find the total number of clusters
max_clust = max(cluster_id_list);
data{max_clust}.stim(num_stim).repeat(num_repeats).spiketimes = []; %Initialize the struct

fprintf('== Getting clusters ==\n');tic;
for c = 1:num_clust
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,num_clust);
    clust_id = cluster_id_list(c);
    spiketimes = clust_info.spikeTimes(clust_info.spikeTimes(:,2)==clust_id,1); %Get the spiketimes that belong to it
    stim_ix = 0; %Initialize the stim ix for every cluster
    for s = 1:num_stim
        
        for r = 1:num_repeats
            stim_ix = stim_ix + 1;
            ix_order = find(grid.randomisedGrid==stim_ix); %Find in which order was this particular stim and rep played
            %Get only the spiketimes that correspond to this stimulus and
            %this repeat
            spiketimes_rep = spiketimes(spiketimes>=start_time_s(ix_order) & spiketimes<(start_time_s(ix_order)+actual_stimlength_s));
            %Normalize relative to the beginning
            data{c}.stim(s).repeat(r).spiketimes = spiketimes_rep - start_time_s(ix_order);
        end
        
    end
    
end
fprintf('== Done! This took %.1fs ==\n',toc);