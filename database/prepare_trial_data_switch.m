function data = prepare_trial_data_switch(clust_info,sorted_dir,cluster_id_list)
%% Define parameters

min_trig_length_s = 39.9; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1; %The minmum inter trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
actual_stimlength_s = 39.9; %The length of each presentation without a noise burst in s
n_reps = 10; 
n_stim = 8; %Number of different stimuli. There are 2 different conditions: R1 start and R2 start. R1 is from small room, R2 is from large room.
%There are 4 different versions of each kind which have different sounds.
%However, the general structure is always:

%R1: |--small 8s--|->|--large 8s--|->|--small 8s--|->|--large 8s--|->|--small 8s--|

%R2: |--large 8s--|->|--small 8s--|->|--large 8s--|->|--small 8s--|->|--large 8s--|

% 1 – switch8_r1_start_1
% 2 – switch8_r1_start_2
% 3 – switch8_r1_start_3
% 4 – switch8_r1_start_4
% 5 – switch8_r2_start_1
% 6 – switch8_r2_start_2
% 7 – switch8_r2_start_3
% 8 – switch8_r2_start_4
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
data{max_clust}.stim(n_stim).repeat(n_reps).spiketimes = []; %Initialize the struct

fprintf('== Getting clusters ==\n');tic;
for c = 1:num_clust
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,num_clust);
    clust_id = cluster_id_list(c);
    spiketimes = clust_info.spikeTimes(clust_info.spikeTimes(:,2)==clust_id,1); %Get the spiketimes that belong to it

    for stim = 1:n_stim
        ix_reps = [];
        ix_reps = find(grid.randomisedGrid==stim); %Find the indices of the repetitions corresponding to this stimulus
        
        for rep = 1:n_reps
            curr_ix = ix_reps(rep);
            %Get only the spiketimes that correspond to this stimulus and
            %this repeat
            spiketimes_rep = spiketimes(spiketimes>=start_time_s(curr_ix) & spiketimes<(start_time_s(curr_ix)+actual_stimlength_s));
            %Normalize relative to the beginning
            data{c}.stim(stim).repeat(rep).spiketimes = spiketimes_rep - start_time_s(curr_ix);
        end
        
    end
    
end
fprintf('== Done! This took %.1fs ==\n',toc);