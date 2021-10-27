function dat_cell = make_ind_func(stim_dir,total_stim_samples,db_target,ramp_ms,fs_desired)

%% Book keeping
cat_dir = dir(stim_dir);
cat_dir = cat_dir(~ismember({cat_dir.name},{'.','..'}));
num_cat = length(cat_dir);
cat_snd_list = cell(num_cat,1);

for c = 1:num_cat
    cat_snd_list{c} = dir(fullfile(cat_dir(c).folder,cat_dir(c).name));
    cat_snd_list{c} = cat_snd_list{c}(~ismember({cat_snd_list{c}.name},{'.','..'}));
end

dat_cell = cell(1);
pass = 0; %Variable to keep track of how many times we are passing trough the data
count = 0; %Counter var for total number of sounds
curr_length_samp = 0; %Counter for total number of samples
flag =0; %Flag to use for breaking out of the loop when ready 
ramp_s = ramp_ms/1000;

while curr_length_samp<total_stim_samples
    pass = pass + 1;
    fprintf('== Iteration #%0.f ==\n',pass);
    for c = 1:num_cat
        count = count + 1;
        data = [];
        ramp = [];
        if pass > length(cat_snd_list{c}) %Check if the number of files in the current category is less than the iteration. If so, just continue ot the next one
            continue
        end
        file_current = fullfile(cat_snd_list{c}(pass).folder,cat_snd_list{c}(pass).name);
        [data,fs]=audioread(file_current);
        %Resample data if different sampling rate
        if fs ~= fs_desired
            data = resample(data,fs_desired,fs);
        end
        %Ramp the beginning and end
        len_s = length(data)/fs_desired;
        ramp = cosrampenv(len_s,ramp_s,fs_desired);
        ramp = ramp(:);
        data = data.*ramp;
        dat_cell{count,1} = data;
        % Check if the total length of the data has reach desired length
        curr_length_samp = curr_length_samp + length(data);
        if curr_length_samp>=total_stim_samples
            flag=1;
            break
        end
    end
    
    if flag == 1
        break
    end
end

temp_data = cell2mat(dat_cell);
adj_coeff = db_adjust(temp_data(:),db_target);
dat_cell = cellfun(@(x) adj_coeff.*x, dat_cell, 'UniformOutput',false);