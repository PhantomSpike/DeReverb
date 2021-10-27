%% Params
dir_name = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_nonspeech_all/Individual_exp_adjusted';
fs_desired = 44100;
total_stim_s = 355; %How long I want the total duration of all the stimuli in s
required_length_samples = round(total_stim_s*fs_desired); %Total duration in samples
ramp_ms = 10; %How much to ramp in ms
start_ms = 10; %How much to cut off from the begining in ms
end_ms = 10; %How much to cut off from the end in ms
%% Make sounds
cat_dir = dir(dir_name);
cat_dir = cat_dir(~ismember({cat_dir.name},{'.','..'}));
num_cat = length(cat_dir);
cat_snd_list = cell(num_cat,1);
snd = [];

for c = 1:num_cat
    cat_snd_list{c} = dir(fullfile(cat_dir(c).folder,cat_dir(c).name));
    cat_snd_list{c} = cat_snd_list{c}(~ismember({cat_snd_list{c}.name},{'.','..'}));
end

pass = 0; %Variable to keep track of how many times we are passing trough the data
ramp_s = ramp_ms/1000;
start_samples = (start_ms/1000)*fs_desired;
end_samples = (end_ms/1000)*fs_desired;
while length(snd)<required_length_samples
    pass = pass + 1;
    fprintf('== Iteration #%0.f ==\n',pass);
    for c = 1:num_cat
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
        %Remove the extra data
        data = data(start_samples:end-end_samples,1);
        %Ramp the beginning and end
        len_s = length(data)/fs_desired;
        ramp = cosrampenv(len_s,ramp_s,fs_desired);
        ramp = ramp(:);
        data = data.*ramp;
        snd = [snd;data];
        if length(snd)>=required_length_samples
            flag=1;
            break
        end
    end
    if flag == 1
        break
    end
end

%Delete extra samples
total_length_samples = length(snd);
clip_samples = total_length_samples - required_length_samples;
snd(1:clip_samples,:) = [];
%Ramp the result
len_s = length(snd)/fs_desired;
ramp = cosrampenv(len_s,ramp_s,fs_desired);
ramp = ramp(:);
snd_final = snd.*ramp;

snd_final = snd_final./max(abs(snd_final(:)));
save_name = fullfile(cat_dir(1).folder,['Non_speech_all.wav']);
audiowrite(save_name,snd_final,fs_desired);