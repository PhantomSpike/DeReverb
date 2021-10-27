%% Params
fname = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_speech/To_use_cropped/Individual_new_dereverb';
fs_desired = 44100;
total_stim_s = 80; %How long I want the total duration of all the stimuli in s
total_stim_samples = round(total_stim_s*fs_desired); %Total duration in samples
remove_s = 0.5; %How many s to clip from beginning and end of each sound 
%to remove recording artefacts when starting/finishing the recording
ramp_ms = 10; %How long to ramp each sound
%% Make sounds
snd_files = dir(fname);
snd_files = snd_files(~ismember({snd_files.name},{'.','..'}));
n_speakers = length(snd_files);
n_samples = round((total_stim_s*fs_desired)/n_speakers); %Divide the total number of samples (data) evenly between the speakers
snd = [];

ramp_s = ramp_ms/1000;
for s = 1:n_speakers
    data = [];
    ramp = [];
    file_current = fullfile(snd_files(s).folder,snd_files(s).name);
    [data,fs]=audioread(file_current);
    %Resample data if different sampling rate
    if fs ~= fs_desired
        data = resample(data,fs_desired,fs);
    end
    remove_samples = round(remove_s*fs_desired);
    %Remove first and last 2s and take a certain chunk of each speaker
    data = data(remove_samples+1:remove_samples+n_samples+1,:);
    %Ramp the beginning and end
    len_s = length(data)/fs_desired;
    ramp = cosrampenv(len_s,ramp_s,fs_desired);
    ramp = ramp(:);
    data = data.*ramp;
    snd = [snd;data];
end

snd = snd./max(abs(snd(:)));
save_name = fullfile(snd_files(1).folder,['Speech_dereverb_all.wav']);
audiowrite(save_name,snd,fs_desired);
