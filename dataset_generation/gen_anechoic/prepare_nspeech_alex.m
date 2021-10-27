input_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/Sounds/Anechoic_nonspeech_alex'; %Name of the folder containing the subfolders with different stimuli
fs_desired = 44100;
remove_s = 1; %How many s to clip from beginning and end of each sound
%to remove recording artefacts when starting/finishing the recording

dir_list = dir(input_dir);
dir_list = dir_list(~ismember({dir_list.name},{'.','..','Cropped'}));

save_folder = fullfile(dir_list(1).folder,'/Cropped');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for j = 1:length(dir_list)
    soundtype_name = fullfile(dir_list(j).folder,dir_list(j).name);
    snd_list = dir(soundtype_name);
    snd_list = snd_list(~ismember({snd_list.name},{'.','..'}));
    
    for s = 1:length(snd_list)
        data = [];
        file_current = fullfile(snd_list(s).folder,snd_list(s).name);
        [data,fs]=audioread(file_current);
        remove_samples = round(remove_s*fs);
        %Remove artifact chunks from beginning and end
        data = data(remove_samples+1:end-remove_samples,:);
        %Resample data if different sampling rate
        if fs ~= fs_desired
            data = resample(data,fs_desired,fs);
        end
        data = data./max(abs(data(:)));
        audiowrite(fullfile(save_folder,snd_list(s).name),data,fs_desired);
    end
    
end
