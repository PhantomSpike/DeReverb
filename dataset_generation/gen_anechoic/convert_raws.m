input_dir = '/home/alex/Desktop/Data/Sounds/RWCP/nospeech/drysrc'; %Directory with the files
save_dir = '/home/alex/Desktop/Data/Sounds/RWCP/nospeech/All_sounds'; %Directory where we want to save
fs = 48000; %Sampling rate of the sound .raw files
fs_new = 44100; %New desired sampling rate

dir_list = dir(input_dir);
dir_list = dir_list(~ismember({dir_list.name},{'.','..','scripts'}));

for f = 1:numel(dir_list)
    fprintf('== Processing directory %s ==\n',dir_list(f).name);tic;
    category_path = fullfile(dir_list(f).folder,dir_list(f).name);
    category_list = dir(category_path);
    category_list = category_list(~ismember({category_list.name},{'.','..'}));
    
    for c = 1:numel(category_list)
        raw_path = fullfile(category_list(c).folder,[category_list(c).name,'/48khz']);
        raw_list = dir(raw_path);
        raw_list = raw_list(~ismember({raw_list.name},{'.','..'}));
        save_folder = fullfile(save_dir,category_list(c).name);
        
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        
        for j = 1:numel(raw_list)
            file_name = fullfile(raw_list(j).folder,raw_list(j).name);
            snd = load_raw(file_name); %Load the raw binary file
            snd = snd./max(abs(snd(:))); %Normalize the sound to be between -1 and 1
            snd = resample(snd,fs_new,fs); %Resample the sound to the new fs
            snd = snd./max(abs(snd(:))); %Normalize the sound to be between -1 and 1
            save_name = fullfile(save_folder,['/',category_list(c).name,'_',raw_list(j).name(1:end-4),'.wav']);
            audiowrite(save_name,snd,fs_new);
        end
        
    end
    fprintf('== Done! Took %.1fs ==\n',toc);
end