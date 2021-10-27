fname = '/mnt/40086D4C086D41D0/Reverb_analysis/Sounds/Anechoic_speech/Stitch';
fs_desired = 44100;
snd_files = dir(fname);
snd_files = snd_files(~ismember({snd_files.name},{'.','..'}));
snd = cell(30,1);
for s = 1:length(snd_files)
    file_current = fullfile(snd_files(s).folder,snd_files(s).name);
    [snd{s},fs]=audioread(file_current);
end

cat_snd = cell2mat(snd);

final_snd = resample(cat_snd,fs_desired,fs);

save_name = fullfile(snd_files(1).folder,['Peter.wav']);
audiowrite(save_name,final_snd,fs_desired);
