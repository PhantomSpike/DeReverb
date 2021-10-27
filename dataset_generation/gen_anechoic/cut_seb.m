fname = '/mnt/40086D4C086D41D0/Reverb_analysis/Sounds/Anechoic_speech/To_use/Sebastian.wav';
save_folder = '/mnt/40086D4C086D41D0/Reverb_analysis/Sounds/Anechoic_speech/To_use';
fs_desired = 44100;
to_cut_s = 103;
[snd,fs] = audioread(fname);
to_cut_samples = to_cut_s*fs;
snd(1:to_cut_samples,:) = [];
snd = resample(snd,fs_desired,fs);
audiowrite(fullfile(save_folder,'Sebastian_short.wav'),snd,fs_desired);
