%% Params
fname = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_dataset/Individual';
jap_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Japanese_nonspeech_anechoic/Individual_refined';
nspeech_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_nonspeech_all/Audacity_new';
speech_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_speech/To_use_cropped/Audacity_new';
fs_desired = 44100;
db_target = 53; %The 3 different classes of stimuli have different intrinsic loudness. It is better to equalize them to the same loudness
jap_stim_s = 30; %Total length of japanese sounds in s
nspeech_stim_s = 360; %Total length of non-speech sounds in s
speech_stim_s = 210; %Total length of speech sounds in s
total_stim_s = 600; %Total length of the stimuli in s
ramp_ms = 10; %How much to ramp in ms
rnd_seed = 7; %Select a random seed number

%% Initialise vars
%Make jap stim
fprintf('== Making japanese sounds ==\n');
total_stim_samples = jap_stim_s*fs_desired;
grande.jap = make_ind_func(jap_dir,total_stim_samples,db_target,ramp_ms,fs_desired);

fprintf('== Making non-speech sounds ==\n');
%Make nspeech stim
total_stim_samples = nspeech_stim_s*fs_desired;
grande.nspeech = make_ind_func(nspeech_dir,total_stim_samples,db_target,ramp_ms,fs_desired);

fprintf('== Making speech sounds ==\n');
%Make speech stim
total_stim_samples = speech_stim_s*fs_desired;
grande.speech = make_ind_func(speech_dir,total_stim_samples,db_target,ramp_ms,fs_desired);

%% Permute the indices
total_snd = [grande.jap;grande.nspeech;grande.speech]; %Put all sounds together
n_clips = length(total_snd); %Find the total number of sound clips
rng(rnd_seed); %Use the selected random seed for reproducible randomness
perm_ix = randperm(n_clips); %Permute the indices of the sound clips
perm_total_snd = total_snd(perm_ix); %Shuffle the sounds
final_snd = cell2mat(perm_total_snd); %Convert to mat files

%% Remove extra samples, ramp and save the results to a .wav file
%Delete extra samples
required_length_samples = total_stim_s*fs_desired;
final_snd(required_length_samples+1:end) = [];
%Ramp the result
ramp_s = ramp_ms/1000;
len_s = length(final_snd)/fs_desired;
ramp = cosrampenv(len_s,ramp_s,fs_desired);
ramp = ramp(:);
final_snd = final_snd.*ramp;
fprintf('== Saving the results ==\n');
save_name = fullfile(fname,'Anechoic_data_all.wav');
audiowrite(save_name,final_snd,fs_desired);




