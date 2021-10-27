%% Params
fname = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_dataset/Individual_dereverb';
japan_name = 'Japanese_dereverb_all.wav';
speech_name = 'Speech_dereverb_all.wav';
ear = 'left'; %Which side to take from the stereo recordings
required_length_s = 100; %The total length I want in seconds
db_target = 53; %The 3 different classes of stimuli have different intrinsic loudness. It is better to equalize them to the same loudness
ramp_ms = 10;
fs_desired = 44100;
%% Make sounds
required_length_samples = required_length_s*fs_desired;
%Load the files
[japan,fs] = audioread(fullfile(fname,japan_name));
[speech,fs] = audioread(fullfile(fname,speech_name));

switch ear
    case 'left'
        side = 1;
    case 'right'
        side = 2;
end

speech = speech(:,side);

%Calculate the adjustments necessary for the loudness of the different sitmuli (dB)
adj_coeff_japan = db_adjust(japan,db_target);
adj_coeff_speech = db_adjust(speech,db_target);

%Make the loudness equal
japan = japan.*adj_coeff_japan;
if max(abs(japan(:)))>1
    warning('Japanese data is clipped!');
end

speech = speech.*adj_coeff_speech;
if max(abs(speech(:)))>1
    warning('Speech data is clipped!');
end

%Make the final file
data = [speech(:,side);japan(:,side)]; %Alex 22/07/2021

%Delete extra samples
total_length_samples = length(data);
clip_samples = total_length_samples - required_length_samples;
data(1:clip_samples,:) = [];
%Ramp the result
ramp_s = ramp_ms/1000;
len_s = length(data)/fs_desired;
ramp = cosrampenv(len_s,ramp_s,fs_desired);
ramp = ramp(:);
data = data.*ramp;
save_name = fullfile(fname,'Anechoic_data_dereverb_all.wav');
audiowrite(save_name,data,fs);