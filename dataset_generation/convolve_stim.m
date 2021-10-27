%Sound data set to be convolved with IR
data_name = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_dataset/Individual_dereverb/Anechoic_data_dereverb_all.wav'; 
save_folder = '/mnt/40086D4C086D41D0/Reverb_normative/Conditions/100s_derverb_new';
%The IRs to be used
ir_anech_name = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Chopped_experimental/experimental/ferret_anech.mat';
ir_small_name = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Chopped_experimental/experimental/ferret_small_chopped.mat';
ir_med_name = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Chopped_experimental/experimental/ferret_med_chopped.mat';
ir_big_name = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Chopped_experimental/experimental/ferret_big_chopped.mat';
ear = 'right'; %The IR from which ear to use for the convolution 
fs_desired = 44100;
db_target = 73; %The target dB for all stimuli. This is the same as the ferret experiments and intentionally lower to prevent clipping
params.ear = ear;
params.ir_anech_name = ir_anech_name;
params.ir_small_name = ir_small_name;
params.ir_med_name = ir_med_name;
params.ir_big_name = ir_big_name;
%Load the data
[snd_data,fs] = audioread(data_name);
anech_ir = load(ir_anech_name);
small_ir = load(ir_small_name);
med_ir = load(ir_med_name);
big_ir = load(ir_big_name);

%Check for mismatch in sampling rate
if fs~=anech_ir.Fs || fs~=small_ir.Fs || fs~=med_ir.Fs || fs~=big_ir.Fs || fs~=fs_desired
    error('Sample rate mismatch between the sound data and IRs OR wrong sampling rate')
end

%Select IR from which ear to convolve
switch ear
    case 'left'
        side = 1;
    case 'right'
        side = 2;
end

%Do the convolution - Be very careful of shifts due to the convolution!!!
options.GPU = false; %Whether to use the GPU for the convolution function
fprintf('== Performing Anech Room Condition Convolution ==\n');tic;
anech_snd_data = conv(snd_data,anech_ir.data(:,side),'full');
fprintf('== Done! This took %.1fs ==\n',toc);

fprintf('== Performing Small Room Condition Convolution ==\n');tic;
small_snd_data = convnfft(snd_data, small_ir.data(:,side), 'full',1,options);
fprintf('== Done! This took %.1fs ==\n',toc);

fprintf('== Performing Medium Room Condition Convolution ==\n');tic;
med_snd_data = convnfft(snd_data, med_ir.data(:,side), 'full',1,options);
fprintf('== Done! This took %.1fs ==\n',toc);

fprintf('== Performing Big Room Condition Convolution ==\n');tic;
big_snd_data = convnfft(snd_data, big_ir.data(:,side), 'full',1,options);
fprintf('== Done! This took %.1fs ==\n',toc);

%Truncate the end from the reverb sounds to make  sure it matches the
%anechoic one
sz_anech = length(anech_snd_data);
small_snd_data = small_snd_data(1:sz_anech);
med_snd_data = med_snd_data(1:sz_anech);
big_snd_data = big_snd_data(1:sz_anech);
%Calculate the adjustments necessary for the loudness of the different sitmuli (dB)
adj_coeff_anech = db_adjust(anech_snd_data,db_target);
adj_coeff_small = db_adjust(small_snd_data,db_target);
adj_coeff_med = db_adjust(med_snd_data,db_target);
adj_coeff_big = db_adjust(big_snd_data,db_target);

%Normalize the data
anech_snd_data = anech_snd_data.*adj_coeff_anech;
small_snd_data = small_snd_data.*adj_coeff_small;
med_snd_data = med_snd_data.*adj_coeff_med;
big_snd_data = big_snd_data.*adj_coeff_big;

%Calculate the final dB
params.db.anech = db_calc(anech_snd_data);
params.db.small = db_calc(small_snd_data);
params.db.med = db_calc(med_snd_data);
params.db.big = db_calc(big_snd_data);

%Save the results as .mat file since .wav can clip them
fprintf('== Saving the results ==\n');tic;
save(fullfile(save_folder,'params.mat'),'params');
data = anech_snd_data;
save(fullfile(save_folder,['anech_room.mat']),'data');
data = small_snd_data;
save(fullfile(save_folder,['small_room.mat']),'data');
data = med_snd_data;
save(fullfile(save_folder,['med_room.mat']),'data');
data = big_snd_data;
save(fullfile(save_folder,['big_room.mat']),'data');

%Normalize between -1 and 1 so they can be saved as .wav files 
anech_snd_data = anech_snd_data./max(abs(anech_snd_data(:)));
if max(abs(anech_snd_data(:)))>1
    warning('Anechoic data is clipped!');
end

small_snd_data = small_snd_data./max(abs(small_snd_data(:)));
if max(abs(small_snd_data(:)))>1
    warning('Small reverb data is clipped!');
end

med_snd_data = med_snd_data./max(abs(med_snd_data(:)));
if max(abs(med_snd_data(:)))>1
    warning('Medium reverb data is clipped!');
end

big_snd_data = big_snd_data./max(abs(big_snd_data(:)));
if max(abs(big_snd_data(:)))>1
    warning('Big reverb data is clipped!');
end

audiowrite(fullfile(save_folder,'Anech_room.wav'),anech_snd_data,fs_desired);
audiowrite(fullfile(save_folder,'Small_room.wav'),small_snd_data,fs_desired);
audiowrite(fullfile(save_folder,'Medium_room.wav'),med_snd_data,fs_desired);
audiowrite(fullfile(save_folder,'Big_room.wav'),big_snd_data,fs_desired);
fprintf('== Done! This took %.0fs ==\n',toc);
