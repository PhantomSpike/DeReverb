%This function will convolve the IRs from the different positions with the
%anechoic stimuli
%% Define params
%Sound data set to be convolved with IR
data_name = '/mnt/40086D4C086D41D0/Reverb_normative/Sounds/Anechoic_dataset/600s/Anechoic_data_all.wav'; 
save_folder = '/mnt/40086D4C086D41D0/Reverb_normative/Conditions/600s/Closed_40pos_600s';
%The IRs to be used
ir_all_name = '/mnt/40086D4C086D41D0/Reverb_normative/IRs/Chopped_experimental/40pos/All/all_irs.mat';
%Bunch of other params
ear = 'right'; %The IR from which ear to use for the convolution 
fs_desired = 44100; %The sampling rate
db_target = 73; %The target dB for all stimuli. This is the same as the ferret experiments 
ramp_s = 0.050; %The length of cosine ramping for each chunk
normalize_pos = 1; %Logical flag whether to normalize each position

%Different rooms
r_name{1} = 'small_room';
r_name{2} = 'med_room';
r_name{3} = 'big_room';
n_rooms = length(r_name);
%Different room types
r_type{1} = 'closed';
r_type{2} = 'anechoic';
n_types = length(r_type);

params.ear = ear;
params.normalize_pos = normalize_pos;
%% Load the data
[snd_data,fs] = audioread(data_name);
load(ir_all_name,'ir');

%Check for mismatch in sampling rate
if fs~=ir.small_room.anechoic{1}.Fs || fs~=ir.small_room.closed{1}.Fs...
        || fs~=ir.med_room.anechoic{1}.Fs || fs~=ir.med_room.closed{1}.Fs...
        || fs~=ir.big_room.anechoic{1}.Fs || fs~=ir.big_room.closed{1}.Fs...
        || fs~=fs_desired 
    error('Sample rate mismatch between the sound data and IRs OR wrong sampling rate')
end

%Select IR from which ear to convolve
switch ear
    case 'left'
        side = 1;
    case 'right'
        side = 2;
end

%% Split the stimulus into n_pos number of equal chunks
l_total = length(snd_data);
n_pos = length(ir.small_room.anechoic);
l_chunk_samp = round(l_total/n_pos);
l_chunk_samp_last = l_total-(n_pos-1)*l_chunk_samp;
start_ix = [1:l_chunk_samp:l_total]; %Get the starting index of each chunk
snd_data_chop = cell(n_pos,1);

for p = 1:n_pos
    ix = start_ix(p);
    if p == n_pos
        snd_data_chop{p} =  snd_data(ix:ix+(l_chunk_samp_last-1));
    else
        snd_data_chop{p} =  snd_data(ix:ix+(l_chunk_samp-1));
    end
    %Ramp every position to make sure there are no clicks
    ramp = cosrampenv(length(snd_data_chop{p})/fs_desired,ramp_s,fs_desired);
    ramp = ramp(:);
    snd_data_chop{p} = snd_data_chop{p}.*ramp;
end
%% Convolve the IRs of the rooms for every position with the chunked data - be aware of shifts in the data due to the convolution!!!
for r = 1:n_rooms
    fprintf('== Performing %s Convolutions ==\n',r_name{r});tic;
    for p = 1:n_pos
        snd_conv.(r_name{r}).(r_type{1}){p,1} = conv(snd_data_chop{p},ir.(r_name{r}).(r_type{1}){p}.data(:,side),'full');
        snd_conv.(r_name{r}).(r_type{2}){p,1} = conv(snd_data_chop{p},ir.(r_name{r}).(r_type{2}){p}.data(:,side),'full');
        params.ir_names.([r_name{r},'_',r_type{1}]){p} = ir.(r_name{r}).(r_type{1}){p}.fname;
        params.ir_names.([r_name{r},'_',r_type{2}]){p} = ir.(r_name{r}).(r_type{2}){p}.fname;
        sz_anech = length(snd_conv.(r_name{r}).(r_type{2}){p}); %Find the length of the anech clip in samples so we can cut the reverb one to same size
        snd_conv.(r_name{r}).(r_type{1}){p} = snd_conv.(r_name{r}).(r_type{1}){p}(1:sz_anech);
    end
    fprintf('== Done! This took %.0fs ==\n',toc);
end
%% Optional normalization
if normalize_pos
    fprintf('== Normalizing the positions for loudness ==\n');tic;
    for r = 1:n_rooms
        for p = 1:n_pos
            adj_coeff = db_adjust(snd_conv.(r_name{r}).(r_type{1}){p},db_target);
            snd_conv.(r_name{r}).(r_type{1}){p} = snd_conv.(r_name{r}).(r_type{1}){p}*adj_coeff;
            adj_coeff = db_adjust(snd_conv.(r_name{r}).(r_type{2}){p},db_target);
            snd_conv.(r_name{r}).(r_type{2}){p} = snd_conv.(r_name{r}).(r_type{2}){p}*adj_coeff;
        end
    end
    fprintf('== Done! This took %.0fs ==\n',toc);
end
%% Concatenate all positions together to normalize the sound levels across anech and reverb conditions
fprintf('== Concatenating and adjusting dB level of stimuli ==\n');tic;
reverb_all = cell(n_rooms,1);
anech_all = cell(n_rooms,1);
for r = 1:n_rooms
    %Adjust reverb
    reverb_all{r} = cell2mat(snd_conv.(r_name{r}).(r_type{1}));
    adj_coeff_reverb = db_adjust(reverb_all{r},db_target);
    reverb_all{r} = reverb_all{r}.*adj_coeff_reverb;

    %Adjust anechoic
    anech_all{r} = cell2mat(snd_conv.(r_name{r}).(r_type{2}));
    adj_coeff_anech = db_adjust(anech_all{r},db_target);
    anech_all{r} = anech_all{r}.*adj_coeff_anech;
end
fprintf('== Done! This took %.0fs ==\n',toc);

%% Double check sound levels


for r = 1:n_rooms
    %Check reverb
    params.db.([r_name{r},'_',r_type{1}]) = db_calc(reverb_all{r});
    
    %Check anechoic
    params.db.([r_name{r},'_',r_type{2}]) = db_calc(anech_all{r});
end

%Normalize the concatenated version before writing to audio files. Don't
%normalize the individual snippets as they can remain as .mat files
for r = 1:n_rooms
    %Find max(abs(x)) value across both stimuli
    max_abs = max(abs([anech_all{r}(:);reverb_all{r}(:)]));
    %normalize reverb
    reverb_all_norm{r} = reverb_all{r}/max_abs;
    
    %normalize anechoic
    anech_all_norm{r} = anech_all{r}/max_abs;
end
%See if there is clipping
for r = 1:n_rooms
    if max(abs(anech_all_norm{r}(:)))>1
        error([r_name{r},' anech data is clipped!']);
    else
        params.clipped.([r_name{r},'_anech']) = 0;
    end
    
    if max(abs(reverb_all_norm{r}(:)))>1
        error([r_name{r},' reverb data is clipped!']);
    else
        params.clipped.([r_name{r},'_reverb']) = 0;
    end
    
end


%% Save the results as wav files
%Transfer to three variables so you don't have to save as -v7.3'
small_room.anechoic = anech_all{1};
med_room.anechoic = anech_all{2};
big_room.anechoic = anech_all{3};

small_room.closed = reverb_all{1};
med_room.closed = reverb_all{2};
big_room.closed = reverb_all{3};

fprintf('== Saving the results ==\n');tic;
save(fullfile(save_folder,'params.mat'),'params');
save(fullfile(save_folder,['small_room.mat']),'small_room');
save(fullfile(save_folder,['med_room.mat']),'med_room');
save(fullfile(save_folder,['big_room.mat']),'big_room');

for r = 1:n_rooms
    audiowrite(fullfile(save_folder,[r_name{r},'_',r_type{1},'_room.wav']),reverb_all_norm{r},fs_desired);
    audiowrite(fullfile(save_folder,[r_name{r},'_',r_type{2},'_room.wav']),anech_all_norm{r},fs_desired);
end
fprintf('== Done! This took %.0fs ==\n',toc);
