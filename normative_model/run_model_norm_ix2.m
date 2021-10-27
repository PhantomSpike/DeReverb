function ker = run_model_norm_ix2(reverb_X_fht,anechoic_X_ft,model,room_sz,dt_s)

val_sz = 0.1; %What percentage to use for CV
fs = 1/dt_s; %Sampling rate
n_t = size(reverb_X_fht,3);
val_chunk_sz = floor(n_t*val_sz); %Size of each val chunk split evenly between stimuli

total_s = 600;
jap_start_s = 0; jap_end_s = 60; val_sz_jap =  floor(val_chunk_sz*((jap_end_s-jap_start_s)./total_s)); %Japanese stim end s
nspeech_start_s = 60; nspeech_end_s = 420; val_sz_nspeech =  floor(val_chunk_sz*((nspeech_end_s-nspeech_start_s)./total_s)); %Non-speech stim end s
speech_start_s = 420; speech_end_s = 600; val_sz_speech =  floor(val_chunk_sz*((speech_end_s-speech_start_s)./total_s)); %Speech stim end s

jap_start_ix = floor(jap_start_s*fs); jap_end_ix = floor(jap_end_s*fs); %Japanese stim start - end ix
nspeech_start_ix = floor(nspeech_start_s*fs); nspeech_end_ix = floor(nspeech_end_s*fs); %Non-speech stim start - end ix
speech_start_ix = floor(speech_start_s*fs); speech_end_ix = floor(speech_end_s*fs); %Speech stim start - end ix

jap_val_start_ix = randi([jap_start_ix jap_end_ix-val_sz_jap],1,1); jap_val_ix = [jap_val_start_ix:jap_val_start_ix+val_sz_jap];
nspeech_val_start_ix = randi([nspeech_start_ix nspeech_end_ix-val_sz_nspeech],1,1); nspeech_val_ix = [nspeech_val_start_ix:nspeech_val_start_ix+val_sz_nspeech];
speech_val_start_ix = randi([speech_start_ix speech_end_ix-val_sz_speech],1,1); speech_val_ix = [speech_val_start_ix:speech_val_start_ix+val_sz_speech];

val_ix = [jap_val_ix(:);nspeech_val_ix(:);speech_val_ix(:)];

if isnumeric(model)
    model_name = 'manual_alpha';
else
    model_name = model;
end
fprintf('== Fitting %s model for %s room ==\n',model_name,room_sz);tic
n_f = size(anechoic_X_ft,1);
ker = cell(n_f,1);
parfor k = 1:n_f
    fprintf('== Fitting frequency band %0.f/%0.f ==\n',k,n_f);
    if ischar(model)
        switch model
            case 'sep'
                ker{k} = sepkerneltensor2(reverb_X_fht, anechoic_X_ft(k,:));
            case {'lasso','ridge','elastic'}
                [ker{k}.main, ker{k}.allKernels, ker{k}.res] = alexnet_fht(reverb_X_fht, anechoic_X_ft(k,:), model,val_ix);
        end
    elseif isnumeric(model)
        [ker{k}.main, ker{k}.allKernels, ker{k}.res] = alexnet_fht(reverb_X_fht, anechoic_X_ft(k,:), model,val_ix);
    else
        error('Unrecognised model type')
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);
end