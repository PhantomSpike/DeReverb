stim_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Sound_stim/Switch_Noah';
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/Switch';
dt_ms = 10;
coch_type = 'specpower';

%% Cochleagram params
if ~exist('coch_type','var') || isempty(coch_type)
    coch_type = 'spechill';
end

if ~exist('dt_ms','var') || isempty(dt_ms)
    dt_ms = 10;
end

fs_desired = 48828; %The sampling rate
ear = 'right'; %Select the ear I want to use for the cochleagram
freq_spacing = 'log'; %The spacing of the frequencies
f_min = 400; %The minimal frequency
f_max = 19000; %The maximal frequency
n_f = 30; %The total number of freqeuncies
actual_stimlength_s = 39.9; %The length of each presentation without a noise burst in s
cor_factor = 10; %If the data comes from Derry/Kilkenny/Cork correct for the loudness to bring to 73dB 
%% Load the data
%R1: |--small 8s--|->|--large 8s--|->|--small 8s--|->|--large 8s--|->|--small 8s--|
%R2: |--large 8s--|->|--small 8s--|->|--large 8s--|->|--small 8s--|->|--large 8s--|

% 1 – switch8_r1_start_1
% 2 – switch8_r1_start_2
% 3 – switch8_r1_start_3
% 4 – switch8_r1_start_4
% 5 – switch8_r2_start_1
% 6 – switch8_r2_start_2
% 7 – switch8_r2_start_3
% 8 – switch8_r2_start_4

fname{1} = 'switch8_r1_start_1';
fname{2} = 'switch8_r1_start_2';
fname{3} = 'switch8_r1_start_3';
fname{4} = 'switch8_r1_start_4';
fname{5} = 'switch8_r2_start_1';
fname{6} = 'switch8_r2_start_2';
fname{7} = 'switch8_r2_start_3';
fname{8} = 'switch8_r2_start_4';
n_stim = numel(fname);

[~,exp_name] = fileparts(stim_dir);

fprintf('== Loading the data ==\n');
store = cell(n_stim,1);
actual_stimlength_samples =  actual_stimlength_s*fs_desired;
for s = 1:n_stim
    [store{s}.data,fs] = audioread(fullfile(stim_dir,[fname{s},'.wav']));
    
    store{s}.data = store{s}.data.*cor_factor; %Correct for the boosting after benware

    if fs ~= fs_desired
        error('Wrong sampling rate');
    end
    
    switch ear
        case 'left'
            store{s}.data = store{s}.data(:,1);
        case 'right'
            store{s}.data = store{s}.data(:,2);
    end
    store{s}.data(actual_stimlength_samples+1:end) = []; %Delete the last 4s which have the noise burst
    end
%% Make the cochleagrams
coch(1).reverb_cond = 'small_1';
coch(2).reverb_cond = 'small_2';
coch(3).reverb_cond = 'small_3';
coch(4).reverb_cond = 'small_4';
coch(5).reverb_cond = 'large_1';
coch(6).reverb_cond = 'large_2';
coch(7).reverb_cond = 'large_3';
coch(8).reverb_cond = 'large_4';
coch(1).type = coch_type;
coch(1).pos = 'single';
coch(1).ear = ear;

fprintf('== Making the cochleagrams ==\n');
for s = 1:n_stim
    fprintf('== Cochleagram %.0f/%0.f ==\n',s,n_stim);tic;
    switch coch_type
        case 'speclog'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_log(store{s}.data, fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'spechill'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_Hill(store{s}.data, fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'specpower'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_spec_power(store{s}.data, fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);     
    end
    coch(s).X_ft = flipud(coch(s).X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
    fprintf('== Done! This took %0.fs ==\n',toc);
end

%% Save the results
exp_name_dir = fullfile(save_dir,exp_name);
if ~exist(exp_name_dir, 'dir')
    mkdir(exp_name_dir);
end

model_dir = fullfile(exp_name_dir,coch_type);
if ~exist(exp_name_dir, 'dir')
    mkdir(model_dir);
end

bin_dir = fullfile(model_dir,[num2str(dt_ms),'ms']);
if ~exist(bin_dir, 'dir')
    mkdir(bin_dir);
end

save_name = fullfile(bin_dir,['coch_all_conditions_',coch_type,'.mat']);
save(save_name,'coch');

%% Plot the full cochleagrams 
freqs = fliplr(coch(1).params.freqs);
num_freq = length(freqs);
skip_f = 3;
freqs = freqs(1:skip_f:end);
n_tlab = 10;
skip_t = round(length(coch(1).t)/n_tlab);
figure('units','normalized','outerposition',[0 0 1 1]);
row = 4;
col = 2;
per = 0.05;
edgel = per; edger = 0.02; edgeh = per; edgeb = 0.06; space_h = 0.06; space_v = 0.1;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(coch(1).t((tm-1)*skip_t +1),'%.1f');
end

for s = 1:size(coch,2)
    subplot('position',pos{s});
    imagesc(coch(s).X_ft);
    colorbar;
    colormap('inferno');
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_t:length(coch(1).t)]);
    xticklabels(x_labels);
    title(coch(s).reverb_cond);
    set(gca,'FontName','Arial','FontSize',17,'FontWeight','Bold');
    xlabel('Time [s]','FontSize',18,'FontWeight','bold');
    ylabel('Freqeuncy [kHz]','FontSize',18,'FontWeight','bold');
end
set(gcf,'color','w');
save_name = fullfile(bin_dir,['All_cochleagrams_',coch_type,'.png']);
export_fig(save_name);
close;

%% Plot the distribution of values in the different freqeuncy bands in the cochleageam
freqs = fliplr(coch(1).params.freqs);
for f = 1:numel(freqs)
    f_labels{f} = num2str(freqs(f)./1000,'%.1f');
end
row = 6;
col = 5;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.1; edgeb = per; space_h = per; space_v = 0.06;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
for s = 1:size(coch,2)
    figure('units','normalized','outerposition',[0 0 1 1]);
    for f = 1:num_freq
        subplot('position',pos{f});
        histogram(coch(s).X_ft(f,:),'Normalization','probability');
        title([f_labels{f},' kHz']);
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    end
    sgtitle([coch(s).reverb_cond,' stim'],'FontSize',18,'Color','red','FontWeight','Bold');
    set(gcf,'color','w');
    save_name = fullfile(bin_dir,[coch(s).reverb_cond,'_fband_value_distribtuion',coch_type,'.png']);
    export_fig(save_name);
    close;
end
