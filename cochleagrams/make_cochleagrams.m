%% Cochleagram params
stim_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Conditions/100s_derverb_new';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/100s_dereverb_new_exp_400_19k_specpower';
r_type = 'closed';
pos = 'single';
chop = 'chopped';
fs_desired = 44100; %The sampling rate
dt_ms = 10; %Time bin in ms
coch_type = 'specpower'; %The type of cochleagram 
freq_spacing = 'log'; %The spacing of the frequencies
f_min = 400; %The minimal frequency
f_max = 19000; %The maximal frequency
n_f = 30; %The total number of freqeuncies
color_map = 'inferno'; %Options are 'inferno', 'magma', 'plasma', 'viridis'
%% Load the data
fname{1} = 'anech_room';
fname{2} = 'small_room';
fname{3} = 'med_room';
fname{4} = 'big_room';
n_stim = numel(fname);

fprintf('== Loading the data ==\n');
store = cell(4,1);
for s = 1:n_stim
    store{s} = load(fullfile(stim_dir,[fname{s},'.mat']));
end
%% Make the cochleagrams
coch(1).reverb_cond = 'anech';
coch(2).reverb_cond = 'small';
coch(3).reverb_cond = 'med';
coch(4).reverb_cond = 'big';
coch(1).type = coch_type;
coch(1).r_type = r_type; 
coch(1).pos = pos;
coch(1).chop = chop;
num_cond = length(coch);

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
        case 'bencoch'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram(store{s}.data, fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'mscpower'
            [coch(s).X_ft, coch(s).t, coch(s).params] = cochleagram_msc_power(store{s}.data, fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
    end
    coch(s).X_ft = flipud(coch(s).X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
    fprintf('== Done! This took %0.fs ==\n',toc);
end

if strcmp(coch_type,'mscpower')
    %Params for the hill function
    n = 1;
    c = 0.0025;
    SAT = 1;
    %Divide by median in every f-band and run hill function
    for s = 2:num_cond
        coch(s).anech_X_ft = coch(1).X_ft./median([coch(1).X_ft,coch(s).X_ft],2);
        revtemp_X_ft{s-1} = coch(s).X_ft./median([coch(1).X_ft,coch(s).X_ft],2);
        coch(s).anech_X_ft=hill_function(coch(s).anech_X_ft,n,c,SAT);
        revtemp_X_ft{s-1}=hill_function(revtemp_X_ft{s-1},n,c,SAT);
    end
    coch = rmfield(coch,'X_ft');
    for s = 1:num_cond-1
        coch(s).X_ft = revtemp_X_ft{s};
    end
end
%% Save the results
room_dir = fullfile(save_dir,r_type);
if ~exist(room_dir, 'dir')
    mkdir(room_dir);
end

pos_dir = fullfile(room_dir,pos);
if ~exist(pos_dir, 'dir')
    mkdir(pos_dir);
end

chop_dir = fullfile(pos_dir,chop);
if ~exist(chop_dir, 'dir')
    mkdir(chop_dir);
end

bin_dir = fullfile(chop_dir,[num2str(dt_ms),'ms']);
if ~exist(bin_dir, 'dir')
    mkdir(bin_dir);
end

coch_dir = fullfile(bin_dir,coch_type);
if ~exist(coch_dir, 'dir')
    mkdir(coch_dir);
end

save_name = fullfile(coch_dir,['coch_all_conditions_',r_type,'_',pos,'pos_',coch_type,'.mat']);
save(save_name,'coch');

%% Plot the full cochleagrams 
freqs = fliplr(coch(1).params.freqs);
num_freq = length(freqs);
skip_f = 2;
freqs = freqs(1:skip_f:end);
n_tlab = 10;
skip_t = round(length(coch(1).t)/n_tlab);
figure('units','normalized','outerposition',[0 0 1 1]);
row = 2;
col = 2;
per = 0.05;
edgel = per; edger = 0.02; edgeh = per; edgeb = 0.06; space_h = 0.06; space_v = 0.1;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

for f = 1:numel(freqs)
    y_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

for tm = 1:n_tlab
    x_labels{tm} = num2str(coch(1).t((tm-1)*skip_t +1),'%.0f');
end

for s = 1:size(coch,2)
    subplot('position',pos{s});
    imagesc(coch(s).X_ft);
    colorbar;
    colormap(color_map);
    if strcmp(coch_type,'spechill')
        caxis([0 1]);
    end
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_t:length(coch(1).t)]);
    xticklabels(x_labels);
    title(coch(s).reverb_cond);
    set(gca,'FontName','Arial','FontSize',15,'FontWeight','Bold');
    xlabel('Time [s]','FontSize',16,'FontWeight','bold');
    ylabel('Freqeuncy [kHz]','FontSize',16,'FontWeight','bold');
end
set(gcf,'color','w');
save_name = fullfile(coch_dir,['All_cochleagrams_',coch_type,'.png']);
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
    sgtitle([coch(s).reverb_cond,' room'],'FontSize',18,'Color','red','FontWeight','Bold');
    set(gcf,'color','w');
    save_name = fullfile(coch_dir,[coch(s).reverb_cond,'_fband_value_distribtuion',coch_type,'.png']);
    export_fig(save_name);
    close;
end


    function y=hill_function(x_in,n,c,SAT)
        % y=hill_function(x,n,c)
        if ~exist('n','var')
            n=1;
        end
        if ~exist('c','var')
            c=0.0025;
        end
        if ~exist('SAT','var')
            SAT=1;
        end
        x_to_n = (c*x_in).^n;
        y=x_to_n./(SAT+x_to_n);
    end
