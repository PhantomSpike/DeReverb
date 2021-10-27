%% Cochleagram params
stim_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Conditions/600s/Closed_40pos_600s';
save_dir = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/600s';
chop = 'chopped';
fs_desired = 44100; %The sampling rate
dt_ms = 10; %Time bin in ms
coch_type = 'spechill'; %The type of cochleagram 
freq_spacing = 'log'; %The spacing of the frequencies
f_min = 200; %The minimal frequency
f_max = 16000; %The maximal frequency
n_f = 30; %The total number of freqeuncies
%% Load the date
%Different rooms
r_name{1} = 'small_room';
r_name{2} = 'med_room';
r_name{3} = 'big_room';
n_rooms = length(r_name);
%Different room types
r_type{1} = 'closed';
r_type{2} = 'anechoic';
n_types = length(r_type);
pos = 'many';
plot_name{1} = 'reverb';
plot_name{2} = 'anechoic';
fprintf('== Loading the data ==\n');
for r = 1:n_rooms
    data{r} = load(fullfile(stim_dir,r_name{r}));
end
%% Make the cochleagrams
coch(1).reverb_cond = 'small';
coch(2).reverb_cond = 'med';
coch(3).reverb_cond = 'big';
coch(1).type = coch_type;
coch(1).r_type = r_type{1} ;
coch(1).pos = pos;
num_cond = length(coch);

fprintf('== Making the cochleagrams ==\n');
for r = 1:n_rooms
    fprintf('== Cochleagrams %.0f/%0.f ==\n',r,n_rooms);tic;
    switch coch_type
        case 'spechill'
            [coch(r).reverb.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_Hill(data{r}.(r_name{r}).(r_type{1}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
            [coch(r).anechoic.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_Hill(data{r}.(r_name{r}).(r_type{2}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'speclog'
            [coch(r).reverb.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_log(data{r}.(r_name{r}).(r_type{1}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
            [coch(r).anechoic.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_log(data{r}.(r_name{r}).(r_type{2}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'specpower'
            [coch(r).reverb.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_power(data{r}.(r_name{r}).(r_type{1}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
            [coch(r).anechoic.X_ft, coch(r).t, coch(r).params] = cochleagram_spec_power(data{r}.(r_name{r}).(r_type{2}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
        case 'mscpower'
            [coch(r).reverb.X_ft, coch(r).t, coch(r).params] = cochleagram_msc_power(data{r}.(r_name{r}).(r_type{1}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
            [coch(r).anechoic.X_ft, coch(r).t, coch(r).params] = cochleagram_msc_power(data{r}.(r_name{r}).(r_type{2}), fs_desired, dt_ms, freq_spacing, f_min, f_max, n_f);
    end
    coch(r).reverb.X_ft = flipud(coch(r).reverb.X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
    coch(r).anechoic.X_ft = flipud(coch(r).anechoic.X_ft); %Keep in mind that cochleagram is output from top -> bottom (low->high) so I flip here
    fprintf('== Done! This took %0.fs ==\n',toc);
end


%% Save the results
room_dir = fullfile(save_dir,r_type{1});
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

save_name = fullfile(bin_dir,['coch_all_conditions_',r_type{1},'_',pos,'pos_',coch_type,'.mat']);
save(save_name,'coch');

%% Plot the full cochleagrams 
freqs = fliplr(coch(1).params.freqs);
num_freq = length(freqs);
skip_f = 2;
freqs = freqs(1:skip_f:end);
n_tlab = 10;
skip_t = round(length(coch(1).t)/n_tlab);
figure('units','normalized','outerposition',[0 0 1 1]);
row = 3;
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

c = 0;
for r = 1:n_rooms
    for t = 1:n_types
        c = c+1;
        subplot('position',pos{c});
        imagesc(coch(r).(plot_name{t}).X_ft);
        colorbar;
        colormap('inferno');
        yticks([1:skip_f:n_f]);
        yticklabels(y_labels);
        xticks([1:skip_t:length(coch(1).t)]);
        xticklabels(x_labels);
        title([r_type{t},' ',r_name{r}]);
        set(gca,'FontName','Arial','FontSize',15,'FontWeight','Bold');
        xlabel('Time [s]','FontSize',16,'FontWeight','bold');
        ylabel('Freqeuncy [kHz]','FontSize',16,'FontWeight','bold');
    end
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
for r = 1:n_rooms
    for t = 1:n_types
        figure('units','normalized','outerposition',[0 0 1 1]);
        for f = 1:num_freq
            subplot('position',pos{f});
            histogram(coch(r).(plot_name{t}).X_ft(f,:),'Normalization','probability');
            title([f_labels{f},' kHz']);
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
        end
        sgtitle([r_type{t},' ',r_name{r}],'FontSize',18,'Color','red','FontWeight','Bold');
        set(gcf,'color','w');
        save_name = fullfile(bin_dir,[r_type{t},' ',r_name{r},'_fband_value_distribtuion',coch_type,'.png']);
        export_fig(save_name);
        close;
    end
end

