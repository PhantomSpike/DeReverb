%% Define vars
snd_path = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Sound_stim/Derry_Kilkenny_Cork/stim1.wav'; %Absolute path to the kernel used for the prediction
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_1';
[data,fs] = audioread(snd_path);

%% Prepare for plotting
start_s = 0;
end_s = 4.01;
dt_s = 1;
n_t = length(data);
t = [1/fs:1/fs:n_t/fs];
dt = mean(diff(t));
dt_samples = round(dt_s/dt);
t(end:2*length(t)) = [t(end):dt:t(end)*2];
[~,ix_start] = min(abs(t - start_s));
[~,ix_end] = min(abs(t - end_s));
t = t(ix_start:ix_end);
t = t - t(1);
n_steps = round(length(t)/dt_samples)+1;
plot_data = data(ix_start:ix_end);
max_val = max(abs(plot_data));
plot_data = plot_data./max_val;

for tm = 1:n_steps
    x_labels{tm} = num2str(t((tm-1)*dt_samples +1),'%.0f');
end

y_labels{1} = num2str(-1,'%.0f'); y_labels{2} = num2str(0,'%.0f'); y_labels{3} = num2str(1,'%.0f');

%% Plot the waveform
sz = 100;
all_font_sz = sz;
x_font_sz = sz;
y_font_sz = sz;
font_type = 'Liberation Sans';
fig1 = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
plot(plot_data,'k','LineWidth',2);
xticks([1:dt_samples:length(t)]);
xticklabels(x_labels);
yticks([-1,0,1]);
yticklabels(y_labels);
% xlabel('Time [s]','FontSize',x_font_sz,'FontWeight','Normal');
% ylabel('Amplitude','FontSize',y_font_sz,'FontWeight','Normal');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
save_name = fullfile(save_dir,'waveform_example.svg');
saveas(gcf,save_name);
close all;

