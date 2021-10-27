function fig = plot_ind_room(kernels,params,room_sz)
%This function will plot all the f-band kernels for a given room size

%% Plotting params
freqs = fliplr(params.freqs);
n_f = length(freqs);
n_h = params.n_h;
dt = params.dt_ms; 
model_type = params.model;
coch_type = params.coch_type; 
normz = params.normz;
h_max_ms = params.h_max_ms;
save_dir = params.save_dir;
h_steps = [0:dt:h_max_ms];
l_width = 2;
axis_sz = 10;
title_sz = 15;
skip_f = 4;
skip_h = 4;

for f = 1:n_f
    f_labels{f} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    y_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for h = 1:skip_h:n_h
    count = count+1;
    h_labels{count} = num2str(h_steps(h),'%.0f');
end

row = 6;
col = 5;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.1; edgeb = per; space_h = per; space_v = 0.05;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

%% STRF - k_fh
figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:n_f
    switch model_type
        case 'sep'
            k_fh = kernels{k}.k_f*flipud(kernels{k}.k_h)'; %Make the full kernel for sep case
        case {'ridge','lasso','manual_alpha','elastic'}
            k_fh = fliplr(kernels{k}.main{11}.k_fh); %Take the full kernel for the elent case
    end
    subplot('position',pos{k});
    maxabsk = max(abs(k_fh(:)));
    k_fh = k_fh./maxabsk;
    imagesc(k_fh,[-1 1]);
    colormap('redblue');
    title([f_labels{k},' kHz']);
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_h:h_max_ms]);
    xticklabels(h_labels);
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gcf,'color','w');
    if k == n_f
        colorbar;
    end
end
sgtitle([coch_type,' ',model_type,' model norm ',normz,' ',room_sz,' room'],'FontSize',title_sz,'Color','k','FontWeight','Bold');
save_name = fullfile(params.save_dir,strjoin({params.coch_type,params.model,params.normz,room_sz,'kfh.png'},'_'));
export_fig(save_name);
close all;
%% Temporal part of kernel - k_h
fig(2) = figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:n_f
    subplot('position',pos{k});
    hold on;
    switch model_type
        case 'sep'
            k_h = flipud(kernels{k}.k_h); %Make the full kernel for sep case
            maxabsk = max(abs(k_h(:)));
            k_h = k_h./maxabsk;
            plot(k_h,'LineWidth',l_width);
        case {'ridge','lasso','manual_alpha','elastic'}
            k_fh = fliplr(kernels{k}.main{11}.k_fh); %Take the full kernel for the elent case
            k_fh_neg = abs(min(k_fh,0));
            k_fh_pos = abs(max(k_fh,0));
            k_h_neg = mean(k_fh_neg);
            k_h_pos = mean(k_fh_pos);
    end
    k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
    k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
    plot(-k_h_neg,'b','LineWidth',l_width);
    plot(k_h_pos,'r','LineWidth',l_width);
    if k == 1
        legend('k_h inhibitory','k_h excitatory');
    end
    
    yline(0,'--','LineWidth',1,'Color','k');
    xticks([1:skip_h:h_max_ms]);
    xticklabels(h_labels);
    title([f_labels{k},' kHz']);
    set(gcf,'color','w');
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    hold off;
end
sgtitle([coch_type,' ',model_type,' model norm ',normz,' ',room_sz,' room'],'FontSize',title_sz,'Color','k','FontWeight','Bold');
save_name = fullfile(params.save_dir,strjoin({params.coch_type,params.model,params.normz,room_sz,'kh.png'},'_'));
export_fig(save_name);
close all;