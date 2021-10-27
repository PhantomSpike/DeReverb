function plot_nb_all_neurons(psth_all, t_edges_ms, save_dir)

transparent = 1;
lw = 2;
axis_sz = 17;
small_col = [0.8902, 0.6353, 0.2510]; 
large_col = 'r';
dt_ms = 10;

if ~exist('clust_per_plot','var') || isempty(clust_per_plot)
    clust_per_plot = 16;
end

n_clust = length(psth_all);
n_groups = ceil(n_clust/clust_per_plot);
last_clust_per_plot = n_clust - clust_per_plot*(n_groups-1);

row = 4;
col = 4;
per = 0.005;
edgel = 0.05; edger = per; edgeh = per; edgeb = 0.08; space_h = 0.012; space_v = 0.025;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
[last_pos]=subplot_pos(row,last_clust_per_plot,edgel,edger,edgeh,edgeb,space_h,space_v);



for g = 1:n_groups
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    counter = (g-1)*clust_per_plot;
    
    if g == n_groups
        pos = last_pos;
        clust_per_plot = last_clust_per_plot;
    end
    
    c = 0;
    for k = 1:clust_per_plot
        c = counter + k;
        
        %Get the psths for this cluster
        psth_small = [psth_all(c).small_1; psth_all(c).small_2];
        psth_large = [psth_all(c).large_1; psth_all(c).large_2];
        psth_small = psth_small*(1000/dt_ms);
        psth_large = psth_large*(1000/dt_ms);
        n_measure = size(psth_small,1);
        
        %Compute mean and std
        psth_small_mean = nanmean(psth_small); psth_large_mean = nanmean(psth_large);
        psth_small_std = nanstd(psth_small); psth_large_std = nanstd(psth_large);
        
        subplot('position',pos{k});
        
        %Plot the mean +/-sem
        shadedErrorBar(t_edges_ms, psth_small_mean, psth_small_std/sqrt(n_measure), {'LineWidth',lw,'Color', small_col}, transparent); hold on;
        shadedErrorBar(t_edges_ms, psth_large_mean, psth_large_std/sqrt(n_measure), {'LineWidth',lw,'Color', large_col}, transparent);
        xline(0,'k-','LineWidth',lw);
        xlim([-20 100])
        hold off;
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        
        if ismember(k,[clust_per_plot-col+1:clust_per_plot])
            xlabel('Time [ms]');
            xticklabels('auto');
        end
        
        if ismember(k,[1:col:clust_per_plot-col+1])
            ylabel('Firing rate [Hz]');
            yticklabels('auto');
        end
        
        set(gca,'TickDir','out');
        set(gca,'FontSize',axis_sz,'FontWeight','Normal');
        set(gcf,'color','w');
    end
    
    
    save_name = fullfile(save_dir,['clusters_group_',num2str(g),'.svg']);
    saveas(gcf, save_name);
    close;
    
end