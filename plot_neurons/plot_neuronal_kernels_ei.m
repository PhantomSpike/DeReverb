function plot_neuronal_kernels_ei(kernels,save_dir,ker_per_plot)

freqs = fliplr(kernels{1}.kernel.freqs);
n_f = length(freqs);
n_h = kernels{1}.kernel.n_h+1;
dt = round(kernels{1}.kernel.dt_ms); 
h_max_ms = kernels{1}.kernel.h_max_ms+dt;
h_steps = [0:dt:h_max_ms];
h_lim_ms = h_max_ms;
h_lim = round(h_lim_ms/dt);
axis_sz = 55;
skip_f = 8;
skip_h = 10;

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    f_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for h = 1:skip_h:n_h+1
    count = count+1;
    h_labels{count} = num2str(h_steps(h),'%.0f');
end

if ~exist('ker_per_plot','var') || isempty(ker_per_plot)
    ker_per_plot = 5;
end

model = kernels{1}.kernel.model;
n_cond = 2;
r_type{1} = 'big';r_type{2} = 'small';
n_clust = length(kernels);
n_groups = ceil(n_clust/ker_per_plot);
last_ker_per_plot = n_clust - ker_per_plot*(n_groups-1);

row = n_cond;
col = ker_per_plot;
per = 0.03;
edgel = 0.097; edger = 0.01; edgeh = 0.096; edgeb = 0.098; space_h = 0.02; space_v = 0.065;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
[last_pos]=subplot_pos(row,last_ker_per_plot,edgel,edger,edgeh,edgeb,space_h,space_v);



for g = 1:n_groups
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    counter = (g-1)*ker_per_plot;
    
    if g == n_groups
        pos = last_pos;
        ker_per_plot = last_ker_per_plot;
    end
    
    
    for r = 1:n_cond
        c = 0;
        for k = 1:ker_per_plot
            c = counter + k;
            ii = (r-1)*ker_per_plot + k;
            subplot('position',pos{ii});
            switch model
                case {'sep','sep_kh'}
                    k_fh = kernels{c}.kernel.(r_type{r}).k_f*fliplr(kernels{c}.kernel.(r_type{r}).k_h');
                case {'lasso','ridge','elastic'}
                    k_fh = fliplr(kernels{c}.kernel.(r_type{r}).main{end}.k_fh);
            end
            k_fh = k_fh./max(abs(k_fh(:)));
            k_fh = [k_fh, zeros(n_f,1)];
            imagesc(k_fh);
            caxis([-1 1]);
            colormap('redblue');
            yticks([1:skip_f:n_f]);
            xticks([0:skip_h:h_lim]);
            set(gca,'Xticklabel',[]);
            set(gca,'Yticklabel',[]);
            set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Normal');
            set(gcf,'color','w');
            set(gca,'TickDir','out');
            
            if r == n_cond && k==1
                xticklabels(h_labels);
                a = get(gca,'XTickLabel');
                set(gca,'XTickLabel',a,'fontsize',axis_sz,'FontWeight','Normal')
                yticklabels(f_labels);
                c = get(gca,'YTickLabel');
                set(gca,'YTickLabel',c,'fontsize',axis_sz,'FontWeight','Normal')
            end
            xlim([0 h_lim]);
        end
    end
    
    save_name = fullfile(save_dir,['kernels_group_',num2str(g),'.svg']);
    saveas(gcf, save_name);
    close all;
end
     
end