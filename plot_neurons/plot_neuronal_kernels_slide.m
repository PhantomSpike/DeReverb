function plot_neuronal_kernels_slide(kernels,save_dir,ker_per_plot)

freqs = fliplr(kernels{1}.kernel{1}.freqs);
n_f = length(freqs);
n_h = kernels{1}.kernel{1}.n_h;
dt = round(kernels{1}.kernel{1}.dt_ms); 
h_max_ms = kernels{1}.kernel{1}.h_max_ms;
h_steps = [0:dt:h_max_ms];
axis_sz = 10;
skip_f = 4;
skip_h = 6;

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    f_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for h = 1:skip_h:n_h
    count = count+1;
    h_labels{count} = num2str(h_steps(h),'%.0f');
end

if ~exist('ker_per_plot','var') || isempty(ker_per_plot)
    ker_per_plot = 10;
end

model = kernels{1}.kernel{1}.model;

title_name{1} = 'SL 1'; title_name{2} = 'SL middle'; title_name{3} = 'SL end';
title_name{4} = 'LS 1'; title_name{5} = 'LS middle'; title_name{6} = 'LS end';
r_type{1} = 'sl'; r_type{2} = 'ls';
n_rooms = length(r_type);
ix = [1,2,3];
n_ix = length(ix);
n_cond = n_rooms*n_ix;
n_clust = length(kernels);
n_groups = ceil(n_clust/ker_per_plot);
last_ker_per_plot = n_clust - ker_per_plot*(n_groups-1);

row = ker_per_plot;
col = n_cond;
per = 0.005;
edgel = 0.035; edger = per; edgeh = 0.07; edgeb = 0.05; space_h = 0.022; space_v = 0.05;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
[last_pos]=subplot_pos(last_ker_per_plot,col,edgel,edger,edgeh,edgeb,space_h,space_v);



for g = 1:n_groups
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    counter = (g-1)*ker_per_plot;
    
    if g == n_groups
        pos = last_pos;
        ker_per_plot = last_ker_per_plot;
    end
    c = 0;
    ii = 0;
    for k = 1:ker_per_plot
        
        c = counter + k;
        for r = 1:n_rooms
            
            for i = 1:n_ix
                
                curr_ix = ix(i);
                ii = ii+1;
                subplot('position',pos{ii});
                switch model
                    case {'sep','sep_kh'}
                        k_fh = kernels{c}.kernel.(r_type{r}).k_f*fliplr(kernels{c}.kernel.(r_type{r}).k_h');
                    case {'lasso','ridge','elastic'}
                        k_fh = fliplr(kernels{c}.kernel{curr_ix}.(r_type{r}).main{end}.k_fh);
                end
                k_fh = k_fh./max(abs(k_fh(:)));
                imagesc(k_fh);
                caxis([-1 1]);
                colormap('redblue');
                xticks([1:skip_h:h_max_ms]);
                xticklabels(h_labels);
                set(gca,'YTickLabel',[]);
                if ismember(ii,[(ker_per_plot-1)*n_cond+1:ker_per_plot*n_cond])
                    xlabel('History [ms]');
                end
                if ismember(ii,[1:n_cond:(ker_per_plot-1)*n_cond+1])
                    ylabel('Freqeuncy [kHz]');
                    yticks([1:skip_f:n_f]);
                    yticklabels(f_labels);
                end
                if ismember(ii,[1:n_cond])
                    title(title_name{ii});
                end
                set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
                set(gcf,'color','w');
            end
        end
    end
    
    save_name = fullfile(save_dir,['kernels_group_',num2str(g),'.png']);
    export_fig(save_name);
    close all;
end
     
end