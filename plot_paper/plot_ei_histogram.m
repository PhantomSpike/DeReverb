function plot_ei_histogram(params)
%This function will plot the histogram comparisons between small and large
%rooms for excitation and inhibition components together

%% Unroll params
fit = params.fit;
units = params.units; 
val_exc_big = params.val_exc_big;
val_exc_small = params.val_exc_small;
val_inh_big = params.val_inh_big;
val_inh_small = params.val_inh_small;

p_val_exc = params.p_val_exc;
p_val_inh = params.p_val_inh;

n_clust = length(val_exc_big);
his_spacing = params.his_spacing;
lim_val = params.lim_val;
hist_type = params.hist_type;
model = params.model;

lw1 = params.lw1;
lw2 = params.lw2;
exc_color = params.exc_color;
inh_color = params.inh_color;
all_font_sz = params.all_font_sz;
font_type = params.font_type;
specific_name = params.specific_name;
calc_method = params.calc_method;
save_dir = params.save_dir;

%% Construct histogram
edges = [-(lim_val+1.5*his_spacing):his_spacing:lim_val+his_spacing/2];
edges(1) = [];

switch units
    case 'normal'
        temp_e = val_exc_big - val_exc_small;
        temp_i = val_inh_big - val_inh_small;
    case 'ratio'
        temp_e = log2(val_exc_big./val_exc_small);
        temp_i = log2(val_inh_big./val_inh_small);
    case 'change ix'
        temp_e = 100*((val_exc_big - val_exc_small)./(val_exc_big + val_exc_small));
        temp_i = 100*((val_inh_big - val_inh_small)./(val_inh_big + val_inh_small));
end

counts_e = histcounts(temp_e,edges);
counts_i = histcounts(temp_i,edges);

switch hist_type
    case 'stairs'
        fc_neg = 'none';
        fc_pos = 'none';
    case 'bar'
        fc_pos = exc_color;
        fc_neg = inh_color;
end
%% Plot
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
hold on;
histogram('BinEdges', edges,'BinCounts', counts_e,'FaceColor',fc_pos,'EdgeColor',exc_color,'DisplayStyle',hist_type,'LineWidth', lw2,'Normalization','probability');
histogram('BinEdges', edges,'BinCounts', counts_i,'FaceColor',fc_neg,'EdgeColor',inh_color,'DisplayStyle',hist_type,'LineWidth', lw2,'Normalization','probability');
% xlabel('Total amount large- Total amount small [AU]','FontSize',x_font_sz,'FontWeight','bold');
% ylabel('Number of neurons','FontSize',y_font_sz,'FontWeight','bold');
% title(['Change of E_{Total} and I_{Total} in AU for',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.63 0.85 0.1 0.1],'String', sprintf('Excitatory^{%s}',p_star(p_val_exc)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.63 0.75 0.1 0.1],'String', sprintf('Inhibitory^{%s}',p_star(p_val_inh)),'LineStyle','none','Color','b','FontSize',all_font_sz,'FontWeight','bold');
xline(0,'k','LineWidth',lw1);
yline(0,'k','LineWidth',lw2);
if ismember('COM',specific_name)
    x_vals = [-lim_val:20:lim_val];
    xticks(x_vals);
    for j = 1:length(x_vals)
        x_labels{j} = num2str(x_vals(j),'%.0f');
    end
    xticklabels(x_labels);
end
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
xlim([-(lim_val+his_spacing/2), lim_val+his_spacing/2]);
hold off;
set(gcf,'color','w');
set(gca,'TickDir','out');
switch fit
    case 'neurons'
        get_neurons = params.get_neurons;
        NPSP_th = params.NPSP_th;
        if isfield(params,'CCnorm_th')
            CCnorm_th = params.CCnorm_th;
            save_name = fullfile(save_dir,['Histogram ',units,' ',hist_type,specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),' ',[model,' '],calc_method,' kernel.svg']);
        else
            save_name = fullfile(save_dir,['Histogram ',units,' ',hist_type,specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',[model,' '],calc_method,' kernel.svg']);
        end
    case 'normative'
        save_name = fullfile(save_dir,['Histogram ',units,' ',hist_type,specific_name,[model,' '],calc_method,' kernel.svg']);
end
saveas(gcf,save_name)
close;