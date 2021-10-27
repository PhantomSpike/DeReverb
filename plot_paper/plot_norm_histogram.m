function plot_norm_histogram(params)
%This function will plot the histogram comparisons between small and large
%rooms for one measure

%% Unroll params
fit = params.fit;
units = params.units; 
big_val = params.big_val;
small_val = params.small_val;

p_val = params.p_val;

n_clust = length(big_val);
his_spacing = params.his_spacing;
lim_val = params.lim_val;
hist_type = params.hist_type;
model = params.model;

lw1 = params.lw1;
lw2 = params.lw2;
his_color = params.his_color;
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
        temp_val = big_val - small_val;
    case 'change ix'
        temp_val = 100*((big_val - small_val)./(big_val + small_val));
end

counts = histcounts(temp_val,edges);

switch hist_type
    case 'stairs'
        fc = 'none';
    case 'bar'
        fc = his_color;
end

%% Plot
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
hold on;
histogram('BinEdges', edges,'BinCounts', counts,'FaceColor',fc,'EdgeColor',his_color,'DisplayStyle',hist_type,'LineWidth', lw2,'Normalization','probability');
% xlabel('Total amount large- Total amount small [AU]','FontSize',x_font_sz,'FontWeight','bold');
% ylabel('Number of neurons','FontSize',y_font_sz,'FontWeight','bold');
% title(['Change of E_{Total} and I_{Total} in AU for',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.72 0.85 0.1 0.1],'String', sprintf('p=%s',p_star(p_val)),'LineStyle','none','Color',his_color,'FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
xline(0,'r','LineWidth',lw1);
set(gcf,'color','w');
set(gca,'TickDir','out');
xlim([-lim_val, lim_val]);
hold off;
switch fit
    case 'neurons'
        get_neurons = params.get_neurons;
        NPSP_th = params.NPSP_th;
        save_name = fullfile(save_dir,['Histogram ',units,' ',hist_type,specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',[model,' '],calc_method,' kernel.svg']);
    case 'normative'
        save_name = fullfile(save_dir,['Histogram ',units,' ',hist_type,specific_name,[model,' '],calc_method,' kernel.svg']);  
end
saveas(gcf,save_name)
close;