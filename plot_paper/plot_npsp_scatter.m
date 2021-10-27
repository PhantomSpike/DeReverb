function plot_npsp_scatter(params)
%This function will plot the scatter plot comparisons between small and large
%rooms for excitation and inhibition components together

%% Unroll params
fit = params.fit;
val_big = params.val_big;
val_small = params.val_small;

p_val = params.p_val;

n_clust = length(val_big);
lim_val = params.lim_val;
model = params.model;

npsp_on = params.npsp_on;
lw = params.lw;
all_font_sz = params.all_font_sz;
font_type = params.font_type;
specific_name = params.specific_name;
sz = params.sz;
save_dir = params.save_dir;

%% Plot
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);

if npsp_on
    c = params.NPSPs;
    scatter(val_small,val_big,sz,c,'filled'); hold on;
    colormap('jet');
    colorbar;
else
    scatter(val_small,val_big,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
end
    


% xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
% ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
% title(['Inhibitory \tau_{big} vs \tau_{small}  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.65 0.2 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.86 0.1 0.1],'String', sprintf('p=%s',p_star(p_val)),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
axis equal;
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
set(gca,'TickDir','out');
switch fit
    case 'neurons'
        l = [0 lim_val];
        get_neurons = params.get_neurons;
        NPSP_th = params.NPSP_th;
        save_name = fullfile(save_dir,[specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',model,' kernel.svg']);
    case 'normative'
        l = [0.7 lim_val];
        save_name = fullfile(save_dir,[specific_name,model,' kernel.svg']);
end
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = lw;
hline.LineStyle = '--';
saveas(gcf,save_name)
close;