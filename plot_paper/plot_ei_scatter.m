function plot_ei_scatter(params)
%This function will plot the scatter plot comparisons between small and large
%rooms for excitation and inhibition components together

%% Unroll params
fit = params.fit;
val_exc_big = params.val_exc_big;
val_exc_small = params.val_exc_small;
val_inh_big = params.val_inh_big;
val_inh_small = params.val_inh_small;

p_val_exc = params.p_val_exc;
p_val_inh = params.p_val_inh;

n_clust = length(val_exc_big);
lim_val = params.lim_val;
model = params.model;

lw = params.lw;
exc_color = params.exc_color;
inh_color = params.inh_color;
all_font_sz = params.all_font_sz;
font_type = params.font_type;
specific_name = params.specific_name;
sz = params.sz;
calc_method = params.calc_method;
save_dir = params.save_dir;

%% Plot
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(val_exc_small, val_exc_big, sz,'filled','MarkerEdgeColor',exc_color,'MarkerFaceColor',exc_color); hold on;
scatter(val_inh_small, val_inh_big, sz,'filled','MarkerEdgeColor',inh_color,'MarkerFaceColor',inh_color);

% xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
% ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
% title(['Inhibitory \tau_{big} vs \tau_{small}  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);

% annotation('textbox',[0.65 0.2 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.77 0.86 0.1 0.1],'String', sprintf('excitatory^{%s}',p_star(p_val_exc)),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.77 0.78 0.1 0.1],'String', sprintf('inhibitory^{%s}',p_star(p_val_inh)),'LineStyle','none','Color','b','FontSize',all_font_sz,'FontWeight','bold');
axis equal;
l = [0 lim_val];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = lw;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Normal');
set(gcf,'color','w');
set(gca,'TickDir','out');
switch fit
    case 'neurons'
        get_neurons = params.get_neurons;
        NPSP_th = params.NPSP_th;
        if isfield(params,'CCnorm_th')
            CCnorm_th = params.CCnorm_th;
            save_name = fullfile(save_dir,[specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),'_CCnorm>',num2str(CCnorm_th),' ',[model,' '],calc_method,' kernel.svg']);
        else
            save_name = fullfile(save_dir,[specific_name,get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',[model,' '],calc_method,' kernel.svg']);
        end
    case 'normative'
        save_name = fullfile(save_dir,[specific_name,[model,' '],calc_method,' kernel.svg']);
end
saveas(gcf,save_name)
close;