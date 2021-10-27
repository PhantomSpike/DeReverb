function plot_tconst(kh_fits,params)
plot_kernels = 1;
plot_comparison = 1;
plot_betas = 1;

%Unroll params
r_type = params.r_type;
save_dir_full = params.save_dir_full;
room_gen_name = params.room_gen_name;
model_name = params.model_name;
int_name = params.int_name;
freqs = fliplr(params.freqs)/1000;
n_rooms = params.n_rooms;
n_ker = params.n_ker;
ub_A_B = params.ub_A_B;
%% Plot the raw k_h, fits on top, neg exp, pos exp and params: A, B, alpha, beta, MSE
if plot_kernels
    %Params for global title textbox
    gtpos_x = 0.3; gtpos_y = 0.97; gtsz_x = 0.45; gtsz_y = 0.03;
    %Params for subplots title textbox
    tpos_x = 0.4; tpos_y = 0.8; tsz_x = 0.3; tsz_y = 0.3;
    %Params for A, alpha textbox
    Apos_x = 0.75; Apos_y = 0.6; Asz_x = 0.28; Asz_y = 0.35;
    %Params for B, beta textbox
    Bpos_x = 0.75; Bpos_y = 0.09; Bsz_x = 0.35; Bsz_y = 0.35;
    %Params for mse title textbox
    mpos_x = 0.3; mpos_y = 0.48; msz_x = 0.5; msz_y = 0.3;
    

    
    title_sz = 15;
    axis_sz = 12;
    mark_sz = 17;
    line_sz = 3;
    
    rows = 6;
    cols = 5;
    per = 0.02;
    edgel = 0.035; edger = per; edgeh = 0.05; edgeb = 0.05; space_h = 0.015; space_v =0.05;
    [pos]=subplot_pos(rows,cols,edgel,edger,edgeh,edgeb,space_h,space_v);
    xdiff = (pos{2}(1) - pos{1}(1)) - space_h;
    ydiff = (pos{1}(2) - pos{cols+1}(2)) - space_v;
    
    h = kh_fits.small{1}.h; %Get the history
    
    
    for r = 1:n_rooms
        fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
        for k = 1:n_ker
            kh_fit = kh_fits.(r_type{r}){k};
            A = kh_fit.A;
            B = kh_fit.B;
            alpha = kh_fit.alpha;
            beta = kh_fit.beta;
            k_h = kh_fit.k_h;
            k_h_hat = kh_fit.k_h_hat;
            exp_pos = kh_fit.exp_pos;
            exp_neg = kh_fit.exp_neg;
            best_Loss = kh_fit.best_Loss;
            subplot('position',pos{k});
            hold on;
            plot(h,k_h,'k.','MarkerSize',mark_sz);
            plot(h,k_h_hat,'k','Linewidth',line_sz);
            plot(h,exp_pos,'r','Linewidth',line_sz-1);
            plot(h,exp_neg,'b','Linewidth',line_sz-1);
            %Title text box
            annotation(fig1,'textbox',pos{k} + [xdiff*tpos_x, ydiff*tpos_y -pos{k}(3)*(1-tsz_x) -pos{k}(4)*(1-tsz_y)],...
                'String', sprintf('%.1f kHz',freqs(k)),'LineStyle','none','FontSize',axis_sz,'FontWeight','bold');
            %A, alpha text box
            annotation(fig1,'textbox',pos{k} + [xdiff*Apos_x, ydiff*Apos_y -pos{k}(3)*(1-Asz_x) -pos{k}(4)*(1-Asz_y)],...
                'String', sprintf('A=%.1f\n\x03b1=%.0fms',A,alpha),'LineStyle','none','FontSize',axis_sz,'FontWeight','bold','Color','r');
            %B, beta text box
            annotation(fig1,'textbox',pos{k} + [xdiff*Bpos_x, ydiff*Bpos_y -pos{k}(3)*(1-Bsz_x) -pos{k}(4)*(1-Bsz_y)],...
                'String', sprintf('B=%.1f\n\x03b2=%.0fms',B,beta),'LineStyle','none','FontSize',axis_sz,'FontWeight','bold','Color','b');
            set(gca,'YTickLabel',[]);
            %mse text box
            annotation(fig1,'textbox',pos{k} + [xdiff*mpos_x, ydiff*mpos_y -pos{k}(3)*(1-msz_x) -pos{k}(4)*(1-msz_y)],...
                'String', sprintf('mse=%.3e',best_Loss),'LineStyle','none','FontSize',axis_sz,'FontWeight','bold');
            %         yline(0,'--k','Linewidth',line_sz);
            ylim([-1 1]);
            xlim([0 max(h)]);
            set(gcf,'color','w');
            hold off;
            if ismember(k,[0:rows-1]*cols + 1)
                ylabel('Weights [AU]');
                yticks('auto');
                yticklabels('auto');
            end
            if ismember(k,[rows*cols-cols+1:rows*cols])
                xlabel('History [ms]');
            end
            set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
        end
        %Global Title text box
        annotation(fig1,'textbox',[gtpos_x, gtpos_y gtsz_x gtsz_y],...
            'String', sprintf('%s room %s %s ub{A,B}=%.0f',r_type{r},room_gen_name,model_name,ub_A_B),'LineStyle','none','FontSize',title_sz,'FontWeight','bold');
        %Save the figure
        save_name = fullfile(save_dir_full,[r_type{r},'_room_',room_gen_name,'_',model_name,'_ub_A_B_',num2str(ub_A_B),int_name,'.png']);
        export_fig(save_name);
        close(fig1);
    end
end
%% Plot the raw k_hs from the different rooms on top of one another
if plot_comparison
    %Params for global title textbox
    gtpos_x = 0.3; gtpos_y = 0.97; gtsz_x = 0.45; gtsz_y = 0.03;
    %Params for subplots title textbox
    tpos_x = 0.4; tpos_y = 0.8; tsz_x = 0.3; tsz_y = 0.3;
    %Params for beta textbox
    betapos_x = 0.7; betapos_y = 0.5; betasz_x = 0.35; betasz_y = 0.65;
    
    title_sz = 15;
    axis_sz = 12;
    beta_sz = 10;
    mark_sz = 17;
    line_sz = 3;
    
    
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    for k = 1:n_ker
        kh_small = kh_fits.small{k}.k_h;
        kh_med = kh_fits.med{k}.k_h;
        kh_big = kh_fits.big{k}.k_h;
        beta_small(k) = kh_fits.small{k}.beta;
        beta_med(k) = kh_fits.med{k}.beta;
        beta_big(k) = kh_fits.big{k}.beta;
        
        subplot('position',pos{k});
        hold on;
        plot(h,kh_small,'b','Linewidth',line_sz);
        plot(h,kh_med,'r','Linewidth',line_sz);
        plot(h,kh_big,'k','Linewidth',line_sz);
        %Title text box
        annotation(fig1,'textbox',pos{k} + [xdiff*tpos_x, ydiff*tpos_y -pos{k}(3)*(1-tsz_x) -pos{k}(4)*(1-tsz_y)],...
            'String', sprintf('%.1f kHz',freqs(k)),'LineStyle','none','FontSize',axis_sz,'FontWeight','bold');
        %beta text box
        annotation(fig1,'textbox',pos{k} + [xdiff*betapos_x, ydiff*betapos_y -pos{k}(3)*(1-betasz_x) -pos{k}(4)*(1-betasz_y)],...
            'String', sprintf('\x03b2_{small}=%.0fms\n\x03b2_{med}=%.0fms\n\x03b2_{big}=%.0fms',beta_small(k),beta_med(k),beta_big(k)),'LineStyle','none','FontSize',beta_sz,'FontWeight','bold','Color','k');
        set(gca,'YTickLabel',[]);
        ylim([min([kh_small;kh_med;kh_big]) 1]);
        xlim([0 round(max(h)/2)]);
        set(gcf,'color','w');
        hold off;
        if ismember(k,[0:rows-1]*cols + 1)
            ylabel('Weights [AU]');
            yticks('auto');
            yticklabels('auto');
        end
        if ismember(k,[rows*cols-cols+1:rows*cols])
            xlabel('History [ms]');
        end
        set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
        if k == n_ker
            legend('small','med','big','FontSize',10,'Location','northwest');
        end
    end
    %Global Title text box
    annotation(fig1,'textbox',[gtpos_x, gtpos_y gtsz_x gtsz_y],...
        'String', sprintf('All rooms \x03b2 %s %s ub{A,B}=%.0f',room_gen_name,model_name,ub_A_B),'LineStyle','none','FontSize',title_sz,'FontWeight','bold');
    %Save the figure
    save_name = fullfile(save_dir_full,['All_rooms_betas_',room_gen_name,'_',model_name,'_ub_A_B_',num2str(ub_A_B),int_name,'.png']);
    export_fig(save_name);
    close(fig1);
end
%% Plot the betas vs frequency for all rooms
if plot_betas
    for k = 1:n_ker
        beta_small(k) = kh_fits.small{k}.beta;
        beta_med(k) = kh_fits.med{k}.beta;
        beta_big(k) = kh_fits.big{k}.beta;
    end
    line_sz = 4;
    axis_sz = 25;
    freqs = fliplr(freqs);
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    plot(freqs,flipud(beta_big(:)),'Linewidth',line_sz,'Color','r');
    plot(freqs,flipud(beta_med(:)),'Linewidth',line_sz,'Color',[0.91, 0.41, 0.17]);
    plot(freqs,flipud(beta_small(:)),'Linewidth',line_sz,'Color','b');
    xlabel('Frequency [kHz]');
    ylabel('\beta values [ms]');
    title('\beta values vs freqeuncy for all rooms');
    set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
    set(gcf,'color','w');
    legend('big','med','small');
%     save_name = fullfile(save_dir_full,['All_rooms_beta_vs_freq',room_gen_name,'_',model_name,'_ub_A_B_',num2str(ub_A_B),int_name,'.png']);
    save_name = fullfile(save_dir_full,['All_rooms_beta_vs_freq_ub_A_B_',num2str(ub_A_B),int_name,'.png']);
    export_fig(save_name);
    close(fig1);
end