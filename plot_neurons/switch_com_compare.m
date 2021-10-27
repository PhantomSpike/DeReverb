%% Load the switch data
switch_com_name = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/stats/switch_com.mat';
switch_info_name = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/stats/switch_info.mat';
swtch_com = load(switch_com_name);
swtch_info = load(switch_info_name);
select_units = 0;

%% Load the normal data
normal_com_name = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/stats/normal_com.mat';
normal_info_name = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/stats/normal_info.mat';
normal_com = load(normal_com_name);
normal_info = load(normal_info_name);

save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/stats';
%% Select the same animals and pens for the normal data using the switch

if select_units
    switch_animals = unique(swtch_info.info.animal_names);
    switch_pens = unique(swtch_info.info.pen_names);
    
    ix_animals = ismember(normal_info.info.animal_names, switch_animals);
    ix_pens = ismember(normal_info.info.pen_names, switch_pens);
    
    ix_select = ix_animals & ix_pens;
    
    com.normal.neg.small = normal_com.com.neg.small(ix_select);
    com.normal.neg.large = normal_com.com.neg.large(ix_select);
    com.normal.pos.small = normal_com.com.pos.small(ix_select);
    com.normal.pos.large = normal_com.com.pos.large(ix_select);
    
else
    
    com.normal.neg.small = normal_com.com.neg.small;
    com.normal.neg.large = normal_com.com.neg.large;
    com.normal.pos.small = normal_com.com.pos.small;
    com.normal.pos.large = normal_com.com.pos.large;
    
end

%Switch Inhibition
com.switch.neg.s1 = swtch_com.com.neg.s1;
com.switch.neg.l1 = swtch_com.com.neg.l1;
com.switch.neg.s2 = swtch_com.com.neg.s2;
com.switch.neg.l2 = swtch_com.com.neg.l2;

%Switch Excitation
com.switch.pos.s1 = swtch_com.com.pos.s1;
com.switch.pos.l1 = swtch_com.com.pos.l1;
com.switch.pos.s2 = swtch_com.com.pos.s2;
com.switch.pos.l2 = swtch_com.com.pos.l2;

%% Compute stats for the differences
display_results = 'off';
mult_corr = 'hsd';
% 'tukey-kramer' or 'hsd'	
%     Tukey's honest significant difference criterion
% 'bonferroni'	
%     Bonferroni method
% 'dunn-sidak'	
%     Dunn and Sidák’s approach
% 'lsd'	
%     Fisher's least significant difference procedure
% 'scheffe'	
%     Scheffé's S procedure


% - Inhibition COM S1, S2 and Small
%Perform Kruskal-Wallis (KW) test to compare the distributions for 
data_small_com_neg = [com.switch.neg.s1'; com.switch.neg.s2'; com.normal.neg.small']; %Make a vector with the data for KW test
group_small_com_neg = [repmat({'S1'},length(com.switch.neg.s1),1); repmat({'S2'},length(com.switch.neg.s2),1); repmat({'Small Normal'},length(com.normal.neg.small),1)]; %Make a vector with the group labels

[pval_kw.com.neg.small.single, pval_kw.com.neg.small.table, pval_kw.com.neg.small.stats] = kruskalwallis(data_small_com_neg, group_small_com_neg, display_results);
pval_kw.com.neg.small.multiple = multcompare(pval_kw.com.neg.small.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Inhibition COM L1, L2 and Large
%Perform Kruskal-Wallis (KW) test to compare the distributions for 
data_large_com_neg = [com.switch.neg.l1'; com.switch.neg.l2'; com.normal.neg.large']; %Make a vector with the data for KW test
group_large_com_neg = [repmat({'L1'},length(com.switch.neg.l1),1); repmat({'L2'},length(com.switch.neg.l2),1); repmat({'Large Normal'},length(com.normal.neg.large),1)]; %Make a vector with the group labels

[pval_kw.com.neg.large.single, pval_kw.com.neg.large.table, pval_kw.com.neg.large.stats] = kruskalwallis(data_large_com_neg, group_large_com_neg, display_results);
pval_kw.com.neg.large.multiple = multcompare(pval_kw.com.neg.large.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Excitation COM S1, S2 and Small
%Perform Kruskal-Wallis (KW) test to compare the distributions for 
data_small_com_pos = [com.switch.pos.s1'; com.switch.pos.s2'; com.normal.pos.small']; %Make a vector with the data for KW test
group_small_com_pos = [repmat({'S1'},length(com.switch.pos.s1),1); repmat({'S2'},length(com.switch.pos.s2),1); repmat({'Small Normal'},length(com.normal.pos.small),1)]; %Make a vector with the group labels

[pval_kw.com.neg.small.single, pval_kw.com.neg.small.table, pval_kw.com.neg.small.stats] = kruskalwallis(data_small_com_pos, group_small_com_pos, display_results);
pval_kw.com.neg.small.multiple = multcompare(pval_kw.com.neg.small.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

% - Excitation COM L1, L2 and Large
%Perform Kruskal-Wallis (KW) test to compare the distributions for 
data_large_pos_pos = [com.switch.neg.l1'; com.switch.neg.l2'; com.normal.neg.large']; %Make a vector with the data for KW test
group_large_pos_pos = [repmat({'L1'},length(com.switch.neg.l1),1); repmat({'L2'},length(com.switch.neg.l1),1); repmat({'Large Normal'},length(com.normal.neg.large),1)]; %Make a vector with the group labels

[pval_kw.com.pos.large.single, pval_kw.com.pos.large.table, pval_kw.com.pos.large.stats] = kruskalwallis(data_large_pos_pos, group_large_pos_pos, display_results);
pval_kw.com.pos.large.multiple = multcompare(pval_kw.com.pos.large.stats,'estimate', 'kruskalwallis', 'CType', mult_corr, 'display', display_results);

%Compute the median
median_val.switch.neg.s1 = nanmedian(com.switch.neg.s1);
median_val.switch.neg.s2 = nanmedian(com.switch.neg.s2);
median_val.switch.neg.l1 = nanmedian(com.switch.neg.l1);
median_val.switch.neg.l2 = nanmedian(com.switch.neg.l2);
median_val.normal.neg.small = nanmedian(com.normal.neg.small);
median_val.normal.neg.large = nanmedian(com.normal.neg.large);

median_val.switch.pos.s1 = nanmedian(com.switch.pos.s1);
median_val.switch.pos.s2 = nanmedian(com.switch.pos.s2);
median_val.switch.pos.l1 = nanmedian(com.switch.pos.l1);
median_val.switch.pos.l2 = nanmedian(com.switch.pos.l2);
median_val.normal.pos.small = nanmedian(com.normal.pos.small);
median_val.normal.pos.large = nanmedian(com.normal.pos.large);

%Compute the mean
median_val.switch.neg.s1 = nanmean(com.switch.neg.s1);
median_val.switch.neg.s2 = nanmean(com.switch.neg.s2);
median_val.switch.neg.l1 = nanmean(com.switch.neg.l1);
median_val.switch.neg.l2 = nanmean(com.switch.neg.l2);
median_val.normal.neg.small = nanmean(com.normal.neg.small);
median_val.normal.neg.large = nanmean(com.normal.neg.large);

median_val.switch.pos.s1 = nanmean(com.switch.pos.s1);
median_val.switch.pos.s2 = nanmean(com.switch.pos.s2);
median_val.switch.pos.l1 = nanmean(com.switch.pos.l1);
median_val.switch.pos.l2 = nanmean(com.switch.pos.l2);
median_val.normal.pos.small = nanmean(com.normal.pos.small);
median_val.normal.pos.large = nanmean(com.normal.pos.large);

%Compute the std
sem_val.switch.neg.s1 = nanstd(com.switch.neg.s1)/sqrt(length(com.switch.neg.s1));
sem_val.switch.neg.s2 = nanstd(com.switch.neg.s2)/sqrt(length(com.switch.neg.s2));
sem_val.switch.neg.l1 = nanstd(com.switch.neg.l1)/sqrt(length(com.switch.neg.l1));
sem_val.switch.neg.l2 = nanstd(com.switch.neg.l2)/sqrt(length(com.switch.neg.l2));
sem_val.normal.neg.small = nanstd(com.normal.neg.small)/sqrt(length(com.normal.neg.small));
sem_val.normal.neg.large = nanstd(com.normal.neg.large)/sqrt(length(com.normal.neg.large));

sem_val.switch.pos.s1 = nanstd(com.switch.pos.s1)/sqrt(length(com.switch.pos.s1));
sem_val.switch.pos.s2 = nanstd(com.switch.pos.s2)/sqrt(length(com.switch.pos.s2));
sem_val.switch.pos.l1 = nanstd(com.switch.pos.l1)/sqrt(length(com.switch.pos.l1));
sem_val.switch.pos.l2 = nanstd(com.switch.pos.l2)/sqrt(length(com.switch.pos.l2));
sem_val.normal.pos.small = nanstd(com.normal.pos.small)/sqrt(length(com.normal.pos.small));
sem_val.normal.pos.large = nanstd(com.normal.pos.large)/sqrt(length(com.normal.pos.large));

%% Violin plot of COM Inhibition S1, S2 and Small Normal
%Select the colors
exc_small_color = [0.9804 0.4196 0.6431];
exc_large_color = [0.8 0 0];

inh_small_color = [0.0745 0.6235 1.0000];
inh_large_color = [0 0 1];

com_neg_small_violin.S1 = com.switch.neg.s1;
com_neg_small_violin.S2 = com.switch.neg.s2;
com_neg_small_violin.Small = com.normal.neg.small;

font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'S1','S2','Small Normal'};
violins = violinplot(com_neg_small_violin, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = inh_small_color;
violins(2).ViolinColor = inh_small_color;
violins(3).ViolinColor = inh_small_color;
title('COM inhibition small');
ylabel(['COM (ms)']);
ylim([30 180]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin_COM_inhibition_small.svg');
saveas(gcf, save_name);
close;

%% Violin plot of COM Inhibition L1, L2 and Large Normal
com_neg_large_violin.L1 = com.switch.neg.l1;
com_neg_large_violin.L2 = com.switch.neg.l2;
com_neg_large_violin.Large = com.normal.neg.large;

font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'L1','L2','Large Normal'};
violins = violinplot(com_neg_large_violin, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = inh_large_color;
violins(2).ViolinColor = inh_large_color;
violins(3).ViolinColor = inh_large_color;
title('COM inhibition large');
ylabel(['COM (ms)']);
ylim([30 180]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin_COM_inhibition_large.svg');
saveas(gcf, save_name);
close;

%% Violin plot of COM Inhibition S1, S2 and Small Normal
com_pos_small_violin.S1 = com.switch.pos.s1;
com_pos_small_violin.S2 = com.switch.pos.s2;
com_pos_small_violin.Small = com.normal.pos.small;

font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'L1','L2','Large Normal'};
violins = violinplot(com_pos_small_violin, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = exc_small_color;
violins(2).ViolinColor = exc_small_color;
violins(3).ViolinColor = exc_small_color;
title('COM excitation small');
ylabel(['COM (ms)']);
ylim([0 100]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin_COM_excitation_small.svg');
saveas(gcf, save_name);
close;

%% Violin plot of COM Inhibition L1, L2 and Large Normal
com_pos_large_violin.L1 = com.switch.pos.l1;
com_pos_large_violin.L2 = com.switch.pos.l2;
com_pos_large_violin.Large = com.normal.pos.large;

font_sz = 55;
figure('units','normalized','outerposition',[0 0 1 1]);
transp = 0.35;
violin_width = 0.45;
int_names = {'L1','L2','Large Normal'};
violins = violinplot(com_pos_large_violin, int_names, 'ShowMean', true, 'ViolinAlpha', transp, 'Width', violin_width, 'ShowData', false, 'ShowNotches', false);
violins(1).ViolinColor = exc_large_color;
violins(2).ViolinColor = exc_large_color;
violins(3).ViolinColor = exc_large_color;
title('COM excitation large');
ylabel(['COM (ms)']);
ylim([0 100]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.25);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Normal');
save_name = fullfile(save_dir, 'Violin_COM_excitation_large.svg');
saveas(gcf, save_name);
close;