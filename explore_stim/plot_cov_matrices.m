coch_file = '/mnt/40086D4C086D41D0/Reverb_normative/Cochleagrams/900s_exp_400_19k_specpower/closed/single/chopped/10ms/coch_all_conditions_closed_singlepos_specpower.mat';
load(coch_file);
rooms{1} = 'anech';
rooms{2} = 'small';
rooms{3} = 'med';
rooms{4} = 'big';
n_rooms = length(rooms);

for r = 1:length(rooms)
    room = rooms{r};
    X_ft.(room) = coch(r).X_ft;
    cov_mat.(room) = cov(X_ft.(room)');
    corr_mat.(room) = corr(X_ft.(room)');
end


row = 2;
col = 4;
per = 0.01;
edgel = per; edger = per; edgeh = per; edgeb = per; space_h = per; space_v = 0.005;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

max_val_cov = max([cov_mat.anech(:);cov_mat.small(:);cov_mat.med(:);cov_mat.big(:);]);
min_val_cov = min([cov_mat.anech(:);cov_mat.small(:);cov_mat.med(:);cov_mat.big(:);]);
for r = 1:n_rooms
    room = rooms{r};
    subplot('position',pos{r});
    imagesc(cov_mat.(room));
    title(['Cov mat ',room]);
    axis equal;
    axis tight;
    caxis([min_val_cov max_val_cov]);
    colorbar;
    subplot('position',pos{col+r});
    imagesc(corr_mat.(room));
    title(['Corr mat ',room]);
    colormap('inferno');
    caxis([-1 1]);
    axis equal;
    axis tight;
    colorbar;
    set(gcf,'color','w');
end
