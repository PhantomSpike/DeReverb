function psth = make_psth(y, t_edges_s, animal_name)

%% Make PSTH for this cluster
n_stim = length(y.stim);
n_rep = length(y.stim(1).repeat);

for s = 1:n_stim
    psth_temp = [];
    for r = 1:n_rep
        psth_temp(r,:) = histc(y.stim(s).repeat(r).spiketimes,t_edges_s);
    end
    y_temp(s,:) = mean(psth_temp);
end

if ismember(animal_name, ["Ronnie","PLP","Cork","Kilkenny","Derry"])
    psth.small_1 = y_temp(3,:);
    psth.small_2 = y_temp(4,:);
elseif ismember(animal_name, ["Noah","Derekah"])
    psth.small_1 = y_temp(1,:);
    psth.small_2 = y_temp(2,:);
else
    error('Unrecognized animal');
end

psth.large_1 = y_temp(5,:);
psth.large_2 = y_temp(6,:);