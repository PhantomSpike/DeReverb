%% Define the primary penetrations
primary_pens.Ronnie = {'P05','P06','P08','P09','P13','P14'}; 
primary_pens.PLP = {'P02'}; 
primary_pens.Derry = {'P01','P05','P06','P08','P09'}; 
primary_pens.Cork = {'P01'}; 
primary_pens.Kilkenny = {'P04','P05'}; 
primary_pens.Noah = {'P01','P02','P03','P04','P06','P08','P09'};
primary_pens.Derekah = {'P01','P04'}; 
animals = fieldnames(primary_pens);

%% Load the info file
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/New_all_neuronal_data_noneuron_norm/perfreq_noneuro/ridge/10ms/200ms';
NPSP_th = 40;
get_neurons = 'all';

load(fullfile(kernel_dir,'info'),'info');

switch get_neurons
    case 'all'
        ix_qualia = ones(length(info.cluster_id),1); %Get all the neurons
    case 'good'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'good'),info.quality,'UniformOutput',false)); %Find all the good units
    case 'mua'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'mua'),info.quality,'UniformOutput',false)); %Find all the mua units
end
ix_npsp = info.NPSP<NPSP_th; %Find all the neurons below certain NPSP
ix = ix_qualia & ix_npsp; %Find the intersection of the two

NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);

%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
qualities = qualities(ix_select);
n_clust = length(cluster_ids);

%% Find proportions of primary neurons in the ferrets

primary_count = 0;
for j = 1:length(animals)
    animal = animals{j};
    ix_animal = ismember(animal_names, animal);
    ix_pens = ismember(pen_names, primary_pens.(animal));
    primary_count = primary_count + sum(ix_animal & ix_pens);
end

per_primary = 100*(primary_count/n_clust);