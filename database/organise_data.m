 %This functions puts all the data in a single folder in the format
%where each file_name  = animal_pen_clsuterid
%Inside there will be the spiketimes organized by stimulus and trial i.e.
%proper cut-and-paste of the date. There will also be other critical info
%like NPSP, good/bad, channel number
%% Define params
NPSP_th = 200; %The threshold for selecting the neurons
type = 'switch'; %The type of stimulus, normal or switch
main_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/sahani_done'; %The folder with all the animals
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Switching_stimuli/For_analysis/All_data'; %Where we will store all the data
%% Load and proccess data
cluster_id_total = [];
NPSP_total = [];
quality_total = [];
animal_name_total = [];
pen_name_total = [];

animal_list = dir(main_dir);
animal_list = animal_list(~ismember({animal_list.name},{'.','..'}));

for animal = 1:numel(animal_list)
    animal_name = animal_list(animal).name;
    fprintf('== Processing animal %s ==\n',animal_name);
    pen_dir = dir([fullfile(animal_list(animal).folder,animal_list(animal).name),'/**/*ap.bin']);
    for pen = 1:length(pen_dir)
        [a,~] = fileparts(pen_dir(pen).folder);
        [~,pen_name] = fileparts(a);
        
        %Load the clust_info
        sorted_dir = pen_dir(pen).folder;
        fprintf('== Penetration %s ==\n',pen_name);
        load(fullfile(sorted_dir,'clust_info'),'clust_info');
        
        NPSP_curr = reach(clust_info.psth_sahani,'NPSP')'; %Get the current NPSP values
        cluster_keep_ix = NPSP_curr<NPSP_th; %Select only clusters that meet this criteria
        
        NPSP_all = NPSP_curr(cluster_keep_ix); %Select only NPSP < NPSP_th
        cluster_id_all = reach(clust_info.psth_sahani,'cluster_id')'; %Get all cluster_ids
        cluster_id_keep = cluster_id_all(cluster_keep_ix); %Select only those with good NPSP
        quality_all = {clust_info.clust_quality{cluster_keep_ix}}';
        n_clust = length(cluster_id_keep);
        
        if ~isempty(cluster_id_keep)
            switch type
                case 'normal'
                    data_all = prepare_trial_data(clust_info,sorted_dir,cluster_id_keep); %Get the spiketimes for every cluster slices according to stim and trial
                case 'switch'
                    data_all = prepare_trial_data_switch(clust_info,sorted_dir,cluster_id_keep); %Get the spiketimes for every cluster slices according to stim and trial
            end
        end
        
        %Save each cluster as a separate file
        for c = 1:n_clust
            data = data_all{c};
            data.params.NPSP = NPSP_all(c);
            data.params.cluster_id = cluster_id_keep(c);
            data.params.quality = quality_all{c};
            data.params.animal_name = animal_name;
            data.params.pen_name = pen_name;
            cluster_name = [animal_name,'_',pen_name,'_',num2str(data.params.cluster_id)];
            save_name = fullfile(save_dir,cluster_name);
            save(save_name,'data');
        end
        
        %Concatenate everything together
        NPSP_total = [NPSP_total;NPSP_all];
        cluster_id_total = [cluster_id_total;cluster_id_keep];
        quality_total = [quality_total;quality_all];
        animal_name_curr = []; [animal_name_curr{1:n_clust,1}] = deal(animal_name);
        pen_name_curr = []; [pen_name_curr{1:n_clust,1}] = deal(pen_name);
        animal_name_total = [animal_name_total;animal_name_curr];
        pen_name_total = [pen_name_total;pen_name_curr];
    end
end

%Save all the necessary info to index the clusters later
info.cluster_id = cluster_id_total;
info.NPSP = NPSP_total;
info.quality = quality_total;
info.animal_name = animal_name_total;
info.pen_name = pen_name_total;
save(fullfile(save_dir,'info'),'info');