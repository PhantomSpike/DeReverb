%This function loads the IRs from the Roomsim output
%% Define params
root_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/IRs/Chopped_experimental/40pos';
%Different rooms
r_name{1} = 'small_room';
r_name{2} = 'med_room';
r_name{3} = 'big_room';
n_rooms = length(r_name);
%Different room types
r_type{1} = 'closed';
r_type{2} = 'anechoic';
n_types = length(r_type);
%% Loop through all IRs and load them

for r = 1:n_rooms
    
    for t = 1:n_types
        
        cur_dir = fullfile(root_dir,[r_name{r},'/',r_type{t}]); %Make dir name
        dirList = dir(cur_dir); %Get all contents
        dirList = dirList(~ismember({dirList.name},{'.','..'}));
        n_pos = length(dirList); %Get number of positions
        ir_names = {dirList.name}';
        ir_names = natsort(ir_names); %Sort them so their order is 1:n_pos
        ir_folder = dirList(1).folder; %Files are all in the same folder so no need to loop over them
        
        for p = 1:n_pos
            ir_name = fullfile(ir_folder,ir_names{p});
            temp = load(ir_name); 
            time = [1/temp.Fs:1/temp.Fs:length(temp.data)/temp.Fs]'; %Make a time vector for latter plotting
            temp.data = temp.data./max(abs(temp.data(:))); %Normalize the IRs to be between -1 and 1
            ir.(r_name{r}).(r_type{t}){p,1} = temp;
            ir.(r_name{r}).(r_type{t}){p,1}.fname = ir_name;
            ir.(r_name{r}).(r_type{t}){p,1}.time = time;
        end
        
    end
    
end

save_dir = fullfile(root_dir,'All');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

save_name = fullfile(save_dir,'all_irs.mat');
save(save_name,'ir');
%% Plot the IRs to make sure they look sensible
font_sz = 10;
title_sz = 18;
row = 8;
col = 5;
per = 0.03;
edgel = per; edger = 0.02; edgeh = 0.08; edgeb = 0.04; space_h = 0.02; space_v = 0.04;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

save_dir = fullfile(root_dir,'All','Plots');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

for r = 1:n_rooms
    for t = 1:n_types
        figure('units','normalized','outerposition',[0 0 1 1]);
        for p = 1:n_pos
            subplot('position',pos{p});
            plot(ir.(r_name{r}).(r_type{t}){p,1}.time,ir.(r_name{r}).(r_type{t}){p}.data);
            xlim([0 ir.(r_name{r}).(r_type{t}){p,1}.time(end)]);
            if ismember(p,[1:col:col*row])
                ylabel('Amplitude');
            end
            if ismember(p,[(col*row-col+1):(col*row)])
                xlabel('Time [s]');
            end
            title(['pos',num2str(p)]);
            set(gca,'FontName','Arial','FontSize',font_sz,'FontWeight','Bold');
        end
        sgtitle([r_type{t},' ',r_name{r}],'FontSize',title_sz,'Color','red','FontWeight','Bold');
        set(gcf,'color','w');
        save_name = fullfile(save_dir,[r_type{t},' ',r_name{r},'_irs.png']);
        export_fig(save_name);
        close;
        
    end
    
end
            
            
