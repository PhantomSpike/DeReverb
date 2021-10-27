%% Define params for run
core_globals = '/home/alex/Desktop/Code/roomsim/Standard_settings/roomsim_core_globals_standard.mat';
run_input = '/home/alex/Desktop/Code/roomsim/Standard_settings/roomsim_run_input_standard.mat';
pos_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/Positions/Closed_manypos'; %Directory with the positions
%Different rooms
room_name{1} = 'small_room';
room_name{2} = 'med_room';
room_name{3} = 'big_room';
%Different room types
room_type{1} = 'closed';
room_type{2} = 'anechoic';
order = [17;50;50];


%% Loop over the different rooms and types

for r = 1:length(room_name)
    
    position_file = fullfile(pos_dir,[room_name{r},'/position_arrays_',room_name{r}]);
    
    for t = 1:length(room_type)
        run_many_positions(room_name{r},room_type{t},position_file,core_globals,run_input,order)
    end
    
end