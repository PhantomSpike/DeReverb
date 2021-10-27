function run_many_positions(room_name,room_type,position_file,core_globals,run_input,order,F_abs)
% run_many_positions(room_name,room_type,position_file,core_globals,run_input,order,F_abs)
%This function gives back the absorption table necessary to run roomsim
%depending on the type of simulation we want to run
%
% input params:
%   room_name -- the name of the room
%   room_type -- the type of room we are making. Currently there are 3 options:
%       'same' - all the walls have the same absorption specified by F_abs
%       'open' - the front (Ax1) adn back (Ax2) walls are completely anechoic
%           ('open') all the other ones are set to F_abs
%       'anechoic' - all the walls have absorption 1 ('completely anechoic')
%   position_file -- path to the position file that contains the necessary
%       coordinates and room size


%
% optional:
%   core_globals -- path to the core_globals parameters
%   run_input -- path to the run_input parameters
%   F_abs -- the absorption values for a single wall. Has dimensions [fx1]
%       where f are the frequencies
%   order -- the maximal order of the reflections calculated in the roomsim
%       simulation. More is better but can be computationally very expensive
%       (impossible) depending on the number. The defualt -1 let's the
%       program choose smth sensible.
%
% e.g.
% run_many_positions('Big_Room','open','folder/positions.mat','folder3/my_globals','folder2/my_input',30,[0.1;0.1;0.1;0.2;0.2;0.3])
%
% Author: Aleksandar Ivanov
% Year: 2019

%% Define default params
if ~exist('core_globals','var')
    core_globals = '/home/alex/Desktop/Code/roomsim/Standard_settings/roomsim_core_globals_standard.mat';
end

if ~exist('run_input','var')
    run_input = '/home/alex/Desktop/Code/roomsim/Standard_settings/roomsim_run_input_standard.mat';
end

if ~exist('order','var')
    order = -1; %This let's roomsim choose some sensible values
end

if ~exist('F_abs','var')
    F_abs = [0.04;0.04;0.04;0.05;0.05;0.05]; %These are the abs coefficients for stone
end

%% Get the positions
load(position_file,'params','receiver_xyz','sensor_xyz','source_xyz','room_size'); %This will load all the necessary positions
positions.room_size = room_size; %Only the room size stays constant. The positions change

%% Loop over the positions
n_pos = params.n_pos;
fprintf('== Running %s %s ==\n',room_type,room_name);tic;
for p = 1:n_pos
    fprintf('== Running position %0.f/%0.f ==\n',p,n_pos);
    %Select the current head, sensors and source positions
    positions.receiver_xyz = receiver_xyz(:,p);
    ix_left = 2*(p-1) + 1;
    ix_right = ix_left + 1;
    positions.sensor_xyz = sensor_xyz(:,ix_left:ix_right);
    positions.source_xyz = source_xyz(:,p);
    %Run roomsim 
    fname = [room_name,'_','pos',num2str(p)];
    run_single_room(fname,core_globals,run_input,room_type,positions,F_abs,order);
end
fprintf('== Done! This took %0.fs ==\n',toc);
end