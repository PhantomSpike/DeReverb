function run_single_room(fname,core_globals,run_input,room_type,positions,F_abs,order)
% A = run_single_room(core_globals)
%This function gives back the absorption table necessary to run roomsim
%depending on the type of simulation we want to run
%
% input params:
%   fname -- the name of your file
%   core_globals -- path to the core_globals parameters
%   run_input -- path to the run_input parameters
%   room_type -- the type of room we are making. Currently there are 3 options:
%       'same' - all the walls have the same absorption specified by F_abs
%       'open' - the front (Ax1) adn back (Ax2) walls are completely anechoic
%           ('open') all the other ones are set to F_abs
%       'anechoic' - all the walls have absorption 1 ('completely anechoic')
%   positions -- a struct consisitng of 3 fields:
%       receiver_xyz - [x;y;z] the Cartesian 3D coordinates of the receiver (head)
%       sensor_xyz - [xl xr; yl yr; zl zr] the Cartesian 3D coordinates of the left and right ears
%       source_xyz - [x;y;z] the Cartesian 3D coordinates of the source (speaker)
%       room_size - the size of the room to be simulated in [Lx Ly Lz] where Lx,
%           Ly and Lz are the x,y and z dimensions of the room in meters
%
% optional:
%   F_abs -- the absorption values for a single wall. Has dimensions [fx1]
%       where f are the frequencies
%   order -- the maximal order of the reflections calculated in the roomsim
%       simulation. More is better but can be computationally very expensive
%       (impossible) depending on the number. The defualt -1 let's the
%       program choose smth sensible.
%
% e.g.
% run_single_room('cool_room',core_globals,run_input,'open',[6 5 2],positions,[0.1;0.1;0.1;0.2;0.2;0.3],30)
%
% Author: Aleksandar Ivanov
% Year: 2019

if ~exist('order','var')
    order = -1; %This let's roomsim choose some sensible values
end

if ~exist('F_abs','var')
    F_abs = [0.04;0.04;0.04;0.05;0.05;0.05]; %These are the abs coefficients for stone
end

roomsim_dir = '/home/alex/Desktop/Code/roomsim/';

%extract parameters
load(core_globals);
receiver_xyz = positions.receiver_xyz;
sensor_xyz = positions.sensor_xyz;
source_xyz = positions.source_xyz;
room_size = positions.room_size;

A = get_absorption(F_abs,room_type); %Get the absorption table based on the room type and abs coefficients

fullname = [room_type,'_',fname];
cd(roomsim_dir); %Change to roomsim directory
tic;
roomsim_run_NH(run_input,fullname,A,order,room_size,receiver_xyz,sensor_xyz,source_xyz); 
fprintf('== Done! This took %0.fs ==\n',toc);

end