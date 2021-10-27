function A = get_absorption(F_abs,room_type)
% A = get_absorption(F_abs,room_type)
%This function gives back the absorption table necessary to run roomsim
%depending on the type of simulation we want to run
%
% input params:
% F_abs -- the absorption values for a single wall. Has dimensions [fx1]
% where f are the frequencies
% room_type -- the type of room we are making. Currently there are 3 options:
% 'same' - all the walls have the same absorption specified by F_abs
% 'open' - the front (Ax1) & back (Ax2) walls are completely anechoic
% ('open'), all the other ones are set to F_abs
% 'anechoic' - all the walls have absorption 1 ('completely anechoic')
%
% output:
% A -- abosrption table of size [fxnumber of walls(6)]
%
% e.g.
% A = get_absorption(F_abs,'open')
%
% Author: Aleksandar Ivanov
% Year: 2019

switch room_type
    case 'closed'
        [Ax1,Ax2,Ay1,Ay2,Az1,Az2] = deal(F_abs); %Make all six surfaces have the same absorption
    case 'open'
        [Ax1,Ax2] = deal(ones(6,1)); %Make the front and back fo the room open (anechoic)
        [Ay1,Ay2,Az1,Az2] = deal(F_abs); %Make all other four surfaces have the same absorption
    case 'anechoic'
        [Ax1,Ax2,Ay1,Ay2,Az1,Az2]  = deal(ones(6,6)); %Make all six surfaces anechoic
end

A = [Ax1 Ax2 Ay1 Ay2 Az1 Az2]; %Combine all six surfaces together
end