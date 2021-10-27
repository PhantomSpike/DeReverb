function make_positions(room)
%This function makes the positions for different runs of roomsim given
%certain room parameters and constraints 
%% Define parameters
pos_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/Positions';
n_pos = 40; %The number of positions
d_ear = 0.014; %The distance between the center of the ferret head and each ear in meters
d_wall = 0.1; %The minimum distance from the wall in meters
d_min = 0.1; %The minimum distance between the source and the ferret in meters
Lx = 3; %The x-dimension of the smallest room in meters
Ly = 0.3; %The y-dimension of the smallest room in meters
Lz = 0.3; %The z-dimension of the smallest room in meters

switch room
    case 'small'
        rng(7);
        scale = 1;
    case 'med'
        rng(8);
        scale = 2.5;
    case 'big'
        rng(9);
        scale = 5;
end

Lx = Lx*scale;
Ly = Ly*scale;
Lz = Lz*scale;

save_dir = fullfile(pos_dir,[room,'_room']);
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% Pick values for the sound source position and the head together making sure they stay off the wells and not closer than d_min
for n = 1:n_pos
    %The speaker
    source_xyz(1,n) = (Lx-2*d_wall)*rand + d_wall;
    source_xyz(2,n) = (Ly-2*d_wall)*rand + d_wall;
    source_xyz(3,n) = (Lz-2*d_wall)*rand + d_wall;
    dis = 0;
    while dis < d_min
        %The ferret head
        receiver_xyz(1,n) = (Lx-2*d_wall)*rand + d_wall;
        receiver_xyz(2,n) = (Ly-2*d_wall)*rand + d_wall;
        receiver_xyz(3,n) = (Lz-2*d_wall)*rand + d_wall;
        dis = norm(source_xyz(:,n) - receiver_xyz(:,n),2);
    end
end

%% Pick values for the two different ears
v = 1:2:(2*n_pos);
for n = 1:n_pos
    theta(n) = 2*pi*rand; %The angle between the ear and the center of the head
    % For the left ear
    sensor_xyz(1,v(n)) = receiver_xyz(1,n) + d_ear*sin(theta(n));
    sensor_xyz(2,v(n)) = receiver_xyz(2,n) + d_ear*cos(theta(n));
    sensor_xyz(3,v(n)) = receiver_xyz(3,n);
    
    % For the right ear
    sensor_xyz(1,v(n)+1) = receiver_xyz(1,n) - d_ear*sin(theta(n));
    sensor_xyz(2,v(n)+1) = receiver_xyz(2,n) - d_ear*cos(theta(n));
    sensor_xyz(3,v(n)+1) = receiver_xyz(3,n);
end
%% Save the positions
params.n_pos = n_pos;
params.d_min = d_min;  
params.d_ear = d_ear; 
params.d_wall = d_wall; 
params.theta = theta;
room_size = [Lx;Ly;Lz];

var_name = ['position_arrays_',room,'_room.mat'];
fullname = fullfile(save_dir,var_name);
save(fullname,'receiver_xyz','sensor_xyz','source_xyz','room_size','params');
%% Plot the resulting positions and distances
axis_sz = 18;
head_sz = 200;
ear_sz = 25;
fig_1 = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
scatter3(receiver_xyz(1,:),receiver_xyz(2,:),receiver_xyz(3,:),head_sz,'k','filled');
scatter3(source_xyz(1,:),source_xyz(2,:),source_xyz(3,:),head_sz,'m','filled');
scatter3(sensor_xyz(1,v),sensor_xyz(2,v),sensor_xyz(3,v),ear_sz,'b','filled');
scatter3(sensor_xyz(1,v+1),sensor_xyz(2,v+1),sensor_xyz(3,v+1),ear_sz,'r','filled');
xlabel('Lx');
ylabel('Ly');
zlabel('Lz');
title(['Speaker and ferret positions in ',room,' room with dimensions = [',num2str(Lx),' ',num2str(Ly),' ',num2str(Lz),']m']);
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
legend('ferret head','speaker','left ear','right ear');
xlim([0 Lx]);
ylim([0 Ly]);
zlim([0 Lz]);
axis equal;
hold off;
fig_name = fullfile(save_dir,[room,'_room_features']);
savefig(fig_1,fig_name);
fig_name = [fig_name,'.png'];
export_fig(fig_name,fig_1);

fig_2 = figure('units','normalized','outerposition',[0 0 1 1]);
abs_dist = vecnorm(receiver_xyz-source_xyz,2);
plot(abs_dist,'Linewidth',2);
title(['Distance between speaker and ferret head in ',room,' room with dimensions = [',num2str(Lx),' ',num2str(Ly),' ',num2str(Lz),']m']);
xlabel('Position #');
ylabel('Distance');
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
yline(d_min,'--','Linewidth',2);
legend('Positions',['min distance = ',num2str(d_min),'m']);
fig_name = fullfile(save_dir,[room,'_room_distances.png']);
export_fig(fig_name,fig_2);
close(fig_2);
end