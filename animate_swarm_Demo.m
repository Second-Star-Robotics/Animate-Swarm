%Execute Cong's script and extract variables
Density_estimation_lastlast

%Assemble Input Data Structure
Input_Data.Density_Plot.T_space = T_space;
Input_Data.Density_Plot.x_esti = x_esti;
Input_Data.Density_Plot.Zeta = Zeta;
Input_Data.t_space = t_space;
Input_Data.Platform(1).z = z(1,:);
Input_Data.Platform(2).z = p_e_max;
Input_Data.Layer.z = p_max;
Input_Data.error = error; %Consider changing the name of this variable as it's a matlab command

%Clean up workspace
close all;
clearvars -except Input_Data

%Set input parameters
platform_stl = '190927 Driftcam.stl';
layer_stl = 'squid.stl';
video_filename = 'test.avi';
speed = 2;
platform_scale = 60;
layer_scale = 75;
smoothing = 15;

animate_result(Input_Data, platform_stl, layer_stl, video_filename, speed, platform_scale, layer_scale, smoothing);

