%function animate_result(input_script, stl_filename, video_filename, speed, model_scale)
%RENDER_SIM_DRIFTCAM Render Video of Driftcam Simulation 
%function render_sim_Driftcam_instagram(Sim_Output, Model_Params, stl_filename, video_filename, speed, model_scale)
%inputs:
%   Sim_Output = Data structure generated by sim_Driftcam.m
%   Model_Params = Model parameters generated by setup_sim_Driftcam.m
%   stl_filename = String with Driftcam 3D stl model
%   video_filename = String with video filename to create filename.avi 
%   speed = simulation time/video time
%   model_scale = rendering size of model stl
%
%outputs:
%   none (saved to filename.avi)
%
%dependencies:


%EXECUTE Script and set vectors
feval(input_script);



t = [0:.01:2*pi]';
z1 = sin(t)*50+50;
z2 = sin(t+pi/2)*50+50;
%stl_filename = '221028-ful-exterior-asm-simplified-low.stl';
stl_filename = 'Drifter.stl';
video_filename = 'test.avi';
speed = 1;
model_scale = 5;
n_platforms = 2;

fps = 30;
dt = (1/fps)*speed;

%Load STL Model
[Vertices_b,faces,normals,name] = stlRead(stl_filename);

%Platform 1 rendering characteristics
platform(1).color = [1 1 0]; %Yellow
platform(1).object.vertices = Vertices_b*model_scale;
platform(1).object.faces = faces;
platform(1).z = z1;

platform(2).color = [1 0.647 0]; %Orange
platform(2).object.vertices = Vertices_b*model_scale;
platform(2).object.faces = faces;
platform(2).z = z2;

%Extract vectors from Sim_State
n = length(t);

%Determine max and min values for the renders and curves
min_t = min(t);
max_t = max(t);
min_z = min(min(platform(:).z))-10;
max_z = max(max(platform(:).z))+10;
range_z = max_z-min_z;
middle_z = mean([min_z max_z]);
min_z_plot = middle_z-(range_z*1.1)/2;
max_z_plot = middle_z+(range_z*1.1)/2;
range_z_plot = max_z_plot - min_z_plot;
x_range = range_z_plot/4;
y_range = x_range;
min_x_plot = 0-x_range/2;
max_x_plot = 0+x_range/2;
min_y_plot = 0-y_range/2;
max_y_plot = 0+y_range/2;


%Determine max and min values for sizing bar plots
%max_thrust = 1.25*max(max(abs(Thrusts)));
%if (max_thrust==0) max_thrust = 1.25; end
% max_tau2 = 1.25*max(max(abs(tau2)));
% if (max_tau2==0) max_tau2 = 1.25; end
% max_Alpha = 1.25*max(max(abs(V2_dot)));
% if (max_Alpha==0) max_Alpha = 1.25; end
% max_V1 = 1.25*max(max(abs(V1)));
% if (max_V1==0) max_V1 = 0.01; end
% max_Nabla_eng = 1.25*max([max([Control_State.Nabla_eng_command]) max([Sim_State.Nabla_eng])]);
% if (max_Nabla_eng == 0) max_Nabla_eng = 0.0001; end
% min_Nabla_eng = 0;
% 
% min_rho = min([rho_Driftcam]);
% max_rho = max([rho_Driftcam]);
% range_rho = max_rho-min_rho;
% middle_rho = mean([min_rho max_rho]);
% if (range_rho==0)
%     range_rho = 1;
% end
% min_rho_plot = middle_rho-(range_rho*1.1)/2;
% max_rho_plot = middle_rho+(range_rho*1.1)/2;

%Initialize Figure and allow user to size
%figure('units','normalized','outerposition',[0 0 1 1]);
%figure('Position',[0 0 1080 1350]);

%figure('Position',[0 0 1920 1080]);
%set(gcf,'color','w');
%disp('Please position and size plot window for rendering');
%input('Press any Enter to continue');
%disp('Rendering...')

h1 = subplot(3,4,[1 2 5 6]);
h2 = subplot(3,4,3);

%Create output video file
current_directory = cd;
clock_time = clock;
date_time = date;
complete_filename = [current_directory '\Outputs\' date_time '-' num2str(clock_time(4)) num2str(clock_time(5)) '-' video_filename];
vidObj = VideoWriter(complete_filename);
open(vidObj);

frame_time = 0;
index = 1;
while(index<=n)
    if (t(index)>frame_time+dt)
        frame_time = frame_time+dt;
 
        %plot Driftcam body
        %Vertices_i = (inv(squeeze(Sim_State.J1_eta2(index,:,:)))*Vertices_b')';
        %Vertices_i = (inv(Sim_State(index).J1_eta2)*Vertices_b')'.*model_scale;
        %eta1 = Sim_State(index).eta1';
        %eta1 = zeros(1,3); %center vehicle in axes and show only attitude change
        delete(h1); %delete the figure to refresh
        
        hold on;
        platform_index = 0;
        while(platform_index<n_platforms)
            platform_index = platform_index + 1;
            h1 = subplot(5,3,[1 4 7 10]);
            %stlPlot(Vertices_i+eta1,faces,['t = ', num2str(frame_time,'%0.2f')]);
            platform_color = platform(platform_index).color;
            platform_object = platform(platform_index).object;
            z = platform(platform_index).z(index);
            platform_object.vertices = platform_object.vertices + [0,0,z];
            patch((platform_object),'FaceColor',       platform_color, ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
            set(gca,'Zdir','reverse','Ydir','reverse')
            camlight('headlight');
            material('shiny');
            axis equal;
            view(45, 35);
            grid on;
            %axis ([-1.5 1.5 -1.5 1.5 -1 1]);
            %axis ([-3 3 -3 3 -12 12]);
            axis([min_x_plot max_x_plot min_y_plot max_y_plot min_z_plot max_z_plot]);
            title_str = ['Depth = ', num2str(z, '%0.2f'), ' m'];
            %xlabel('x');
            %ylabel('y');
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            zlabel('Depth [m]');
            %title(title_str);
        end

        hold off;

        
%         %plot Depth Track Setpoint
%         %CS_Vertices_i = (inv(squeeze(Sim_State.J1_eta2_setpoint(index,:,:)))*CS_Vertices_b')';
%         subplot(5,3,[2 3 5 6 8 9 11 12]);
%         %keyboard;
%         %plot(t(1:index),eta1(1:index,3),'-b',t(index),eta1(index,3),'dk', 'color', Actual_Color, t(1:index), z_SP(1:index), '--','color', Setpoint_Color,'MarkerSize',20,'LineWidth', 3);
%         plot(t(1:index)/3600,eta1(1:index,3),'-b',t(index),eta1(index,3),'dk', 'color', Platform1_Color, 'LineWidth', 3);
%         hold on;
%         plot(t(index)/3600,eta1(index,3),'dk', 'color', Platform1_Color, 'MarkerSize',20,'LineWidth', 3);
%         plot(t(1:index)/3600, z_SP(1:index), '-','color', Platform2_Color,'LineWidth', 2);
%         hold off;
%         %axis([min_t max_t min_z_plot max_z_plot]);
%         axis([min_t/3600 max_t/3600 0 280]);
%         grid;
%         set(gca,'Ydir','reverse')
%         xlabel('Time [hr]');
%         ylabel('Depth [m]');
%         t_str = ['t = ', num2str(frame_time,'%0.0f'), ' s'];
%         title(t_str);
%         
%         hold on
%         I = imread('Diel_Migration_With_Predator_Pray_Interaction.jpg'); 
%         h = image(xlim,ylim,I); 
%         uistack(h,'bottom')
%         hold off
% 
%         
%         %plot Velocity
%         subplot(5,3,13);
%         V1 = Sim_State(index).V1;
%         V_dot = Sim_State(index).V_dot;
%         accel = V_dot(3);
%         Control_State_V1 = [Control_State.V1];
%         V1_command = Control_State_V1(index).command';
%         b = bar([accel*1e4 V1(3)*100]);
%         b(1).FaceColor = Platform1_Color;
%         %b(2).FaceColor = Setpoint_Color;  
%         title('Dynamics');
%         axis([0.5, 2.5, -max_V1*100, max_V1*100]);
%         grid on;
%         set(gca,'Ydir','reverse');
%         u_tick = ['u = ',num2str(V1(2)), ' [m/s]'];
%         v_tick = ['v = ',num2str(V1(2)), ' [m/s]'];
%         accel_tick = ['accel = ', num2str(accel*1e6, '%0.1f'), ' \mum/s^2'];
%         w_tick = ['vel = ',num2str(V1(3)*100, '%0.2f'), ' cm/s'];
%         %xticklabels({u_tick v_tick w_tick});
%         xticklabels({accel_tick w_tick});
% 
%         %plot Engine Volume
%         subplot(5,3,14);
%         Nabla_eng = Sim_State(index).Nabla_eng; %[L]
%         sps = Sim_State(index).sps/12812.3;
%         Control_State_Nabla_eng = [Control_State.Nabla_eng];
%         %Nabla_eng_command = Control_State_Nabla_eng(index).command; %[L]
%         b = bar([abs(sps)*0.1 Nabla_eng]);
%         b.FaceColor = [0 0 0];       
%         title('Buoyancy Engine');
%         axis([0.5, 2.5, 0, 1.5]);
%         grid on;
%         xtick_volume = ['Volume = ',num2str(Nabla_eng, '%0.2f'),' L'];
%         xtick_pumpspeed = ['Speed = ',num2str(abs(sps), '%0.2f'),' rpm'];
%         xticklabels({xtick_pumpspeed xtick_volume}); 
%               
%         %plot Environment
%         subplot(5,3,15);
%         b = bar([rho(index) Ta(index)]);
%         b.FaceColor = [0 0 0];        
%         title('Environment');
%         axis([0.5, 2.5, 0, 30]);
%         grid on;
%         xtick_seawater = ['Density = ', num2str(rho(index), '%0.3f'), ' kg/L'];
%         xtick_temp = ['T = ', num2str(Ta(index), '%0.1f'), ' ', char(176), 'C'];
%         xticklabels({xtick_seawater, xtick_temp});

        drawnow;

        writeVideo(vidObj,getframe(gcf)); %Write frame to video file
        hold off
    end
       
    index = index+1;
end

close(vidObj); %close video file
close; %close the rendering figure

%UNUSED PLOTS

        %plot forward vector components
        %subplot(3,4,11);
        %bar(Sim_Output.Forward_i(index,:));
        %title('Forward Unit Vector');
        %axis([0.5, 3.5, -1.25 1.25]);
        %grid on;
        %xticklabels({'x' 'y' 'z'});

        %plot forward vector
        %subplot(3,4,12);
        %FW = [Sim_Output.Forward_i(index,:); [0,0,0]];
        %plot3(FW(:,1), FW(:,2), FW(:,3),'-o');
        %axis ([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
        %grid on;
        %xlabel('X');
        %ylabel('Y');
        %zlabel('Z');
        %title('Forward Vector');
        %grid on;
        
        %plot rotational acceleration vector components
        %subplot(3,4,9);
        %bar(Sim_Output.Alpha(index,:));
        %title('Angular Accel [rad/s^2]');
        %title('9');
        %axis([0.5, 3.5, -max_Alpha, max_Alpha]);
        %grid on;
        %xticklabels({'\alphax' '\alphay' '\alphaz'});
        
        %plot time derivative of phi, theta, and psi angles
        %subplot(3,4,11);
        %bar([d_phi_dt, d_theta_dt, d_psi_dt]);
        %bar(Sim_Output.d_EA_dt(index,:));
        %title('dEulerAngle/dt [rad/s]');
        %title('11');
        %axis([0.5, 3.5, -max_d_EA_dt max_d_EA_dt]);
        %grid on;
        %xticklabels({'d\phi/dt' 'd\theta/dt' 'd\psi/dt'});

        
        
