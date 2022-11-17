%Execute Cong's script and extract variables
Density_estimation_lastlast

Input_Data.T_space = T_space;
Input_Data.x_esti = x_esti;
Input_Data.Zeta = Zeta;
Input_Data.t_space = t_space;
Input_Data.platform(1).z = z(1,:);
Input_Data.platform(2).z = p_e_max;
Input_Data.layer.z = p_max;

%MOVE THIS DOWN AND EXTRACT INPUT VARIABLES FIRST

figure(1)
plot(t_space,tr_C,'LineWidth',2);
xlabel('Time (hr)','FontSize',20)
ylabel('trace(C)','FontSize', 40,'Interpreter','latex')
title('Measurement updated every 10 min','FontSize',20)
figure(2)
% plot(t_space,mean(Rho_error,1),'LineWidth',2);
plot(t_space(1:k-1),error(1:k-1),'LineWidth',2);
xlabel('Time (hr)','FontSize',20)
ylabel('Error','FontSize', 40,'Interpreter','latex')
title('Measurement updated every 10 min','FontSize',20)
keyboard;
%

close all;
clearvars -except Input_Data

%Extract Data
T_space = Input_Data.T_space;
x_esti = Input_Data.x_esti;
Zeta = Input_Data.Zeta;
t_space = Input_Data.t_space;
z1 = Input_Data.platform(1).z;
z2 = Input_Data.platform(2).z;
z_layer = Input_Data.layer.z;

index = 1000;

%Plot Density
hold on
h = surf(T_space,Zeta,(x_esti));
view(0,90)
shading interp
colorbar
co = colorbar;
co.Label.String = 'Density';
%co.FontSize = 20;

yticks([-800 -700 -600 -500 -400 -300 -200 -100 0]);
yticklabels({'800' '700' '600' '500' '400' '300' '200' '100' '0'})
ylabel('Depth [m]');
xlabel('Time [hrs]');

p2 = plot3(t_space(1:index),z1(1:index),250*ones(1,length(t_space(1:index))),'LineStyle','--','Color','r');
p1 = plot3(t_space(1:index),z_layer(1:index),250*ones(1,length(t_space(1:index))),'k','LineWidth',2,LineStyle='--');
p3 = plot3(t_space(1:index),z2(1:index),250*ones(1,length(t_space(1:index))),'b','LineWidth',2);