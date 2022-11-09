% tot_Time = 48;             % total time 24 Hours
% freq = 60;                   % time increment between two flow measurements and position updates
% dt = 1/freq;
% 
% iter = tot_Time*freq;       % total number of iterations
% 
% t_space = 0:dt:tot_Time;          % time vector 
% L = 100;
% % error_1 = sqrt(1/L*diag(rho_error_1'*rho_error_1));
% % error_2 = sqrt(1/L*diag(rho_error_2'*rho_error_2));
% % error_2_C = sqrt(1/L*diag(rho_error_2_C'*rho_error_2_C));
% figure(1)
% 
% p0 = plot(t_space,error_1./error_1,'LineWidth',2);
% hold on
% p1 = plot(t_space,error_2./error_1,'LineWidth',2);
% p2 = plot(t_space,error_2_C./error_1,'LineWidth',2);
% xlabel('Time (hr)','FontSize',20)
% ylabel('$\%\sqrt{1/n_o\sum_{n_o}(\bar{\hat{x}}_i-x_i)^2}$','FontSize', 40,'Interpreter','latex')
% % title('Estimation Error','FontSize',40)
% legend([p0 p1 p2],{'Single Dirftcam','Two Driftcam','Two Driftcam, Closed-loop'},'FontSize',15)
% ax = gca; % current axes
% ax.FontSize = 20;
% axis([0 48 -0.05 1.5])
% grid on;
% 
% figure(2)
% q0 = plot(t_space,tr_C_1./tr_C_1,'LineWidth',2);
% hold on
% q1 = plot(t_space,tr_C_2./tr_C_1,'LineWidth',2);
% q2 = plot(t_space,tr_C_2_C./tr_C_1,'LineWidth',2);
% xlabel('Time (hr)','FontSize',20)
% ylabel('%trace(C)','FontSize', 40,'Interpreter','latex')
% % ylabel('$tr[\hat{e}\hat{e}^T]$','FontSize', 40,'Interpreter','latex')
% axis([0 48 0 1.2])
% legend([q0 q1 q2],{'Single Dirftcam','Two Driftcam','Two Driftcam, Closed-loop'},'FontSize',15)
% % title('Estimation Error covariance trace','FontSize',20)
% ax = gca; % current axes
% ax.FontSize = 20;
% grid on;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% T_space_outer1 = t_space.*ones(length(zeta_outer1(:,1)),length(t_space));
% T_space_outer2 = t_space.*ones(length(zeta_outer2(:,1)),length(t_space));
% x_hat_mean = reshape(mean(x_hat,2),[L,length(t_space)]);
% [M3,I3] = max(x_hat_mean,[],1);
% for i1 = 1:length(T_space)
% p_e_max(i1) = zeta(I3(i1),i1);
% end
% 
% t_m = t_space.*ones(M,length(T_space));
% k_m = find(z_measure);
% t_m = t_m(k_m);
figure(3)
h = surf(T_space,Zeta,(x_esti));
% h = meshz(T_space,Zeta,(x_esti));
% h = surf(T_space,Zeta,(rho));
% h = surf(T_space,Zeta,rho);
% set(h,'LineStyle','none')
view(0,90)
shading interp
colorbar
co = colorbar;
co.Label.String = 'Density';
co.FontSize = 20;
hold on

% h_outer1 = surf(T_space_outer1,zeta_outer1,rho_outer1);
% h_outer2 = surf(T_space_outer2,zeta_outer2,rho_outer2);
% set(h_outer1,'LineStyle','none')
% set(h_outer2,'LineStyle','none')
% p0 = scatter3(t_m,z_measure(k_m),25*ones(1,length(k_m)),'red','filled');
p0 = scatter3(t_space(kk),z(1,kk),250*ones(1,length(kk)),'r','filled');

% scatter3(t_space(kk),z(2,kk),250*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(3,kk),25*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(4,kk),25*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(5,kk),25*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(6,kk),25*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(7,kk),25*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(8,kk),25*ones(1,length(kk)),'r','filled');

p2 = plot3(t_space,z(1,:),250*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(2,:),250*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(3,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(4,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(5,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(6,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(7,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
% plot3(t_space,z(8,:),25*ones(1,length(t_space)),'LineStyle','--','Color','r');
p1 = plot3(t_space,p_max,250*ones(1,length(T_space)),'w','LineWidth',2,LineStyle='--');
p3 = plot3(t_space,p_e_max,250*ones(1,length(T_space)),'w','LineWidth',2);
xlabel('Time (hr)','FontSize',20)
ylabel('Depth (m)','FontSize',20)
axis([0 48 -800 0])
ax = gca; % current axes
ax.FontSize = 20;
lgnd = legend([p0 p1 p2 p3],{'Measurement points','True density peak','Driftcam trajectory','Estimated density peak'},'FontSize',15);
set(lgnd,'color','none');
