clc
clear all

tot_Time = 48;             % total time 24 Hours
freq = 60;                   % time increment between two flow measurements and position updates
dt = 1/freq;

iter = tot_Time*freq;       % total number of iterations

t_space = 0:dt:tot_Time;          % time vector 
k = 1;                      % number of iterations of flow measurements
%% number of ensemble
N = 100;

%% state scattering layer

L = 100; % discretization of Ocean depth
omega = 1/12*pi; % frequency of scattering layer: 24 Hours

zeta = linspace(-800,-1,L)'; % full depth of ocean

%% coefficients of initial density distribution
A_pdf = 1000;
B_pdf = 0.005;
mu_z = -450;
sigma = 200;


%% PF operator related coefficients
D = zeros(length(t_space),1);
Delta = zeros(L,length(t_space));
%% information about Driftcam
beta = 1; % motion frequency of Driftcam
M = 2; % number of Driftcams
phi = 0;
% phi = 0;
z(1,:) = -400+400*sin(beta*t_space+phi); % position of Driftcam 1 depth
z(2,:) = -400+400*cos(beta/2*t_space+pi/4); % position of Driftcam 2 depth
%% initial density at fixed depth
rho = zeros(L,length(t_space));
rho(:,1) = B_pdf+A_pdf*normpdf(zeta,mu_z,sigma);
%% information about noise
a_rho = 0.00001; % process noise
a_y = 0.2; % observation noise

H = zeros(M,L); % output matrix
%%estimation
y = zeros(M,length(t_space));
x_tilde = zeros(L,N,length(t_space)); % perdict 
x_hat = zeros(L,N,length(t_space)); % update

x_hat(:,:,1) =  10*A_pdf*normpdf(zeta,randi([-500,500],1,N),randi([100,300],1,N))+0.0*randn(L,N);

E = mean(x_hat(:,:,1),2); 

A = x_hat(:,:,1)-E;
C(:,:,1) = A*A'/(N-1); % initial emsemble covariance

tr_C = zeros(1,length(t_space)); % trace of C
det_C = zeros(1,length(t_space));
non_zero_eigenC = zeros(1,length(t_space));
tr_C(1) = trace(C(:,:,1)); 
det_C(1) = det(C(:,:,1));
non_zero_eigenC(1) = numel(find(~eig(C(:,:,1))));
rho_error = zeros(L,length(t_space));
error = zeros(1,length(t_space));
%% peak points
p_max = zeros(1,length(t_space));
p_e_max = zeros(1,length(t_space)); % log of the depth of density peak
%% control switch
flag_1 = 2*ones(1,length(t_space));
flag_2 = zeros(1,length(t_space));

while k <= iter

t = k*dt; % time update
rho_error(:,k) = rho(:,k)-mean(x_hat(:,:,k),2); % error of the true density and its estimation (average for ensembles)
error(k) = 1/sqrt(L)*sqrt(rho_error(:,k)'*rho_error(:,k)); % root-squre-mean-error
[M1,I1] = max(rho(:,k));
[M3,I3] = max(mean(x_hat(:,:,k),2)); % find the peak
p_max(k) = zeta(I1);
p_e_max(k) = zeta(I3);


%% feedback controller
if tr_C(k) < 20 && M == 3 && z(1,k)-(p_e_max(k)+40*sin(pi*t)) < 0 && t>10
        z(1,k+1) = z(1,k)+2;
        flag_1(k+1) = 1;
elseif tr_C(k) < 20 && M == 3 && z(1,k)-(p_e_max(k)+40*sin(pi*t)) > 0 && t>10
        z(1,k+1) = z(1,k)-2;
        flag_1(k+1) = 1;
elseif tr_C(k) >= 20 && M == 3
        flag_1(k+1) = 0;
     if flag_1(k+1) ~= flag_1(k)
         k_1 = k;
     end
        z(1,k+1) = -400+400*sin(beta*(k-k_1)*dt+asin(z(1,k_1)/400+1));
end
%% propagation of density at fixed depth

D(k+1) = (0.55+0.45*cos(omega*t))/(0.55+0.45*cos(omega*(t+dt))); % matrix determinant part of PF operator 
z_1 = zeta./(0.55+0.45*cos(omega*(t+dt)));
z_2 = zeta./(0.55+0.45*cos(omega*t));
Delta(:,k+1) = D(k+1)*(B_pdf+A_pdf*normpdf(z_1,mu_z,sigma))./(B_pdf+A_pdf*normpdf(z_2,mu_z,sigma));
rho(:,k+1) = Delta(:,k+1).*rho(:,k);


%% EnKF forecast
   w_f = a_rho*randn(L,N);
   w_fs = zeros(size(w_f));
   for j_f = 1:N
    w_fs(:,j_f) = conv(w_f(:,j_f),ones(1,90)/10,'same'); % smoothing
   end
   x_tilde(:,:,k+1) = Delta(:,k+1).*x_hat(:,:,k)+w_f;

%% EnKF update

[M2,I2] = min(abs(z(:,k)'-zeta)); % measurement is collected at the nearest depth from that of each Driftcam
SBM = eye(L);
for i = 1:M
    
        H(i,:) = SBM(I2(i),:);

end
            E = mean(x_tilde(:,:,k+1),2); 
            A = x_tilde(:,:,k+1)-E;
            C(:,:,k+1) = 1^2*(A*A')/(N-1);
            tr_C(k+1) = trace(C(:,:,k+1));
            det_C(k+1) = det(C(:,:,k+1));
            non_zero_eigenC(k+1) = numel(find(eig(C(:,:,k+1))<0.0001));
        if mod(k,30) == 0 
            

            v = a_y*randn(M,1);            
            y(:,k+1) = H*rho(:,k+1)+v;
            v1 = a_y*randn(M,N);
            K = C(:,:,k+1)*H'*(H*C(:,:,k+1)*H'+a_y^2*eye(M))^(-1);


            x_hat(:,:,k+1) = x_tilde(:,:,k+1)+K*(y(:,k+1)+v1-H*x_tilde(:,:,k+1));
        else
        x_hat(:,:,k+1) = x_tilde(:,:,k+1);

        end
k = k+1;

end
T_space = t_space.*ones(L,length(t_space));
Zeta =  zeta.*ones(L,length(t_space));
x_esti = reshape(mean(x_hat,2),[L,length(t_space)]);
kk = find(mod(1:k,20)==0);

figure(1)
% meshz(T_space,Zeta,x_esti);
h = surf(T_space,Zeta,x_esti);
% h = surf(T_space,Zeta,rho);
set(h,'LineStyle','none')
view(0,90)
shading interp
colorbar
hold on
% q3 = plot3(t_space,z(2,:),100*ones(1,length(t_space)),'LineStyle','--','Color','k');
 plot3(t_space,z(1,:),100*ones(1,length(t_space)),'LineStyle','--','Color','k');
scatter3(t_space(kk),z(1,kk),100*ones(1,length(kk)),'r','filled');
% scatter3(t_space(kk),z(2,kk),100*ones(1,length(kk)),'r','filled');
q1 = plot3(t_space,p_max,50*ones(1,length(T_space)),'w','LineWidth',2);
xlabel('Time (hr)','FontSize',20)
ylabel('Depth (m)','FontSize',20)
axis([0 48 -800 0])
ax = gca; % current axes
ax.FontSize = 20;
