%clear all

% addpath('C:\Program Files\MATLAB\casadi-windows-matlabR2016a-v3.4.5');
% import casadi.*;

%addpath('C:\Users\Friends\Documents\MATLAB\PMP on phevs');
%import tables.*;

%lookup tables
% ngen= transpose([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%                 0,0,77,77,78,78,79,79,80,80,80,80,81,82,83,84,85,86,87,88,0;
%                 0,0,83,83,84,84,85,85,86,86,87,87,88,89,89,90,91,91,91,89,0;
%                  0,0,84,85,86,86,87,88,88,89,89,90,90,91,91,92,92,92,92,90,0;
%                  0,0,87,88,88,89,89,90,90,90,91,91,92,92,92,93,93,92,92,90,0;
%                  0,0,88,89,90,90,91,91,91,91,92,92,92,93,93,93,93,93,92,90,0;
%                  0,0,0,90,91,91,92,92,92,92,93,93,93,93,93,94,93,93,92,90,0;
%                  0,0,0,91,91,92,92,92,93,93,93,93,93,94,94,94,94,94,93,91,0;
%                  0,0,0,0,0,91,92,93,93,93,94,94,94,94,94,94,94,94,93,90,0;
%                  0,0,0,0,0,0,0,92,93,93,94,94,94,94,94,94,94,94,93,90,0;
%                  0,0,0,0,0,0,0,0,92,92,93,94,94,94,94,94,94,94,93,90,0;
%                  0,0,0,0,0,0,0,0,0,92,92,93,94,94,94,94,94,94,93,90,0;
%                  0,0,0,0,0,0,0,0,0,0,92,92,93,94,94,94,94,94,93,90,0]);
% tgen= [-200,-190,-180,-170,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0];            
% wgen= [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000];
% 
% lut_ngen = interpolant('LUT','bspline',{tgen, wgen},ngen(:));
% %lut_ngen([-180 500])=77
% 
% nmot= transpose([69,69,69,69,69,69,69,69,69,69,69;
%                  69,70,80,80,80,80,80,80,80,80,80;
%                  69,80,84,84,84,84,84,84,84,84,84;
%                  69,80,86,88,88,86,86,84,0,0,0;
%                  69,80,84,86,86,84,80,0,0,0,0;
%                  69,70,80,84,80,0,0,0,0,0,0;
%                  69,70,70,70,0,0,0,0,0,0,0;
%                  69,70,70,70,0,0,0,0,0,0,0]);
% wmot= [0,250,500,750,1000,1250,1500,1750,2000,2250,2500];
% tmot= [0,25,50,75,100,125,150,175];
% 
% lut_nmot = interpolant('LUT','bspline',{wmot, tmot},nmot(:));
% %lut_nmot([250 50])=80


% Variable Definition
N = 1000; % number of control intervals

opti = casadi.Opti(); % Optimization problem

X = opti.variable(2,N+1); % state trajectory
xi   = X(1,:);
v   = X(2,:);

U = opti.variable(2,N+1);   % control trajectory 
wm   = U(1,:);
twh  = U(2,:);

T = zeros(1,N+1);      % time
tf = 60*60;             %final time secs
dt = tf/N;            %length of a control interval
for k=1:N+1
    T(1,k)=k*dt;
end

v_des = 20*sin(0.5*pi*T/tf); %max velocity 72km/hr

%Objective Function
opti.minimize(-X(1,N+1)); % minimize final soc

%Constants

r_wh = 33e-2; %m
c0 = 105.95; %N
c1 = 0.01; %Ns/m
c2 = 0.4340; %Ns^2/m^2
Jv = 207; %kg m^2
wb = 2685e-3; %m
cogh = 550e-3;%m
g = 9.81; %m/s^2
mv = 1812; %kg

rho = 2.24;
Rd = 2.16;

% Functions

%elec machines-
tm = @(twh) twh/(Rd*(rho+1));
%effm = @(wm,twh) lut_nmot([wm,tm(twh)]); %lookuptable

% effm = opti.variable(1,N+1);
% for k=1:N+1
%     effm(1,k)= @(wm,twh) lut_nmot([wm(k),twh(k)]);
% end

wg = @(v,wm) ((rho+1).*Rd.*v./r_wh-wm)/rho;
tg = @(twh) twh*rho/(Rd.*(rho+1));
%effg = @(v,wm,twh) lut_ngen([-tg(twh),wg(v,wm)]); %lookuptable

% effg = opti.variable(1,N+1);
% for k=1:N+1
%     effg(1,k) = @(v,wm,twh) lut_ngen([-tg(twh(k)),wg(v(k),wm(k))]);
% end

%battery-
% Pbatt = @(v,wm,twh) wm.*tm(twh)./effm + wg(v,wm).*tg(twh)./effg ;
Pbatt = @(v,wm,twh) wm.*tm(twh) + wg(v,wm).*tg(twh);

n_c=1;
Voc = @(xi) (-0.5911*(xi.*xi) + 1.1892*(xi) + 3.242)*96; %volts
Q = 160000; %As
Ro = @(xi) 32/1000*(-0.3501*(xi.*xi.*xi)+0.952*(xi.*xi)-0.9387*(xi)+2.2806); %ohms

I_bat = @(xi,v,wm,twh) (Voc(xi)-sqrt(Voc(xi)^2 - 4*Pbatt(v,wm,twh).*Ro(xi)))./(2*Ro(xi));
p_ech = @(xi,v,wm,twh) Voc(xi).*I_bat(xi,v,wm,twh);

xi_dot= @(xi,v,wm,twh) -n_c*I_bat(xi,v,wm,twh)/Q;
v_dot= @(v,twh) (r_wh/Jv)*(twh - r_wh*(c0 + c1.*v + c2*v.*v));


% Dynamics
f = @(xi,v,wm,twh) [xi_dot(xi,v,wm,twh); v_dot(v,twh)]; % dx/dt = f(x,u)

% Integrator

for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(xi(:,k), v(:,k), wm(:,k), twh(:,k));
   k2 = f(xi(:,k)+dt/2*k1(1), v(:,k)+dt/2*k1(2), wm(:,k), twh(:,k));
   k3 = f(xi(:,k)+dt/2*k2(1), v(:,k)+dt/2*k2(2), wm(:,k), twh(:,k));
   k4 = f(xi(:,k)+dt*k3(1), v(:,k)+dt*k3(2), wm(:,k), twh(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% Boundary Conditions
opti.subject_to(xi(1)==1);  
opti.subject_to(v==v_des);


% Constraints
opti.subject_to(-60000<=Pbatt(v,wm,twh)<=160000);   

%solver
opti.solver('ipopt'); % set numerical backend
p_opts= struct('expand',true);
s_opts= struct('max_iter',1000);
opti.solver('ipopt',p_opts,s_opts);
sol = opti.solve();   % actual solve

%plots

% figure
% hold on
% plot(sol.value(v));
% plot(sol.value(xi));
% plot(sol.value(wm));
% plot(sol.value(twh));
% legend('v','xi','wm','twh')

% plot(T/60,wg(sol.value(v),sol.value(wm)));
% plot(T/60,sol.value(wm));
% plot(T/60,wg(sol.value(v),sol.value(wm))*60/(2*3.14));
% plot(T/60,Pbatt(sol.value(v),sol.value(wm),sol.value(twh)));

