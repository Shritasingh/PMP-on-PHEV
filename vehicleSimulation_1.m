%From Source
r_wh = 33e-2; %m
c0 = 105.95; %N
c1 = 0.01; %Ns/m
c2 = 0.4340; %Ns^2/m^2
Jv = 207; %kg m^2
wb = 2685e-3; %m
cogh = 550e-3;%m
g = 9.81; %m/s^2
mv = 1812; %kg

%Set by me
alpha = 0; %no slope (let)
twheel = 100; %Nm
tbrake = 0; %matlab?

n = 1;
dt = 0.1;
v = [0];
a = [fun(v(n))];
t = 0;
tf = 1000; %s (Let)
x = [0];

while( t < tf)
    v_new = v(n) + a(n)*dt;
    v = [v v_new];
    t = t + dt;
    n = n + 1;
    x = [x t];
    a = [a fun(v(n))];
end

figure(1)
p1 = plot(x, v);
figure(2)
p2 = plot(x, a, 'g');

function v_dot = fun(v)
    r_wh = 33e-2; %m
    c0 = 105.95; %N
    c1 = 0.01; %Ns/m
    c2 = 0.4340; %Ns^2/m^2
    Jv = 207; %kg m^2
    wb = 2685e-3; %m
    cogh = 550e-3;%m
    g = 9.81; %m/s^2
    mv = 1812; %kg
    
    %Set by me
    alpha = 0; %no slope (let)
    twheel = 100; %Nm
    tbrake = 0; %matlab?

    v_dot = (r_wh/Jv)*(tbrake + twheel - mv*r_wh*g*sin(deg2rad(alpha)) - r_wh*(c0 + c1*v + c2*v*v));
end   
