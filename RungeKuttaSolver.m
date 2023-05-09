function[psi,eta] = RungeKuttaSolver(alpha, yc, betacor, U, initialy, endy, initialv, initialvdiff)
% It calculates ODE using Runge-Kutta 4th order method
% Author Ido Schwartz
% Originally available form: http://www.mathworks.com/matlabcentral/fileexchange/29851-runge-kutta-4th-order-ode/content/Runge_Kutta_4.m
% Edited by Amin A. Mohammed, for 2 ODEs(April 2016)
% Edited 31/07/18 to include rotation
global dy
Uyy = diff(U,2);
c = U(yc);

%dy=0.05;  % step size
h = sign(endy-initialy)*dy;
y = initialy:h:endy;                                         % Calculates upto y(1)
psi = zeros(1,length(y)); 
eta = zeros(1,length(y)); 
psi(1) = initialv(initialy,alpha,yc);                                          % initial condition
eta(1) = initialvdiff(initialy,alpha,yc);                                          % initial condition
% F_xy = @(t,r) 3.*exp(-t)-0.4*r;                  % change the function as you desire
F_xyz = @(y,psi,eta) eta;                                  % change the function as you desire
G_xyz = @(y,psi,eta) (alpha^2 + (Uyy(y)-betacor)/(U(y)-c))*psi;

tic
for i=1:(length(y)-1)  % calculation loop
    i
    k_1 = F_xyz(y(i),psi(i),eta(i));
    L_1 = G_xyz(y(i),psi(i),eta(i));
    k_2 = F_xyz(y(i)+0.5*h,psi(i)+0.5*h*k_1,eta(i)+0.5*h*L_1);
    L_2 = G_xyz(y(i)+0.5*h,psi(i)+0.5*h*k_1,eta(i)+0.5*h*L_1);
    k_3 = F_xyz((y(i)+0.5*h),(psi(i)+0.5*h*k_2),(eta(i)+0.5*h*L_2));
    L_3 = G_xyz((y(i)+0.5*h),(psi(i)+0.5*h*k_2),(eta(i)+0.5*h*L_2));
    k_4 = F_xyz((y(i)+h),(psi(i)+k_3*h),(eta(i)+L_3*h)); % Corrected        
    L_4 = G_xyz((y(i)+h),(psi(i)+k_3*h),(eta(i)+L_3*h));

    psi(i+1) = psi(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    eta(i+1) = eta(i) + (1/6)*(L_1+2*L_2+2*L_3+L_4)*h;  % main equation

end
toc

psi = psi(length(y));
eta = eta(length(y));