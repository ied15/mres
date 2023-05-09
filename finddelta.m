function [Delta] = finddelta(alpha,yc, v1, v2,Yinf, U, psiinf)
%finddelta: script used to search for neutral curve. Solves above/below
%critical layer using RungeKutta method. Also uses jump condition from
%inside critical layer, and forces determinant of resulting matrix to be
%zero to get non-trivial solution.
global Re J1 dy
syms y

%need to establish yc and dpos,dneg fall on grid points
if rem(yc,dy) == 0
    dneg = dy;%/alpha;
    dpos = dy;%/alpha;
else
    dneg = dy + rem(yc,dy);
    dpos = dy + (dy-rem(yc,dy));
end

Uy(y) = diff(U);
Uyy(y) = diff(U,2);
v1y = diff(v1,y);
v2y = diff(v2,y);
psiinfy = diff(psiinf,y);
c = U(yc);

% % Solve below critical layer using Runge Kutta
% [psib1,etab1] = RungeKuttaSolver(alpha,yc, U, yc-dneg, 0, v1, v1y);
% [psib2,etab2] = RungeKuttaSolver(alpha,yc, U, yc-dneg, 0, v2, v2y);
% % Boundary condition below critical layer
% f1 = psib1 + (-1i*alpha*c*Re)^(-1/2)*etab1;
% f2 = psib2 + (-1i*alpha*c*Re)^(-1/2)*etab2;

% Solve below critical layer using Runge Kutta
[Phib1] = RungeKutta(alpha, yc, U, v1, v1y,...
yc-dneg, 0, 1);
[Phib2] = RungeKutta(alpha, yc, U, v2, v2y,...
yc-dneg, 0, 1);
% Boundary condition below critical layer
f1 = Phib1(1) + (-1i*alpha*c*Re)^(-1/2)*Phib1(2);
f2 = Phib2(1) + (-1i*alpha*c*Re)^(-1/2)*Phib2(2);

% %Solve above critical layer using Runge Kutta
% [psia,etaa] = RungeKuttaSolver(alpha,yc, U, Yinf, yc+dpos, psiinf, psiinfy);
% % Boundary condition above critical layer
% zeta = etaa/psia;
% g1 = v1(yc+dpos,alpha,yc)*zeta-v1y(yc+dpos,alpha,yc);
% g2 = v2(yc+dpos,alpha,yc)*zeta-v2y(yc+dpos,alpha,yc);

%Solve above critical layer using Runge Kutta
[Phia] = RungeKutta(alpha, yc, U, psiinf, psiinfy,...
Yinf, yc+dpos, 1);
% Boundary condition above critical layer
zeta = Phia(2)/Phia(1);
g1 = v1(yc+dpos,alpha,yc)*zeta-v1y(yc+dpos,alpha,yc);
g2 = v2(yc+dpos,alpha,yc)*zeta-v2y(yc+dpos,alpha,yc);

matr = [0 f1 f2; g1 0 g2; Uy(yc)/Uyy(yc) (-Uy(yc)/Uyy(yc)) (-2*J1)];
Delta = double(det(matr));

end

