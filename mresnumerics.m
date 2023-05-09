%----------------------------------
%MRes Numerics
%----------------------------------
clear all
close all
syms y alpha yc

dy = 0.03; %step size
d = 2*dy; % distance above and below critical layer
Re = 10e4; %typical 10^3 or 10^4

U(y) = tanh(y);%1-exp(-y);
Uy(y) = diff(U);
Uyy(y) = diff(U,2);
Uyyy(y) = diff(U,3);
psiinf = symfun( exp(-alpha*y), [y, alpha,yc]);
psiinfy = diff(psiinf,y);

% Calculate jump for particular value of lambda and N = 1
N = 1;
% lambdac = 0.1;
% Jn = CriticalLayer(lambdac, N);
J1 = 0;
%J1 = 1i*imag(Jn(N+1,1));
% phi = 2i*J1;
v1 = symfun(((y-yc) + (Uyy(yc)/(2*Uy(yc)))*(y-yc)^2 + (alpha^2/6 + Uyyy(yc)/(6*Uy(yc)))*(y-yc)^3),[y, alpha, yc]);
v2 = symfun((1 + (alpha^2/2 + Uyyy(yc)/(2*Uy(yc))-(Uyy(yc)/Uy(yc))^2)*(y-yc)^2 + (Uyy(yc)/Uy(yc))*v1(y,alpha,yc)*log(abs(y-yc))), [y,alpha, yc]);
v1y = diff(v1,y);
v2y = diff(v2,y);

for yc = 0.2:0.1:1 %make sure this is bigger than dneg
    if yc/dy ~= floor(yc/dy)
        Nc = floor(yc/dy);
        Nd = 1;
        dneg = yc - (Nc-Nd)*dy;
        dpos = (Nc+1+Nd)*dy-yc;
    else
        dneg = d;
        dpos = d;
    end
    dneg = d;
    dpos = d;
    c = U(yc);
    for alpha = 0.8:0.1:2
        if alpha < 0.3
            Yinf = 30;
        else
            Yinf = 10;
        end
        tic
        % Below critical layer
        [psib1,etab1] = RungeKuttaSolver(alpha,yc, U, yc-dneg, 0, v1, v1y);
        [psib2,etab2] = RungeKuttaSolver(alpha,yc, U, yc-dneg, 0, v2, v2y);
        f1 = psib1 - (1i*alpha*c*Re)^(-1/2)*etab1;
        f2 = psib2 - (1i*alpha*c*Re)^(-1/2)*etab2;
%         f1 = psib1 - (-1i*alpha*c*Re)^(-1/2)*etab1;
%         f2 = psib2 - (-1i*alpha*c*Re)^(-1/2)*etab2;
        
        % Above critical layer
        [psia,etaa] = RungeKuttaSolver(alpha,yc, U, yc+Yinf, yc+dpos, psiinf, psiinfy);
        g1 = (psia/etaa)*v1y(yc+dpos,alpha,yc) - v1(yc+dpos,alpha,yc);
        g2 = (psia/etaa)*v2y(yc+dpos,alpha,yc) - v2(yc+dpos,alpha,yc);
        
        %matr = [Uy(yc)/Uyy(yc) (-Uy(yc)/Uyy(yc)) (-2*J1);0 f1 f2; g1 0 g2];
        %Delta = double(det(matr))
        
        fprintf('alpha = %8.2f \t yc = %8.3f \t c = %8.5f \n',alpha,yc,c);
        Delta = double((-Uy(yc)/Uyy(yc))*(g2/g1-f2/f1)-2*J1)
        toc
    end
end