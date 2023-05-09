%----------------------------------
%MRes Numerics Rossby
%----------------------------------
%A wave impinging on a critical layer is considered with reflection
%coefficient R. Outer solution is governed by the Rayleigh equation, solved
%using a Runge-Kutta method. Jump condition across critical layer is found
%using a finite difference method.

%181114 EDITED so code runs with alpha, omega and not c, before
% it ran with omega, c and not alpha.

clear all
close all
y = sym('y');
alpha = sym('alpha');
yc = sym('yc');
global J1 dy betacor
%Step size
dy = 0.01;
delta = dy;

%Coriolis parameter
betacor = 0.06;

%Velocity Profile
U(y) = 1/2*(1+tanh(y));%tanh(y);%1-exp(-y);
Uy(y) = diff(U);
Uyy(y) = diff(U,2);
Uyyy(y) = diff(U,3);

%far field condition
k = symfun(sqrt(betacor/(1-U(yc))-alpha^2), [alpha,yc]);
psiinf1 = symfun(exp(1i*k(alpha,yc)*y), [y, alpha,yc]);
psiinf2 = symfun(exp(-1i*k(alpha,yc)*y), [y, alpha,yc]);

% Calculate jump for particular value of lambda and N = 1
N = 1;
%lambdac = 0.1;
%Jn = CriticalLayer(lambdac, N);
%J1 = 1i*imag(Jn(N+1,1));
J1 = 1i*pi/2; %viscous

% Tollmien solutions close to critical layer
v1 = symfun(((y-yc) + ((Uyy(yc)-betacor)/(2*Uy(yc)))*(y-yc)^2 ...
    + 1/6*(alpha^2 + Uyyy(yc)/(Uy(yc)) - betacor*(Uyy(yc)-betacor)...
    /(2*Uy(yc)^2))*(y-yc)^3),[y, alpha, yc]);
v2 = symfun((1 + (alpha^2/2 + Uyyy(yc)/(2*Uy(yc))-(4*Uyy(yc)-3*betacor)...
    *(Uyy(yc)-betacor)/(4*Uy(yc)^2))*(y-yc)^2 + ((Uyy(yc)- betacor)...
    /Uy(yc))*v1(y,alpha,yc)*log(abs(y-yc))), [y,alpha, yc]);

% Vector of values of yc to loop over
%ycvec = 0.1:0.1:1; %0.2:0.01:0.3;
omegavec = 0.057:0.00025:0.0625;
alphavec = 0.26:0.0025:0.28;
Rmat = zeros(length(omegavec)*length(alphavec), 4);
realsurfmat = zeros(length(alphavec),length(omegavec));
imagsurfmat = zeros(length(alphavec),length(omegavec));
cmat = zeros(length(alphavec),length(omegavec));
omegamat = zeros(length(alphavec),length(omegavec));
alphamat = zeros(length(alphavec),length(omegavec));

Yinf = 10; 

for j = 1:length(alphavec)
    alpha = alphavec(j);
    for l = 1:length(omegavec)
        omega = omegavec(l);
        c = omega/alpha;
        if omega < alpha % this avoids cases where c greater than 1
        yc =  solve(U(y)==c,y);
            if alpha^2 < betacor/(1-c)
                keval = k(alpha,yc);
                theta = atan(keval/alpha);
                fprintf('omega = %8.8f \t alpha = %8.8f \t yc = %8.8f \t',...
                    omega, alpha,yc)
                fprintf('c = %8.8f \t k = %8.8f \n theta = %8.8f \n',c,...
                    double(keval),theta);
                tic
                [Aplus, Aminus, B, R] = findR(alpha,yc, v1, v2,Yinf, U,...
                    psiinf1, psiinf2,delta) %edited to go from -Yinf
                toc
                Rmat((j-1)*length(alphavec) + l, :) = [c omega theta R];
                realsurfmat(j,l) = real(R);
                imagsurfmat(j,l) = imag(R);
                cmat(j,l) = c;
                omegamat(j,l) = omega;
                alphamat(j,l) = alpha;
            else
                fprintf('omega = %8.2f \t alpha = %8.2f \t yc = %8.3f \t',...
                    omega, alpha,yc);
                fprintf('c = %8.5f \t k = %8.5f \n',c,double(k(alpha,yc)) );
                %disp("alpha > sqrt(betacor/(1-U(yc)))")
                Rmat((j-1)*length(alphavec) + l, :) = [c omega nan nan];
                realsurfmat(j,l) = nan;
                imagsurfmat(j,l) = nan;
                cmat(j,l) = c;
                omegamat(j,l) = omega;
            end
        else 
           realsurfmat(j,l) = nan;
            imagsurfmat(j,l) = nan;
            cmat(j,l) = c;
            omegamat(j,l) = omega;
            alphamat(j,l) = alpha;
        end
    end
end
% figure(1)
% plot(Rmat(:,2), real(Rmat(:,3)),'-g', 'LineWidth',2); hold on
% plot(Rmat(:,2), imag(Rmat(:,3)),'-b','LineWidth',2);
% plot(Rmat(:,2), abs(Rmat(:,3)),'-m','LineWidth',2);
% legend('real(R)','imag(R)','abs(R)');
% xlabel('\theta')
% title(['real(R), imag(R) and abs(R) plotted against incident angle \theta for c = ' num2str(c) ', Yinf = ' num2str(Yinf) ', dy = ' num2str(dy)])
% hold off

%figure(1)
%surf(cmat,omegamat,realsurfmat)

csvwrite('cmat.txt',cmat)
csvwrite('omegamat.txt',omegamat)
csvwrite('alphamat.txt',alphamat)
csvwrite('realsurfmat.txt',realsurfmat)
csvwrite('imagsurfmat.txt',imagsurfmat)