%----------------------------------
%Master script to find neutral curves
%----------------------------------
%DESCRIPTION GOES HERE
clear all
close all
syms y alpha yc
global J1 Re dy

%step size - one step also indicates dpos and dneg
dy = 0.05; %make sure yc - dneg falls on grid point

%velocity profile
U(y) = 1-exp(-y);
Uy(y) = diff(U);
Uyy(y) = diff(U,2);
Uyyy(y) = diff(U,3);

%far field condition
psiinf = symfun( exp(-alpha*y), [y, alpha,yc]);
psiinfy = diff(psiinf,y);

%Reynolds number
Re = 10e3;

% Calculate jump for particular value of lambda and N = 1
N = 1;
lambdac = 0.1;
Jn = CriticalLayer(lambdac, N);
%J1 = 0;
J1 = 1i*imag(Jn(N+1,1));

% Tollmien solutions close to critical layer
v1 = symfun(((y-yc) + (Uyy(yc)/(2*Uy(yc)))*(y-yc)^2 + (alpha^2/6 + Uyyy(yc)/(6*Uy(yc)))*(y-yc)^3),[y, alpha, yc]);
v2 = symfun((1 + (alpha^2/2 + Uyyy(yc)/(2*Uy(yc))-(Uyy(yc)/Uy(yc))^2)*(y-yc)^2 + (Uyy(yc)/Uy(yc))*v1(y,alpha,yc)*log(abs(y-yc))), [y,alpha, yc]);
v1y = diff(v1,y);
v2y = diff(v2,y);

% Vector of values of yc and alpha to loop over
ycvec = 0.1:0.1:1; %0.2:0.01:0.3;
alphavec = 0.1:0.1:1;
rrootmat = [ycvec.' zeros(length(ycvec),1)];
irootmat = [ycvec.' zeros(length(ycvec),1)];

% Shooting method for finding specific value of c and alpha
for k = 1:length(ycvec) %make sure this is bigger than dneg
    yc = ycvec(k);
    c = U(yc);
    rDeltaold = 0;
    iDeltaold = 0;
    for a = 1:length(alphavec)
        alpha = alphavec(a);
        fprintf('alpha = %8.4f \t yc = %8.4f \t c = %8.5f \n',alpha,yc,c);
        %Yinf = 5;
        if alpha < 0.3
            Yinf = 10;
        else
            Yinf = 5;
        end
        tic
        [Delta] = finddelta(alpha,yc, v1, v2,Yinf, U, psiinf)
        rDeltanew = real(Delta);
        iDeltanew = imag(Delta);
        if sign(rDeltanew*rDeltaold)==-1
            guess = 1/2*(alphavec(a)+alphavec(a-1));
            rroot = muller1d(guess,yc, v1, v2,Yinf, U, "real",psiinf)
            if rrootmat(k,2) == 0 %this loop ensures no root is overwritten when
                                  % there are two or more roots for one value of yc
                rrootmat(k,2) = rroot;
            else
                rrootmat = [rrootmat zeros(length(ycvec),1)];
                rrootmat(k,end) = rroot;
            end
        end
        if sign(iDeltanew*iDeltaold)==-1
            guess = 1/2*(alphavec(a)+alphavec(a-1));
            iroot = muller1d(guess,yc, v1, v2,Yinf, U, "imag",psiinf)
            if irootmat(k,2) == 0
                irootmat(k,2) = iroot;
            else
                irootmat = [irootmat zeros(length(ycvec),1)];
                irootmat(k,end) = iroot;
            end
        end
        rDeltaold = rDeltanew;
        iDeltaold = iDeltanew;
        toc
    end
end