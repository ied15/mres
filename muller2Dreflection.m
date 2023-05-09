function [root] = muller2Dreflection(guessa,guessyc, v1, v2,Yinf, U,...
    psiinf1, psiinf2, delta)
% muller2D finds the root of delta for alpha and yc, given the roots of 
% im(delta) and re(delta), using a 2D muller method.
% Failure due to A having zero determinant may be due to the rows being
% linearly dependent, in which case adjust the intial three guesses
% for alpha and yc.

%create initial three guesses for alpha
%subtraction and addition of number should be decided by precision of
%knowledge of location of root.
%it needs to be small enough so that the algorithm doesn't diverge
alpha1 = guessa-0.01
alpha2 = guessa+0.01
alpha3 = guessa

%create initial three guesses for yc
yc1 = guessyc-0.01
yc2 = guessyc+0.01
yc3 = guessyc

%create delta values of root guesses
[Aplus1, Aminus1, B1, R1] = findR(alpha1,yc1, v1, v2, Yinf, U,...
    psiinf1, psiinf2, delta)
[Aplus2, Aminus2, B2, R2] = findR(alpha2,yc2, v1, v2, Yinf, U,...
    psiinf1, psiinf2, delta)
[Aplus3, Aminus3, B3, R3] = findR(alpha3,yc3, v1, v2, Yinf, U,...
    psiinf1, psiinf2, delta)

R1hat = 1/R1
R2hat = 1/R2
R3hat = 1/R3
% precision condition could either rely on f3 close enough to zero or
% alpha and yc close enough to previous values of alpha and yc
while abs(real(R3hat)) > 1e-10 || abs(imag(R3hat)) > 1e-10 
    % construct plane passing through three guesses for real(f)
    matA = [alpha1 yc1 1;alpha2 yc2 1;alpha3 yc3 1];
    matB = [real(f1); real(f2); real(f3)];
    
    if det(matA) == 0
        disp('ERROR: rows are linearly dependent')
        break
    end
    
    % obtain coefficients X(1),X(2),X(3) of the plane
    % z = X(1)*alpha + X(2)*yc + X(3)
    X = double(linsolve(matA,matB)) 
    
    % intersect with plane z = 0 and solve for yc
    yc_alpha =   (-X(1)*alpha3-X(3))/X(2);
    
    % use 1d muller method on linearlised function imag(f) of one variable
    disp('beginning 1d muller for alpha')
    tic
    alphanew = muller1d(alpha3,yc_alpha, v1, v2,Yinf, U, "imag", psiinf);
    toc
    disp('ending 1d muller for alpha')
    
    % update values of alpha
    alpha1 = alpha2
    alpha2 = alpha3
    alpha3 = alphanew
    
    % update values of yc using relation found in intersection of z=0
    % plane
    yc1 = yc2
    yc2 = yc3
    yc3 =  (-X(1)*alpha3-X(3))/X(2)
    
    % update values of f
    R1hat = R2hat
    R2hat = R3hat
    [Aplus3, Aminus3, B3, R3] = findR(alpha3,yc3, v1, v2, Yinf, U,...
    psiinf1, psiinf2, delta)
    R3hat = 1/R3
end

root = [alpha3; yc3];
end