function [root] = muller1d(guess,yc, v1, v2,Yinf, U, realimag,psiinf)
%MULLER METHOD: this script finds the root of Delta using the 
%muller method for a fixed value of yc. Precision condition must be entered
%by the user.

%create initial three guesses
x1 = guess*1.01; %xk-3
x2 = guess*0.99; %xk-2
x3 = guess; %xk-1

%find delta values for root guesses
if realimag == "real"
    f1 = real(finddelta(x1,yc, v1, v2,Yinf, U,psiinf));
    f2 = real(finddelta(x2,yc, v1, v2,Yinf, U,psiinf));
    f3 = real(finddelta(x3,yc, v1, v2,Yinf, U,psiinf));
end
if realimag == "imag"
    f1 = imag(finddelta(x1,yc, v1, v2,Yinf, U,psiinf));
    f2 = imag(finddelta(x2,yc, v1, v2,Yinf, U,psiinf));
    f3 = imag(finddelta(x3,yc, v1, v2,Yinf, U,psiinf));
end


%run muller method using precision condition of 1e-6
% precision condition could either rely on f3 close enough to zero or
% alpha and yc close enough to previous values of alpha and yc
while abs(f3) > 1e-6
    f32 = (f3-f2)/(x3-x2);
    f31 = (f3-f1)/(x3-x1);
    f21 = (f2-f1)/(x2-x1);
    f321 = (f3-f2)/((x3-x2)*(x3-x1))-(f2-f1)/((x2-x1)*(x3-x1));
    omega = f32 + f31 - f21;
    denom = [omega+sqrt(omega^2-4*f3*f321), omega-sqrt(omega^2-4*f3*f321)];
    [maxdenom, indexdenom] = max(abs(denom));
    maxdenom = maxdenom*sign(denom(indexdenom));
    
    x4 = x3-(2*f3)/maxdenom;
    
    %update values of three guesses
    x1 = x2;
    x2 = x3
    x3 = x4
    
    %find new values of f
    f1 = f2;
    f2 = f3;
    if realimag == "real"
        f3 = real(finddelta(x3,yc, v1, v2,Yinf, U,psiinf))
    end
    if realimag == "imag"
        f3 = imag(finddelta(x3,yc, v1, v2,Yinf, U,psiinf))
    end
        
end
root = x3;
end