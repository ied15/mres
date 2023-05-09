%--------------------------------------
% Solve for Omega_N, eliminate n = 0
%--------------------------------------
clear all
close all
format long

% Grid setup inside crit layer
Y = 40; %we choose Y from 40 to 80
dyc = 0.08; %we choose a uniform step size from 0.02 to 0.1
y = -Y:dyc:Y;
yc = (y(1:end-1)+y(2:end))/2; %cell centered approximation
M = length(yc);


%...truncate omega
N = 10;
Z1 = -N:1:-1;
Z2 = 1:1:N;
Z = [Z1 Z2]; %creates vector excluding N = 0...
P = 2*N;

lmin = 0;
lmax = 15;
lstep = 0.05;
lvec = lmin:lstep:lmax;
lambdacmat = zeros(length(lvec),1);


%vector to be populated with jump for each lambda
jump = zeros(length(lvec),1);

for l = 1:length(lvec)
    lambda = lvec(l);
    lambdacmat(l,1) = lambda;
    %...first term
    A = zeros(M*P,M*P);
    for n = 1:P
        for j = 1:M
            A((n-1)*M + j, (n-1)*M + j) = 1i*Z(n)*yc(j);
        end
    end
    A = sparse(A);
    
    
    
    %...second term
    BB = (-lambda)*spdiags(ones(M,1)*[1 -2 1]/dyc/dyc...
        ,-1:1,M,M);
    BB(1,1) = -3*(-lambda)/dyc^2;
    
    BB(M,M) = -3*(-lambda)/dyc^2;
    B = kron(speye(P),BB);
    
    
    
    %...third term
    CC1 = (1i/(4*dyc))*spdiags(ones(M,1)*[1 0 -1]...
        ,-1:1,M,M);
    CC1(1,1) = -1i/(4*dyc);
    CC1(M,M) = 1i/(4*dyc);
    d1 = spdiags(ones(P,1),-1,P,P);
    
    for n = 1:P
        if Z(n) == 1
            d1(n,:) = 0;
        end
    end
    
    C1 = kron(d1, CC1);
    for n = 1:P %special case for Omega_1
        for j = 1:M
            if Z(n) == 1
                C1((n-1)*M+j,(n-2)*M+j) = -1/(4*lambda);
                C1((n-1)*M+j,(n-1)*M+j) = 1/(4*lambda);
            end
        end
    end
    
    
    
    %...fourth term
    CC2 = (1i/(4*dyc))*spdiags(ones(M,1)*[-1 0 1]...
        ,-1:1,M,M);
    CC2(1,1) = 1i/(4*dyc);
    CC2(M,M) = -1i/(4*dyc);
    d2 = spdiags(ones(P,1),1,P,P);
    
    
    for n = 1:P
        if Z(n) == -1
            d2(n,:) = 0;
        end
    end
    C2 = kron(d2, CC2);
    for n = 1:P %special case for Omega_{-1}
        for j = 1:M
            if Z(n) == -1
                C2((n-1)*M+j,n*M+j) = -1/(4*lambda);
                C2((n-1)*M+j,(n-1)*M+j) = 1/(4*lambda);
            end
        end
    end
    %...rhs delta matrix
    deltamat = zeros(M*P,1);
    for n = 1:P
        for j = 1:M
            if Z(n) == -1
                deltamat((n-1)*M + j,1) = -1i/2;
            end
            if Z(n) == 1
                deltamat((n-1)*M + j,1) = 1i/2;
            end
        end
    end
    deltamat = sparse(deltamat);
    %...rhs boundary conditions
    BCs = zeros(M*P,1);
    
    
    for n = 1:P
        for j = 1:M
            if Z(n) == -3
                BCs((n-1)*M + 1,1) = (1i/(2*dyc))*...
                    (-1/(8*yc(1)^3));
                BCs((n-1)*M + M,1) = -(1i/(2*dyc))*...
                    (-1/(8*yc(M)^3));
            end
            if Z(n) == 3
                BCs((n-1)*M + 1,1) = -(1i/(2*dyc))*...
                    (-1/(8*yc(1)^3));
                BCs((n-1)*M + M,1) = (1i/(2*dyc))*...
                    (-1/(8*yc(M)^3));
            end
            if Z(n) == -2
                BCs((n-1)*M + 1,1) = -lambda/...
                    (4*yc(1)^3*dyc^2)...
                    + (1i/(2*dyc))*(1/(2*yc(1)));
                BCs((n-1)*M + M,1) = -lambda/...
                    (4*yc(M)^3*dyc^2)...
                    - (1i/(2*dyc))*(1/(2*yc(M)));
            end
            if Z(n) == 2
                BCs((n-1)*M + 1,1) = -lambda/...
                    (4*yc(1)^3*dyc^2)...
                    - (1i/(2*dyc))*(1/(2*yc(1)));
                BCs((n-1)*M + M,1) = -lambda/...
                    (4*yc(M)^3*dyc^2)...
                    + (1i/(2*dyc))*(1/(2*yc(M)));
            end
            if Z(n) == -1
                BCs((n-1)*M + 1,1) = lambda/...
                    (yc(1)*dyc^2)...
                    - (1i/(2*dyc))*(-1/(8*yc(1)^3));
                BCs((n-1)*M + M,1) = lambda/...
                    (yc(M)*dyc^2)...
                    + (1i/(2*dyc))*(-1/(8*yc(M)^3));
            end
            if Z(n) == 1
                BCs((n-1)*M + 1,1) = lambda/...
                    (yc(1)*dyc^2)...
                    + (1i/(2*dyc))*(-1/(8*yc(1)^3));
                BCs((n-1)*M + M,1) = lambda/...
                    (yc(M)*dyc^2)...
                    - (1i/(2*dyc))*(-1/(8*yc(M)^3));
            end
        end
    end
    BCs = sparse(BCs);
    %...full equation
    LHS = A + B + C1 + C2;
    RHS = deltamat + BCs;
    OmegaVector = LHS\RHS;
    %... we integrate using Simpson's rule
    for n = 1:P
        integral = 0;
        for j = 1:M
            if mod(j,2) == 0
                if j == M
                    integral = integral...
                        + OmegaVector((n-1)*M+j,1);
                else
                    integral = integral...
                        + 4*OmegaVector((n-1)*M+j,1);
                end
            else
                if j == 1
                    integral = integral...
                        + OmegaVector((n-1)*M+j,1);
                elseif j == M
                    integral = integral...
                        + OmegaVector((n-1)*M+j,1);
                else
                    integral = integral...
                        + 2*OmegaVector((n-1)*M+j,1);
                end
            end
        end
        jump(l,n) = dyc*integral/3;
    end
    
    
    
    %... output
    fprintf('lambdac = %8.3f \t N = %8d \n',lambda,N);
    fprintf('Y = %8d \t dyc = %8.3f \n',Y,dyc);
    omega1jump = jump(l,P/2+1);
    fprintf('phi = %10.8f \n',2*imag(omega1jump));
end

figure
plot(lambdacmat(:,1), imag(jump(:,P/2+1)),'-');
xlabel('lambda')
ylabel('-phi/2')
title('The lambdac-phase shift relation')
figure
integral = integral...
    + 2*OmegaVector((n-1)*M+j,1);


plot(lambdacmat(:,1), real(jump(:,P/2+1)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+1)),'-');
xlabel('lambda')
ylabel('J_1')
legend('real(jump)','imag(jump)')
figure
plot(lambdacmat(:,1), real(jump(:,P/2+2)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+2)),'-');
xlabel('lambda')
ylabel('J_2')
legend('real(jump)','imag(jump)')
figure
plot(lambdacmat(:,1), real(jump(:,P/2+3)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+3)),'-');
xlabel('lambda')
ylabel('J_3')
legend('real(jump)','imag(jump)')
figure
plot(lambdacmat(:,1), real(jump(:,P/2+4)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+4)),'-');
xlabel('lambda')
ylabel('J_4')
legend('real(jump)','imag(jump)')
figure
plot(lambdacmat(:,1), real(jump(:,P/2+1)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+1)),'-',...
    lambdacmat(:,1), real(jump(:,P/2+2)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+2)),'-',...
    lambdacmat(:,1), real(jump(:,P/2+3)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+3)),'-',...
    lambdacmat(:,1), real(jump(:,P/2+4)),'-',...
    lambdacmat(:,1), imag(jump(:,P/2+4)),'-');
xlabel('lambda')
ylabel('J_n')
legend('real(jump1)','imag(jump1)',...
    'real(jump2)','imag(jump2)',...
    'real(jump3)','imag(jump3)',...
    'real(jump4)','imag(jump4)')