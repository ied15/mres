function[Jn] = CriticalLayer2(lambdac, N,Bmat)

%Still need to implement boundary conditions
%--------------------------------------
% Solve for Omega_N, eliminate n = 0
%--------------------------------------
format long

%...truncate omega
Z1 = -N:1:-1;
Z2 = 1:1:N;
Z = [Z1 Z2]; %creates vector excluding N = 0...
P = 2*N;

% Grid setup inside crit layer
Y = 40; %we choose Y from 40 to 80
dyc = 0.05; %we choose a uniform step size from 0.02 to 0.1
y = -Y:dyc:Y;
yc = (y(1:end-1)+y(2:end))/2; %cell centered approximation
M = length(yc);

%vector to be populated with jump for each lambda
jump = zeros(P,1);
    
    
%...first term
A = zeros(M*P,M*P);
for n = 1:P
    for j = 1:M
        A((n-1)*M + j, (n-1)*M + j) = 1i*Z(n)*yc(j);
    end
end
A = sparse(A);



%...second term
BB = (-lambdac)*spdiags(ones(M,1)*[1 -2 1]/dyc/dyc...
    ,-1:1,M,M);
BB(1,1) = -3*(-lambdac)/dyc^2;

BB(M,M) = -3*(-lambdac)/dyc^2;
B = kron(speye(P),BB);




%...third term
CC1 = (1/(2*dyc))*spdiags(ones(M,1)*[-1 0 1]...
    ,-1:1,M,M);
CC1(1,1) = 1/(2*dyc);
CC1(M,M) = -1/(2*dyc);
dd1 = -1i/2*rot90(hankel([Bmat(2*N-1:-1:N+1);0]));
dd2 = 1i/2*rot90(hankel([Bmat(2:N);0]),3);
d1 = [dd1 zeros(N,N); zeros(N,N) dd2];

% d1 = spdiags(ones(P,1),-1,P,P);
% 
% for n = 1:P
%     if Z(n) == 1
%         d1(n,:) = 0;
%     end
% end
tic
C1 = kron(d1, CC1);
toc
% for n = 1:P %special case for Omega_1
%     for j = 1:M
%         if Z(n) == 1
%             C1((n-1)*M+j,(n-2)*M+j) = -1/(4*lambdac);
%             C1((n-1)*M+j,(n-1)*M+j) = 1/(4*lambdac);
%         end
%     end
% end



%...fourth term
CC2 = (1/(2*dyc))*spdiags(ones(M,1)*[-1 0 1]...
    ,-1:1,M,M);
CC2(1,1) = 1i/(4*dyc);
CC2(M,M) = -1i/(4*dyc);
f1 = 1i/2*rot90( hankel([Bmat(2:N);0]),3);
f2 = 1i/2*rot90(hankel([Bmat(N-1:-1:1);0]));
f3 = -1i/2*rot90(hankel([Bmat(N+2:2*N);0]),3);
f4 = -1i/2*rot90(hankel([Bmat(2*N-1:-1:N+1);0]));
d2 = [f1 f2; f3 f4];

%d2 = spdiags(ones(P,1),1,P,P);


% for n = 1:P
%     if Z(n) == -1
%         d2(n,:) = 0;
%     end
% end

C2 = kron(d2, CC2);
% for n = 1:P %special case for Omega_{-1}
%     for j = 1:M
%         if Z(n) == -1
%             C2((n-1)*M+j,n*M+j) = -1/(4*lambdac);
%             C2((n-1)*M+j,(n-1)*M+j) = 1/(4*lambdac);
%         end
%     end
% end

DD2 = ones(M,M);
D2 = repmat(Bmat,1,2*N)';
for n = 1:N
    D2(n,:) = 1i/(2*lambdac)*Bmat(n)*D2(n,:);
    D2(N+n,:) = -1i/(2*lambdac)*Bmat(N+n)*D2(N+n,:);
end

tic
F2 = kron(D2, DD2);
toc
%...rhs delta matrix
deltamat = zeros(M*P,1);
for n = 1:P
    for j = 1:M
        deltamat((n-1)*M + j,1) = -1i/2*Bmat(n);
        deltamat((n-1)*M + j,1) = 1i/2*Bmat(N+n);
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
            BCs((n-1)*M + 1,1) = -lambdac/...
                (4*yc(1)^3*dyc^2)...
                + (1i/(2*dyc))*(1/(2*yc(1)));
            BCs((n-1)*M + M,1) = -lambdac/...
                (4*yc(M)^3*dyc^2)...
                - (1i/(2*dyc))*(1/(2*yc(M)));
        end
        if Z(n) == 2
            BCs((n-1)*M + 1,1) = -lambdac/...
                (4*yc(1)^3*dyc^2)...
                - (1i/(2*dyc))*(1/(2*yc(1)));
            BCs((n-1)*M + M,1) = -lambdac/...
                (4*yc(M)^3*dyc^2)...
                + (1i/(2*dyc))*(1/(2*yc(M)));
        end
        if Z(n) == -1
            BCs((n-1)*M + 1,1) = lambdac/...
                (yc(1)*dyc^2)...
                - (1i/(2*dyc))*(-1/(8*yc(1)^3));
            BCs((n-1)*M + M,1) = lambdac/...
                (yc(M)*dyc^2)...
                + (1i/(2*dyc))*(-1/(8*yc(M)^3));
        end
        if Z(n) == 1
            BCs((n-1)*M + 1,1) = lambdac/...
                (yc(1)*dyc^2)...
                + (1i/(2*dyc))*(-1/(8*yc(1)^3));
            BCs((n-1)*M + M,1) = lambdac/...
                (yc(M)*dyc^2)...
                - (1i/(2*dyc))*(-1/(8*yc(M)^3));
        end
    end
end
BCs = sparse(BCs);
%...full equation
LHS = A + B + C1 + C2 + F2;
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
    jump(n,1) = dyc*integral/3;
end



%... output
fprintf('lambdac = %8.3f \t N = %8d \n',lambdac,N);
fprintf('Y = %8d \t dyc = %8.3f \n',Y,dyc);
Jn = jump;
%fprintf('phi/2 = %10.8f \n',imag(Jn));
