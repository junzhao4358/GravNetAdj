function [XX,Pend] = Helmert(A,L,P,P1,P2)
% L = [L1' L2'] = [A1'; A2';]'*X + e e~N(0,[P1 0;0 P2];
%This program is about the variance component estimation (VCE) for two type
%observations
%L:observation
%A: coefficient matrix
%L: observation in observation equation
[n1,m1] = size(P1);[n2,m2] = size(P2);
[n,m] = size(A);
L1 = L(1:n1);
L2 = L(n1+1:end);
A1 = A(1:n1,:); A2 = A(n1+1:end,:);
iter_lim = 40;                  %the maximun number of iteration
vc = [1 1;zeros(iter_lim,2)];   %variance component for each iteration
n = max(size(A));
iter = 0;
for i = 1:iter_lim                 
    XX = WLS(A,L,P);
    N = A'*P*A;V1 = L1 - A1*XX;V2 = L2 - A2*XX;N1 = A1'*P1*A1;N2 = A2'*P2*A2;
    N11 =N\N1;N22 = N\N2;
    NX = [n1-2*trace(N11)+trace(N11*N11) trace(N11*N22);trace(N11*N22) n2-2*trace(N22)+trace(N22*N22)];  %normal matrix
    lX = [V1'*P1*V1 V2'*P2*V2]';         % right-hand term for normal equation
    vc(i+1,:) = NX\lX;                    % Estimator of variance components
    P1 = P1;P2 = vc(i+1,1)/vc(i+1,2)*P2;  % Updating the weighted matrix
    P(1:n1,1:n1) = P1;P(n1+1:end,n1+1:end) = P2;
    iter = i;                             % The number of iteration
    if norm(vc(i+1,:)-vc(i,:),inf)<=0.01  % The condition for ending of iteration
       break;
    end
end
X_vc = WLS(A,L,P);            %LS estimation
Pend = P;

%***********************************
%Attention: the most improtant is that the parameter you input is right.


        
      
      
        
        
        