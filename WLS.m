%----------------------The weight least square--------------------------
function  X = WLS(A,L,P)
%--------------weighted least squares fit----------------
 chol_P = chol(P);
 [row col] = size(A);
 chol_PL = chol_P*L;
 chol_PA = chol_P*A;
 [QQ,RR] = qr(chol_PA,0); 
 X = RR\(QQ'*chol_PL);
% X1 = (A'*P*A)\(A'*P*L);