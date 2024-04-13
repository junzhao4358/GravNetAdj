function  [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e)
format long
[n,t] = size(A);
X1 = (A'*P*A)\(A'*P*L);               %给定初值：初始权的最小二乘估计
X2 = 99999*ones(t,1);                 %给出初始比较值
V = L - A*X1;                         %计算初始残差
sigma0 = sqrt((V'*P*V)/(n-t));         %计算初始验后单位权方差
% PV = inv(P)- A*inv(A'*P*A)*A';  %残差的权逆
PV = diag(1./diag(P)) - A*((A'*P*A)\A');
PVV = sigma0*sqrt(diag(PV));          %计算初始残差的单位权方差
iter = 0;
while max(abs(X2-X1))>e/1000 && iter <15
      %选取两个相邻迭代估值最大差值分量小于一个小数或迭代次数超过10
      iter = iter + 1;
      X2 = X1;
      V = A*X2 - L;                     %计算新的残差
      abVV = abs(V./PVV);               %计算标准化残差
      tao =  down_t(abVV,k1,k2,numJ);   %计算降权因子
      down_P = diag(tao.*diag(P));      %更新权阵
      X1 = (A'*down_P*A)\(A'*down_P*L); %最小二乘估计
end
X_robust = X1;        %输出最终的抗差估值
Pend = down_P;        %输出最终的权阵