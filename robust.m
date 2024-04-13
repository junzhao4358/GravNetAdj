function  [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e)
format long
[n,t] = size(A);
X1 = (A'*P*A)\(A'*P*L);               %������ֵ����ʼȨ����С���˹���
X2 = 99999*ones(t,1);                 %������ʼ�Ƚ�ֵ
V = L - A*X1;                         %�����ʼ�в�
sigma0 = sqrt((V'*P*V)/(n-t));         %�����ʼ���λȨ����
% PV = inv(P)- A*inv(A'*P*A)*A';  %�в��Ȩ��
PV = diag(1./diag(P)) - A*((A'*P*A)\A');
PVV = sigma0*sqrt(diag(PV));          %�����ʼ�в�ĵ�λȨ����
iter = 0;
while max(abs(X2-X1))>e/1000 && iter <15
      %ѡȡ�������ڵ�����ֵ����ֵ����С��һ��С���������������10
      iter = iter + 1;
      X2 = X1;
      V = A*X2 - L;                     %�����µĲв�
      abVV = abs(V./PVV);               %�����׼���в�
      tao =  down_t(abVV,k1,k2,numJ);   %���㽵Ȩ����
      down_P = diag(tao.*diag(P));      %����Ȩ��
      X1 = (A'*down_P*A)\(A'*down_P*L); %��С���˹���
end
X_robust = X1;        %������յĿ����ֵ
Pend = down_P;        %������յ�Ȩ��