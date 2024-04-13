function [AA,LL,PP] = Quasi_Detection(A,L,P,num_jd)
[n,m] = size(A);
R = eye(n) - A*((A'*P*A)\A')*P;
W = zeros(n);
aa = zeros(n,1);
PRL = P*R*L;
PR = P*R;
k = 1.5; %其中k一般取为[1.5,2.5],可以自行调整
for i = 1:n
    w(i) = abs(PRL(i)/PR(i,i));
end
for i = 1:n
    if w(i) < k*1.4826*sqrt(median(w.^2)) %其中k一般取为[1.5,2.5]
       W(i,i) = 1;
    else
       aa(i) = i;
    end
end
aa(find(aa == 0)) = [];
W([aa],:) = [];
G = A'*W'*((W*inv(P)*W')\W);
invpr = inv(P*R + G'*G);
delt = -invpr*(P*R*L);
figure(4)
plot(num_jd+1:n,delt(num_jd+1:n),'.')     %给出相应的真误差图形
Qt = invpr*PR*invpr;
posi_det = zeros(n,1);
for i = num_jd+1:n
    delt_st(i-num_jd) = abs(delt(i)/Qt(i,i));
    if delt_st(i) >= 2.5
        posi_det(i) = i;
    end
end
posi_det(find(posi_det == 0)) = [];
figure(6)
plot(num_jd+1:n,delt_st,'.',num_jd+1:n,2.5*ones(n-num_jd,1))
A(posi_det,:) = [];L(posi_det) =[];
P(posi_det,:) = [];P(:,posi_det)= [];
AA = A;LL=L;PP = P;

