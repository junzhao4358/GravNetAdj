%-------------------------Downweight function-----------------------------
function tao = down_t(abVV,k1,k2,numJ)
%仅考虑独立等权条件下的抗差估计
%P:weight matrix
%k1,k2:为IGGIII方案的参数
%A:表示系数矩阵
%L:误差方程右端项
%down_p：为最终的权
%----------------------确定降权函数-------------------------
n = length(abVV);
tao = ones(n,1);
for i = numJ+1:n
       if abVV(i) <= k1
          tao(i,1) = 1;
       elseif abVV(i) >k1 && abVV(i) <=k2
          tao(i,1) = k1/abVV(i)*((k2-abVV(i))/(k2-k1))^2;
       else
          tao(i,1) = 10^-10;
       end
end
%***********************************************
% [n,m] = size(A);
% N = A'*P*A; 
% S = A*(N\A');
% R = eye(n) - S*P;         %平差因子阵
% V = R*L;                  %残差
% PV = (diag(1./diag(P))-S);%残差的权逆阵
% num = 0;
%**********************************************
%计算标准化残差和统计权阵为零的个数
% for i = 1:n
%       VV(i) = V(i)/sqrt(PV(i,i));   %计算标准化残差，这里没有处理单位权方差
% %     if P(i,i) <= 10^-8
% %         num = num +1;             %统计权为零的个数
% %     end
% end
%**********************************************
%标准化残差除以sigma的方案：共三种
% mad_sigma = sqrt((V'*P*V)/(n-m-num));   %利用验后单位权方差确定单位权中误差，剔除权为零的观测量，自由度需减去相权为零的个数
% % mad_sigma = mad(abs(VV),1)*1.4826;     %利用中位数法确定单位权中误差
% % mad_sigma =sqrt(median(VV.^2))*1.4826; %利用抗差LMS(抗差最下二乘中位数估计);
% %**********************************************
% % VV = VV/mad_sigma;                       %得到最终的标准化残差

%**********************************************
% 依据降权IGGIII方案，计算降权因子
% tao = zeros(n,1);
% for i = 1:n
%        abVV = abs(VV(i));
%        if abVV  <= k1
%           tao(i,1) = 1;
%        elseif abVV >k1 && abVV <=k2
%           tao(i,1) = k1/abVV*((k2-abVV)/(k2-k1))^2;
%        else
%           tao(i,1) = 10^-8;
%        end
% end
%***************************************
%解算方案：确定在397开始的时候才变权阵且在飞机固定权时
%***************************************
%方案1：将已知重力基准值的权确定为无限大，故在开展抗差估计时，其权无需改变
%这里已知重力基准点共8个
% tao(2011:end) = 1;                  %确定相对观测中前396个观测的权不变:程序2选择
%***************************************
%方案2：将已知重力基准值当做已知点进行解算且不参与平差计算，故在开展抗差估计时，无需考虑对应的权。
% tao(numJ+1:numJ+396) = 1;           %确定相对观测中前396个观测的权不变
% tao(ts) = 1;                        %396以后为飞机的权不变
% % tao(1:85,1) = ones(1,85);
% tao(1:numJ,1) = ones(1,numJ);
% for i = 86:n
%     if P(i,i) <= 10^-5
%         tao(i) = 1;
%     end
%     if tao(i) == 0
%         P(i,i) = 10^-8;
%         tao(i) = 1;
%     end
% end
% for i = 145:n
%     if P(i,i) <= 10^-5
%         tao(i) = 1;
%     end
%     if tao(i) == 0
%         P(i,i) = 10^-8;
%         tao(i) = 1;
%     end
% end
% tao(2023+85:end,1) = ones(1,1114);
% down_p = diag(tao.*diag(P));        %得到最终的权矩阵
% num0 = 0;
% for i = 1:n 
%     if down_p(i,i) == 10^-8
%        num0 = num0 + 1;
%     end
% end
% num0
%***********************************************
%        for i = 1:n
%            tao(i) = V(i)/sqrt(PV(i,i));    
%        end
% %      mad_sigma = median(abs(tao-ones(n,1)*median(tao)))*1.4826;
%        mad_sigma = mad(abs(tao),1)*1.4826;
%        tao = tao/mad_sigma;
%        w = min(1,tune./abs(tao));           %这里采用Huber权函数 
%        ts = numJ + weiN;
%        w(ts) = 1;
%        down_p = diag(w'.*diag(P)); 
%mad_sigma = mad(abs(tao),1)*1.4826;
%tao = tao/mad_sigma;
%w = min(1,tune./abs(tao));                   %这里采用Huber权函数 
%w = 1./w;
%down_Q = diag([w])*Q*diag([w]);
%***********************************************