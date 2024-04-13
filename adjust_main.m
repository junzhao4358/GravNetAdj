function results = adjust_main(yiparaN,yipara,unObsN,jdObsN,jdObs,numJ,numX,xdObs,xdObsN,app_RobusteParameter_Value,app_K1Parameter_Value,app_K2Parameter_Value,app_Adjustment_MethodButtonGroup_SelectedObject_Text,app_InstrumentParameter_Value,app_Datum_PointCheckBox_Value,app_Basic_PointCheckBox_Value,app_PointofreferenceCheckBox_Value,app_ShortBaselineCheckBox_Value,app_residualCheckBox_Value,app_instrumentCheckBox_Value)             
feature('DefaultCharacterSet','GBK'); 
slCharacterEncoding('GBK');
f = waitbar(0.1,'Computation is processing');
numyp = length(yiparaN);                 %仪器观测文件中的观测个数:总仪器个数
%确定参与平差计算的仪器参数个数stN
stN = zeros(numyp,2);
for i = 1:numyp
    for j = 1:9
        if j < 3
           if yipara(i,j) == 1
              stN(i,1) = stN(i,1) +1;
           end
        else
           if yipara(i,j)  == 1
              stN(i,2) = stN(i,2) +2;
           end
        end
    end
end
waitbar(0.2,f)
%************************************************
% unknowN为平差参数总个数：gravN为平差重力点个数，sum(sum(stN))为周期项平差个数
gravN = length(unObsN);                 %确定参与平差计算的重力点个数
unknowN = gravN + sum(sum(stN));        %平差总的未知参数个数
%************************************************
%判断绝对点是否在网中
% jdObss = jdObs;
% jdObsn = jdObsN;
ss = 0;
for i = 1:numJ
    for j = 1:length(unObsN)
        if strcmp(jdObsN(i,1),unObsN(j)) == 1
           ss = ss + 1;
           jdObsn(ss,1) = jdObsN(i,1);
           jdObss(ss,1:2) = jdObs(i,:);
         end
    end
end
numJ = ss;
%************************************************
%绝对误差观测方程建立******************************                           
%************************************************
waitbar(0.3,f);
Ajd = zeros(numJ,unknowN);             % 给出初始的绝对观测系数矩阵
Pjd = diag(jdObss(:,2));               % 绝对观测方程的权阵
for i =1:numJ
    pp = find(strcmp(jdObsn(i,1),unObsN) == 1);
    Ajd(i,pp) = 1;
end
init = 9.78*10^8;conS = 10^3;
Ljd = (jdObss(:,1)-init)/conS;         %绝对误差方程：观测值向量
%************************************************
%相对误差观测方程建立******************************                           
%************************************************
Axd = zeros(numX,unknowN);              %确定相对观测方程维数
Lxd = zeros(numX,1);
Pxd =  diag(xdObs(:,8));                %取单位阵：diag(xdObs(:,12));或者读取观测数据中的权
T = [220 110 73.33 36.67 7.33 3.67 1];                                %需要输入的参数
              
%***********************************************************************
for i = 1:numX
    %*********************************************
    %*********************************************
    %得到观测值中观测点在平差参数中的位置
    %没有给出已知重力点时的观测点在系数矩阵中的位置
    poinS = find(strcmp(xdObsN(i,2),unObsN) == 1);
    poinE = find(strcmp(xdObsN(i,3),unObsN) == 1);
    Axd(i,poinS) = -1;
    Axd(i,poinE) = 1;
    %确定仪器格值对应的系数矩阵
    %寻找相对观测值对应的仪器参数
    ge = find(strcmp(xdObsN(i,1),yiparaN(:,1)) == 1);
    for j = 1:2
        %得到C1项和C2项平差参数的系数
        if yipara(ge,j) == 1
           if j == 1
              if ge >1
                 Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+j) = (xdObs(i,3) - xdObs(i,4))/conS;
              else
                 Axd(i,gravN+j) =  (xdObs(i,3) - xdObs(i,4))/conS; 
              end
           else
              if ge >1
                 Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+j) = (xdObs(i,3)^2 - xdObs(i,4)^2)/conS^2;  
              else
                 Axd(i,gravN+j)  =   (xdObs(i,3)^2- xdObs(i,4)^2)/conS^2; 
              end
           end
        end
   end
   st = 0;
   %********************************************
   % 确定周期项平差参数的系数
   for j = 1:7
       if yipara(ge,2+j) == 1
          st = st + 1;
          if st > 0
             if ge >1
                Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+stN(ge,1)+2*st-1) = (cos(xdObs(i,1)*2*pi/T(j))-cos(xdObs(i,2)*2*pi/T(j)));
                Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+stN(ge,1)+2*st) = (sin(xdObs(i,1)*2*pi/T(j))-sin(xdObs(i,2)*2*pi/T(j)));
             else
                Axd(i,gravN+stN(ge,1)+2*st-1) = (cos(xdObs(i,1)*2*pi/T(j))-cos(xdObs(i,2)*2*pi/T(j)));
                Axd(i,gravN+stN(ge,1)+2*st) =(sin(xdObs(i,1)*2*pi/T(j))-sin(xdObs(i,2)*2*pi/T(j)));
             end
          end
       end
   end
end
waitbar(0.4,f);
A = [Ajd;Axd];
L = [Ljd' Lxd']';
P = diag([diag(Pjd)' diag(Pxd)'])';
save('data.mat','A','L','P');
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'最小二乘')
   NN = A'*P*A;NL = A'*P*L;
   gravX = NN\NL; 
   Vjd = -(Ljd - Ajd*gravX)*conS;                   %绝对观测残差
   Vxd = -(Lxd - Axd*gravX)*conS;                   %相对观测残差
   V = [Vjd' Vxd']';
   %****************************************
   gravX(1:gravN) = gravX(1:gravN)*conS+ init;      %利用初始值，计算最终参数估值
   rms = Vjd'*Pjd*Vjd + Vxd'*Pxd*Vxd;               %残差平方和
   freedom = numJ + numX - unknowN;                 %自由度
   sigma = sqrt(rms/freedom);                       %计算验后单位权中误差
   QQ  = inv(NN);                                   %未知参数的协方差矩阵
   varX = sigma*sqrt(QQ);                           %计算未知参数的中误
   varM = sigma*sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN); %计算平差重力点的平均中误差 
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'抗差估计')
    e = app_RobusteParameter_Value;
    k1 =app_K1Parameter_Value ;k2= app_K2Parameter_Value;
    [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e);
    gravX = X_robust;
    Vxd = -(Lxd - Axd*gravX)*conS;    %相对观测残差
    Vjd = -(Ljd - Ajd*gravX)*conS;    %绝对观测残差
    V = [Vjd' Vxd']';                 %总的观测残差
    %**************************************************
    %**************************************************
    gravX(1:gravN)= X_robust(1:gravN)*conS+ init;
    num0 = 0;
    for i = numJ+1:numJ+numX 
        if Pend(i,i) == 10^-10
           num0 = num0 + 1;
        end
    end
    Pend_point = find(diag(Pend) == 10^-10) - numJ;
    sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));  %计算验后单位权中误差
    QQ  = inv(A'*Pend*A);                                      %未知参数的协方差矩阵
    varX = sigma*sqrt(QQ);                                     %计算未知参数的中误差
    varM = sigma*(sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN)); %计算平差重力点的平均中误差
    app_iter_number_Value = num2str(iter);
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'方差分量估计')
          S = [numJ,numX];
          [XX,Pend] = Helmert(A,L,P,Pjd,Pxd);
          gravX = XX;
          Vxd = -(Lxd - Axd*gravX)*conS;    %相对观测残差
          Vjd = -(Ljd - Ajd*gravX)*conS;    %绝对观测残差
          V = [Vjd' Vxd']';                 %总的观测残差
          sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));         %计算验后单位权中误差
          QQ  = inv(A'*Pend*A);                                     %未知参数的协方差矩阵
          varX = sigma*sqrt(QQ);                                     %计算未知参数的中误差
          varM = sigma*(sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN)); %计算平差重力点的平均中误差
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'抗差方差分量估计') 
       s =3;
       [XX,Q] = RLSVCE_bjw(A,L,P,numJ,numX,s);
       Pend = inv(Q);
       gravX = XX;
       NN = A'*Pend*A;
       gravX = NN\(A'*Pend*L);
       Vjd = -(L(1:numJ) - A(1:numJ,:)*gravX)*conS;                   %绝对观测残差
       Vxd = -(L(numJ+1:end) - A(numJ+1:end,:)*gravX)*conS;                   %相对观测残差
       V = [Vjd' Vxd']';
       %*****************************************************************
       gravX(1:gravN) = gravX(1:gravN)*conS+ init;      %利用初始值，计算最终参数估值
%        rms = Vjd'*Pend(1:numJ,1:numJ)*Vjd + Vxd'*Pend(numJ+1:end,numJ+1:end)*Vxd;  
%       % sigma = sqrt(rms/freedom);                       %计算验后单位权中误差
%        num0 = 0;
%        for i = numJ+1:numJ+numX 
%            if Pend(i,i) == 10^-10
%               num0 = num0 + 1;
%            end
%        end
%        Pend_point = find(diag(Pend) == 10^-10) - numJ;
       sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));  %计算验后单位权中误差
       QQ  = inv(NN);                                   %未知参数的协方差矩阵
       varX = sigma*sqrt(QQ);                           %计算未知参数的中误
       varM = sigma*sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN); %计算平差重力点的平均中误差 
end
waitbar(0.5,f);
%**********************************************************************
%残差统计
%统计小于-3倍中误差或大于3倍中误差
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'最小二乘')
   PP = P;
else
   PP = Pend;
end
PV = diag(1./diag(PP)) - A*((A'*PP*A)\A');
PVV = sigma*sqrt(diag(PV)); 
%对各项结果进行精度统计分析 
Cpath = strsplit(app_InstrumentParameter_Value,'\'); %存储路径
if app_residualCheckBox_Value
   [stat_chao,numVV1,numVV2,stat_V,stat_VJ,stat_VX] = stat_res(V,PVV,Vjd,Vxd,numJ);
end
if app_instrumentCheckBox_Value
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'仪器参类型精度统计.txt'),'w');
   AS = xdObsN;
   num_AS = length(AS(:,1));
   for i = 1:num_AS
       AS(i,4) = num2cell(V(numJ+i));
   end
   AS = sortrows(AS,1);
   num_AS = length(AS(:,1));
   aa ={};
   aa{1} = AS{1,1}(1);
   %************************************************
   %找出仪器的种类
   for i = 2:num_AS
       aa{i} = AS{i,1}(1);
   end
   leii = unique(aa);  
   %***********************************************
   for j = 1:length(leii)
       ASJ = [];
       posi_yiqi = find(strcmp(aa,leii{j}) == 1);
       ASJ = cell2mat(AS(posi_yiqi,4));
       num_lei = length(posi_yiqi);
       yiqi_pre = sqrt(sum((ASJ - mean(ASJ)).^2)/num_lei/(num_lei-1));
       if j == 1
          fprintf(fid,' 仪器类型     测量次数     精度\n');
       end
       fprintf(fid, '%5s     %5d    %7.3f\n',strcat(leii{j},'型仪器'),num_lei, yiqi_pre);
       figure(j+6)
       histogram(ASJ(intersect(find(ASJ<150), find(ASJ>-150))));
       title(strcat(leii{j},'型仪器残差图'));xlabel('残差(uGal)');ylabel('个数');
   end
   msgbox('仪器参数统计完成:Operation Successfully');   
   fclose(fid);
end                 
waitbar(0.6,f);
%*******************************************************************
[yiE,QyiE,An,alf,QAn,Qalf]= stat_yiqi(yipara,gravX,varX,stN,conS,gravN,numyp);
%***************************************************
%输出剔除的观测值
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'抗差估计') ==1   
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'抗差估计剔除数据统计.txt'),'w');
   for i = 1:num0
        fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{Pend_point(i),1:3},xdObs(Pend_point(i),:));
   end
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'粗差探测') ==1 
       fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'拟准粗差探测剔除数据统计.txt'),'w');
       for i = 1:num0
           fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{posi_det(i),1:3},xdObs(posi_det(i),:));
       end
       fclose(fid)  
end
waitbar(0.7,f);
%*****************************************************
%输出重力点参数结果
fid=fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'自由网平差结果.txt'),'w','n','UTF-8');
fprintf(fid,'%7s   %7d      %7s   %4d\n','1.平差观测量:  ',numX+numJ,' 2.重力点平差参数和仪器参数个数:  ',unknowN);
fprintf(fid,'%7s  %7f    %7s   %7f\n\n','3.单位权中误差:  ',sigma,'4.重力点的平均中误差:          ',varM);
fprintf(fid,'%6s %6s    %7s     %7s\n','5.重力点平差结果:  ','点号   ','重力值','精度');
results = {};
for i = 1:gravN
    results(i,1:3) = {unObsN{i}  gravX(i) varX(i,i)};
end
results = sortrows(results,1);
for i =1:gravN
    if i == gravN
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n\n',results{i,:}); 
    else
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n',results{i,:});     %输出各个重力点的精度
    end
end
%if app.East_westCheckBox.Value   %统计中西部精度
%   fprintf(fid,'重力点类型    个数   精度');
%   fpritnf(fid,' 基本点   %4d     %7.3f',num_Eastwest);
% end
C = firstalp(results);
if app_Datum_PointCheckBox_Value       %统计基准点精度                    
   CC1 = find(strcmp(C,'A') == 1);
   num_Datumpoint = length(CC1);
   abs_pre = sqrt(sum(cell2mat(results(CC1,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.重力点精度统计:');
   fprintf(fid,'重力点类型     个数      精度\n');
   fprintf(fid,'                   基准点      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_Basic_PointCheckBox_Value      %统计基本点精度
    CC2 = find(strcmp(C,'B') == 1);
    abs_pre = sqrt(sum(cell2mat(results(CC2,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.重力点精度统计:');
    fprintf(fid,'重力点类型     个数      精度\n');
    fprintf(fid,'                   基本点      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if app_PointofreferenceCheckBox_Value   %统计引点精度
   CC3 = find(strcmp(C,'C') == 1);
   num_Datumpoint = length(CC3);
   abs_pre = sqrt(sum(cell2mat(results(CC3,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.重力点精度统计:');
   fprintf(fid,'重力点类型     个数      精度\n');
   fprintf(fid,'                   引点      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_ShortBaselineCheckBox_Value     %统计短基线精度
    CC4 = find(strcmp(C4,'S') == 1);
    num_Datumpoint = length(CC4);
    abs_pre = sqrt(sum(cell2mat(results(CC4,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.重力点精度统计:');
    fprintf(fid,'重力点类型     个数      精度\n');
    fprintf(fid,'                   短基线      %4d     %7.3f\n',num_Datumpoint,abs_pre);   
end
waitbar(0.8,f);
%***************************************************************
%输入仪器参数平差结果和精度
fprintf(fid,'\n');
fprintf(fid,'%12s  %4d  %4s %4s   %9s\n','7.格值和周期误差参数:',unknowN-gravN,'C1','C2','Xi,Yi(i=1:7)');
for i= 1:numyp
    if i == numyp
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n\n',yiparaN{i,1},yiE(i,:));
    else
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n',yiparaN{i,1},yiE(i,:));
    end
end
fprintf(fid,'%12s\n','8.格值和周期误差参数的中误差:');
for i= 1:numyp
    if i == numyp
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n\n',yiparaN{i,1},QyiE(i,:));
    else
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n',yiparaN{i,1},QyiE(i,:));
    end     
end
waitbar(0.9,f);
%*****************************************************
for i=1:numyp
    fprintf(fid,'      %6s       %13s    %10s      %10s     %9s    %10s    %11s     %11s     %11s\n',yiparaN{i,1},'x(C1)','mx','y(C2)','my','A','mA','alp','malp');
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f\n','一二次因子:',yiE(i,1),QyiE(i,1),yiE(i,2),QyiE(i,2));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项1:',yiE(i,3),QyiE(i,3),yiE(i,4),QyiE(i,4),An(i,1),QAn(i,1),alf(i,1),Qalf(i,1));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项2:',yiE(i,5),QyiE(i,5),yiE(i,6),QyiE(i,6),An(i,2),QAn(i,2),alf(i,2),Qalf(i,2));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项3:',yiE(i,7),QyiE(i,7),yiE(i,8),QyiE(i,8),An(i,3),QAn(i,3),alf(i,3),Qalf(i,3));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项4:',yiE(i,9),QyiE(i,9),yiE(i,10),QyiE(i,10),An(i,4),QAn(i,4),alf(i,4),Qalf(i,4));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项5:',yiE(i,11),QyiE(i,11),yiE(i,12),QyiE(i,12),An(i,5),QAn(i,5),alf(i,5),Qalf(i,5));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','仪器周期项6:',yiE(i,13),QyiE(i,13),yiE(i,14),QyiE(i,14),An(i,6),QAn(i,6),alf(i,6),Qalf(i,6));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n\n','仪器周期项7:',yiE(i,15),QyiE(i,15),yiE(i,16),QyiE(i,16),An(i,7),QAn(i,7),alf(i,7),Qalf(i,7));
end   
fclose(fid)
if app_residualCheckBox_Value
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'自由网平差残差结果.txt'),'w');
   fprintf(fid,'%8s  %8s\n','总残差大于零个数','总残差小于零的个数');
   fprintf(fid,'%    5d          %5d\n\n',numVV1,numVV2);
   fprintf(fid,'%6s   %6s  %6s\n','总残差最大值','总残差最小值','总残差平均值');
   fprintf(fid,'%10.3f   %10.3f   %10.3f \n',stat_V);
   fprintf(fid,'%9s     %9s     %9s\n','绝对重力残差最大值','绝对重力残差最小值','绝对重力残差平均值');
   fprintf(fid,'%10.3f     %10.3f       %10.3f\n\n',stat_VJ);
   fprintf(fid,'%9s  %9s  %9s\n','相对重力残差最大值','相对重力残差最小值','相对重力残差平均值');
   fprintf(fid,'%10.3f       %10.3f       %10.3f\n\n',stat_VX);
   fprintf(fid,'-3*m0=<V<-2*m0:');fprintf(fid,'%5d\n',stat_chao(1));
   fprintf(fid,'-3*m0=<V<-2*m0:');fprintf(fid,'%5d\n',stat_chao(2));
   fprintf(fid,'-2*m0=<V<-1*m0:');fprintf(fid,'%5d\n',stat_chao(3));
   fprintf(fid,'-1*m0=<V<0*m0:');fprintf(fid,'%5d\n',stat_chao(4));
   fprintf(fid,'0*m0=<V<1*m0:');fprintf(fid,'%5d\n',stat_chao(5));
   fprintf(fid,'1*m0=<V<2*m0:');fprintf(fid,'%5d\n',stat_chao(6));
   fprintf(fid,'2*m0=<V<3*m0:');fprintf(fid,'%5d\n',stat_chao(7));
   fprintf(fid,'3*m0=<V:');fprintf(fid,'%5d\n',stat_chao(8));
   fprintf(fid,'%5s    %5s\n', '残差','中误差');
   for i = 1:numJ + numX
       fprintf(fid,'%12f  %12f \n',V(i),PVV(i));
   end
   fclose all
end
waitbar(1,f,{'The processing is ending' 'Operation Successfully'});
results = 0;