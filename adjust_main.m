function results = adjust_main(yiparaN,yipara,unObsN,jdObsN,jdObs,numJ,numX,xdObs,xdObsN,app_RobusteParameter_Value,app_K1Parameter_Value,app_K2Parameter_Value,app_Adjustment_MethodButtonGroup_SelectedObject_Text,app_InstrumentParameter_Value,app_Datum_PointCheckBox_Value,app_Basic_PointCheckBox_Value,app_PointofreferenceCheckBox_Value,app_ShortBaselineCheckBox_Value,app_residualCheckBox_Value,app_instrumentCheckBox_Value)             
feature('DefaultCharacterSet','GBK'); 
slCharacterEncoding('GBK');
f = waitbar(0.1,'Computation is processing');
numyp = length(yiparaN);                 %�����۲��ļ��еĹ۲����:����������
%ȷ������ƽ������������������stN
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
% unknowNΪƽ������ܸ�����gravNΪƽ�������������sum(sum(stN))Ϊ������ƽ�����
gravN = length(unObsN);                 %ȷ������ƽ���������������
unknowN = gravN + sum(sum(stN));        %ƽ���ܵ�δ֪��������
%************************************************
%�жϾ��Ե��Ƿ�������
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
%�������۲ⷽ�̽���******************************                           
%************************************************
waitbar(0.3,f);
Ajd = zeros(numJ,unknowN);             % ������ʼ�ľ��Թ۲�ϵ������
Pjd = diag(jdObss(:,2));               % ���Թ۲ⷽ�̵�Ȩ��
for i =1:numJ
    pp = find(strcmp(jdObsn(i,1),unObsN) == 1);
    Ajd(i,pp) = 1;
end
init = 9.78*10^8;conS = 10^3;
Ljd = (jdObss(:,1)-init)/conS;         %�������̣��۲�ֵ����
%************************************************
%������۲ⷽ�̽���******************************                           
%************************************************
Axd = zeros(numX,unknowN);              %ȷ����Թ۲ⷽ��ά��
Lxd = zeros(numX,1);
Pxd =  diag(xdObs(:,8));                %ȡ��λ��diag(xdObs(:,12));���߶�ȡ�۲������е�Ȩ
T = [220 110 73.33 36.67 7.33 3.67 1];                                %��Ҫ����Ĳ���
              
%***********************************************************************
for i = 1:numX
    %*********************************************
    %*********************************************
    %�õ��۲�ֵ�й۲����ƽ������е�λ��
    %û�и�����֪������ʱ�Ĺ۲����ϵ�������е�λ��
    poinS = find(strcmp(xdObsN(i,2),unObsN) == 1);
    poinE = find(strcmp(xdObsN(i,3),unObsN) == 1);
    Axd(i,poinS) = -1;
    Axd(i,poinE) = 1;
    %ȷ��������ֵ��Ӧ��ϵ������
    %Ѱ����Թ۲�ֵ��Ӧ����������
    ge = find(strcmp(xdObsN(i,1),yiparaN(:,1)) == 1);
    for j = 1:2
        %�õ�C1���C2��ƽ�������ϵ��
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
   % ȷ��������ƽ�������ϵ��
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
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'��С����')
   NN = A'*P*A;NL = A'*P*L;
   gravX = NN\NL; 
   Vjd = -(Ljd - Ajd*gravX)*conS;                   %���Թ۲�в�
   Vxd = -(Lxd - Axd*gravX)*conS;                   %��Թ۲�в�
   V = [Vjd' Vxd']';
   %****************************************
   gravX(1:gravN) = gravX(1:gravN)*conS+ init;      %���ó�ʼֵ���������ղ�����ֵ
   rms = Vjd'*Pjd*Vjd + Vxd'*Pxd*Vxd;               %�в�ƽ����
   freedom = numJ + numX - unknowN;                 %���ɶ�
   sigma = sqrt(rms/freedom);                       %�������λȨ�����
   QQ  = inv(NN);                                   %δ֪������Э�������
   varX = sigma*sqrt(QQ);                           %����δ֪����������
   varM = sigma*sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN); %����ƽ���������ƽ������� 
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'�������')
    e = app_RobusteParameter_Value;
    k1 =app_K1Parameter_Value ;k2= app_K2Parameter_Value;
    [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e);
    gravX = X_robust;
    Vxd = -(Lxd - Axd*gravX)*conS;    %��Թ۲�в�
    Vjd = -(Ljd - Ajd*gravX)*conS;    %���Թ۲�в�
    V = [Vjd' Vxd']';                 %�ܵĹ۲�в�
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
    sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));  %�������λȨ�����
    QQ  = inv(A'*Pend*A);                                      %δ֪������Э�������
    varX = sigma*sqrt(QQ);                                     %����δ֪�����������
    varM = sigma*(sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN)); %����ƽ���������ƽ�������
    app_iter_number_Value = num2str(iter);
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'�����������')
          S = [numJ,numX];
          [XX,Pend] = Helmert(A,L,P,Pjd,Pxd);
          gravX = XX;
          Vxd = -(Lxd - Axd*gravX)*conS;    %��Թ۲�в�
          Vjd = -(Ljd - Ajd*gravX)*conS;    %���Թ۲�в�
          V = [Vjd' Vxd']';                 %�ܵĹ۲�в�
          sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));         %�������λȨ�����
          QQ  = inv(A'*Pend*A);                                     %δ֪������Э�������
          varX = sigma*sqrt(QQ);                                     %����δ֪�����������
          varM = sigma*(sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN)); %����ƽ���������ƽ�������
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'������������') 
       s =3;
       [XX,Q] = RLSVCE_bjw(A,L,P,numJ,numX,s);
       Pend = inv(Q);
       gravX = XX;
       NN = A'*Pend*A;
       gravX = NN\(A'*Pend*L);
       Vjd = -(L(1:numJ) - A(1:numJ,:)*gravX)*conS;                   %���Թ۲�в�
       Vxd = -(L(numJ+1:end) - A(numJ+1:end,:)*gravX)*conS;                   %��Թ۲�в�
       V = [Vjd' Vxd']';
       %*****************************************************************
       gravX(1:gravN) = gravX(1:gravN)*conS+ init;      %���ó�ʼֵ���������ղ�����ֵ
%        rms = Vjd'*Pend(1:numJ,1:numJ)*Vjd + Vxd'*Pend(numJ+1:end,numJ+1:end)*Vxd;  
%       % sigma = sqrt(rms/freedom);                       %�������λȨ�����
%        num0 = 0;
%        for i = numJ+1:numJ+numX 
%            if Pend(i,i) == 10^-10
%               num0 = num0 + 1;
%            end
%        end
%        Pend_point = find(diag(Pend) == 10^-10) - numJ;
       sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));  %�������λȨ�����
       QQ  = inv(NN);                                   %δ֪������Э�������
       varX = sigma*sqrt(QQ);                           %����δ֪����������
       varM = sigma*sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN); %����ƽ���������ƽ������� 
end
waitbar(0.5,f);
%**********************************************************************
%�в�ͳ��
%ͳ��С��-3�����������3�������
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'��С����')
   PP = P;
else
   PP = Pend;
end
PV = diag(1./diag(PP)) - A*((A'*PP*A)\A');
PVV = sigma*sqrt(diag(PV)); 
%�Ը��������о���ͳ�Ʒ��� 
Cpath = strsplit(app_InstrumentParameter_Value,'\'); %�洢·��
if app_residualCheckBox_Value
   [stat_chao,numVV1,numVV2,stat_V,stat_VJ,stat_VX] = stat_res(V,PVV,Vjd,Vxd,numJ);
end
if app_instrumentCheckBox_Value
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'���������;���ͳ��.txt'),'w');
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
   %�ҳ�����������
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
          fprintf(fid,' ��������     ��������     ����\n');
       end
       fprintf(fid, '%5s     %5d    %7.3f\n',strcat(leii{j},'������'),num_lei, yiqi_pre);
       figure(j+6)
       histogram(ASJ(intersect(find(ASJ<150), find(ASJ>-150))));
       title(strcat(leii{j},'�������в�ͼ'));xlabel('�в�(uGal)');ylabel('����');
   end
   msgbox('��������ͳ�����:Operation Successfully');   
   fclose(fid);
end                 
waitbar(0.6,f);
%*******************************************************************
[yiE,QyiE,An,alf,QAn,Qalf]= stat_yiqi(yipara,gravX,varX,stN,conS,gravN,numyp);
%***************************************************
%����޳��Ĺ۲�ֵ
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'�������') ==1   
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'��������޳�����ͳ��.txt'),'w');
   for i = 1:num0
        fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{Pend_point(i),1:3},xdObs(Pend_point(i),:));
   end
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'�ֲ�̽��') ==1 
       fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'��׼�ֲ�̽���޳�����ͳ��.txt'),'w');
       for i = 1:num0
           fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{posi_det(i),1:3},xdObs(posi_det(i),:));
       end
       fclose(fid)  
end
waitbar(0.7,f);
%*****************************************************
%���������������
fid=fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'������ƽ����.txt'),'w','n','UTF-8');
fprintf(fid,'%7s   %7d      %7s   %4d\n','1.ƽ��۲���:  ',numX+numJ,' 2.������ƽ�������������������:  ',unknowN);
fprintf(fid,'%7s  %7f    %7s   %7f\n\n','3.��λȨ�����:  ',sigma,'4.�������ƽ�������:          ',varM);
fprintf(fid,'%6s %6s    %7s     %7s\n','5.������ƽ����:  ','���   ','����ֵ','����');
results = {};
for i = 1:gravN
    results(i,1:3) = {unObsN{i}  gravX(i) varX(i,i)};
end
results = sortrows(results,1);
for i =1:gravN
    if i == gravN
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n\n',results{i,:}); 
    else
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n',results{i,:});     %�������������ľ���
    end
end
%if app.East_westCheckBox.Value   %ͳ������������
%   fprintf(fid,'����������    ����   ����');
%   fpritnf(fid,' ������   %4d     %7.3f',num_Eastwest);
% end
C = firstalp(results);
if app_Datum_PointCheckBox_Value       %ͳ�ƻ�׼�㾫��                    
   CC1 = find(strcmp(C,'A') == 1);
   num_Datumpoint = length(CC1);
   abs_pre = sqrt(sum(cell2mat(results(CC1,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.�����㾫��ͳ��:');
   fprintf(fid,'����������     ����      ����\n');
   fprintf(fid,'                   ��׼��      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_Basic_PointCheckBox_Value      %ͳ�ƻ����㾫��
    CC2 = find(strcmp(C,'B') == 1);
    abs_pre = sqrt(sum(cell2mat(results(CC2,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.�����㾫��ͳ��:');
    fprintf(fid,'����������     ����      ����\n');
    fprintf(fid,'                   ������      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if app_PointofreferenceCheckBox_Value   %ͳ�����㾫��
   CC3 = find(strcmp(C,'C') == 1);
   num_Datumpoint = length(CC3);
   abs_pre = sqrt(sum(cell2mat(results(CC3,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.�����㾫��ͳ��:');
   fprintf(fid,'����������     ����      ����\n');
   fprintf(fid,'                   ����      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_ShortBaselineCheckBox_Value     %ͳ�ƶ̻��߾���
    CC4 = find(strcmp(C4,'S') == 1);
    num_Datumpoint = length(CC4);
    abs_pre = sqrt(sum(cell2mat(results(CC4,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.�����㾫��ͳ��:');
    fprintf(fid,'����������     ����      ����\n');
    fprintf(fid,'                   �̻���      %4d     %7.3f\n',num_Datumpoint,abs_pre);   
end
waitbar(0.8,f);
%***************************************************************
%������������ƽ�����;���
fprintf(fid,'\n');
fprintf(fid,'%12s  %4d  %4s %4s   %9s\n','7.��ֵ������������:',unknowN-gravN,'C1','C2','Xi,Yi(i=1:7)');
for i= 1:numyp
    if i == numyp
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n\n',yiparaN{i,1},yiE(i,:));
    else
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n',yiparaN{i,1},yiE(i,:));
    end
end
fprintf(fid,'%12s\n','8.��ֵ�������������������:');
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
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f\n','һ��������:',yiE(i,1),QyiE(i,1),yiE(i,2),QyiE(i,2));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������1:',yiE(i,3),QyiE(i,3),yiE(i,4),QyiE(i,4),An(i,1),QAn(i,1),alf(i,1),Qalf(i,1));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������2:',yiE(i,5),QyiE(i,5),yiE(i,6),QyiE(i,6),An(i,2),QAn(i,2),alf(i,2),Qalf(i,2));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������3:',yiE(i,7),QyiE(i,7),yiE(i,8),QyiE(i,8),An(i,3),QAn(i,3),alf(i,3),Qalf(i,3));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������4:',yiE(i,9),QyiE(i,9),yiE(i,10),QyiE(i,10),An(i,4),QAn(i,4),alf(i,4),Qalf(i,4));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������5:',yiE(i,11),QyiE(i,11),yiE(i,12),QyiE(i,12),An(i,5),QAn(i,5),alf(i,5),Qalf(i,5));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n','����������6:',yiE(i,13),QyiE(i,13),yiE(i,14),QyiE(i,14),An(i,6),QAn(i,6),alf(i,6),Qalf(i,6));
    fprintf(fid,'    %8s    %13.10f   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n\n','����������7:',yiE(i,15),QyiE(i,15),yiE(i,16),QyiE(i,16),An(i,7),QAn(i,7),alf(i,7),Qalf(i,7));
end   
fclose(fid)
if app_residualCheckBox_Value
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'������ƽ��в���.txt'),'w');
   fprintf(fid,'%8s  %8s\n','�ܲв���������','�ܲв�С����ĸ���');
   fprintf(fid,'%    5d          %5d\n\n',numVV1,numVV2);
   fprintf(fid,'%6s   %6s  %6s\n','�ܲв����ֵ','�ܲв���Сֵ','�ܲв�ƽ��ֵ');
   fprintf(fid,'%10.3f   %10.3f   %10.3f \n',stat_V);
   fprintf(fid,'%9s     %9s     %9s\n','���������в����ֵ','���������в���Сֵ','���������в�ƽ��ֵ');
   fprintf(fid,'%10.3f     %10.3f       %10.3f\n\n',stat_VJ);
   fprintf(fid,'%9s  %9s  %9s\n','��������в����ֵ','��������в���Сֵ','��������в�ƽ��ֵ');
   fprintf(fid,'%10.3f       %10.3f       %10.3f\n\n',stat_VX);
   fprintf(fid,'-3*m0=<V<-2*m0:');fprintf(fid,'%5d\n',stat_chao(1));
   fprintf(fid,'-3*m0=<V<-2*m0:');fprintf(fid,'%5d\n',stat_chao(2));
   fprintf(fid,'-2*m0=<V<-1*m0:');fprintf(fid,'%5d\n',stat_chao(3));
   fprintf(fid,'-1*m0=<V<0*m0:');fprintf(fid,'%5d\n',stat_chao(4));
   fprintf(fid,'0*m0=<V<1*m0:');fprintf(fid,'%5d\n',stat_chao(5));
   fprintf(fid,'1*m0=<V<2*m0:');fprintf(fid,'%5d\n',stat_chao(6));
   fprintf(fid,'2*m0=<V<3*m0:');fprintf(fid,'%5d\n',stat_chao(7));
   fprintf(fid,'3*m0=<V:');fprintf(fid,'%5d\n',stat_chao(8));
   fprintf(fid,'%5s    %5s\n', '�в�','�����');
   for i = 1:numJ + numX
       fprintf(fid,'%12f  %12f \n',V(i),PVV(i));
   end
   fclose all
end
waitbar(1,f,{'The processing is ending' 'Operation Successfully'});
results = 0;