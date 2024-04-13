function adj = Gravity_Lice(numtObs,numtgezhi)
%******************************************
format long 
% numtgezhi = dir('./gezhi/*.txt');         %�������ȡtxt�ļ�
numfiles= length(numtObs);                %�õ��ļ��ĸ���
%feature('DefaultCharacterSet','UTF8');  %�޸ı����ʽ����ֹ��ȡ����ʱ��������
%feature('DefaultCharacterSet','GBK');
k2 = 0;
%��ȡȫ����������������
fid = fopen(strcat(numtObs.folder,'\',numtObs.name),'r');
i = 0;
while (~feof(fid))
          line = fgetl(fid);
          if isempty(deblank(line))
              continue;
          end
          i = i + 1;
          if i == 1
             numobs = textscan(line,'%f'); 
          end
          if length(line) > 5 && length(line) < 60
             aa = textscan(line,'%s %s %f %f');
          end
          if length(line) > 60
             k2 = k2 + 1;
             liance = textscan(line,'%s %s %f %f %f %f %f %f %f %f %f %s %s %s  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f'); 
             xdobs(k2,1) = strcat(aa{1},aa{2});
             xdobs(k2,2:41) = liance;
          end          
end
j = 1;
while j <length(xdobs(:,1))
    if j>1 && strcmp(xdobs{j-1,2},xdobs{j,2}) && strcmp(xdobs{j-1,1},xdobs{j,1})
        xdobs(j-1,:) = [];
    end
    j = j + 1;
end
fclose(fid);
%���������������ݶ�ȡ���
%******************************************
%��ʼ����������
yiqi_name = xdobs{1,1};
%******************************************
%���峱����
for i = 1:length(xdobs(:,1))
    %���������ĸ�ֵת��
    if strcmp(yiqi_name,xdobs{i,1})~=0
    else
       yiqi_name = xdobs{i,1};
    end
    fid = fopen(strcat(numtgezhi(i).folder,'\',strcat(yiqi_name ,'.txt')),'r');
    gezhi = textscan(fid,'%f %f %f');
    fclose(fid);
    posi = find((floor(xdobs{i,8}/100)*100) == gezhi{1,1});
    %********************************************
    %�����ֵת����Ĺ۲����(΢��ֵ)
    getr(i) = gezhi{1,2}(posi) + (xdobs{i,8}- (floor(xdobs{i,8}/100))*100)*gezhi{1,3}(posi);
    %���峱����(΢��ֵ)
    deg(i) = floor(xdobs{i,end})+(xdobs{i,end}-floor(xdobs{i,end}))/0.6;
%   if deg(i)-xdobs{i,7}>0.001
%       error('�ȷ���ת��С�����ڴ�������ԭʼ����');
%   end
    geodeticfai = dms2degrees([xdobs{i,22:24}]);
    geodeticlamda = dms2degrees([xdobs{i,19:21}]);
    gutu(i) = guitichao([xdobs{i,4:6} deg(i)],geodeticfai,geodeticlamda); 
    %�����߸���(΢��ֵ)
    yiqigai(i) = 308.6*xdobs{i,10}/1000/1000;  %�����ߵ�λΪ������%ע����ʵ��������ݶ�ֵ,���ĳ�΢��ֵ
    %��ѹ����(΢��ֵ)
    qiya(i) = 0.3*(xdobs{i,9}-1.01325*10^3*(1-0.0065*xdobs{i,12}/288.15)^5.2559)/1000;
    %��������  
end
%*************************************************
%��Ư�ʼ������Ư����
yiqi_name = xdobs(1,1);
yiqi_num = 0;st_num = 0;
num_obs = max(1:aa{3}:length(xdobs(:,1)));
for i = 1:aa{3}:length(xdobs(:,1))
   %��ͬ����ģʽ����Ư���������ͬ
   if aa{3} == 3
       %��Ư�ʼ���:����A-B-A����ģʽ
       gRC1 =  (getr(i)  +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));    %������������۲�ֵ
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC2 =  (getr(i+1) +  gutu(i+1) /1000+ yiqigai(i+1) +  qiya(i+1));
       k = -(gRC3-gRC1)/(deg(i+2)- deg(i));           %������Ư��
       gRCend1 = (gRC1+ k*(deg(i)- deg(i)))*10^3;     %ת��Ϊ����ֵ
       gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3; %ת��Ϊ����ֵ
       gRCend3 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3; %ת��Ϊ����ֵ
       %********************%
       %�β��Ƿ�����˸�ֵϵ��
%      duancha = gezhi{1,1}(1)*((gRCend1 - gRCend2) + (gRCend3 -gRCend2))/2;
       %********************%
       duancha = ((gRCend2 - gRCend1) + (gRCend2 -gRCend3))/2;  %����β�
       %�洢����
       adj((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
   elseif aa{3} == 4
       %��Ư�ʼ���:����A-B-B-A����ģʽ
       %������������۲�ֵ
       gRC1 =  (getr(i)  +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));
       gRC2 =  (getr(i+1) +  gutu(i+1)/1000+ yiqigai(i+1) +  qiya(i+1));
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC4 =  (getr(i+3) +  gutu(i+3)/1000 + yiqigai(i+3) +  qiya(i+3));
       k = -(gRC4-gRC1-(gRC3-gRC2))/(deg(i+3)- deg(i)-(deg(i+2)- deg(i+1))); %������Ư��
       gRCend1 = (gRC1+ k*(deg(i)- deg(i)))*10^3;     %ת��Ϊ����ֵ
       gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3; %ת��Ϊ����ֵ
       gRCend3 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3; %ת��Ϊ����ֵ
       gRCend4 = (gRC4 + k*(deg(i+3) - deg(i)))*10^3; %ת��Ϊ����ֵ
       %�β��Ƿ�����˸�ֵϵ��
%      duancha = gezhi{1,1}(1)*((gRCend1 - gRCend2) + (gRCend3 -gRCend2))/2;
       %********************%
       duancha = ((gRCend2 - gRCend1) + (gRCend3 -gRCend4))/2; %����β�
       %�洢����
       adj((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
   elseif aa{3} == 5
       %��Ư�ʼ���:����A-B-C-B-A����ģʽ
       %������������۲�ֵ
       gRC1 =  (getr(i)   +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));
       gRC2 =  (getr(i+1) +  gutu(i+1)/1000 + yiqigai(i+1) +  qiya(i+1));
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC4 =  (getr(i+3) +  gutu(i+3)/1000 + yiqigai(i+3) +  qiya(i+3));
       gRC5 =  (getr(i+4) +  gutu(i+4)/1000 + yiqigai(i+4) +  qiya(i+4));
       %*****************************************************************
       %A-B��
       k1 = -(gRC5-gRC1 )/(deg(i+4)- deg(i));           %������Ư��
       gRCend1 = (gRC1 + k1*(deg(i)- deg(i)))*10^3;     %ת��Ϊ����ֵ
       gRCend2 = (gRC2 + k1*(deg(i+1) - deg(i)))*10^3;  %ת��Ϊ����ֵ
       duancha1 = gRCend2 - gRCend1;                    %����β�
       %*****************************************************************
       %B-C��
       k2 = -(gRC4-gRC2 )/(deg(i+3)- deg(i+1));  
       gRCend3 = (gRC2 + k1*(deg(i+1)- deg(i)))*10^3;  %ת��Ϊ����ֵ
       gRCend4 = (gRC3 + k1*(deg(i+2) - deg(i)))*10^3;   %ת��Ϊ����ֵ
       duancha2 = gRCend4 - gRCend3;                     %����β�
       adj1((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha1, k1,xdobs{i,26},1];
       adj2((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i+1,2},xdobs{i+2,2},xdobs{i+1,8},xdobs{i+2,8},gRCend3,gRCend4,duancha2, k1,xdobs{i,26},1];
       adj((i-1)/aa{3}+1,:) =  adj1((i-1)/aa{3}+1,:) ;
       adj((num_obs-1)/aa{3}+1+(i-1)/aa{3}+1,:) =  adj2((i-1)/aa{3}+1,:);
   end
end
for i = 1:length(adj(:,1))
    if length(num2str(adj{i,3})) == 1
      adj{i,3} = strcat('0',num2str(adj{i,3}));
    else
      adj{i,3} = num2str(adj{i,3});
    end
    if length(num2str(adj{i,4})) == 1
        adj{i,4} = strcat('0',num2str(adj{i,4}));
    else
        adj{i,3} = num2str(adj{i,4});
    end
end


