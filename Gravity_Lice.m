function adj = Gravity_Lice(numtObs,numtgezhi)
%******************************************
format long 
% numtgezhi = dir('./gezhi/*.txt');         %批处理读取txt文件
numfiles= length(numtObs);                %得到文件的个数
%feature('DefaultCharacterSet','UTF8');  %修改编码格式，防止读取中文时出现乱码
%feature('DefaultCharacterSet','GBK');
k2 = 0;
%读取全部的重力联测数据
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
%所有重力联测数据读取完毕
%******************************************
%开始计算各项改正
yiqi_name = xdobs{1,1};
%******************************************
%固体潮改正
for i = 1:length(xdobs(:,1))
    %仪器读数的格值转换
    if strcmp(yiqi_name,xdobs{i,1})~=0
    else
       yiqi_name = xdobs{i,1};
    end
    fid = fopen(strcat(numtgezhi(i).folder,'\',strcat(yiqi_name ,'.txt')),'r');
    gezhi = textscan(fid,'%f %f %f');
    fclose(fid);
    posi = find((floor(xdobs{i,8}/100)*100) == gezhi{1,1});
    %********************************************
    %计算格值转换后的观测读数(微加值)
    getr(i) = gezhi{1,2}(posi) + (xdobs{i,8}- (floor(xdobs{i,8}/100))*100)*gezhi{1,3}(posi);
    %固体潮改正(微加值)
    deg(i) = floor(xdobs{i,end})+(xdobs{i,end}-floor(xdobs{i,end}))/0.6;
%   if deg(i)-xdobs{i,7}>0.001
%       error('度分秒转度小数存在错误，请检查原始数据');
%   end
    geodeticfai = dms2degrees([xdobs{i,22:24}]);
    geodeticlamda = dms2degrees([xdobs{i,19:21}]);
    gutu(i) = guitichao([xdobs{i,4:6} deg(i)],geodeticfai,geodeticlamda); 
    %仪器高改正(微加值)
    yiqigai(i) = 308.6*xdobs{i,10}/1000/1000;  %仪器高单位为：毫米%注意用实测的重力梯度值,均改成微加值
    %气压改正(微加值)
    qiya(i) = 0.3*(xdobs{i,9}-1.01325*10^3*(1-0.0065*xdobs{i,12}/288.15)^5.2559)/1000;
    %海潮改正  
end
%*************************************************
%零漂率计算和零漂改正
yiqi_name = xdobs(1,1);
yiqi_num = 0;st_num = 0;
num_obs = max(1:aa{3}:length(xdobs(:,1)));
for i = 1:aa{3}:length(xdobs(:,1))
   %不同测量模式下零漂计算改正不同
   if aa{3} == 3
       %零漂率计算:采用A-B-A测量模式
       gRC1 =  (getr(i)  +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));    %计算各项改正后观测值
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC2 =  (getr(i+1) +  gutu(i+1) /1000+ yiqigai(i+1) +  qiya(i+1));
       k = -(gRC3-gRC1)/(deg(i+2)- deg(i));           %计算零漂率
       gRCend1 = (gRC1+ k*(deg(i)- deg(i)))*10^3;     %转换为毫珈值
       gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3; %转换为毫珈值
       gRCend3 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3; %转换为毫珈值
       %********************%
       %段差是否乘以了格值系数
%      duancha = gezhi{1,1}(1)*((gRCend1 - gRCend2) + (gRCend3 -gRCend2))/2;
       %********************%
       duancha = ((gRCend2 - gRCend1) + (gRCend2 -gRCend3))/2;  %计算段差
       %存储数据
       adj((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
   elseif aa{3} == 4
       %零漂率计算:采用A-B-B-A测量模式
       %计算各项改正后观测值
       gRC1 =  (getr(i)  +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));
       gRC2 =  (getr(i+1) +  gutu(i+1)/1000+ yiqigai(i+1) +  qiya(i+1));
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC4 =  (getr(i+3) +  gutu(i+3)/1000 + yiqigai(i+3) +  qiya(i+3));
       k = -(gRC4-gRC1-(gRC3-gRC2))/(deg(i+3)- deg(i)-(deg(i+2)- deg(i+1))); %计算零漂率
       gRCend1 = (gRC1+ k*(deg(i)- deg(i)))*10^3;     %转换为毫珈值
       gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3; %转换为毫珈值
       gRCend3 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3; %转换为毫珈值
       gRCend4 = (gRC4 + k*(deg(i+3) - deg(i)))*10^3; %转换为毫珈值
       %段差是否乘以了格值系数
%      duancha = gezhi{1,1}(1)*((gRCend1 - gRCend2) + (gRCend3 -gRCend2))/2;
       %********************%
       duancha = ((gRCend2 - gRCend1) + (gRCend3 -gRCend4))/2; %计算段差
       %存储数据
       adj((i-1)/aa{3}+1,:) = [xdobs{i,1},xdobs{i,4:6},xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
   elseif aa{3} == 5
       %零漂率计算:采用A-B-C-B-A测量模式
       %计算各项改正后观测值
       gRC1 =  (getr(i)   +  gutu(i)/1000   + yiqigai(i)   +  qiya(i));
       gRC2 =  (getr(i+1) +  gutu(i+1)/1000 + yiqigai(i+1) +  qiya(i+1));
       gRC3 =  (getr(i+2) +  gutu(i+2)/1000 + yiqigai(i+2) +  qiya(i+2));
       gRC4 =  (getr(i+3) +  gutu(i+3)/1000 + yiqigai(i+3) +  qiya(i+3));
       gRC5 =  (getr(i+4) +  gutu(i+4)/1000 + yiqigai(i+4) +  qiya(i+4));
       %*****************************************************************
       %A-B段
       k1 = -(gRC5-gRC1 )/(deg(i+4)- deg(i));           %计算零漂率
       gRCend1 = (gRC1 + k1*(deg(i)- deg(i)))*10^3;     %转换为毫珈值
       gRCend2 = (gRC2 + k1*(deg(i+1) - deg(i)))*10^3;  %转换为毫珈值
       duancha1 = gRCend2 - gRCend1;                    %计算段差
       %*****************************************************************
       %B-C段
       k2 = -(gRC4-gRC2 )/(deg(i+3)- deg(i+1));  
       gRCend3 = (gRC2 + k1*(deg(i+1)- deg(i)))*10^3;  %转换为毫珈值
       gRCend4 = (gRC3 + k1*(deg(i+2) - deg(i)))*10^3;   %转换为毫珈值
       duancha2 = gRCend4 - gRCend3;                     %计算段差
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


