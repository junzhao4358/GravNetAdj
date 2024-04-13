function [xdobs,deg,getr,gutu,yiqigai,qiya,ocean_load,kk,sst] = Gravity_correction(numtObs,numtgezhi,Coord,gradient,app.OceanTideCheckBox.Value)
%******************************************
format long 
% numtgezhi = dir('./gezhi/*.txt');           %�������ȡtxt�ļ�
% numfiles= length(numtObs);                  %�õ��ļ��ĸ���
% feature('DefaultCharacterSet','UTF-8');      %�޸ı����ʽ����ֹ��ȡ����ʱ��������
% feature('DefaultCharacterSet','GBK');
k2 = 0;
k1 = 0;
member = ['A','B','C','D','S'];
%��ȡȫ����������������
fid = fopen(strcat(numtObs.folder,'\',numtObs.name),'r', 'n','UTF-8');
[filename,~,~,encoding] = fopen(fid);
if strcmp(encoding,'GBK')
   feature('DefaultCharacterSet','GBK');
else
   feature('DefaultCharacterSet','UTF-8');
end
i = 0;
sst = zeros(100,1);
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
         k1 = k1 + 1;
         aa(k1,:) = textscan(line,'%s %s %f %f');
      end
      if length(line) > 60
         k2 = k2 + 1;
         sst(k1) = sst(k1) + 1;
         liance = textscan(line,'%s %s %f %f %f %f %f %f %f %f %f %s %s %s  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f'); 
         if  double(liance{37}) == 1
             k2 = k2 -1;
             sst(k1) = sst(k1) - 1;
         else      
            xdobs(k2,1) = liance{12};
            xdobs(k2,2:41) = liance;
            if double(ismember(member,liance{1}{1}(1)))== zeros(4,1)
                cc = liance{1}{1};
                cc(1) = member(str2num(cc(1)));
                xdobs{k2,2} = cc;
            end
        end
     end          
end
sst(find(sst == 0)) = [];
% j = 1;
% % while j <length(xdobs(:,1))
% %     if j>1 && strcmp(xdobs{j-1,2},xdobs{j,2}) && strcmp(xdobs{j-1,1},xdobs{j,1})
% %         xdobs(j-1,:) = [];
% %     end
% %     j = j + 1;
% % end
fclose(fid);
numk = unique([aa{:,3}]);
for i = 1:length(numk)
    abk(i) = length(find((numk(i) == [aa{:,3}]) == 1));
end
kk = numk(find((max(abk) == abk) == 1));
if k2/k1 ~=kk
    kk = k2/k1;
end
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
    fid = fopen(strcat(numtgezhi(1).folder,'\',strcat(yiqi_name ,'.txt')),'r');
    gezhi = textscan(fid,'%f %f %f');
    fclose(fid);
    posi = find((floor(xdobs{i,8}/100)*100) == gezhi{1,1});
    %********************************************
    %�����ֵת����Ĺ۲����(΢��ֵ)
    getr(i) = gezhi{1,2}(posi) + (xdobs{i,8}- (floor(xdobs{i,8}/100))*100)*gezhi{1,3}(posi);
    %���峱����(΢��ֵ)
    time_st(i) = floor(xdobs{i,end})+(xdobs{i,end}-floor(xdobs{i,end}))/0.6;
    deg(i) = Julian(xdobs{i,4},xdobs{i,5},xdobs{i,6})*24 + time_st(i);
%   if deg(i)-xdobs{i,7}>0.001
%       error('�ȷ���ת��С�����ڴ�������ԭʼ����');
%   end
    if isempty(Coord)
       %���峱����
       if xdobs{i,23}>=60 || xdobs{i,24}>=60
          geodeticfai = xdobs{i,22} + xdobs{i,23}/100 + xdobs{i,24}/10000;
       else
          geodeticfai = dms2degrees([xdobs{i,22:24}]);
       end
       if xdobs{i,20}>=60 || xdobs{i,21}>=60
          geodeticlamda = xdobs{i,19} + xdobs{i,20}/100 + xdobs{i,21}/10000;
       else
          geodeticlamda = dms2degrees([xdobs{i,19:21}]);
       end
       %��ѹ����(΢��ֵ)
       qiya(i) = 0.3*(xdobs{i,9}-1.01326*10^3*(1-0.0065*xdobs{i,12}/288.15)^5.2559)/1000;
    else
        possi = find(double(strcmp(Coord{1},xdobs{i,2})) == 1);
        if isempty(possi)
           if xdobs{i,23}>=60 || xdobs{i,24}>=60
               geodeticfai = xdobs{i,22} + xdobs{i,23}/100 + xdobs{i,24}/10000;
           else
               geodeticfai = dms2degrees([xdobs{i,22:24}]);
           end
            if xdobs{i,20}>=60 || xdobs{i,21}>=60
               geodeticlamda = xdobs{i,19} + xdobs{i,20}/100 + xdobs{i,21}/10000;
            else
               geodeticlamda = dms2degrees([xdobs{i,19:21}]);
            end
           %��ѹ����(΢��ֵ)
            qiya(i) = 0.3*(xdobs{i,9}-1.01326*10^3*(1-0.0065*xdobs{i,12}/288.15)^5.2559)/1000;
        else
            possi = find(double(strcmp(Coord{1},xdobs{i,2})) == 1);
            geodeticfai = dmstodeg(Coord{2}(possi));
            geodeticlamda = dmstodeg(Coord{3}(possi));
            qiya(i) = 0.3*(xdobs{i,9}-1.01326*10^3*(1-0.0065*Coord{4}(possi)/288.15)^5.2559)/1000;
        end
    end
    gutu(i) = guitichao([xdobs{i,4:6} time_st(i)],geodeticfai,geodeticlamda); 
    %�����߸���(΢��ֵ)
    if isempty(gradient)
       yiqigai(i) = 308.6*xdobs{i,10}/1000/1000;  %�����ߵ�λΪ������%ע����ʵ��������ݶ�ֵ,���ĳ�΢��ֵ
    else
       weizhi = find(double(strcmp(gradient{1},xdobs{i,2}) == 1));
       if isempty(weizhi)
           yiqigai(i) = 308.6*xdobs{i,10}/1000/1000; 
       else
           yiqigai(i) = xdobs{i,10}/1000*(-0.1)*gradient{2}(weizhi);
       end
    end
    if app.OceanTideCheckBox.Value   %��������  
        file_raod = strcat(numtObs.folder,'\',spotl,'\');
        if isempty(Coord)
           str_road = writetoocean(file_raod,geodeticfai,geodeticlamda,xdobs{i,12},xdobs{i,4},xdobs{i,5},xdobs{i,6},xdobs{i,7});
           [status,cmdout] = system(str_road);
           ocean_load(i) = load(strcat(file_raod ,'haichao.txt')); 
        else
           possi = find(double(strcmp(Coord{1},xdobs{i,2})) == 1);
           str_road = writetoocean(file,raod,geodeticfai,geodeticlamda,Coord{4}(possi),xdobs{i,4},xdobs{i,5},xdobs{i,6},xdobs{i,7});
           [status,cmdout]  = system(str_road);
           ocean_load(i) = load(strcat(file_raod ,'haichao.txt')); 
        end
    else
        ocean_load(i)  =0;
    end  
end
%*************************************************