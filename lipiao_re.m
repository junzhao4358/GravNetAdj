function  asd = lipiao_re(xdobs,deg,getr,gutu,yiqigai,qiya,kk,app_SolidTide,app_AtmosphericPressure,app_InstrumentHight)         
%零漂率计算和零漂改正
yiqi_name = xdobs(1,1);
yiqi_num = 0; st_num = 0;
num_Obss = length(xdobs(:,1));
num_obs = max(1:kk:num_Obss);
adj = {};adj1 ={};adj2 ={};
for i = 1:kk:num_Obss
%不同测量模式下零漂计算改正不同
   if kk == 3
      gRC1 = getr(i) ; gRC3  = getr(i+2); gRC2 = getr(i+1);
   %零漂率计算:采用A-B-A测量模式
%     if app.SolidTideCheckBox.Value
      if app_SolidTide
         gRC1  = gRC1  + gutu(i)/1000;
         gRC3  = gRC3  + gutu(i+2)/1000;
         gRC2  = gRC2  + gutu(i+1)/1000;
      end
%     if app.AtmosphericPressureCheckBox.Value
      if app_AtmosphericPressure
         gRC1  = gRC1  + qiya(i);
         gRC3  = gRC3  + qiya(i+2);
         gRC2  = gRC2  + qiya(i+1);
      end
%     if app.InstrumentHightCheckBox.Value
      if app_InstrumentHight
         gRC1  = gRC1  + yiqigai(i) ;
         gRC3  = gRC3  + yiqigai(i+2);
         gRC2  = gRC2  + yiqigai(i+1);
      end
      k = -(gRC3-gRC1)/(deg(i+2)- deg(i));           %计算零漂率
      gRCend1 = (gRC1+ k*(deg(i)- deg(i)))*10^3;     %转换为毫珈值
      gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3; %转换为毫珈值
      gRCend3 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3; %转换为毫珈值
      %********************%
      %段差是否乘以了格值系数
      %duancha = gezhi{1,1}(1)*((gRCend1 - gRCend2) + (gRCend3 -gRCend2))/2;
      %********************%
      duancha = ((gRCend2 - gRCend1) + (gRCend2 -gRCend3))/2;  %计算段差
      %存储数据
      adj((i-1)/kk+1,:) = [xdobs{i,1},xdobs(i,4:6),xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
    elseif kk == 4
           %零漂率计算:采用A-B-B-A测量模式
           %计算各项改正后观测值
           gRC1=getr(i);gRC2 =getr(i+1);gRC3=getr(i+2); gRC4=getr(i+3);
           if app_SolidTide
              gRC1  = gRC1  + gutu(i)/1000;
              gRC3  = gRC3  + gutu(i+2)/1000;
              gRC2  = gRC2  + gutu(i+1)/1000;
              gRC4  = gRC4  + gutu(i+3)/1000;
           end
           if app_AtmosphericPressure
              gRC1  = gRC1  + qiya(i);
              gRC3  = gRC3  + qiya(i+2);
              gRC2  = gRC2  + qiya(i+1);
              gRC4  = gRC4  + qiya(i+3);
           end
           if app_InstrumentHight
              gRC1  = gRC1  + yiqigai(i) ;
              gRC3  = gRC3  + yiqigai(i+2);
              gRC2  = gRC2  + yiqigai(i+1);
              gRC4  = gRC4  + yiqigai(i+3);
            end
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
            adj((i-1)/kk+1,:) = [xdobs{i,1},xdobs(i,4:6),xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1];
     elseif kk == 5
            %零漂率计算:采用A-B-C-B-A测量模式
            %计算各项改正后观测值
            gRC1 = getr(i) ; gRC2  = getr(i+1); gRC3=getr(i+2); gRC4=getr(i+3);gRC5=getr(i+4);
            if app_SolidTide
               gRC1  = gRC1  + gutu(i)/1000;
               gRC3  = gRC3  + gutu(i+2)/1000;
               gRC2  = gRC2  + gutu(i+1)/1000;
               gRC4  = gRC4  + gutu(i+3)/1000;
               gRC5  = gRC5  + gutu(i+4)/1000;
             end
             if app_AtmosphericPressure
                gRC1  = gRC1  + qiya(i);
                gRC3  = gRC3  + qiya(i+2);
                gRC2  = gRC2  + qiya(i+1);
                gRC4  = gRC4  + qiya(i+3);
                gRC5  = gRC5  + qiya(i+4);
              end
              if app_InstrumentHight
                 gRC1  = gRC1  + yiqigai(i) ;
                 gRC3  = gRC3  + yiqigai(i+2);
                 gRC2  = gRC2  + yiqigai(i+1);
                 gRC4  = gRC4  + yiqigai(i+3);
                 gRC5  = gRC5  + yiqigai(i+4);
              end
               %*****************************************************************
               %A-B段
               k = -(gRC5-gRC1 )/(deg(i+4)- deg(i));           %计算零漂率
               gRCend1 = (gRC1 + k*(deg(i)- deg(i)))*10^3;     %转换为毫珈值
               gRCend2 = (gRC2 + k*(deg(i+1) - deg(i)))*10^3;  %转换为毫珈值
               duancha1 = gRCend2 - gRCend1;                   %计算段差
               %*****************************************************************
               %B-C段
               gRCend3 = (gRC2 + k*(deg(i+1)- deg(i)))*10^3;   %转换为毫珈值
               gRCend4 = (gRC3 + k*(deg(i+2) - deg(i)))*10^3;  %转换为毫珈值
               duancha2 = gRCend4 - gRCend3;                   %计算段差
               adj1((i-1)/kk+1,:) = [xdobs{i,1},xdobs(i,4:6),xdobs{i,2},xdobs{i+1,2},xdobs{i,8},xdobs{i+1,8},gRCend1,gRCend2,duancha1, k,xdobs{i,26},1];
               adj2((i-1)/kk+1,:) = [xdobs{i,1},xdobs(i,4:6),xdobs{i+1,2},xdobs{i+2,2},xdobs{i+1,8},xdobs{i+2,8},gRCend3,gRCend4,duancha2, k,xdobs{i,26},1];
               adj((i-1)/kk+1,:) =  adj1((i-1)/kk+1,:);
               adj((num_obs-1)/kk+1+(i-1)/kk+1,:) =  adj2((i-1)/kk+1,:);
   else     
            gRC = getr;
            if app_SolidTide
               gRC = gRC + gutu/1000;
            end
            if app_AtmosphericPressure
               gRC  = gRC +  qiya;
            end
            if app_InstrumentHight
               gRC  = gRC +  yiqigai;
            end
            if mod(kk,2) == 1
               k = -(gRC(i) -  gRC(i+kk-1))/(deg(i)- deg(i+kk-1));  %零漂率
               num_ducha = (kk-1)/2;
               for j = 1:num_ducha 
                   gRCend1 = (gRC((i-1)+j) + k*(deg((i-1)+j) - deg((i-1)+1)))*10^3;     %转换为毫珈值
                   gRCend2 = (gRC((i-1)+j+1) + k*(deg((i-1)+j+1) - deg((i-1)+1)))*10^3;  %转换为毫珈值
                   duancha = gRCend2 - gRCend1;
                   adj = [adj;[xdobs{i,1},xdobs(i,4:6),xdobs{(i-1)+j,2},xdobs{(i-1)+j+1,2},xdobs{(i-1)+j,8},xdobs{(i-1)+j+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1]];
               end
            else
                for j = 1:kk-1
                    if strcmp(xdobs{i+j-1,2},xdobs{i+j,2})
                        break
                    end
                end
                k = -(gRC(i+kk-1) - gRC(i)-(gRC(i+j) - gRC(i+j-1)))/(deg(i+kk-1)-deg(i)- (deg(i+j) - deg(i+j-1)));  %零漂率;
                num_ducha = (kk-1)/2;
                for j = 1:num_ducha 
                   gRCend1 = (gRC((i-1)+j) + k*(deg((i-1)+j) - deg((i-1)+1)))*10^3;     %转换为毫珈值
                   gRCend2 = (gRC((i-1)+j+1) + k*(deg((i-1)+j+1) - deg((i-1)+1)))*10^3;  %转换为毫珈值
                   duancha = gRCend2 - gRCend1;
                   adj = [adj;[xdobs{i,1},xdobs{i,4:6},xdobs{(i-1)+j,2},xdobs{(i-1)+j+1,2},xdobs{(i-1)+j,8},xdobs{(i-1)+j+1,8},gRCend1,gRCend2,duancha, k,xdobs{i,26},1]];
                end
            end
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
       adj{i,4} = num2str(adj{i,4});
    end
    asd{i,1} = adj{i,1}; 
    asd{i,2} = strcat(num2str(adj{i,2}),adj{i,3},adj{i,4});
    asd(i,3:12) = adj(i,5:14);
 end