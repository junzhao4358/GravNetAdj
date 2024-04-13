function   Axd = obsequ_xd(xdObs,xdObsN,unObsN,yipara,yiparaN,stN,conS,numX,gravN,T)
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