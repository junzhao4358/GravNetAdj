function success = yiqi_stat(file_raod,yiObs,xdObsN,AA,xdObs)
fid = fopen(file_raod,'w');
yiObs = sort(yiObs);
for i = 1:length(yiObs)
    posi_yiqi = find(double(strcmp(yiObs(i),xdObsN(:,1))) == 1);
    dx = 0;
    numyiqi = length(posi_yiqi);
    for j = 1:numyiqi 
        dd = find(double(strcmp(strcat(xdObsN(posi_yiqi(j),2),xdObsN(posi_yiqi(j),3)),strcat(AA(:,1),AA(:,2)))) == 1);
        if isempty(dd)
           dd = find(double(strcmp(strcat(xdObsN(posi_yiqi(j),2),xdObsN(posi_yiqi(j),3)),strcat(AA(:,2),AA(:,1)))) == 1);
        end  
        yiqivalue(j) = AA{dd,4};
        dx = dx + (xdObs(posi_yiqi(j),5) - yiqivalue(j))^2;
     end
     stat_yiqi(i,1) = numyiqi;
     stat_yiqi(i,2) = sqrt(dx/(stat_yiqi(i,1))/(stat_yiqi(i,1)-1));
     if  i==1
         fprintf(fid,'仪器名称   观测次数   平均精度\n')
     end
     fprintf(fid,'%5s  %5d %7.2f\n',yiObs{i},stat_yiqi(i,:));
end
f = msgbox('仪器参数统计完成:Operation Successfully');
success = 1;