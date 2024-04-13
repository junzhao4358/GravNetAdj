%相对重力数据预处理程序
%************************************
%仪器读数R转换成毫伽
R1 = fix(R/100)*100;
gr = F1 +(R- R1)*df1
%*************************************
%气压改正值
deltga = qiya(P,H);
%固体潮改正
delt = gutichao(t,geodeticfai,geodeticlamda,T0);
%仪器高改正
deltgh = yiqigao(h,xita);
gRC = gR + delt + deltga + delgh;
%**************************************
%计算零漂改正
k = -(gRC2 - gRC1)/(t2 - t1); %计算零漂率
deltgz = k*di_time;
%**************************************
%计算段差
