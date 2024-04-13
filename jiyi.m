%该函数实现极移改正
function deltgp = jiyi(lamda,fai)
MI = pi/180;
x = ; y = ; %地极坐标，采用IERS公布的数值
omga = 7292115*10^-11; %地球自转速度
a =  6378136; %地球长半轴
deltgp = -1.164*10^8*omga^2*a*sin(2*fai*MI)*(x*cos(lamda*MI) -y *sin(lamda*MI));
