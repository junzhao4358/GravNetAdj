%�ú���ʵ�ּ��Ƹ���
function deltgp = jiyi(lamda,fai)
MI = pi/180;
x = ; y = ; %�ؼ����꣬����IERS��������ֵ
omga = 7292115*10^-11; %������ת�ٶ�
a =  6378136; %���򳤰���
deltgp = -1.164*10^8*omga^2*a*sin(2*fai*MI)*(x*cos(lamda*MI) -y *sin(lamda*MI));
