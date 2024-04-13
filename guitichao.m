function delt = guitichao(ymdh,geodeticfai,geodeticlamda)
%绝对重力数据处理各项改正函数
%固体潮改正
%***************************************
%计算时刻的儒略世纪数
y = ymdh(1);
m = ymdh(2);
d = ymdh(3);
t = ymdh(4);
Mi = pi/180;
T0 = Julian(y,m,d);
T = (T0  - 2415020.0 + (t-8)/24)/36525;
geodeticlamda = geodeticlamda*Mi;
geodeticfai =  geodeticfai *Mi;
%****************************************
%计算六个天文引数
s = 270.43659 + 481267.89057*T + 0.00198*T^2 + 0.000002*T^3;
h = 279.69668 + 36000.76892*T + 0.0003*T^2;
p = 334.32956 + 4069.03403*T - 0.01032*T^2 - 0.00001*T^3;
N = 259.18328 - 1934.14201*T + 0.00208*T^2 + 0.000002*T^3;
ps = 281.22083 + 1.71902*T + 0.00045*T^2 + 0.000003*T^3;
eslon = 23.45229 - 0.01301*T - 0.000002*T^2;
s = s*Mi;
h = h*Mi;
p = p*Mi;
N = N*Mi;
ps = ps*Mi;
eslon = eslon*Mi;
%********************************************************
%计算月亮的c/r及cosZ
c_r = 1 + 0.0545*cos(s-p ) + 0.0030*cos(2*(s-p))+0.01*cos(s-2*h+p) + 0.0082*cos(2*(s-h) )+0.0006*cos(2*s-3*h+ps)+0.0009*cos(3*s -2*h -p );
lamda = s - 0.0032*sin(h-ps) - 0.001*sin(2*h-2*p) + 0.001*sin(s-3*h+p+ps )+0.0222*sin(s-2*h+p)+0.0007*sin(s-h-p+ps )-0.0006*sin(s-h) + 0.1098*sin(s-p) - 0.0005*sin(s+h-p-ps)+0.0008*sin(2*s-3*h+ps) +0.0115*sin(2*s-2*h) + 0.0037*sin(2*s-2*p)-0.0020*sin(2*s-2*N) +0.0009*sin(3*s-2*h-p);
beta = -0.0048*sin(p-N )-0.0008*sin(2*h-p-N) +0.003*sin(s-2*h+N) + 0.0895*sin(s-N) + 0.001*sin(2*s-2*h+p-N) + 0.0049*sin(2*s-p-N)+0.0006*sin(3*s-2*h-N);
sin_delt = sin(eslon) *sin(lamda)*cos(beta) + cos(eslon)*sin(beta);
hexingshi =(t-8)*15*Mi + h + geodeticlamda- pi;
sphericalfai = geodeticfai/Mi - 0.193296*sin(2*geodeticfai);
cos_delt_H = cos(beta)*cos(lamda)*cos(hexingshi) + sin(hexingshi)*(cos(eslon)*cos(beta)*sin(lamda) - sin(eslon)*sin(beta));
cos_Z = sin(sphericalfai*Mi)*sin_delt + cos(sphericalfai*Mi)*cos_delt_H;
%***********************************************
%求太阳的Cs/rs cos(zs)
Cs_rs =  1+0.0168*cos(h-ps)+0.0003*cos(2*h-2*ps);
lamda_s = h + 0.0335*sin(h-ps)+0.0004*sin(2*(h-ps));
cos_zs = sin(sphericalfai*Mi)*sin(eslon)*sin(lamda_s) + cos(sphericalfai*Mi)*(cos(lamda_s)*cos(hexingshi) + sin(hexingshi)*cos(eslon)*sin(lamda_s));
%**********************************************
%固体潮改正计算
F_fai = 0.998327 + 0.00167*cos(2*geodeticfai);
deth = 1.16 ;      %潮汐因子平均值
defc = -4.83 + 15.73 *(sin(sphericalfai*Mi))^2 - 1.59*(sin(sphericalfai*Mi))^4;
Gt  = - 165.17*F_fai *c_r^3*(cos_Z^2 - 1/3) - 1.37*F_fai ^2*c_r^4*cos_Z*(5*cos_Z^2-3) - 76.08 *F_fai*Cs_rs ^3*(cos_zs^2 - 1/3);
delt = - (deth *Gt - defc);






