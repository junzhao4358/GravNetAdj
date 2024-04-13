function deltag = gravity_earthtide(ymdh,Longtitude,Latitude)
%理论固体潮函数
%output
MI = pi/180;
y = ymdh(1);
m = ymdh(2);
d = ymdh(3);
t = ymdh(4);
Phi = Latitude;
%换算
Phip = Phi - 0.192424* sin(2*Phi*MI);
Longtitude = Longtitude*MI;
Latitude  = Latitude*MI;
Phip = Phip*MI;
To = Julian(y,m,d);
T = (To - 2415020 + (t -8)/24)/36525;
%天文参数
S = 270.43659 + 481267.89057*T - 0.00198*T^2 + 0.000002*T^3;
H = 279.69668 + 36000.76892*T + 0.0003*T^2;
P = 334.32956 + 4069.03403*T -0.01032*T^2 - 0.00001*T^3;
N = 259.18328 - 1934.14201*T + 0.00208*T^2 +0.000002*T^3;
Ps = 281.22083 + 1.71902*T + 0.00045*T^2 + 0.000003*T^3;
sigma = 23.45229 - 0.01301*T - 0.000002*T^2;
S = S*MI;
h = H*MI;
p = P*MI;
N = N*MI;
ps = Ps*MI;
e = sigma*MI;
crm = 1 + 0.0545*cos(S-p) + 0.0030*cos(2*(S -p)) + 0.01*cos(S- 2*h+p)+0.0082*cos(2*(S-h))+0.0006*cos(2*S - 3*h +ps) + 0.0009*cos(3*S-2*h-p);
Lambdam = S -0.0032*sin(h-ps)-0.001*sin(2*h-2*p)+ 0.001*sin(S-3*h+p+ps)+ 0.0222*sin(S-2*h+p)+0.0007*sin(S-h-p+ps)-0.0006*sin(S-h)+ 0.1098*sin(S-p) - 0.0005*sin(S+h -p -ps) + 0.0008*sin(2*S - 3*h +ps)+ 0.0115*sin(2*S-2*h)+ 0.0037*sin(2*S-2*p) -0.0020*sin(2*S-2*N) + 0.0009*sin(3*S - 2*h  -p);
Beltam = -0.0048*sin(p-N)-0.0008*sin(2*h -p -N) +0.003*sin(S -2*h+N)+0.0895*sin(S-N) + 0.001*sin(2*S - 2*h + p -N)+0.0049*sin(2*S -p -N)+0.0006*sin(2*S -2*h-N);
Theta = (t -8)*15*MI + h + Longtitude - pi;
Delta = sin(e)*sin(Lambdam)*cos(Beltam) + cos(e)*sin(Beltam);
H = cos(Beltam)*cos(Lambdam)*cos(Theta) +sin(Theta)*(cos(Beltam)*sin(Lambdam) - sin(e)*sin(Beltam));
Zm = sin(Phip)*Delta + cos(Phip)*H;
crs = 1 + 0.0168*cos(h-ps)+0.0003*cos(2*h-2*ps);
Lambdas = h +0.0335*sin(h-ps)+0.0004*sin(2*h-2*ps);
Beltas = 0;
Zs = sin(Phip)*sin(e)*sin(Lambdas) +cos(Phip)*(cos(Lambdas)*cos(Theta) + sin(Theta)*cos(e)*sin(Lambdas));
F = 0.998327 + 0.00167*cos(2*Phi);
Gt = -165.17 *F*crm*crm*crm*(Zm*Zm - 1/3) - 1.37*F*F*crm^4*Zm*(5*Zm*Zm-3) -76.08*F*crs^3*(Zs^2 - 1/3);
deltfc = -4.83 + 15.73*sin(Phi)^2 - 1.59*sin(Phi)^4;
deltag = -Gt*1.16+deltfc;
end


