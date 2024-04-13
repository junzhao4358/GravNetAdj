%该函数实现仪器高改正
function deltgh = yiqigao(h,xita)
%h:重力仪面板高度
%xita:重力垂直梯度
%g0:墩面重力值
%gp:重力观测值
deltgh = xita*h;