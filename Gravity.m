clear
clc
%相对重力数据处理程序
% feature('DefaultCharacterSet','UTF8');  %修改编码格式，防止读取中文时出现乱码
%feature('DefaultCharacterSet','GBK');
numtObs = dir('./liance/*.txt');          %批处理读取txt文件
numtgezhi = dir('./gezhi/*.txt');         %批处理读取txt文件
fid  = fopen('平差格式数据.txt','w');
for i = 1:length(numtObs)
   adj = Gravity_Lice(numtObs(i));
   for j = 1:length(adj(:,1))
       fprintf(fid,'%5s  %d%s%s %5s %5s %10.3f %10.3f %14.3f %14.3f %15.3f %10.5f %3d %3d\n',adj{j,:});
   end
end
fclose(fid)

