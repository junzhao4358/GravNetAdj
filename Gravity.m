clear
clc
%����������ݴ������
% feature('DefaultCharacterSet','UTF8');  %�޸ı����ʽ����ֹ��ȡ����ʱ��������
%feature('DefaultCharacterSet','GBK');
numtObs = dir('./liance/*.txt');          %�������ȡtxt�ļ�
numtgezhi = dir('./gezhi/*.txt');         %�������ȡtxt�ļ�
fid  = fopen('ƽ���ʽ����.txt','w');
for i = 1:length(numtObs)
   adj = Gravity_Lice(numtObs(i));
   for j = 1:length(adj(:,1))
       fprintf(fid,'%5s  %d%s%s %5s %5s %10.3f %10.3f %14.3f %14.3f %15.3f %10.5f %3d %3d\n',adj{j,:});
   end
end
fclose(fid)

