function  results = plot_distr(file)
f = waitbar(0,'画图正在进行中！请稍后!');
coord_file = importdata(file);
rode = strsplit(file,'\');
rode1 = '';
for i = 1:length(rode)-1
    rode1 =  strcat(rode1,rode(i),'\');
end
waitbar(0.2,f);
fid1 = fopen(strcat(rode1,'data1.txt'),'w');   %基准点
fid2 = fopen(strcat(rode1,'data2.txt','w');    %基本点
fid3 = fopen(strcat(rode1,'data3.txt','w');    %引点
fid4 = fopen(strcat(rode1,'data4.txt','w');    %岛礁点 
fid5 = fopen(strcat(rode1,'data5.txt','w');    %短基线点
for i = 1:length(coord_file(:,1))
    if strcmp(coord_file{i,1}(1),'A')
       fprintf(fid1,'%15.13f %15.13f',str2num(coord_file{i,3})*100,str2num(cood_file{i,2})*100));
    elseif strcmp(coord_file{i,1}(1),'B')
       fprintf(fid2,'%15.13f %15.13f',str2num(coord_file{i,3})*100,str2num(cood_file{i,2})*100));
    elseif strcmp(coord_file{i,1}(1),'C')
       fprintf(fid3,'%15.13f %15.13f',str2num(coord_file{i,3})*100,str2num(cood_file{i,2})*100));
    elseif strcmp(coord_file{i,1}(1),'D')
       fprintf(fid4,'%15.13f %15.13f',str2num(coord_file{i,3})*100,str2num(cood_file{i,2})*100));
    elseif strcmp(coord_file{i,2}(1),'S')
       fprintf(fid5,'%15.13f %15.13f',str2num(coord_file{i,3})*100,str2num(cood_file{i,2})*100));
    end
end  
waitbar(0.4,f);
fid = fopen(strcat(rode1,'dis_point.bat'),'w');
fprintf(fid,'gmt begin');
fprintf(fid,'    gmt figure global_taizhan png  A+m0.5c  E700');
fprintf(fid,'    gmt set FORMAT_GEO_MAP ddd:mm:ssF MAP_FRAME_WIDTH 5p');
fprintf(fid,'    gmt set FONT_ANNOT_PRIMARY 8p  MAP_TITLE_OFFSET 8p');
fprintf(fid,'    gmt set FONT_TITLE  12p,Helvetica,black  FONT_LABEL 10p,Helvetica,black');
fprintf(fid,'    gmt coast -Rd  -JN7i -Bxa30 -Bya30 -BWSEN  -G');
fprintf(fid,'    gmt grdimage earth_relief_05m.grd  -JN7i  -I+d');
fprintf(fid,'    gmt plot  data1.txt -Sc0.2c -Gred');
fprintf(fid,'    gmt plot  data2.txt -Sa0.2c -Gblue');
fpritnf(fid,'    gmt plot  data3.txt -St0.2c -Ggreen');
fprintf(fid,'    gmt plot  data4.txt -Sx0.2c -Gyellow');
fprintf(fid,'    gmt plot  data5.txt -Si02.c -Gblack');
fprintf(fid,'    echo S 0.3i f0.5+t+l 0.4i/0.3c 2/138/210 2.0p,red  0.7i 基准点>legend.dat');
fprintf(fid,'    echo S 0.3i - 0.5i 2/138/210 1.0p,blue 0.7i 基本点>>legend.dat');
fprintf(fid,'    echo S 0.3i - 0.44i - 1.p,green, -0.7i 引点>>legend.dat');
fprintf(fid,'    echo S 0.3i - 0.50i yellow 1.0p,yellow 0. 7i 岛礁点>>legend.dat');
fprintf(fid,'    echo S 0.3i - 0.50i black 1.0p,black 0. 7i 短基线>>legend.dat');
fprintf(fid,'    gmt legend legend.dat -DjTL+w1.5i+jTl+o0.1c/0.1c -F+gwhite+p0.5p --FONE_ANNOT_PRIMARY = 8p,37');
fprintf(fid,'    gmt colorbar -Crelief -Dx13c/-1.2c+w8c/0.25c+jMR+h  -Bx2000f2000 -By+l"m" ');
fprintf(fid,'gmt end');
waitbar(0.6,f);
[status,cmdout] = dos('dis_point.bat') ;
waitbar(0.8,f);
fclose all 
waitbar(1,f);
delete(f);
msgbox('       　画图完毕!请关闭！')
results = 1;