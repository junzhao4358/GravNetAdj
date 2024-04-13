function [stat_chao,numVV1,numVV2,stat_V,stat_VJ,stat_VX] = stat_res(V,PVV,Vjd,Vxd,numJ)
% feature('DefaultCharacterSet','UTF-8');
slCharacterEncoding('UTF-8')
s1 = 0;s2 = 0;s3 =0;s4 =0;s5=0;s6=0;s7 = 0;s8 = 0;
              for i= 1:length(V)
                  if V(i)<-3*PVV(i)
                     s1 = s1 + 1;
                  elseif V(i)>=-3*PVV(i) && V(i) <-2*PVV(i)
                     s2 = s2 + 1;
                  elseif V(i)>=-2*PVV(i) && V(i) <-1*PVV(i)
                     s3 = s3 + 1;
                  elseif V(i)>=-1*PVV(i) && V(i) <0
                     s4 = s4 + 1;
                  elseif V(i)>=0*PVV(i) && V(i) <1*PVV(i)
                     s5 = s5 + 1;
                  elseif V(i)>=1*PVV(i) && V(i) <2*PVV(i)
                     s6 = s6 + 1;
                  elseif  V(i)>=2*PVV(i) && V(i) <3*PVV(i)
                     s7 = s7 + 1;
                 elseif  V(i) >3*PVV(i)
                    s8 = s8 + 1;
                  end
              end
             stat_chao = [s1 s2 s3 s4 s5 s6 s7 s8]; 
             figure(1)
             %���Ʋв�ı�ͼ
             for i = 1:length(stat_chao)
                 anal(i) = categorical({num2str(stat_chao(i))});
             end
             labels = {'С��-3�������','����-3�������С��-2�������','����-2�������С��-1�������','����-1�������С��0','����0С��1�������','����1�������С��2�������','����2�������С��3�������','����3�������'};                         
             t = tiledlayout(1,1);
             ax1 = nexttile;
             pie(ax1,stat_chao)
             legend(labels)
             
             title('�в�������ͳ��ͼ')
             %***********************************
             %ͳ�Ʋв�(�ܲв�����������������)�������Ÿ���
             numVV1 = length(find(V>0) == 1);
             numVV2 = length(find(V<0) == 1);
             %*************************************
%              numjiajian = 0;
%              for i = 1:length(V)
%                 if V(i) < 150 && V(i) > -150
%                    numjiajian = numjiajian + 1;
%                    sv(numjiajian,1) = V(i);
%                 end
%              end
             figure(2)
             %���Ʋв�ķֲ�ͼ
%            histogram(sv);
             histogram(V(intersect(find(V<150),find(V>-150))));
             title('ȫ���۲�в�ͳ��ͼ');xlabel('�в�(uGal)');ylabel('����');
             figure(3)
             %���ƾ��������в�ֲ�ͼ
             histogram(V(1:numJ));
             title('���Թ۲�в�ͳ��ͼ');xlabel('�в�(uGal)');ylabel('����');
             figure(4) 
             vvv = V(numJ+1:end);
%              numjiajian1 = 0;
%              for i = 1:length(vvv)
%                 if vvv(i) < 200 && vvv(i) > -200
%                    numjiajian1 = numjiajian1 + 1;
%                    svv(numjiajian,1) = vvv(i);
%                 end
%              end
             histogram(vvv(intersect(find(vvv<150),find(vvv>-150))));
             title('��Թ۲�в�ͳ��ͼ');xlabel('�в�(uGal)');ylabel('����');
             %******************************************
             %�в���ͳ�ƣ����ֵ����Сֵ��ƽ��ֵ��
             stat_V = [max(V) min(V) mean(V)];
             stat_VJ = [max(Vjd) min(Vjd) mean(Vjd)];
             stat_VX = [max(Vxd) min(Vxd) mean(Vxd)];