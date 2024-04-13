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
             %绘制残差的饼图
             for i = 1:length(stat_chao)
                 anal(i) = categorical({num2str(stat_chao(i))});
             end
             labels = {'小于-3倍中误差','大于-3倍中误差小于-2倍中误差','大于-2倍中误差小于-1倍中误差','大于-1倍中误差小于0','大于0小于1倍中误差','大于1倍中误差小于2倍中误差','大于2倍中误差小于3倍中误差','大于3倍中误差'};                         
             t = tiledlayout(1,1);
             ax1 = nexttile;
             pie(ax1,stat_chao)
             legend(labels)
             
             title('残差正负号统计图')
             %***********************************
             %统计残差(总残差，相对重力，绝对重力)的正负号个数
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
             %绘制残差的分布图
%            histogram(sv);
             histogram(V(intersect(find(V<150),find(V>-150))));
             title('全部观测残差统计图');xlabel('残差(uGal)');ylabel('个数');
             figure(3)
             %绘制绝对重力残差分布图
             histogram(V(1:numJ));
             title('绝对观测残差统计图');xlabel('残差(uGal)');ylabel('个数');
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
             title('相对观测残差统计图');xlabel('残差(uGal)');ylabel('个数');
             %******************************************
             %残差结果统计（最大值、最小值、平均值）
             stat_V = [max(V) min(V) mean(V)];
             stat_VJ = [max(Vjd) min(Vjd) mean(Vjd)];
             stat_VX = [max(Vxd) min(Vxd) mean(Vxd)];