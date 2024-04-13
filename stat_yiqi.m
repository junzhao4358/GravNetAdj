function   [yiE,QyiE,An,alf,QAn,Qalf]= stat_yiqi(yipara,gravX,varX,stN,conS,gravN,numyp)
yiE = zeros(numyp,2+2*7);
              QyiE =  zeros(numyp,2+2*7);
              numY = zeros(numyp,9);
              for i = 1:numyp
                  for j = 1:2
                     if yipara(i,j) == 1
                        if i == 1
                            if j == 1
                               yiE(i,j) = gravX(gravN+j);
                            else
                               yiE(i,j) = gravX(gravN+j)/conS;
                            end
                               QyiE(i,j) = varX(gravN+j,gravN+j);
                        else
                            if j == 1
                               yiE(i,j) = gravX(gravN+sum(sum(stN(1:i-1,:)))+j);
                            else
                               yiE(i,j) = gravX(gravN+sum(sum(stN(1:i-1,:)))+j)/conS;
                            end
                              QyiE(i,j) = varX(gravN+sum(sum(stN(1:i-1,:)))+j,gravN+sum(sum(stN(1:i-1,:)))+j);
                        end
                     end
                  end
                  st = 0;
                  for j = 1:7
                     if yipara(i,2+j) == 1
                        st = st + 1;
                        if i == 1
                           if st >0
                              poN = gravN+stN(1,1);
                              yiE(i,(j-1)*2+1+2:2*j+2) = gravX(poN+(st-1)*2+1:poN+st*2)*conS;
                              QyiE(i,(j-1)*2+1+2:2*j+2) = [varX(poN+(st-1)*2+1,poN+(st-1)*2+1) varX(poN+st*2,poN+st*2)];
                           end
                        else
                           if st >0
                              poN = gravN+sum(sum(stN(1:i-1,:))) + stN(i,1);
                              yiE(i,(j-1)*2+2+1:2*j+2) = gravX(poN+(st-1)*2+1:poN+st*2)*conS; 
                              QyiE(i,(j-1)*2+2+1:2*j+2) = [varX(poN+(st-1)*2+1,poN+(st-1)*2+1)  varX(poN+(st-1)*2+2,poN+(st-1)*2+2)];
                           end
                        end
                     end
                  end
               end
               % ***********************************************
               % 获取仪器参数周期函数的振幅和相位以及相应精度
               An = zeros(numyp,7);
               alf  = zeros(numyp,7);
               QAn  = zeros(numyp,7);
               Qalf = zeros(numyp,7);
               for i = 1:length(yipara(:,1))
                   for j = 1:7
                       if yiE(i,j*2+1) ~=0
                          An(i,j) = sqrt(yiE(i,j*2+1)^2+yiE(i,j*2+2)^2);
                          alf(i,j) = atan(yiE(i,j*2+2)/yiE(i,j*2+1));
                          QAn(i,j) = 1/An(i,j)*sqrt(yiE(i,j*2+1)^2* QyiE(i,j*2+1)^2+ yiE(i,j*2+2)^2* QyiE(i,j*2+2)^2);
                          Qalf(i,j) = 1/An(i,j)^2*sqrt(yiE(i,j*2+1)^2* QyiE(i,j*2+2)^2+ yiE(i,j*2+2)^2* QyiE(i,j*2+1)^2);
                       end
                   end
               end