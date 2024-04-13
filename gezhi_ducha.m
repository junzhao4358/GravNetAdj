function AA = gezhi_ducha(xdObsN,xdObs,gezhii)
%利用仪器格值系数，得到更新的段差
for i =1:length(xdObsN(:,1))
    xdObs(i,5) = xdObs(i,5)*gezhii{2}(find(strcmp(gezhii{1},xdObsN{i,1})==1));
end
    
cedui(1,1) = xdObsN(1,2);cedui(1,2) = xdObsN(1,3);
unObs(1) = xdObsN(1,1);
kk = 1;ss = 0;
st = 0;
num_posi = 0;
AA = {};
delete_po = [];
for i = 1:length(xdObsN(:,1))
    if double(strcmp(xdObsN(i,2),cedui(1,1))) ==1 && double(strcmp(xdObsN(i,3),cedui(1,2))) == 1
        st = st + 1;
    else
        ss  = ss + 1;
        AA(ss,1) = cedui(1,1);  %起始点号
        AA(ss,2) = cedui(1,2);  %结束点号
        AA{ss,3} = st;
        st = st - 1;
        AA{ss,4} = mean(xdObs(kk:kk+st,5));
        AA{ss,5} = xdObs(kk+st,7);
        %注意算法不一样，计算结果不一样
        if strcmp(strcat(AA(ss,2),AA(ss,1)),strcat(AA(:,1),AA(:,2))) == zeros(length(AA(:,1)),1)
        else
        posi_XD = find(double(strcmp(strcat(AA(ss,2),AA(ss,1)),strcat(AA(:,1),AA(:,2)))) == 1);
        sAA = find(double(strcmp(strcat(xdObsN(:,2),xdObsN(:,3)),strcat(AA(ss,2),AA(ss,1)))) == 1);
        AA{posi_XD,3} = AA{posi_XD,3} + AA{ss,3};
        AA{posi_XD,4} = (sum(xdObs(kk:kk+st,5)) - sum(xdObs(sAA,5)))/AA{posi_XD,3};
        AA{posi_XD,6} = sqrt(sum(([xdObs(kk:kk+st,5)' -xdObs(sAA,5)'] - AA{posi_XD,4}).^2)/AA{posi_XD,3}/(AA{posi_XD,3}-1));
        num_posi = num_posi + 1;
        delete_po(num_posi) = ss;
        end
        AA{ss,6} = sqrt(sum((xdObs(kk:kk+st,5) - AA{ss,4}).^2)/AA{ss,3}/(AA{ss,3}-1));
        kk =  kk + st + 1;
        cedui(1,1) = xdObsN(i,2);
        cedui(1,2) = xdObsN(i,3);
        st = 1;
    end
    if i == length(xdObsN(:,1))
       ss = ss + 1;
       AA(ss,1) = cedui(1,1);  %起始点号
       AA(ss,2) = cedui(1,2);  %结束点号
       AA{ss,3} = st;
       st = st - 1;
       AA{ss,4} = mean(xdObs(kk:kk+st,5));
       AA{ss,5} = xdObs(kk+st,7);
       AA{ss,6} = sqrt(sum((xdObs(kk:kk+st,5) - AA{ss,4}).^2)/AA{ss,3}/(AA{ss,3}-1));
     end
end
if isempty(delete_po)
else
AA(delete_po,:) = [];
end