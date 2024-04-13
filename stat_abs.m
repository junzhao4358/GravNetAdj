function [stat_judui,tt] = stat_abs(jdObsN,jdObs)
%**************************************
%绝对重力数据统计
stat_judui(1,1) = jdObsN(1,1);
st1 = 1;
for i = 1:length(jdObsN(:,1))
     if strcmp(stat_judui,jdObsN(i,1)) == zeros(length(stat_judui),1)
        st1 = st1 + 1;
        stat_judui(st1,1) = jdObsN(i,1);
      end
end
for i = 1:length(stat_judui)
    posi_juidui = find(double(strcmp(stat_judui(i,1),jdObsN)) == 1);
    tt(i,1) = length(posi_juidui);
    tt(i,2) = mean(jdObs(posi_juidui,1));
    if tt(i,1) == 1
       tt(i,3) = 0;
    else
       tt(i,3) = sqrt(sum((jdObs(posi_juidui,1) - tt(i,2)).^2)/tt(i)/(tt(i)-1));
    end
 end