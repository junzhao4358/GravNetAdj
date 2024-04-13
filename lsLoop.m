function success = lsLoop(AAA,m_Sigma,unObsN,file_raod)
fid = fopen(file_raod,'w');
f = waitbar(0,'环闭合差正在计算中');
m_Lnumber = length(AAA(:,1));
PP = ones(m_Lnumber,1);
m_Pnumber = length(unObsN);
for i = 1:m_Lnumber
    for j = 1:2
        if j == 1
           StartP(i) = find(double(strcmp(AAA{i,j},unObsN)) == 1);
        else
            EndP(i) = find(double(strcmp(AAA{i,j},unObsN))== 1);
         end
     end
end
%***********************************************************
S = zeros(m_Pnumber,1);
diff = zeros(m_Pnumber,1);
neighbor = zeros(m_Pnumber,1);
used = zeros(m_Lnumber,1);
NumInOneLoop = 1;
LoopClosureNum = 0;
for i = 1:m_Lnumber 
        k1 = StartP(i);
        k2 = EndP(i);
        if used(i)
            continue;
        end
        [S,diff,neighbor] = FindShortPath(k2,i,m_Lnumber,m_Pnumber,StartP,EndP,AAA,PP);
        if neighbor(k1) > 0
           used(i) = 1;
           LoopClosureNum = LoopClosureNum  + 1;
           p1 = k1;
           while 1
                 p2 = neighbor(p1);
                 lsLoopClosure(NumInOneLoop) = p1;
                 NumInOneLoop = NumInOneLoop + 1;
                 for r  = 1:m_Lnumber 
                     if StartP(r) == p1 && EndP(r) == p2
                         used(r) = 1;
                         break;
                     elseif StartP(r) == p2 && EndP(r) == p1
                         used(r) = 1;
                         break;
                     end
                 end
                 if p2 == k2 
                    break;
                 else
                    p1 = p2;
                 end
             end
             lsLoopClosure(NumInOneLoop) = k2;
             lsLoopClosure(NumInOneLoop+1) = k1;
             %NumInOneLoop = NumInOneLoop + 1;
             W  = AAA{i,4} + diff(k1);
             SS = S(k1) + 1/PP(i);
             fprintf(fid,'环闭合差%3d:',LoopClosureNum);
             fprintf(fid,'W =%9.4f   边数 = %2d  限差 = %8.4f (2*sqrt(n)*m0)', -W, NumInOneLoop,2*sqrt(NumInOneLoop*1)*m_Sigma);
             if  abs(W) > 2*sqrt(NumInOneLoop*1)*m_Sigma
                 fprintf(fid,'%s   ','超限');
             else
                 fprintf(fid,'%s','    ');
             end
             for j = 1:NumInOneLoop
                 fprintf(fid,'%s - ',unObsN{lsLoopClosure(j)});
             end
             fprintf(fid,'%s\n',unObsN{lsLoopClosure(j+1)});
             NumInOneLoop = 1;
        end
        if i == m_Lnumber
           waitbar(i/m_Lnumber,f,{'环闭合差检验完毕' 'Operation Successfully'});
        else
           waitbar(i/m_Lnumber,f);
        end
end 
success = 1;