function [S,diff,neighbor] = FindShortPath(p,exclude,m_Lnumber,m_Pnumber,StartP,EndP,AA,PP)
%m_Lnumber:表示高差总数;
%m_Pnumber:总点数;
%global m_Lnumber m_Pnumber StartP EndP AA PP
S = zeros(m_Pnumber,1);
diff = zeros(m_Pnumber,1);
neighbor = zeros(m_Pnumber,1);
for i = 1:m_Pnumber
    neighbor(i) = -1;
    S(i) = 1e30;
end
S(p) = 0;
diff(p) = 0;
neighbor(p) = p;
for i = 1:1000000
    unchanged = true;
    for j = 1:m_Lnumber
        if j == exclude
            continue;
        end
        p1 = StartP(j);
        p2 = EndP(j);
        S12 = 1/PP(j);
        if neighbor(p1)<0 && neighbor(p2)<0
            continue;
        end
        if S(p2) > S(p1) + S12
            neighbor(p2) = p1;
            S(p2) = S(p1) + S12;
            diff(p2) = diff(p1) + AA{j,4};
            unchanged  = false;
        elseif S(p1) > S(p2) + S12
            neighbor(p1) = p2;
            S(p1) = S(p2) + S12;
            diff(p1) = diff(p2) - AA{j,4};
            unchanged = false;
        end
    end
    if unchanged
        break;
    end
end

            
        
        