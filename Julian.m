function R = Julian(y,m,d)
%ÈåÂÔÊı
yy = y -1900;
mm = m -1;
dd = d;
w =  floor(yy/4);
if (y == 4*w & mm<2)
    dd = dd -1;
end
if mm == 0
    d1 = yy*365 + w - 0.5 + mm +dd;
else
    if (mm == 1)
        mm = 31;
        d1 = yy*365 + w - 0.5 + mm + dd;
    else
        mm = floor(mm*365/12)-floor(10/(4+mm));
        d1 = yy *365 + w -0.5 +mm + dd;
    end
end
R = d1 +2415020.0;
end