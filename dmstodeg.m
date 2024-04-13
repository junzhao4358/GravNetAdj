function  deg = dmstodeg(dms)
dms = dms/10000;
DD = floor(dms);
MM = floor((dms - floor(dms))*100);
SS = (dms - DD - MM/100)*10000;
deg = dms2degrees([DD MM SS]);
