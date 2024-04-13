function str_road = writetoocean(file_raod,B,L,H,Y,M,D,HH)
      str_road = strcat(file_raod,'haichao.bat');
      fid =  fopen(str_road,'w');
      fprintf(fid,'nloadf qing %5f %5f  %7f  k1.fes.2004.asc gr.gbaver.wef.p01.ce l >   out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  k2.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  m2.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  m4.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  mf.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  mm.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  n2.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  o1.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  p1.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  q1.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  s1.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n',B,L,H);
      fprintf(fid,'nloadf qing %5f %5f  %7f  s2.fes.2004.asc gr.gbaver.wef.p01.ce l >>  out_fes2004.txt\n\n',B,L,H);
      fprintf(fid,'type out_fes2004.txt | harprp g>out1.txt\n');
      fprintf(fid,'type out1.txt | hartid %5d %2d %2d  %10f %2d  %2d %2d %2d  > haichao.txt',Y,M,D,HH,0,0,1,1);
      fclose(fid)
