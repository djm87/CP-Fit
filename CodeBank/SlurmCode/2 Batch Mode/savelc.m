function [] = savelc(lc)
%savelc Saves lc with a time stamp
if lc.isMatlab  
  t=char(datetime('now'));
  t=strrep(t,' ', '_');
  t=strrep(t,'-', '_');
  t=strrep(t,':', '_');
  save(['lc_' t '.mat'],'lc');
else
  t= strftime ("%d_%b_%y_%T", localtime (time ()));
  t=strrep(t,':', '_');
  save(['lc_' t '.mat'],'lc','-mat7-binary');
end

