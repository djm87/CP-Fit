function [BP] = checkBP(BP) 
  if isunix
    if strcmp(BP(1:3),'///')
      BP=strrep(BP,'///','//');
    end
    BP=strrep(BP,'//','/');  
    BP=['/' BP];
  elseif ispc
    if strcmp(BP(1:3),'\\\')
      BP=strrep(BP,'\\\','\\');
    end
    BP=strrep(BP,'\\','\');  
    BP=['\' BP];
  end
end