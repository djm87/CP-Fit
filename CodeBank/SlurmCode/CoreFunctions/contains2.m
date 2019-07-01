function [output] = contains2(strings2Search,string2find)
%output is a logical arrayfun
  %find where string occurs
%   if ~iscell(string2find) 
%     string2find={string2find};
%   end
  if ~iscell(strings2Search) 
    strings2Search={strings2Search};
  end  
  
  loc=strfind(strings2Search,string2find);
  output=~cellfun(@isempty,loc);
end