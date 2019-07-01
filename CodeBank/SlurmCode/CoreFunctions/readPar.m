function [c] = readPar(parName,option)
%[] = readPar(parName,varargin) 
%   Description:readPar reads in a MAUD parameter (.par) file into a cell
%   structure wherein each cell is string containing a line of the
%   parameter file
%
%   Input:  (1) name if local or full path to par file
%   Output: (1) a cell list of strings with each line of the parameter file.
    if ispc
      loc=strfind(parName,'\');
    elseif isunix
      loc=strfind(parName,'/');
    end
    tmpPar = [tempname,'.par'];
    [~,~]=system(sprintf('sed -n "/#custom_object_intensity_data/,/#end_custom_object_intensity_data/!p" <%s >%s',parName,tmpPar));
    c=textread(tmpPar,'%s','delimiter','\n');
    if and(~isempty(c),strcmp(option,'write back'))
      copyfile(tmpPar,parName);
      printMessage('Removing intensity data at import and rewriting par without intensity data for future processing');
    else
      printMessage('Removing intensity data at import');
    end


end
