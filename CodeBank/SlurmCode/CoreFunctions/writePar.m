function [] = writePar(c,parNameOut,varargin)
%[] = writePar(c,parNameOut,varargin) 
%   Description: readPar writes a MAUD parameter (.par) file from a cell 
%   structure wherein each cell is string containing a line of the 
%   parameter file
%
%   Input:  (1) a cell list of strings with each line of the parameter file,
%           (2)name of local or full path to par file
%   Option: (1) set flag = 'keep intensity' to keep intensity data on write in 
%   Note: Removing intensity data will improve write time
    if ( length(varargin) == 0) %default if no option specified
        flag='remove intensity';
    elseif strcmp( varargin{1},'keep intensity') %Keep intensity if specified
        flag=varargin{1};
    elseif strcmp( varargin{1},'no check') %Keep intensity if specified
        flag=varargin{1};
    else %If something specified but not handled tell user default option is being used
        disp('Undefined option, setting default "remove intensity"')
        flag='remove intensity';
    end
        
    %no check addded because octave contains is very slow.
    if ~strcmp(flag, 'no check')
      %See if intensity data exists
      StartInt=contains2(c,'#custom_object_intensity_data');
      EndInt=contains2(c,'#end_custom_object_intensity_data');

      StartIntPos=find(StartInt==true);

    
      if isempty(StartIntPos)
%         disp('No intensity data found in MAUD parameter file')
      else
        if strcmp(flag,'keep intensity')
            disp('Keeping intensity data')
        else
            %Remove the intensity data from c
            EndIntPos=find(EndInt==true);

            for i =1:numel(StartIntPos)
                if i==1
                    toRemove=StartIntPos(i):EndIntPos(i);
                else
                    toRemove=horzcat(toRemove,StartIntPos(i):EndIntPos(i));
                end
            end
            c(toRemove)=[];       
        end
      end
    end

    %Write c to parameter file
    fid=fopen(parNameOut,'w');
    for i=1:length(c)
      fprintf(fid,'%s\n',char(c(i)));
    end
    fclose(fid);

end
