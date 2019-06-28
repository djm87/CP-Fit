function [phase] = ExtractPhaseFromPar(parName)

% import E-WIMV odf from MAUD par file (tested for MAUD version 2.8)
%
% Syntax
%   odf = loadODF_MAUD(fname)
%
% Input
%  fname - filename
%
% Output
%  phase - a struct of phases containing: symmetries, resolution, odf, etc..
%
    c = readPar(parName,'');
    keywords={'_pd_phase_name','#subordinateObject_E-WIMV','_symmetry_space_group_name_H-M',...
        '_cell_length_a','_cell_length_b','_cell_length_c',...
        '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',...
        '_rita_generate_symmetry','_rita_wimv_odf_resolution',...
        '_rita_wimv_odf_values'};

    % Find the number of phases in par file 
    loctmp=find(~cellfun(@isempty,strfind(c,keywords{1})));
    phase.num=numel(loctmp);

    %Get information for each phase
    phase.c=cell(phase.num,1);
    for i=1:phase.num
        %Get the phase names and the beginning and end of the phase objects
        text=strrep(c{loctmp(i)},[keywords{1},' '],''); 
        phase.name{i}=text(text~='''');
        tmp=find(~cellfun(@isempty,...
        strfind(c,['#subordinateObject_',phase.name{i}])));
        phase.begin(i)=tmp(1);
        tmp=find(~cellfun(@isempty,...
        strfind(c,['#end_subordinateObject_',phase.name{i}])));
        phase.end(i)=tmp(end);

        %Create subsets of c for each phase
        phase.c{i}={c{phase.begin(i):phase.end(i)}};

        %find if phase.c contains an E-WIMV odf 
        phase.hasODF(i)=~isempty(find(~cellfun(@isempty,...
            strfind(phase.c{i},keywords{2}))));  
    end
    
    %loop over each phase and extract cell and odf parameters
    phase.axes=zeros(phase.num,3);
    phase.angle=zeros(phase.num,3);
    phase.odfBounds=zeros(phase.num,3);
    for i=1:phase.num
        
            %Get the space group for the phase
            loctmp=find(~cellfun(@isempty,strfind(phase.c{i},keywords{3})));
            text=strrep(phase.c{i}{loctmp},[keywords{3},' '],'');
            
            if and(i==2,any(text==':')) 
                disp('Warning: MAUD cell origin option not handled, removing :# from space group')
                text=text(1:end-2); %Removes the differientiation of cell origins
            end
            phase.sym{i}=text;

            %Get the Cell lattice length for the phase
            for j=1:3
                loctmp=find(~cellfun(@isempty,strfind(phase.c{i},keywords{j+3})));
                text=regexp(phase.c{i}{loctmp}, ' ', 'split');
                text=text{2};
                if any(text=='(')
                    phase.axes(i,j)=str2num(text(1:find(text=='(')-1));
                else
                    phase.axes(i,j)=str2num(text);
                end
            end

            %Get the Cell lattice angles for the phase
            for j=1:3
                loctmp=find(~cellfun(@isempty,strfind(phase.c{i},keywords{j+6})));
                text=regexp(phase.c{i}{loctmp}, ' ', 'split');
                text=text{2};
                if any(text=='(')
                    phase.angle(i,j)=str2num(text(1:find(text=='(')-1));
                else
                    phase.angle(i,j)=str2num(text);
                end
            end

            %Get crystal symmetry 
%             phase.cs{i}=crystalSymmetry(phase.sym{i},phase.axes(i,:),...
%                 phase.angle(i,:)*degree,'mineral',phase.name{i},'X||a','Z||c');
        if phase.hasODF(i)==true
            %Get odf sample symmetry
            loctmp=find(~cellfun(@isempty,strfind(phase.c{i},keywords{10})));
            text=regexp(phase.c{i}{loctmp}, ' ', 'split');
            if strcmp(text{2},'none')
                phase.ss{i}='-1';
            elseif strcmp(text{2},'orthorhombic')
                phase.ss{i}='222';
            else
               error('Unhandled specimen symmetry used in ODF');
            end

            %Get odf resolution
            loctmp=find(~cellfun(@isempty,strfind(phase.c{i},keywords{11})));
            text=regexp(phase.c{i}{loctmp}, ' ', 'split');
            phase.res(i)=str2num(text{2});


            %Get odf bounds
%             [a,b,g] = fundamentalRegionEuler(phase.cs{i},phase.ss{i});
%             phase.odfBounds(i,:)=[a,b,g]./degree;
            
            %Find odf start location
            phase.odfStartLoc(i)=find(~cellfun(@isempty,...
                strfind(phase.c{i},keywords{12})));
            
            phase=GetODFBounds(phase,i);
            
        end
    end
    %% Extract ODF
    for w=1:phase.num
        if phase.hasODF(w)==true
            res=phase.res(w);

            alphaBlk=int16(phase.odfBounds(w,1)/res+1); %across a section    
            betaBlk=int16(phase.odfBounds(w,2)/res+1); %down a section
            gammaBlk=int16(phase.odfBounds(w,3)/res+1); %Number of sections
            odf=zeros(gammaBlk,betaBlk,alphaBlk);
            odflin=zeros(gammaBlk*betaBlk*alphaBlk,4);

            cnt=phase.odfStartLoc(w);
            cnt2=0;
            for i=1:double(gammaBlk)
                for j=1:double(betaBlk)
                    cnt=cnt+1;
                    %Since ZYZ convention 10,0,0 must be the same as 0,0,10,
                    %20,0,0 and 0,0,20 must be the same etc..
                    if j==1 && i>1
                        odf(i,j,1:alphaBlk-1)=circshift(odf(1,1,1:alphaBlk-1),-i+1);
    %                     odf(i,j,1:alphaBlk-1)=NaN;
                    else
                        odf(i,j,:)=str2num(phase.c{w}{cnt});
                    end
                    %Set ODF so alpha 0 is the same as 360
                    odf(i,j,end)=odf(i,j,1);
    %                 odf(i,j,end)=NaN;

                    for k=1:double(alphaBlk)
                        cnt2=cnt2+1;
                        odflin(cnt2,:)=[(k-1)*res,(j-1)*res,(i-1)*res,odf(i,j,k)];
                    end
                end
                cnt=cnt+1; %To skip the blank lines
            end
            phase.odf{w}=odflin;
        end
    end
end
function phase = GetODFBounds(phase,i)
%GetODFBounds counts the ODF format to extract the odf bounds
    cnt=phase.odfStartLoc(i);
    recorda=true;
    recordb=true;
    cntb=0; 
    cntg=0; 
    while ~strcmp(phase.c{i}{cnt+2},'#end_custom_object_odf')
        cnt=cnt+1;
        cntb=0;
        while ~isempty(phase.c{i}{cnt})
            if recorda
                cnta=length(str2num(phase.c{i}{cnt}));
                phase.odfBounds(i,1)=(cnta-1)*phase.res(i);%*(pi/180);
            end
            cnt=cnt+1;
            cntb=cntb+1;
        end
        if recordb
            phase.odfBounds(i,2)=(cntb-1)*phase.res(i);%*(pi/180);
        end
        cntg=cntg+1;
    end
    phase.odfBounds(i,3)=(cntg-1)*phase.res(i);%*(pi/180);
end