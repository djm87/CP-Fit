function [lc] = DistributeRun2CPU(lc)
%DistributeRun2CPU This distributes cases accross the cores
    lc.runtime=cell(1,2,lc.numBatchStep);
    
    %Set the concurrent environement for calling MAUD                   
    lc=SetCPUs(lc); 
    
    %Make a 1D struct for each core and fill each with a list from CaseNums
    while true
        for i=1:ceil(lc.ncases{lc.BS}/lc.options{8,2})
            for j=1:lc.options{8,2}
                cnt=lc.options{8,2}*(i-1)+j;
                lc.caseID{j}(i)=cnt; %The indices of the appended folders
                if cnt==lc.ncases{lc.BS}; break; end
            end 
            if cnt==lc.ncases{lc.BS}; break; end
        end
        if cnt==lc.ncases{lc.BS}; break; end
    end
end
