function [lc] = extractMAUDPhases(lc)
%extractMAUDPhases pulls all the phase information from MAUD
    phase=cell(lc.ncases{lc.BS} ,1);    
    if and(lc.isMatlab,~lc.options{6,2})
        parfor i=1:lc.ncases{lc.BS} 
            if isfile(lc.OutputPar.FullPath{i,lc.BS})
                phase{i}=ExtractPhaseFromPar(lc.OutputPar.FullPath{i,lc.BS});
            end
        end
    else %Octave
        for i=1:lc.ncases{lc.BS} 
            if isfile(lc.OutputPar.FullPath{i,lc.BS})
                phase{i}=ExtractPhaseFromPar(lc.OutputPar.FullPath{i,lc.BS});
            end
        end        
    end
    lc.phase{lc.BS}=phase;
end
