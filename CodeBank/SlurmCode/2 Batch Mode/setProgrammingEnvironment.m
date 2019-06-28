function [lc] = setProgrammingEnvironment(lc)
%getProgrammingEnvironment determines if using matlab or octave and it
%determines if you have the parallel toolbox installed
    lc.isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if lc.isOctave
        lc.isMatlab=false;
    else
        lc.isMatlab=true;
    end
    
    %Turn off the parallel compute toolbox if specified
    if and(lc.options{5,2},lc.isMatlab)
       ps = parallel.Settings;
       ps.Pool.AutoCreate = true; 
       p = gcp('nocreate');
       if isempty(p)
           parpool;
       end
    elseif lc.isMatlab
       ps = parallel.Settings;
       ps.Pool.AutoCreate = false; 
       poolobj = gcp('nocreate');
       delete(poolobj);
    end
end

