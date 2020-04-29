function [euler,weights] = loadVPSC(fname)
%loadVPSC loads the vpsc file and tries to delimiter options that are
%messing up the normal imputs
    try
        tmp = importdata(fname, ' ', 4);
        tmp = tmp.data;
    catch
        tmp = importdata(fname, '\t', 4);
        tmp = tmp.data;
    end
    weights=tmp(:,4);
    euler=tmp(:,1:3);
end

