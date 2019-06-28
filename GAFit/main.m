function main()

caseFname = 'Inp_CaseList.in';
infoFname = 'Inp_Info.in';
GAFname = 'Inp_GACtl.in';

cases = readCaseFile(caseFname);

info = readInfoFile(infoFname);

GAopt = readGAControlFile(GAFname);

optimizationStart();

end