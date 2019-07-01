function [] = printMessage(message)
  fprintf('%s\n',message)
  if and(exist('OCTAVE_VERSION', 'builtin') ~= 0,isunix)
    fflush(stdout);
  end
end