function recipe = readRecipeFile(recipeFname,totalParams)

file = fileread(recipeFname);
fileText = regexp(file, '\r\n|\r|\n', 'split')';

nrecipe = str2double(fileText{2});

if (fileText{3+nrecipe} == ' ')
    recipe = 1:totalParams;
else
    recipe = str2num(fileText{3+nrecipe});
end

end