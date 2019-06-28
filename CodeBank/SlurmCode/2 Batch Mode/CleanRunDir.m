function [] = CleanRunDir
%CleanRunDir Removes .ins and .bat with MAUD_Batch_input as part of name
delete('Maud_Batch_input_*.*')
end

