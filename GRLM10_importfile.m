function [TPJO, GRLM10_Altimetry] = GRLM10_importfile3(filename)%,startRow, endRow)
 %function [TPJO,VarName1,VarName2,VarName3] = importfile(filename, startRow, endRow)

%% Initialize variables.
% if nargin<=2
% %     startRow = 44;
startRow = 1;
    endRow = inf;
% end

formatSpec = '%5s%5s%9s%4s%3s%7s%7s%7s%[^\n\r]';
%formatSpec = '%5s%4s%9s%4s%3s%7s%7s%7s%4s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Remove white space around all cell columns.
dataArray{1} = strtrim(dataArray{1});
dataArray{2} = strtrim(dataArray{2});
dataArray{3} = strtrim(dataArray{3});
dataArray{4} = strtrim(dataArray{4});
dataArray{5} = strtrim(dataArray{5});
dataArray{6} = strtrim(dataArray{6});
dataArray{7} = strtrim(dataArray{7});
dataArray{8} = strtrim(dataArray{8});

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
TPJO = dataArray{:, 1};
VarName1 = dataArray{:, 3};
VarName2 = dataArray{:, 6};
VarName3 = dataArray{:, 7};

VarName= [VarName1 VarName2 VarName3];  

VarName=cellfun(@str2double, VarName);  
I=find((VarName(:,1)<99999990)&(VarName(:,1)>19800000));  

GRLM10_Altimetry=VarName(I,:); 
TPJO=TPJO(I); 

% GRLM10_Altimetry = [VarName1 VarName2 VarName3]; 



