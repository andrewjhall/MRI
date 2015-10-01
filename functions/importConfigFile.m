function [id,date1,folderName,goal,contents1,pathNiftis,pathMasks,nScans,scan001,scan002,scan003,scan004,scan005,scan006,scan007,scan008,scan010,scan011,scan012,scan013,scan014,scan015,scan016,scan017] = importConfigFile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [ID,DATE1,FOLDERNAME,GOAL,CONTENTS1,PATHNIFTIS,PATHMASKS,NSCANS,SCAN001,SCAN002,SCAN003,SCAN004,SCAN005,SCAN006,SCAN007,SCAN008,SCAN010,SCAN011,SCAN012,SCAN013,SCAN014,SCAN015,SCAN016,SCAN017]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [ID,DATE1,FOLDERNAME,GOAL,CONTENTS1,PATHNIFTIS,PATHMASKS,NSCANS,SCAN001,SCAN002,SCAN003,SCAN004,SCAN005,SCAN006,SCAN007,SCAN008,SCAN010,SCAN011,SCAN012,SCAN013,SCAN014,SCAN015,SCAN016,SCAN017]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [id,date1,folderName,goal,contents1,pathNiftis,pathMasks,nScans,scan001,scan002,scan003,scan004,scan005,scan006,scan007,scan008,scan010,scan011,scan012,scan013,scan014,scan015,scan016,scan017]
%   = importfile('config.csv',2, 3);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/09/29 15:42:50

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,8]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [1,2,8]);
rawCellColumns = raw(:, [3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
id = cell2mat(rawNumericColumns(:, 1));
date1 = cell2mat(rawNumericColumns(:, 2));
folderName = rawCellColumns(:, 1);
goal = rawCellColumns(:, 2);
contents1 = rawCellColumns(:, 3);
pathNiftis = rawCellColumns(:, 4);
pathMasks = rawCellColumns(:, 5);
nScans = cell2mat(rawNumericColumns(:, 3));
scan001 = rawCellColumns(:, 6);
scan002 = rawCellColumns(:, 7);
scan003 = rawCellColumns(:, 8);
scan004 = rawCellColumns(:, 9);
scan005 = rawCellColumns(:, 10);
scan006 = rawCellColumns(:, 11);
scan007 = rawCellColumns(:, 12);
scan008 = rawCellColumns(:, 13);
scan010 = rawCellColumns(:, 14);
scan011 = rawCellColumns(:, 15);
scan012 = rawCellColumns(:, 16);
scan013 = rawCellColumns(:, 17);
scan014 = rawCellColumns(:, 18);
scan015 = rawCellColumns(:, 19);
scan016 = rawCellColumns(:, 20);
scan017 = rawCellColumns(:, 21);

