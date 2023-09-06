function data = import_SRTT_MRI(filename, dataLines)
%IMPORTFILE Import data from a text file
%  P0011SRTTVERSION12022MAR281445 = IMPORTFILE(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the data as a
%  table.
%
%  P0011SRTTVERSION12022MAR281445 = IMPORTFILE(FILE, DATALINES) reads
%  data for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  P0011SRTTversion12022Mar281445 = importfile("/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/P001_1_SRTT_version_1_2022_Mar_28_1445.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 30-Mar-2022 12:56:51

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 54, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["colourZ", "colourX", "colourC", "colourV", "corrButton", "corrAns", "trial", "blocksthisRepN", "blocksthisTrialN", "blocksthisN", "blocksthisIndex", "blockoftrialsthisRepN", "blockoftrialsthisTrialN", "blockoftrialsthisN", "blockoftrialsthisIndex", "keyPressLoopthisRepN", "keyPressLoopthisTrialN", "keyPressLoopthisN", "keyPressLoopthisIndex", "textBeginstarted", "textBeginstopped", "wait_for_triggerkeys", "wait_for_triggerrt", "wait_for_triggerstarted", "wait_for_triggerstopped", "stimulus_start_time", "circleZstarted", "circleZstopped", "circleXstarted", "circleXstopped", "circleCstarted", "circleCstopped", "circleVstarted", "circleVstopped", "key_resp_1keys", "key_resp_1corr", "key_resp_1rt", "key_resp_1started", "key_resp_1stopped", "key_resp_2keys", "key_resp_2started", "key_resp_2stopped", "reststarted", "reststopped", "key_resp_2rt", "endTextstarted", "endTextstopped", "participant", "session", "date", "expName", "psychopyVersion", "frameRate", "VarName54"];
opts.VariableTypes = ["categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "string", "double", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "categorical", "double", "double", "double", "categorical", "categorical", "double", "categorical", "double", "string", "double", "string", "string", "categorical", "double", "categorical", "double", "categorical", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["textBeginstopped", "wait_for_triggerstopped", "reststopped", "endTextstarted", "endTextstopped", "VarName54"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["colourZ", "colourX", "colourC", "colourV", "textBeginstopped", "wait_for_triggerstopped", "circleZstopped", "circleXstopped", "circleCstopped", "circleVstopped", "key_resp_1keys", "key_resp_1stopped", "key_resp_2keys", "key_resp_2stopped", "reststopped", "endTextstarted", "endTextstopped", "participant", "date", "psychopyVersion", "VarName54"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["corrButton", "corrAns", "expName"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["corrButton", "corrAns", "expName"], "ThousandsSeparator", ",");

% Import the data
data = readtable(filename, opts);

end