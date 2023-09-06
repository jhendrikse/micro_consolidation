function [block,BlockAccuracy,blockAccuracy,blockNcorrect,meanCorrKeyPressPerSecNoOutlier,meankeypressperblock,MicroOnline,MicroOffline,TotalMicroOnline,TotalMicroOffline,TotalTaskLearning, deltaperformance, TotalLearning] = srtt_work_MRI(filename)

% This function can be used to output a 'learning score' from the serial
% reaction time task (SRTT)


% import data

%filename = '/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/P001_1_SRTT_version_1_2022_Mar_28_1445.csv'
%filename = '/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/S34_RY/S34_RY_r1_SRTT.csv'
%filename = '/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/S30_WD/S30_WD_r1_SRTT.csv'


data = import_SRTT_MRI(filename); % import individual data file using out import_data function (requires an input variable 'filename' which contains a string of full path to data file


%% calculate overall performance on the task

% First, we will calculate the total number of hits and errors made across the entire task

% Create logical variable to index all trials where a response was made (i.e. when a key was
% pressed, regardless of whether it was correct or not)

logResp = ((data.key_resp_1keys ~= 'None') & (~isundefined(data.key_resp_1keys))) ; % identify trials where values of either a,b,c,or d are contained in key_resp_1keys. This line will return an array of values with 1s for trials with response, 0s for trials without response
data.logResp = logResp ; % add logResp to original data table
trialResp = (data.key_resp_1keys(logResp == 1)) ; % Extract all trials where a response was made
totalTrialResp = length(trialResp) ; % determine the total number of trials where a response was made

% identify all 'correct' keypresses (when the correct key was the first key pressed).

% For trials that have a response in key_resp_1keys, determine correct responses from key_resp_2keys (i.e. if 'None' is listed in key_resp_2keys
... participant has responded correctly on trial and pressed correct button in first instance)
    
% for trials that have a response, look for a 'None' value in key_resp_2keys. This will create another logical variable with 1s for correct trials and 0s for incorrect trials
logCorr = (strcmp('None',string(data.key_resp_2keys(logResp == 1))) & (~isundefined(data.key_resp_2keys(logResp == 1)))) ; % index trials in key_resp_2keys where no second button was pressed (i.e. trials where correct button is pressed in first instance)

% determine the number of correct trials
TrialCorr = trialResp(logCorr == 1) ; % use the logical variable logCorr to index correct trials from all trial responses
TrialIncorr = trialResp(logCorr == 0) ; % use the logical variable logCorr to index incorrect trials from all trial responses
totalTrialIncorr = length(TrialIncorr) ; % find the length of this variable (gives us the total number of items) to determine total no. of incorrect responses

% calculate mean and SD reaction time for correct responses
trialRespRTs = data.key_resp_1rt(logResp) ; % extract all RT's for trials where a responses was made
trialRTcorr = trialRespRTs(logCorr == 1) ; % extract RT's corresponding to correct response
meanRTcorr = mean(trialRTcorr, 'omitnan') ; % calculate mean of RT's corresponding to correct response
stdevRTcorr = std(trialRTcorr, 'omitnan') ; % calculate sd of RT's corresponding to correct response

% calculate mean and SD of reaction times for incorrect responses
trialRTincorr = trialRespRTs(logCorr == 0) ; % extract RT's corresponding to correct response
meanRTincorr = mean(trialRTincorr) ; % calculate mean of RT's corresponding to incorrect response
stdevRTincorr = std(trialRTincorr) ; % calculate sd of RT's corresponding to incorrect response

logCorrAlldata = (strcmp('None',string(data.key_resp_2keys)) & (~isundefined(data.key_resp_2keys)) & (logResp == 1)) ; % adding a log index of correct responses to put into original data variable
... (this needs to be the same size as the original dataset i.e. contain ~2353 rows, so the code is a little different to other logCorr variable which was only contains trials w/ a response)
    
data.logCorrAlldata = logCorrAlldata ; % add a log index (1s,0s) of correct responses into original data table
data.blocksthisRepN_add1 = data.blocksthisRepN + 1 ; % in the original table 'data', the variable which indicates the current block number count upwards from 0. This can be annoying to index in matlab code.
... adding +1 to values in this variable to refer to the first block as 'block 1' instead of 'block 0'.
    
%% reduce dataset to key variables

% Our original datafile now contains two additional columns denoting trials
% with a response and trials with correct responses, but it also contains a
% lot of variables that aren't essential for our analysis.

% create another table variable 'primaryData' containing our main variables
% of interest to make life easier

primaryData = data(:,{'trial','blocksthisRepN','blocksthisRepN_add1','corrButton','corrAns','logResp','logCorrAlldata','key_resp_1keys','key_resp_1rt','key_resp_2keys','key_resp_2rt','key_resp_1started'}) ;


%% calculate general performance in each block

% calculate mean RTs and number of correct responses and errors per block

block = struct() ; % initialise the structure variable 'block' - helps the code run faster if you create the variable before adding contents to it in for loop

for i = 1:(max(primaryData.blocksthisRepN_add1)) % this is the start of a 'for loop' which we will use to extract data separately for each block in the task.
    % here, 'i' is used as a counter to indicate how many times we want the following section of code to continuously loop/repeat. You can read the line above as: run
    % this loop 'i' times, where i = the values of 1 through to the total
    % number of blocks in the task (i.e 24)
    
    % we will save this data in a 'structure' - a hierarchical way to
    % organise data whereby variables can have many fields
    
    block.RTallResponse{i} = primaryData.key_resp_1rt(primaryData.logResp & (primaryData.blocksthisRepN_add1 == (i))) ; % extract all response RT's in each block (regardless of correctness)
    block.RTcorrResponse{i} = primaryData.key_resp_1rt(primaryData.logCorrAlldata & (primaryData.blocksthisRepN_add1 == (i))) ; % extract all correct RT's in each block
    block.RTincorrResponse{i} = primaryData.key_resp_1rt(primaryData.logCorrAlldata == 0 & primaryData.logResp & (primaryData.blocksthisRepN_add1 == (i))) ; % extract all incorrect RT's in each block
    block.Nresponse{i} = length(block.RTallResponse{i}) ; % calculate number of responses made in each block (based on the number of items contained in block.RTallResponse)
    block.Ncorrect{i} = length(block.RTcorrResponse{i}) ; % calculate number of correct responses made in each block
    block.Nerror{i} = length(block.RTincorrResponse{i}) ; % calculate number of incorrect responses made in each block
    block.accuracy{i} = 1 - (block.Nerror{i} ./ block.Ncorrect{i}) ; % calculate accuracy according to Bonstrup paper (1 - number of erroneous relative to correct key presses in each block)
    
end

% extract values from cell and calculate mean and sd for RT's in each block - we will save these into separate arrays

blockMeanRTallResponse = cellfun(@mean,block.RTallResponse(:)) ; % array of mean RTs for all responses
blockSDallResponse = cellfun(@std,block.RTallResponse(:)) ; % SD's of RTs from all responses
blockMeanRTcorrResponse = cellfun(@mean,block.RTcorrResponse(:)) ; % array of mean RTs for correct responses
blockSDcorrResponse = cellfun(@std,block.RTcorrResponse(:)) ; % SD's of RTs from correct responses
blockMeanRTincorrResponse = cellfun(@mean,block.RTincorrResponse(:)) ; % array of mean RTs for incorrect responses
blockSDincorrResponse = cellfun(@std,block.RTincorrResponse(:)) ; % SD's of RTs from incorrect responses
blockNresponse = cell2mat(block.Nresponse(:)) ; % number of responses made in each block
blockNcorrect = cell2mat(block.Ncorrect(:)) ; % number of correct responses made in each block
blockNerror = cell2mat(block.Nerror(:)) ; % number of incorrect responses made in each block
blockAccuracy = cell2mat(block.accuracy(:)) ; % block accuracy, quantified according to Bonstrup paper
blockNsummary = [blockNcorrect,blockNerror,blockNresponse,blockAccuracy] ; % correct, incorrect, and total number of responses made in each block, and overall accuracy
BlockAccuracy = blockAccuracy.' ;

% Up to this point, calculations have been without screening for outliers.
% Use the mean and SDs values to identify outliers (i.e. RT's outside mean +/- 2.7 SDs).

% index all correct responses across task and determine outliers based on overall mean/sd
logRTcorrOutlier = ((trialRTcorr < (meanRTcorr - (2.7.* stdevRTcorr))) | (trialRTcorr > (meanRTcorr + (2.7.* stdevRTcorr)))) ; % detect reaction times for correct response trials that fall outslide M +/- 2.7 SD
allOutlier = trialRTcorr(logRTcorrOutlier) ; % create an array of outlier reaction times (for correct trials)
totalOutlier = length(allOutlier) ; % determine the total number of outliers (based on overall mean and sd)
meanRTcorrOutlierCorr = mean(trialRTcorr(~logRTcorrOutlier)) ; %

% create a log index of outlier trials to put back into our primaryData table
logRTcorrOutlierAlldata = (primaryData.logCorrAlldata) & (~isnan(primaryData.key_resp_1rt)) & ((primaryData.key_resp_1rt < (meanRTcorr - (2.7.* stdevRTcorr))) | (primaryData.key_resp_1rt > (meanRTcorr + (2.7.* stdevRTcorr)))) ;
primaryData.logRTCorrOutlier = logRTcorrOutlierAlldata ;

% create logical variable that denotes all correct responses that are not
% outliers
logCorrNoOutlier = primaryData.logCorrAlldata & ~logRTcorrOutlierAlldata ;
primaryData.logCorrNoOutlier = logCorrNoOutlier ; % place this log variable into primaryData table

% run another for loop to extract correct responses from each block,
% excluding outliers

for j = 1:(max(primaryData.blocksthisRepN_add1))
    
    % save this data in the same structure 'block'
    block.RTcorrResponseNoOutlier{j} = primaryData.key_resp_1rt(primaryData.logCorrNoOutlier & (primaryData.blocksthisRepN_add1 == (j))) ; % extract all correct RT's in each block, excluding outliers
    block.NcorrNoOutlier{j} = length(block.RTcorrResponseNoOutlier{j}) ; % calculate number of correct responses made in each block, excluding outliers
    
end

%calculate mean, SD, N of the data per block, outliers removed
blockMeanRTcorrNoOutlier = cellfun(@mean,block.RTcorrResponseNoOutlier(:)) ; % array of mean RTs for correct responses, transpose into column vector
blockSDcorrNoOutlier = cellfun(@std,block.RTcorrResponseNoOutlier(:)) ; % SD's of RTs from correct responses, excluding outliers, transpose into column vector
blockNcorrNoOutlier = cell2mat(block.NcorrNoOutlier(:)) ; % extract number of correct responses in each block, excluding outliers, transpose into column vector

%updated block summary array
blockNsummary = [blockNsummary,blockNcorrNoOutlier,blockMeanRTcorrResponse,blockSDcorrResponse,blockMeanRTcorrNoOutlier,blockSDcorrNoOutlier] ; % add this information to the block summary array (& mean RTs from each block without outlier removal)
tableBlockNsummary = table(blockNcorrect,blockNerror,blockAccuracy,blockNresponse,blockNcorrNoOutlier,blockMeanRTcorrResponse,blockSDcorrResponse,blockMeanRTcorrNoOutlier,blockSDcorrNoOutlier) ; % create a table with this information, which includes headers

%% calculate inter-tap intervals

% to this point, we have calculated basic performance measures for each block (e.g. number correct, errors, mean RT, etc.).
% Now, we want to calculate 1/RT (keypresses per second) as a more specific performance measure

%first, we will calculate 1/RT from mean RTs of each block
meanCorrKeypressPerSec = 1./ blockMeanRTcorrResponse ;
meanCorrKeyPressPerSecNoOutlier = 1./ blockMeanRTcorrNoOutlier ;

% next, calculate inter-tap interval using the timestamps in key_resp_1_started variable

% create a table that just contains data from valid keypresses (i.e removes the empty 'rest' sequences)
responseData = table(primaryData.trial(primaryData.logResp),primaryData.blocksthisRepN_add1(primaryData.logResp),primaryData.key_resp_1keys(primaryData.logResp),primaryData.key_resp_1started(primaryData.logResp),primaryData.key_resp_1rt(primaryData.logResp),primaryData.logCorrAlldata(primaryData.logResp),primaryData.logRTCorrOutlier(primaryData.logResp),primaryData.logCorrNoOutlier(primaryData.logResp)) ;
responseData.Properties.VariableNames = {'trial','blocks','key_resp_1keys','key_resp_1started','key_resp_1rt','logCorrAlldata','logRTCorrOutlier','logCorrNoOutlier'} ;

interTapInterval = responseData.key_resp_1started(2:end) - responseData.key_resp_1started(1:end-1) ; % subtract the timestamp of key_resp_1started from the previous trial key_resp_1started timestamp (excluding first entry which is not a true inter-tap interval)
responseData.interTapInterval = nan(length(responseData.trial),1) ; % initialise a new varaible filled with NaNs in responseData table
responseData.interTapInterval(2:end,1) = interTapInterval ; % enter interTapInterval data into table

% run another for loop to extract inter-tap intervals (1. for all responses, and 2. for correct responses from each block, excluding outliers)
for k = 1:(max(responseData.blocks))
    
    %save this data in the same structure 'block'
    block.interTapIntervalAllresponse{k} = responseData.interTapInterval((responseData.blocks == (k)) & (responseData.trial ~= 1) & ((~isnan(responseData.interTapInterval)))) ; % extract all inter-tap intervals where a response was made (excluding first entry in each block, which is not a true inter-tap interval)
    block.interTapIntervalCorr{k} = responseData.interTapInterval((responseData.blocks == (k)) & (responseData.logCorrNoOutlier) & (responseData.trial ~= 1) & ((~isnan(responseData.interTapInterval)))) ; % extract all correct inter-tap intervals in each block, excluding outliers
    
end

% next, we will calculate a running avearge of keypresses per second, over
% 5 presses (aka the Bonstrup paper) and over 4 presses (works well with
% our 12-item sequence)

% inside the following for loops, the variables will change size as more values are saved into them with every iteration of the loop. Matlab doesn't like this very much.
% % in the line below, we have initialised the cell variable blockInterTap with a pre-allocated
% % size (1 row, 24 columns) so matlab stays calm and carries on.

ave4 = 4; % average across 4 presses
ave5 = 5; % average across 5 presses

block.interTapAve4 = cell(1,24) ; % initialise variables
block.interTapAve5 = cell(1,24) ;

block.interTapSD4 = cell(1,24) ;
block.interTapSD5 = cell(1,24) ;

block.TapPerSecAve4 = cell(1,24) ;
block.TapPerSecAve5 = cell(1,24) ;

block.TapPerSecSD4 = cell(1,24) ;
block.TapPerSecSD5 = cell(1,24) ;


for z = 1:length(block.interTapIntervalCorr)
    
    InterTapSize = numel(block.interTapIntervalCorr{z}) ; % each block has a different number of values (i.e. responses). This number will not always divide evenly by 4 or 5.
    % 'numel' returns the number of elements in each block of inter-tap interval values
    
    % There is a fair bit in the following lines of code.
    % We are taking one block of inter-tap intervals and dividing it by ave4 (i.e. 4 presses). If it does not divide evenly, we will add the remainder in NaN values to the array to contain a divisible total number of values ...
    ... (i.e. if block contains 14 values, and we wish to divide by 4, we would have a remainder of 2, add 2 NaN values to array to bring array size to divisible size of 16).
        ... We are then using function 'reshape' to restructure our single column array of inter-tap intervals into an array with 'ave4' number of columns (i.e. 4) and automatically determined number of rows '[]'.
        ... This array configuration allows us to then use nanmean to calculate the average of each column of values in the array (i.e. average across 4 inter-tap intervals). These averages are saved back into the structure 'block'.
        
    if mod(InterTapSize,ave4) == 0 % if inter-tap interval block divides by 4, proceed with following line.
        block.interTapAve4{z} = mean(reshape(block.interTapIntervalCorr{z},ave4,[])) ; % use function 'reshape' to restructure our single column array of inter-tap intervals into an array with 'ave4' number of columns (i.e. 4) and automatically determined number of rows '[]'.
        ... This array configuration allows us to then calculate the average of each column of values in the array (i.e. average across 4 inter-tap intervals). These averages are saved back into the structure 'block'.
            block.interTapSD4{z} = std(reshape(block.interTapIntervalCorr{z},ave4,[])) ;  % standard deviation
        
    elseif mod(InterTapSize,ave4) ~= 0 % else if inter-tap interval block does not divide by 4, do the following..
        %Here, we run the same calculation as above, except we remove the remainder (i.e. incomplete sequence)
        block.interTapAve4{z} = mean(reshape(block.interTapIntervalCorr{z}(1:end-(mod(InterTapSize,ave4))),ave4,[])) ;
        block.interTapSD4{z} = std(reshape(block.interTapIntervalCorr{z}(1:end-(mod(InterTapSize,ave4))),ave4,[])) ;
        
    end
    
    %% add in moving average
    moveavelength = 4; %set the moving average window length. This should be either 4 or 5.. (any other number may have implications for the code below)
    moveaveTruncate = floor(moveavelength / 2); %identify how many need to be removed (round down to identify how many trials need to be removed from inter-tap intervals)
    
    %for loop...
    for z = 1:length(block.interTapIntervalCorr) %loop will run over each block
        
        block.movemeanintertap{z} = movmean(block.interTapIntervalCorr{1,z},moveavelength); % calculates the moving mean (as previously set by length) of the intertap intervals of correct trials
        block.movemeanintertap{z} = block.movemeanintertap{z}(moveaveTruncate+1:length(block.interTapIntervalCorr{1,z})-moveaveTruncate); %removes the necessary number of values from start and end of the column that are not averaged over the previosuly set number
        block.movemeankeypresss{z} = 1 ./ block.movemeanintertap{z}; % converts intertap intervals into key presses per second
        
        
        block.movemeanrests{z} = nan(10,1); %creates 1 column with 10 rows of NaNs
        %for plotting - need to insert gaps to seperate blocks and rests
        
    end
    %%
 if filename == '/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/S30_WD/S30_WD_r1_SRTT.csv'
     block.interTapIntervalCorr(1,1)={nan(8,1)} %only for S30_WD, remove for others 
 end 

    %%

    block.meankeypressperblock = cellfun(@mean,block.movemeankeypresss(:)); %mean of the moving means for each block
    
    meankeypressperblock = block.meankeypressperblock' ;%save variable outside of block
    
    allmovemeanintertap = vertcat(block.movemeanintertap{1,:}); % transforms intertap intervals into a single column
    allmovemeankeypresss = vertcat(block.movemeankeypresss{1,:}); % transforms key presses per second into a single column
    %%
    
    %Bonstrup 'the tapping speed of incomplete sequences was averaged with
    %the previous complete sequence' - at present script removes incomplete
    %4/5 item sequences.
    
    % repeat this process as above, but averaging across 5 inter-tap intervals to align with Bonstrup et al. 2019
    
    %     if mod(InterTapSize,ave5) == 0
    %
    %         block.interTapAve5{z} = mean(reshape(block.interTapIntervalCorr{z},ave5,[])) ; % averages
    %         block.interTapSD5{z} = std(reshape(block.interTapIntervalCorr{z},ave5,[])) ;  % standard deviation
    %
    %     elseif mod(InterTapSize,ave5) ~= 0
    %
    %         block.interTapAve5{z} = mean(reshape(block.interTapIntervalCorr{z}(1:end-(mod(InterTapSize,ave5))),ave5,[])) ; % averages
    %         block.interTapSD5{z} = std(reshape(block.interTapIntervalCorr{z}(1:end-(mod(InterTapSize,ave5))),ave5,[])) ; % standard deviation
    %
    %     end
    %
    %     % calculate 1/inter-tap intervals (presses per second)
    %
    %     block.TapPerSecAve4{z} = 1./ block.interTapAve4{z} ; %Tap per second (5 sec running average)
    %     block.TapPerSecAve5{z} = 1./ block.interTapAve5{z} ; %Tap per second (4 sec running average)
    %
    %     block.TapPerSecSD4{z} = 1./ block.interTapSD4{z} ; %Tap per second SD (5 sec running average)
    %     block.TapPerSecSD5{z} = 1./ block.interTapSD5{z} ; %Tap per second SD (4 sec running average)
    
end



%% Micro online & offline learning

% still need to consider how/if we model total learning to identify plateau
% (e.g. 95% of total learning akin to Bonstrup paper)

% micro online learning - difference in tapping speed between the first and
% last correct sequence of a block

for g = 1:length(block.movemeankeypresss) %loop over each block
    %take the last movemeankeypresss value in the block and subtract the
    %first value from the same block - save as new variable (MicroOnline)
    if isempty(block.movemeankeypresss{1,g}) %if the movemeankeypresss column is empty (i.e participant did not complete enough trials in the block for moving mean to to generate a number - possibly error/interruption in the task)
        block.MicroOnline{g} = nan; %fill empty cell with nan
    else
        block.MicroOnline{g} = block.movemeankeypresss{1,g}(end,1) - block.movemeankeypresss{1,g}(1,1);
    end
end
MicroOnline = cell2mat(block.MicroOnline);%save variable outside of block

% sum all values in MicroOnline for TotalMicroOnline
block.TotalMicroOnline = sum((cell2mat(block.MicroOnline)), 'omitnan');
TotalMicroOnline = block.TotalMicroOnline; %save variable outside of block

% micro offline learning - difference in tapping speed of the last correct
% sequence of a block and first correct sequence of the next block
for h = 1:length(block.movemeankeypresss)-1%loop over each block
    %take the first movemeankeypresss value in block n+1 and subtract the
    %last value from the block n - save as new variable (MicroOffline)
    if isempty(block.movemeankeypresss{1,h}) | isempty(block.movemeankeypresss{1,h+1})%if the movemeankeypresss column is empty (i.e participant did not complete enough trials in the block for moving mean to to generate a number - possibly error/interruption in the task)
        block.MicroOffline{h} = nan; %fill empty cell with nan
    else
        block.MicroOffline{h} = block.movemeankeypresss{1,h+1}(1,1) - block.movemeankeypresss{1,h}(end,1);
    end
end
MicroOffline = cell2mat(block.MicroOffline);%save variable outside of block

% sum all values in MicroOffline for TotalMicroOffline
block.TotalMicroOffline = sum((cell2mat(block.MicroOffline)), 'omitnan');
TotalMicroOffline = block.TotalMicroOffline; %save variable outside of block

% total task learning - very last 4 item chunk - very first 4 item chunk
block.TotalTaskLearning = block.movemeankeypresss{1,end}(end,1) - block.movemeankeypresss{1,1}(1,1);
TotalTaskLearning = block.TotalTaskLearning;%save variable outside of block

% Total learning (used for total early learning) - "sum of single-trial
% performance changes" (bonstrup)
for w = 1:length(meankeypressperblock)-1
    deltaperformance{w} = meankeypressperblock(1,w+1)- meankeypressperblock(1,w);% calculates change in mean key presses per second for each block
end 
deltaperformance = cell2mat(deltaperformance);
TotalLearning = sum(deltaperformance, 'omitnan');


% let's plot some outputs
%
% .. similar figures could be adapted and used in your theses (i.e. averaged
% across subjects with error bars)

%% Uncomment section below to generate plots

% figure('color','w') ; % opens a plain figure with white background
% 
% subplot(2,2,1) % subplots enable you to plot multiple axes on same figure, read this as 'create a plot with 2 rows and 2 columns, and assign the following axes into position 1'
% plot(meanCorrKeypressPerSec,'*-k') % plot keypress per second, without outlier removal
% hold on ; % 'hold on' keeps the current plot open to allow you to run more code without overriding (e.g. add axis labels, title, etc.)
% title('Correct keypresses per second, per block') ;
% xlabel('Block') ;
% ylabel('Correct keypresses per second') ;
% 
% subplot(2,2,2) % plot second set of axes in position 2 of the 2 x 2 grid
% plot(meanCorrKeyPressPerSecNoOutlier,'*-k') % plot keypress per second, with outlier removal
% hold on ;
% title('Correct keypresses per second, per block, outliers removed') ;
% xlabel('Block') ;
% ylabel('Correct keypresses per second') ;
% 
% subplot(2,2,3) % plot third set of axes in position 3 of the 2 x 2 grid
% plot(blockNcorrNoOutlier,'*-k') % plot number of correct responses per block
% hold on ;
% title('Number of correct responses per block (s), outliers removed') ;
% xlabel('Block') ;
% ylabel('Number of correct responses') ;
% 
% subplot(2,2,4) % plot fourth set of axes in position 4 of the 2 x 2 grid
% plot(blockNresponse,'*-k') % plot number of responses per block
% hold on ;
% title('Number of responses per block (s)') ;
% xlabel('Block') ;
% ylabel('Number of responses') ;
% 
% figure ; % plot errors in separate figure (can append to subplot if needed)
% plot(blockNerror,'*-k') % plot number of errors per block
% hold on ;
% title('Number of errors per block (s)') ;
% xlabel('Block') ;
% ylabel('Number of errors') ;
% 
% figure ; % plot errors in separate figure (can append to subplot if needed)
% plot(blockAccuracy,'*-k') % plot accuracy per block
% hold on ;
% title('Accuracy per block (s)') ;
% xlabel('Block') ;
% ylabel('Accuracy') ;

%end - uncomment with finalised function