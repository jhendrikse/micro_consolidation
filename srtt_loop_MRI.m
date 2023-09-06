clear; close all; clc;

% This script can be used to return the main outcome measures for the
% serial reaction time task across the entire sample 

%ID = {'S02_CD', 'S03_DC', 'S04_VC', 'S05_IS', 'S06_TL', 'S07_YX', 'S08_TP', 'S09_LA', 'S10_EC','S11_JA','S12_MH','S13_AA','S14_DP','S15_HF','S16_LA','S17_CR','S18_SC','S19_CC','S20_TF','S21_EB','S22_BB','S23_DS','S24_GM','S25_DR','S26_SG','S27_MM', 'S28_JW', 'S29_IH','S30_WD', 'S31_CW', 'S32_TT', 'S33_RF', 'S34_RY', 'S35_RG', 'S36_CP', 'S37_ED', 'S38_EF', 'S39_JM'}; % All particpants 

%ID = {'S02_CD', 'S03_DC', 'S04_VC', 'S05_IS', 'S06_TL', 'S07_YX', 'S08_TP', 'S09_LA', 'S12_MH','S14_DP','S19_CC','S22_BB','S23_DS','S24_GM','S25_DR', 'S27_MM','S31_CW', 'S35_RG','S38_EF'}; % HIIT 

%ID = {'S10_EC','S11_JA','S13_AA','S15_HF','S16_LA','S17_CR','S18_SC','S20_TF','S21_EB', 'S26_SG', 'S28_JW', 'S29_IH', 'S30_WD', 'S32_TT', 'S33_RF', 'S34_RY',  'S36_CP', 'S37_ED',  'S39_JM'}; % Active control

ID = {'S02_CD', 'S04_VC', 'S05_IS', 'S06_TL', 'S07_YX', 'S08_TP', 'S09_LA', 'S10_EC','S11_JA','S12_MH','S13_AA','S14_DP','S16_LA','S17_CR','S18_SC','S20_TF','S22_BB','S23_DS','S24_GM','S25_DR','S26_SG','S27_MM', 'S28_JW', 'S29_IH','S30_WD', 'S31_CW', 'S32_TT', 'S33_RF', 'S34_RY', 'S35_RG', 'S36_CP', 'S37_ED', 'S38_EF', 'S39_JM'}; % Removed awareness 

%ID = {'S10_EC','S11_JA','S13_AA','S16_LA','S17_CR','S18_SC','S20_TF','S26_SG','S28_JW', 'S29_IH', 'S30_WD', 'S32_TT', 'S33_RF', 'S34_RY',  'S36_CP','S37_ED',  'S39_JM'}; % Active control (awareness removed)

%ID = {'S02_CD', 'S04_VC', 'S05_IS', 'S06_TL', 'S07_YX','S08_TP', 'S09_LA','S12_MH','S14_DP','S22_BB','S23_DS','S24_GM','S25_DR','S27_MM','S31_CW', 'S35_RG','S38_EF'}; % HIIT (awareness removed)

%timePoint = {'pre';'post'};

%timePoint = {'pre'} ;

%timePoint = {'post'};

%timePoint_number = {'1';'2'} ;

%timePoint_number = {'1'} ;

%timePoint_number = {'2'} ;

%timePoint = {'r1';'r2'};

timePoint = {'r1'} ;

%timePoint = {'r2'};

pathIn = '/Users/student/Desktop/PhD_Research/SRTT_data_analysis/SRTT_MRI_data/';

for x = 1:length(ID)
    
    for y = 1:length(timePoint)
        
        %filename = [char(pathIn),char(ID(x,1)),'/',char(ID(x,1)),'_',char(timePoint(y,1)),'_SRTT.csv'] ;
        subDir = ID{x};
        filename = [pathIn, subDir, '/' , subDir, '_', timePoint{y}, '_SRTT.csv'] ;

        [block{x,y},BlockAccuracy{x,y},blockAccuracy{x,y},blockNcorrect{x,y},meanCorrKeyPressPerSecNoOutlier{x,y},meankeypressperblock{x,y},MicroOnline{x,y},MicroOffline{x,y},TotalMicroOnline{x,y},TotalMicroOffline{x,y},TotalTaskLearning{x,y}, deltaperformance{x,y}, TotalLearning{x,y}] = srtt_work_MRI(filename);
        
    end
end
   
%% Generate group averages 

% Group average key press per block
groupavekeypressperblock = mean((cell2mat(meankeypressperblock)), 'omitnan');
stdavekeypressperblock = std((cell2mat(meankeypressperblock)), 'omitnan');%standard devation

% Group average micro-online effects
groupmicroonline = mean((cell2mat(MicroOnline)), 'omitnan'); 
stdmicroonline = std((cell2mat(MicroOnline)), 'omitnan'); %standard devation

% Group average micro-offline effects
groupmicrooffline = mean((cell2mat(MicroOffline)), 'omitnan'); 
stdmicrooffline = std((cell2mat(MicroOffline)), 'omitnan'); %standard devation

% Group average total micro-online effects
grouptotalmicroonline = mean((cell2mat(TotalMicroOnline)), 'omitnan');
stdtotalmicroonline = std((cell2mat(TotalMicroOnline)), 'omitnan');%standard devation

% Group average total micro-offline effects
grouptotalmicrooffline = mean((cell2mat(TotalMicroOffline)), 'omitnan');
stdtotalmicrooffline = std((cell2mat(TotalMicroOffline)), 'omitnan');%standard devation

% Group average total task learning effects
groupTotalTaskLearning = mean((cell2mat(TotalTaskLearning)), 'omitnan');
stdTotalTaskLearning = std((cell2mat(TotalTaskLearning)), 'omitnan');%standard devation

% Group average total early learning (for now over entire task, will
% truncate later)
groupdeltaperformance = mean((cell2mat(deltaperformance)), 'omitnan');
stddeltaperformance = std((cell2mat(deltaperformance)), 'omitnan');%standard devation

%% Early learning (i.e 95%)

%Make a call to create fit to model the group average key presses per block
%[fitresult, gof] = createFit(meankeypressperblock);
[fitresult, gof] = createFit(groupavekeypressperblock);

%use a for loop to store the modelled fit data for keypressespersecond for
%each block. e.g. this put in row 2 of groupavekeypressperblock
blocks = 16; % enter in number of blocks 
for e = 1:blocks
    modelkp{e} = fitresult.a+(fitresult.b/(1+exp(-fitresult.c*e)));
end
%obsVSmodel = [cell2mat(meankeypressperblock); cell2mat(modelkp)]; % puts obsevered values in row 1 and model values in row 2
obsVSmodel = [vertcat(meankeypressperblock); cell2mat(modelkp)]; 

%once you have this, you can calculate the block where the group modelled
%learning curve reaches 95% of total learning.(according to Bonstrup equ)
fitmax = fitresult.a+(fitresult.b/(1+exp(-fitresult.c*600))); %600 is just a big enough value to sub in for infinity
fitmin = fitresult.a+(fitresult.b/(1+exp(-fitresult.c*1)));
NFlearning = (0.95*(fitmax-fitmin))+fitmin ; 

%find the block where the modelled fit data is first above this 95% value. we need this value to
%truncate the micro on and off totals.
%NFblocknum = find(cell2mat(modelkp) > NFlearning,1);
NFblocknum = 7; % to override 95% uncomment this line and comment line above

%make a nice plot and save the figure as a .pdf or a .tiff


% truncated variables (i.e early learning)(means)
NFGroupMicroOnline = groupmicroonline(1,[1:NFblocknum]);
NFGroupMicroOffline = groupmicrooffline(1,[1:NFblocknum-1]);
NFTotalMicroOnline = sum(NFGroupMicroOnline);
NFTotalMicroOffline = sum(NFGroupMicroOffline);
EarlyLearning = groupdeltaperformance(1,[1:NFblocknum]);
TotalEarlyLearning = sum(EarlyLearning);
%standard deviations 
stdNFGroupMicroOnline = stdmicroonline(1,[1:NFblocknum]);
stdNFGroupMicroOffline = stdmicrooffline(1,[1:NFblocknum-1]);
stdEarlyLearning = stddeltaperformance(1,[1:NFblocknum]);

% participant level variables for plotting 
MOff = ((cell2mat(MicroOffline)));
NFMOff = MOff(:,[1:NFblocknum-1]);
MOn = ((cell2mat(MicroOnline)));
NFMOn = MOn(:,[1:NFblocknum]);
EL = (cell2mat(deltaperformance));
NFEL = EL(:,[1:NFblocknum]);
Accuracy = ((cell2mat(BlockAccuracy)));
NFAccuracy = Accuracy(:,[1:NFblocknum-1]);

%particiapnt averages for analysis 
AVENFMOff = mean(NFMOff,2, 'omitnan');
AVENFMOn = mean(NFMOn,2, 'omitnan');
AVENFEL = mean(NFEL,2, 'omitnan');

%particiapnt variables summed for analysis 
SUMNFMOff = sum(NFMOff,2, 'omitnan');
SUMNFMOn = sum(NFMOn,2, 'omitnan');
SUMNFEL = sum(NFEL,2, 'omitnan');

% Summary table
% % 95% learning summary data
% SummaryData = struct();
% SummaryData.NFblocknumber = NFblocknum % indicates up to which block analyses is focused on
% SummaryData.NFGroupMicroOnline = NFGroupMicroOnline; % mean micro-online changes for each block across participants 
% SummaryData.NFGroupMicroOffline = NFGroupMicroOffline;% mean micro-offline changes for each block across participants 
% SummaryData.EarlyLearning = EarlyLearning; % mean total learning
% SummaryData.NFTotalMicroOnline = NFTotalMicroOnline; % summed micro-online effects
% SummaryData.NFTotalMicroOffline = NFTotalMicroOffline; % summed micro-offline effects 
% SummaryData.TotalEarlyLearning = TotalEarlyLearning; %summed  totoal learning 
% save SummaryData.csv
% 
% particiapntlevelsummary = table(NFMOn, NFMOff, NFEL, 'VariableNames',{'MicroOnline','MicroOffline','EarlyLearning'}); % table for 95% micro-on/offline/total effects for each particiapnt 
% MeansSummaryTable = table(NFGroupMicroOnline', NFGroupMicroOffline', EarlyLearning', 'VariableNames',{'MicroOnline','MicroOffline','EarlyLearning'}); % mean 95% micro-on/offline/total effects
% TotalsSummaryTable = table(NFTotalMicroOnline, NFTotalMicroOffline, TotalEarlyLearning, 'VariableNames',{'TotalMicroOnline','TotalMicroOffline','TotalEarlyLearning'}); %summed effects 
% writetable(particiapntlevelsummary,'particiapntlevelsummary');
% writetable(MeansSummaryTable,'MeansSummaryTable');
% writetable(TotalsSummaryTable,'TotalsSummaryTable');
% 


% %% Plots 
% 
% % all micros 
% figure(2) ;
% plot(groupmicrooffline, 'r');
% hold on
% plot(groupmicroonline, 'b');
% plot(groupdeltaperformance, 'k');
% title('Block-wise learning', 'FontSize', 15);
% xlabel('Block', 'FontSize', 13);
% ylabel('Delta', 'FontSize', 13);
% legend('Micro-offline', 'Micro-online', 'Total', 'FontSize', 11);
% box off 
% 
% saveas(figure(2),'Block-wiseLearning.tiff')
% 
% MOff = ((cell2mat(MicroOffline)));
% NFMOff = MOff(:,[1:NFblocknum]);
% MOn = ((cell2mat(MicroOnline)));
% NFMOn = MOn(:,[1:NFblocknum]);
% EL = (cell2mat(deltaperformance));
% NFEL = EL(:,[1:NFblocknum]);
% 
% % truncated micro-on/offline & total (i.e Early learning)
% figure(3) ;
% plot(NFGroupMicroOffline, 'r', 'LineWidth',2);
% %shadedErrorBar([],NFMOff, {@mean,@std}, 'lineProps',{'r','markerfacecolor','r'});
% shadedErrorBar([],NFMOff,{@mean,@(NFMOff) std(NFMOff)./sqrt(size(NFMOff,1))},'lineProps',{'r','markerfacecolor','r'});
% hold on
% plot(NFGroupMicroOnline, 'b', 'LineWidth',2);
% %shadedErrorBar([],NFMOn, {@mean,@std}, 'lineProps',{'b','markerfacecolor','b'});
% shadedErrorBar([],NFMOn,{@mean,@(NFMOn) std(NFMOn)./sqrt(size(NFMOn,1))},'lineProps',{'b','markerfacecolor','b'});
% plot(EarlyLearning, 'k', 'LineWidth',2);
% %shadedErrorBar([],NFEL, {@mean,@std}, 'lineProps',{'k','markerfacecolor','k'});
% shadedErrorBar([],NFEL,{@mean,@(NFEL) std(NFEL)./sqrt(size(NFEL,1))},'lineProps',{'k','markerfacecolor','k'});
% title('Block-wsie early learning', 'FontSize', 15);
% xlabel('Block', 'FontSize', 13);
% ylabel('Delta', 'FontSize', 13);
% %legend('Micro-offline', 'Micro-online', 'Total', 'FontSize', 12);
% box off 
% 
% saveas(figure(3),'Block-wiseEarlyLearning.tiff')
% 
% % % Subplots(seperate)
% % %micro-offline 
% % figure(4) ;
% % plot(NFGroupMicroOffline, 'r', 'LineWidth',2);
% % %shadedErrorBar([],NFMOff, {@mean,@std}, 'lineProps',{'r','markerfacecolor','r'});
% % shadedErrorBar([],NFMOff,{@mean,@(NFMOff) std(NFMOff)./sqrt(size(NFMOff,1))},'lineProps',{'r','markerfacecolor','r'});
% % title('Micro-offline effect', 'FontSize', 15);
% % xlabel('Block', 'FontSize', 13);
% % ylabel('Delta', 'FontSize', 13);
% % legend('Micro-offline', 'FontSize', 12);
% % box off 
% % saveas(figure(4),'Micro-Offline.tiff')
% 
% % %micro-online 
% % figure(5) ;
% % plot(NFGroupMicroOnline, 'b', 'LineWidth',2);
% % %shadedErrorBar([],NFMOn, {@mean,@std}, 'lineProps',{'b','markerfacecolor','b'});
% % shadedErrorBar([],NFMOn,{@mean,@(NFMOn) std(NFMOn)./sqrt(size(NFMOn,1))},'lineProps',{'b','markerfacecolor','b'});
% % title('Micro-online effect', 'FontSize', 15);
% % xlabel('Block', 'FontSize', 13);
% % ylabel('Delta', 'FontSize', 13);
% % legend('Micro-online', 'FontSize', 12);
% % box off 
% % saveas(figure(5),'Micro-Online.tiff')
% 
% % %total early learning 
% % figure(6) ;
% % plot(EarlyLearning, 'k', 'LineWidth',2);
% % %shadedErrorBar([],NFEL, {@mean,@std}, 'lineProps',{'k','markerfacecolor','k'});
% % shadedErrorBar([],NFEL,{@mean,@(NFEL) std(NFEL)./sqrt(size(NFEL,1))},'lineProps',{'k','markerfacecolor','k'});
% % title('Total early learning', 'FontSize', 15);
% % xlabel('Block', 'FontSize', 13);
% % ylabel('Delta', 'FontSize', 13);
% % legend('Total', 'FontSize', 12);
% % box off 
% %saveas(figure(6),'TotalEarlyLearning.tiff')
% 
% % Subplots
% %micro-offline 
% figure(7);
% subplot(1,3,1) ;
% plot(NFGroupMicroOffline, 'r', 'LineWidth',2);
% %shadedErrorBar([],NFMOff, {@mean,@std}, 'lineProps',{'r','markerfacecolor','r'});
% shadedErrorBar([],NFMOff,{@mean,@(NFMOff) std(NFMOff)./sqrt(size(NFMOff,1))},'lineProps',{'r','markerfacecolor','r'});
% title('Micro-offline effect', 'FontSize', 15);
% xlabel('Block', 'FontSize', 13);
% ylabel('Delta', 'FontSize', 13);
% legend('Micro-offline', 'FontSize', 12);
% axis([0,NFblocknum,-1,1]); % edit with final numbers 
% box off 
% %micro-online 
% subplot(1,3,2) ;
% plot(NFGroupMicroOnline, 'b', 'LineWidth',2);
% %shadedErrorBar([],NFMOn, {@mean,@std}, 'lineProps',{'b','markerfacecolor','b'});
% shadedErrorBar([],NFMOn,{@mean,@(NFMOn) std(NFMOn)./sqrt(size(NFMOn,1))},'lineProps',{'b','markerfacecolor','b'});
% title('Micro-online effect', 'FontSize', 15);
% xlabel('Block', 'FontSize', 13);
% ylabel('Delta', 'FontSize', 13);
% legend('Micro-online', 'FontSize', 12);
% axis([0,NFblocknum,-1,1]);% edit with final numbers 
% box off 
% %total early learning 
% subplot(1,3,3) ;
% plot(EarlyLearning, 'k', 'LineWidth',2);
% %shadedErrorBar([],NFEL, {@mean,@std}, 'lineProps',{'k','markerfacecolor','k'});
% shadedErrorBar([],NFEL,{@mean,@(NFEL) std(NFEL)./sqrt(size(NFEL,1))},'lineProps',{'k','markerfacecolor','k'});
% title('Total early learning', 'FontSize', 15);
% xlabel('Block', 'FontSize', 13);
% ylabel('Delta', 'FontSize', 13);
% legend('Total', 'FontSize', 12);
% axis([0,NFblocknum,-1,1]); %edit with final numbers 
% box off 
% 
% saveas(figure(7),'subplots.tiff')
% 
% save srtt_loop
% 


 %% Group analyses 
% 
% [h_MON,p_MON,ci_MON,stats_MON] = ttest(NFGroupMicroOnline)
% 
% [h_MOFF,p_MOFF,ci_MOFF,stats_MOFF] = ttest(NFGroupMicroOffline)
% 
% [h_TOTAL,p_TOTAL,ci_TOTAL,stats_TOTAL] = ttest(EarlyLearning)
% 

%% Accuracy 

% Ave accuracy across blocks - enter appropriate data set 
AveBlockAccuracy = mean(cell2mat(BlockAccuracy));



