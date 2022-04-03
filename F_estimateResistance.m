function [cellStruct] = F_estimateResistance(data,sampling_interval,header_parameters,cellStruct, cellCount, I_rec_cnt, flag, cell_name)

%% Initiate constant values

avgPeriod = 30; %in ms %for averaging baseline
constantSamps = 5; %constant samples deducted to be sure to onlz take the baseline we want
peak_period = 200; %ms %peak during which negative peak is determined
smoothing_period = 1; %in ms
upPercent = 0.2; 
downPercent = 0.8;
% percentage of trace between baseline and negative peak taken for fit
% during F_estimateConductance
rate = 1000/sampling_interval;  
avgPeriod_samps = rate*avgPeriod;
peak_period_samps = peak_period*rate;
smoothing_samps = 1*rate;

%% Outlier names for cells that need a different up- & down percent

outlier1 = '2018_02_22 Mice #3&4 Cell 2';
outlier2 = '2018_01_25 Mouse #1 Cell 2';
outlier3 = '2018_06_25 Mouse #56 Cell 6';
outlier4 = '2018_10_09 Mouse #61 Cell 5';
outlier5 = '2020_08_16 Mouse #101 Cell 2';

%% Calculate resistance for cell

% Arrange data
mV_channel = find(strcmp(header_parameters.recChNames,'Vm_prime'));
mV_data = squeeze(data(:,mV_channel,:));
I_channel = find(strcmp(header_parameters.recChNames,'Im_sec'));
I_data = squeeze(data(:,I_channel,:)); 

%Get start and end point 
[squareBeg, squareEnd] = squarePulseDetect(I_data(:,:));
trigPoints = sort([squareBeg,squareEnd]);
trigPoints = [trigPoints(1),trigPoints(2)];
onLengths = length(trigPoints(1):trigPoints(2));

%Identify least hyperpolarizing sweep
baseline = mean(I_data(trigPoints(1)-constantSamps-avgPeriod:trigPoints(1)-constantSamps,:),1);
stimulus = mean(I_data(trigPoints(1)+constantSamps:trigPoints(1)+constantSamps+avgPeriod,:),1);
sweepNo_Comp_Neg = find(stimulus < baseline-1);
[min_val min_ind] = min(abs(stimulus(sweepNo_Comp_Neg)));
mV_data_Neg = mV_data(:,min_ind); 
I_data_Neg = I_data(:,min_ind);

%Get change in voltage for hyperpolarization
mV_bl_samps = [(trigPoints(1)-constantSamps-avgPeriod_samps):(trigPoints(1)-constantSamps)];
peak_search_samps = [trigPoints(1)+constantSamps:trigPoints(1)+constantSamps+peak_period_samps];
[min_val min_ind] = min(mV_data_Neg(peak_search_samps));
mV_ss_samps_Neg = trigPoints(1)+constantSamps+min_ind-1;
smoothing_period_samps = mV_ss_samps_Neg-(smoothing_samps/2):mV_ss_samps_Neg+(smoothing_samps/2);
min_val_smooth = mean(mV_data_Neg(smoothing_period_samps,1));
mV_bl_Neg = mean(mV_data_Neg(mV_bl_samps),1);
mV_ss_Neg = min_val_smooth;
deltamV_Neg = mV_ss_Neg-mV_bl_Neg;

%Get change in current for hyperpolarization
I_bl_Neg = mean(I_data_Neg(mV_bl_samps),1);
I_ss_Neg = mean(I_data_Neg(mV_ss_samps_Neg),1);
deltaI_Neg = I_ss_Neg-I_bl_Neg;

%Calculate resistance and time of recording
membrane_resistance_Neg = deltamV_Neg/deltaI_Neg;
time_for_cell = header_parameters.uFileStartTimeMS/1000/60/60; %time of day in hours

%Adapt extend of fitting trace for outlier cells
if strcmp(outlier1,cell_name) || strcmp(outlier2,cell_name) || strcmp(outlier4,cell_name) 
    upPercent = 0.1;
    downPercent = 0.9;
elseif strcmp(outlier3,cell_name) || strcmp(outlier5,cell_name)
    upPercent = 0.025;
    downPercent = 0.975;
end

%Get fit samps for hyperolarization
%Get whole trace for fitting
sampsStart_Neg = trigPoints(1)+constantSamps;
sampsEnd_Neg = mV_ss_samps_Neg;
allSamps_Neg = [sampsStart_Neg:sampsEnd_Neg];
%Choose samples from 20% to 80% of fitting trace
fitStart_Neg = allSamps_Neg(round(upPercent*length(allSamps_Neg)));
fitEnd_Neg = allSamps_Neg(round(downPercent*length(allSamps_Neg)));
fitStEnd_Neg = [fitStart_Neg fitEnd_Neg];
fitStEnd_Neg_Length = fitEnd_Neg-fitStart_Neg;

%% Save in struct
cellStruct.Resistance.membrane(cellCount, I_rec_cnt) = membrane_resistance_Neg;
cellStruct.Resistance.time(cellCount, I_rec_cnt) = time_for_cell;
cellStruct.Resistance.deltamV_Neg(cellCount, I_rec_cnt) = deltamV_Neg;
cellStruct.Resistance.deltaI_Neg(cellCount, I_rec_cnt) = deltaI_Neg;
cellStruct.Resistance.fitSamps_Neg(cellCount, I_rec_cnt) = {fitStEnd_Neg};

%% Sanity check
if flag == 1
    t = (0:size(data,1)-1)/rate;
    figure(cellCount*10+1)
    plot(t, [mV_data_Neg,I_data_Neg])
    hold on
    plot([fitStart:fitEnd]/rate, mV_data_Neg([fitStart:fitEnd]))
    hold on
    plot(sampsStart/rate, mV_data_Neg(sampsStart),'ro')
    hold on
    plot(sampsEnd/rate, mV_data_Neg(sampsEnd),'ro')
    hold on
    plot(mV_bl_samps/rate, mV_data_Neg(mV_bl_samps))
    hold on
    plot(mV_ss_samps/rate, mV_data_Neg(mV_ss_samps))
    hold on
    plot(mV_bl_samps/rate, I_data_Neg(mV_bl_samps))
    hold on
    plot(mV_ss_samps/rate, I_data_Neg(mV_ss_samps))
    xlabel('Time in ms')
    ylabel('Activity in mV')
    sgtitle(strrep(['Resistance Sanity check',cell_name],'_',' '))
end


end