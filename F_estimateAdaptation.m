function [cellStruct] = F_estimateAdaptation(data,sampling_interval,header_parameters,cellStruct, cellCount, flag, cell_name)

%% Initiate constant values 
constantSamps = 5; %constant samples deducted from pulse start and end to account
%for inconsistencies in squarePulseDetect (sometimes doesn't find the exact start and end sample of the pulse)
avgPeriod = 50; %average period for baseline and stimulus period %in samples
rate = 1000/sampling_interval;  % sampling rate, in kHz % check out what the relationship between kHz and mV is

%% Define outliers - i.e. cells for which we must adapt parameters to
% extract their features

outlier_string4 = ['2018_06_25 Mouse #56 Cell 6'];
outlier_string3 = ['2017_12_20 Mouse #33 Cell 3'];
outlier_string2 = ['2017_12_20 Mouse #33 Cell 2'];
outlier_string1 = ['2018_01_18 Mouse #38 Cell 3'];

%% Calculate adaptation for cell

%Arrange data
mV_channel = find(strcmp(header_parameters.recChNames,'Vm_prime'));
mV_data = squeeze(data(:,mV_channel,:));
I_channel = find(strcmp(header_parameters.recChNames,'Im_sec'));
I_data = squeeze(data(:,I_channel,:)); 
I_trace = I_data(:,:);

%Get start and end of pulse
[squareBeg, squareEnd] = squarePulseDetect(I_data(:,:));
trigPoints = sort([squareBeg,squareEnd]);
trigPoints = [trigPoints(1),trigPoints(2)];
interval=(trigPoints(2)-trigPoints(1)+1)/rate;

%get deltaI for each sweep
baseline = mean(I_data(trigPoints(1)-constantSamps-avgPeriod:trigPoints(1)-constantSamps,:),1);
stimulus = mean(I_data(trigPoints(2)-constantSamps-avgPeriod:trigPoints(2)-constantSamps,:),1);
deltaI = stimulus-baseline;

% Detect spikes
all_spike_ind = cell(size(mV_data,2),1);
for j = 1:size(mV_data,2)    
    mV_trace = mV_data(:,j); 
    spike_ind = SpikeDetect_2015_09_24(mV_trace,rate);
    %Throw out edge artifacts that were idetnfied as spikes
    if ~isempty(traceInds) & ismember(traceInds(2),spike_ind) 
        spike_ind(find(spike_ind==traceInds(2))) = [];
    end
    %Get firing rate
    if spike_ind == 0
        firing_rate(j) = 0;
    else
        firing_rate(j) = size(spike_ind,1)/interval;   
        if strcmp(outlier_string3,cell_name) & j == 13
            spike_ind = spike_ind(1:end-1);
            firing_rate(j) = size(spike_ind,1)/interval;  
        end
    end
all_spike_ind(j) = {spike_ind}; 
end

%Get spike train that reaches peak spike frequency first - we call this
%variable norm_spike_train
norm_spike_train = find(firing_rate==max(firing_rate),1,'first');

%Find spike train with frequency closest to 80% of norm_spike_train - this
%will be our peak_spike_train
firing_rate_percent = firing_rate/firing_rate(norm_spike_train)*100;    
[peak_spike_train_v peak_spike_train_i] = min(abs(firing_rate_percent(1:norm_spike_train)-80));
peak_spike_train = peak_spike_train_i;

% Correct outlier cells - either add undetected spike or choose different
% spike train
if strcmp(outlier_string1,cell_name)
    %change to previous train
    peak_spike_train = peak_spike_train-1;
elseif strcmp(outlier_string4,cell_name)
    %change to previous train
    peak_spike_train = peak_spike_train-1;
elseif strcmp(outlier_string2,cell_name)
    %add spike
    firing_rate(norm_spike_train) = (firing_rate(norm_spike_train)*interval+1)/interval;
elseif strcmp(outlier_string3,cell_name)
    %add spike
    firing_rate(peak_spike_train) = (firing_rate(peak_spike_train)*interval+1)/interval; 
    peak_spike_train = peak_spike_train-1; 
end

%ratio100
spike_times = cell2mat(all_spike_ind(norm_spike_train))/rate;
spike_intervals = spike_times(1+1:end)-spike_times(1:end-1); 
first2second_ratio_100 = spike_intervals(2)/spike_intervals(1);
first2last_ratio_100 = spike_intervals(end)/spike_intervals(1);
second2last_ratio_100 = spike_intervals(end)/spike_intervals(2);
third2last_ratio_100 = spike_intervals(end)/spike_intervals(3);

%ratio 80
spike_times = cell2mat(all_spike_ind(peak_spike_train))/rate;
spike_intervals = spike_times(1+1:end)-spike_times(1:end-1); 
first2second_ratio_80 = spike_intervals(2)/spike_intervals(1);
first2last_ratio_80 = spike_intervals(end)/spike_intervals(1);
second2last_ratio_80 = spike_intervals(end)/spike_intervals(2);
third2last_ratio_80 = spike_intervals(end)/spike_intervals(3);

%% save relevant data in struct
%for I/O graph
cellStruct.Adaptation.firing_rate(cellCount,:) = {firing_rate};
cellStruct.Adaptation.spike_intervals(cellCount,:) = {spike_intervals'};
cellStruct.Adaptation.all_spike_ind(cellCount,:) = {all_spike_ind'};
cellStruct.Adaptation.peak_spike_train(cellCount,:) = norm_spike_train; %take highest spike_train
cellStruct.Adaptation.deltaIs(cellCount,:) = {deltaI};
cellStruct.Adaptation.first2second_ratio_80(cellCount,:) = first2second_ratio_80;
cellStruct.Adaptation.first2second_ratio_100(cellCount,:) = first2second_ratio_100;
cellStruct.Adaptation.second2last_ratio_100(cellCount,:) = second2last_ratio_100;
cellStruct.Adaptation.second2last_ratio_80(cellCount,:) = second2last_ratio_80;

%% Sanity check - check 80% and 100% spike train to see that all spikes were
% detected
plot_name_100 = 'peak spike train';
plot_name_80 = '80% of peak spike train';
if flag == 1 
    %Prep data for plot
    t = (0:size(data,1)-1)/rate;
    figure(cellCount*10+3)
    plotCnt = 1;
    for m = [norm_spike_train peak_spike_train]
        if plotCnt == 1
            plotInd = [1:2];
            name = [plot_name_100, ' firing rate: ', num2str(firing_rate(norm_spike_train))];
            plotInd2 = 5;
            curr_spikes = cell2mat(all_spike_ind(m));
            dispTrain_end = curr_spikes(3)+((curr_spikes(4)-curr_spikes(3))/2);
            axis_lim = [traceInds(1)/rate dispTrain_end/rate];
            ratio = first2second_ratio_100;
        elseif plotCnt == 2
            plotInd = [3:4];
            name = [plot_name_80, ' firing rate: ', num2str(firing_rate(peak_spike_train))];
            plotInd2 = 6;
            curr_spikes = cell2mat(all_spike_ind(m));
            dispTrain_end = curr_spikes(3)+((curr_spikes(4)-curr_spikes(3))/2);
            axis_lim = [traceInds(1)/rate dispTrain_end/rate];
            ratio = first2second_ratio_80;
        end
        %plot whole spike train
        subplot(3,2,plotInd)
        plot(t, mV_data(:,m))
        hold on
        plot(cell2mat(all_spike_ind(m))/rate, mV_data(cell2mat(all_spike_ind(m)), m), 'ro')
        text(cell2mat(all_spike_ind(m))/rate,mV_data(cell2mat(all_spike_ind(m)), m),[cellstr(string(1:length(cell2mat(all_spike_ind(m)))))])    
        xlabel('Time in ms')
        ylabel('Activity in mV')
        title(name)
        plotCnt = plotCnt + 1;
        xlim([traceInds(1) traceInds(2)]/rate)
        subplot(3,2,plotInd2)
        plot(t, mV_data(:,m))
        xlim(axis_lim)
        text(curr_spikes(2)/rate, mV_data(curr_spikes(2),m), num2str(ratio))
        title(name)
    end
    sgtitle(strrep(['Adaptation ',cell_name],'_',' '))
end



