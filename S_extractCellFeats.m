%% extractCellFeats
%This script extracts electrophysiological features from I-V data and saves
%it in a struct. 

%% Intialize analysis

clear all
folder = 'C:\Users\maxab\OneDrive\Desktop\Yizhar Matlab\Files for ephys properties\';
%Get folder names (one for each cell)
d = dir(folder);
dfolders = d([d(:).isdir]); 
%remove first two entries
cellNames = {dfolders(3:end).name};
cellCount = 1;

%Turn on/off graphs for sanity checks
%1. estimate resistance
    %a. plot resistance values for each sweep for each cell (first measurement)
scp1_flag = 0;
run1_flag = 1;
%2. estimate conductance
    %a. plot fits for each curve showing return to steady-stae
scp2_flag = 0;
run2_flag = 1;
%3. adapatation
    %a. plot peak_spike_train to make sure we have peak one
scp3_flag = 0;
run3_flag = 1;
%4. H current
    %a. check fits for H current
scp4_flag = 0;
run4_flag = 1;
%5. AP duration
    %a. check if start and end of duration are in the right place
scp5_flag = 0;
run5_flag = 1;

%% Run feature extraction functions on data

cellStruct = [];
for h = 1:size(cellNames,2) 
    %Get current cell name & create list of files for cell folder
    clear listOfFilesM
    current_cell = cell2mat(cellNames(h));
    cellNum = str2num(current_cell(end));
    current_folder = [folder, current_cell,'\'];
    DirList = dir(fullfile(current_folder, '*.abf'));
    listOfFiles = {DirList.name};
    for i = 1:length(listOfFiles)
        listOfFilesM(i,:) = cell2mat(listOfFiles(i));
    end 
    % Initiate iteration variables for loop
    mV_rec_cnt = 1;
    I_rec_cnt = 1;
    Fig_count = 1;
    file_exist = 0;
    % Skip folder if it contains no files
    if isempty(listOfFiles)
       cellStruct.Names(cellCount) = {current_cell};
       cellStruct.Date(cellCount) = NaN;
       cellStruct.cellNum(cellCount) = cellNum;
       cellCount = cellCount+1;
       continue
    end
    % Extract cell features from data
    for j = 1:size(listOfFilesM,1)
            [data,sampling_interval,header_parameters] = abfload([current_folder,listOfFilesM(j,:)]);
        % Skip Im_prime files
        if strcmp(cell2mat(header_parameters.recChNames(1)), 'Im_prime') == 1 || strcmp(cell2mat(header_parameters.recChNames(2)), 'Im_prime') == 1
            continue; 
        else
            file_exist = 1;
            % Add cell name, number, and date
            if mV_rec_cnt == 1
            cellStruct.Names(cellCount) = {current_cell};
            cellStruct.Date(cellCount) = header_parameters.uFileStartDate;
            cellStruct.cellNum(cellCount) = cellNum;
            end
            % Extract resistance measures
            if run1_flag == 1
            [cellStruct] = F_estimateResistance(data,sampling_interval,header_parameters,cellStruct, cellCount, mV_rec_cnt, scp1_flag, current_cell); 
            end
            % Extract conductance measures
            if run2_flag == 1
            [cellStruct] = F_estimateConductance(data,sampling_interval,header_parameters,cellStruct, cellCount, mV_rec_cnt, scp2_flag, current_cell);
            end
            % Extract adaptation measures
            if run3_flag == 1 & mV_rec_cnt == 1
            [cellStruct] = F_estimateAdaptation(data,sampling_interval,header_parameters, cellStruct, cellCount, scp3_flag, current_cell); 
            end
            % Extract hcurrent measures
            if run4_flag == 1 & mV_rec_cnt == 1
            [cellStruct] = F_estimateHcurrent(data,sampling_interval,header_parameters,cellStruct, cellCount,scp4_flag, current_cell);
            end
            % Extract AP duration measures
            if run5_flag == 1 & mV_rec_cnt == 1
            [cellStruct] = F_estimateAPDuration(data,sampling_interval,header_parameters,cellStruct, cellCount,scp5_flag, current_cell);
            end
            mV_rec_cnt = mV_rec_cnt+1;  
        end 
    end
    % If the folder is empty, add folder name and number anyway
    if file_exist == 0
       cellStruct.cellNum(cellCount) = cellNum;
       cellStruct.Date(cellCount) = NaN;
    end
    cellCount = cellCount+1;
end

%% Substitute 0s with NaNs in cell struct for easier analysis

%Switch zeros to NaN
    %Header parameters
cellStruct.cellNum(find(cellStruct.cellNum == 0)) = NaN;
    %AP Duration
cellStruct.APDuration.spike_duration(find(cellStruct.APDuration.spike_duration == 0)) = NaN;
cellStruct.APDuration.spike_ind(find(cellStruct.APDuration.spike_ind == 0)) = NaN;
cellStruct.APDuration.sweepNo(find(cellStruct.APDuration.sweepNo == 0)) = NaN;
cellStruct.APDuration.spike_threshold(find(cellStruct.APDuration.spike_threshold == 0)) = NaN;
cellStruct.APDuration.spike_amplitude(find(cellStruct.APDuration.spike_amplitude == 0)) = NaN;
    %Voltage difference for hCurrent
cellStruct.hCurrent.voltage_diff(find(cellStruct.hCurrent.voltage_diff == 0)) = NaN;
cellStruct.hCurrent.voltage_sag_rate(find(cellStruct.hCurrent.voltage_sag_rate == 0)) = NaN;
cellStruct.hCurrent.voltage_diff_unsmooth(find(cellStruct.hCurrent.voltage_diff_unsmooth == 0)) = NaN;
    %Adaptation
cellStruct.Adaptation.first2second_ratio_80(find(cellStruct.Adaptation.first2second_ratio_80 == 0)) = NaN;
cellStruct.Adaptation.first2second_ratio_100(find(cellStruct.Adaptation.first2second_ratio_100 == 0)) = NaN;
cellStruct.Adaptation.second2last_ratio_80(find(cellStruct.Adaptation.second2last_ratio_80 == 0)) = NaN;
cellStruct.Adaptation.second2last_ratio_100(find(cellStruct.Adaptation.second2last_ratio_100 == 0)) = NaN;
cellStruct.Adaptation.peak_spike_train(find(cellStruct.Adaptation.peak_spike_train == 0)) = NaN;
for i = 1:size(cellStruct.Adaptation.deltaIs,2)
    for j = 1:size(cellStruct.Adaptation.deltaIs,1)
        if isempty(cell2mat(cellStruct.Adaptation.deltaIs(j,i)))
            cellStruct.Adaptation.deltaIs(j,i) = {NaN};
        end
    end
end
for i = 1:size(cellStruct.Adaptation.firing_rate,2)
    for j = 1:size(cellStruct.Adaptation.firing_rate,1)
        if isempty(cell2mat(cellStruct.Adaptation.firing_rate(j,i)))
            cellStruct.Adaptation.firing_rate(j,i) = {NaN};
        end
    end
end
for i = 1:size(cellStruct.Adaptation.spike_intervals,2)
    for j = 1:size(cellStruct.Adaptation.spike_intervals,1)
        if isempty(cell2mat(cellStruct.Adaptation.spike_intervals(j,i)))
            cellStruct.Adaptation.spike_intervals(j,i) = {NaN};
        end
    end
end
    %Capacitance
cellStruct.Capacitance.membrane(find(cellStruct.Capacitance.membrane == 0)) = NaN;
cellStruct.Capacitance.time(find(cellStruct.Capacitance.time == 0)) = NaN;
    %Resistance
cellStruct.Resistance.membrane(find(cellStruct.Resistance.membrane == 0)) = NaN;
cellStruct.Resistance.time(find(cellStruct.Resistance.time == 0)) = NaN;

%% Save struct

varname = 'featureStruct';
save(varname,'cellStruct')

