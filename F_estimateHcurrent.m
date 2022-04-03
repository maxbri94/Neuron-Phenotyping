function [cellStruct] = F_estimateHcurrent(data,sampling_interval,header_parameters,cellStruct, cellCount, flag, cell_name)

%% Initiate constant values

mVPeriod = 1000; %in samples (i.e. 50ms --> 1000 samples)
constantSamps = 5; %constant samples deducted to be sure to only take the baseline we want
smoothing_period = 1; %in ms, smoothing around peak chosen for h current
avgPeriod = 50; %period for determining baseline and steady state
hC_period = 0.2; %hC must occur during first 20% of period during pulse
rate = 1000/sampling_interval;  % sampling rate, in kHz 
smoothing_samps = smoothing_period*rate;

%% Calculate hCurrent for cell

% Arrange data
mV_channel = find(strcmp(header_parameters.recChNames,'Vm_prime'));
mV_data = squeeze(data(:,mV_channel,:));
I_channel = find(strcmp(header_parameters.recChNames,'Im_sec'));
I_data = squeeze(data(:,I_channel,:));
onLengths = header_parameters.sweepLengthInPts/2; %get length of sweep in samples %change this

%Get start and end of pulse
[squareBeg, squareEnd] = squarePulseDetect(I_data(:,:));
trigPoints = sort([squareBeg,squareEnd]);
trigPoints = [trigPoints(1),trigPoints(2)];

% Determine which sweeps are hyperpolarizing and take first sweep
baseline = mean(I_data(trigPoints(1)-constantSamps-avgPeriod:trigPoints(1)-constantSamps,:),1);
stimulus = mean(I_data(trigPoints(2)-constantSamps-avgPeriod:trigPoints(2)-constantSamps,:),1);
sweepNo_Comp = find(stimulus < baseline-1);
hyperSweep = sweepNo_Comp(1);

% Find minimum index and baseline, steady-state samples
half_mV_trace = mV_data(1:onLengths,hyperSweep);
[min_val min_ind] = min(half_mV_trace);
ind_during_pulse = min_ind-trigPoints(1);
hC_cutOff = hC_period*onLengths;
ss_samps = [trigPoints(2)-constantSamps-avgPeriod*rate:trigPoints(2)-constantSamps];
bl_samps = [trigPoints(1)-constantSamps-avgPeriod*rate:trigPoints(1)-constantSamps];

%determine h current ratio
%determine h-current only if it occurs during first 
if ind_during_pulse > hC_cutOff
    cellStruct.hCurrent.voltage_diff(cellCount) = NaN;
    cellStruct.hCurrent.voltage_diff_unsmooth(cellCount) = NaN;
    cellStruct.hCurrent.voltage_sag_rate(cellCount) = NaN;
else
    %get smoothed min_val
    smoothing_period_samps = min_ind-(smoothing_samps/2):min_ind+(smoothing_samps/2);
    min_val_smooth = mean(half_mV_trace(smoothing_period_samps,1));
    %get differences between min_value, steady state and baseline
    steady_state = mean(mV_data(ss_samps,hyperSweep));
    baseline_of_interest = mean(mV_data(bl_samps,hyperSweep));
    diff_hC_ss = min_val_smooth - steady_state;
    diff_hC_bl = min_val_smooth - baseline_of_interest;

    %% Save in struct
    cellStruct.hCurrent.voltage_diff(cellCount) = diff_hC_ss;
    cellStruct.hCurrent.voltage_diff_unsmooth(cellCount) = min_val - steady_state;
    cellStruct.hCurrent.voltage_sag_rate(cellCount) = diff_hC_ss/diff_hC_bl;
end

%% Sanity check with flagged graph
if flag == 1
    t = (0:size(data,1)-1)/rate;
    cutOff_samps = [trigPoints(1):trigPoints(1)+hC_cutOff];
    figure(cellCount*10+4)
    plot(t, mV_data(:,hyperSweep))
    y_bottom = min(mV_data(:,hyperSweep))-5;
    y_top = max(mV_data(:,hyperSweep))+5;
    ylim([y_bottom y_top])
    hold on
    xline(min_ind/rate)
    if ind_during_pulse <= hC_cutOff
        text(min_ind/rate, mV_data(min_ind,hyperSweep), num2str(cellStruct.hCurrent.voltage_sag_rate(cellCount)))
    end
    hold on
    ss_start = ss_samps(1);
    ss_dis = ss_samps(end)-ss_start;
    mv_dis = y_top-y_bottom;
    rectangle('Position',[ss_start/rate,y_bottom,ss_dis/rate,mv_dis],'FaceColor',[0.3010 0.7450 0.9330])
    title(strrep(['h-current ',cell_name],'_',' '))
end

end