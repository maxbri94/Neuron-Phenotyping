function [cellStruct] = F_estimateAPDuration(data,sampling_interval,header_parameters,cellStruct, cellCount, flag, cell_name)

%% Initiate constant values 
rate = 1000/sampling_interval;  % sampling rate, in kHz 
velocity_threshold = 5; % hard threshold for determining AP start, in ms/s
sec_per_samp = 1/rate;
segTime = 20; %in ms, segment around first spike (i.e. 20ms prior to and following AP)
baselines_samps = 10;
find_exact_threshold_samples = 5;

%% Calculate action potential duration for cell

%Arrange data
mV_channel = find(strcmp(header_parameters.recChNames,'Vm_prime'));
mV_data = squeeze(data(:,mV_channel,:));
I_channel = find(strcmp(header_parameters.recChNames,'Im_sec'));
I_data = squeeze(data(:,I_channel,:)); 

%Get puls start, end, and interval
[squareBeg, squareEnd] = squarePulseDetect(I_data(:,:));
trigPoints = sort([squareBeg,squareEnd]);
traceInds = [trigPoints(1),trigPoints(2)];
interval=(trigPoints(2)-trigPoints(1)+1)/rate;

%Detect spikes
all_spike_ind = cell(size(mV_data,2),1);
for j = 1:size(mV_data,2)    
    mV_trace = mV_data(:,j); 
    spike_ind = SpikeDetect_2015_09_24(mV_trace,rate);
    %Throw out spikes that are edge artifacts
    if ~isempty(traceInds) & ismember(traceInds(2),spike_ind) 
        spike_ind(find(spike_ind==traceInds(2))) = [];
    end
    %Calculate firing rate
    if spike_ind == 0
        firing_rate(j) = 0;
    else
        firing_rate(j) = size(spike_ind,1)/interval;         
    end
all_spike_ind(j) = {spike_ind};
end

%Get spike train with lowest firing rate & first spike in that train
min_spike_train = find(firing_rate>0,1,'first');
min_spikes = cell2mat(all_spike_ind(min_spike_train));
first_spike = min_spikes(1);

%define segment around peak and extract AP onset from it 
%scale segTime according to number of spikes (i.e. higher spike frequency, shorter segment, to avoid 
%analyzing more than the first spike)
clear spike_segs
if size(min_spikes,1) > 4 & size(min_spikes,1) < 8
    segTime = segTime/2;
elseif size(min_spikes,1) > 7
    segTime = segTime/4;
elseif size(min_spikes,1) > 10
    segTime = segTime/8;
end

%get segment around spike in samples
segSamps = segTime*rate;
spike_segs = [first_spike-segSamps:first_spike+segSamps];
seg_original_length = size(spike_segs,2);
if length(min_spikes) > 1 && any(find(min_spikes(2)==spike_segs))
    second_spike = min_spikes(2);
    [cutOffVal,cutOffInd] = min(mV_data([first_spike:second_spike],min_spike_train));
    cutOff = find(spike_segs == first_spike+cutOffInd, 1, 'first');
    spike_segs = spike_segs(1:cutOff);
end

%get velocity and acceleration from spike of interest    
t = (0:size(data,1)-1)/rate;
t_seg = (0:size(spike_segs,2)-1)/rate;
seg_trace = mV_data(spike_segs,min_spike_train);
seg_whole = mV_data(:,min_spike_train);
dt = diff(t_seg); 
dt = dt(1); 
d_mV = diff(seg_trace);
v = d_mV./dt;
dv = diff(v);
a = dv./dt;
%identify first point of velocity > 5 m/s
ind_v = find(v > velocity_threshold, 1, 'first'); 
half_a = a(1:round((seg_original_length-2)/2));
%identify lowest acceleration in window comprising 5 samples prior to
%velocity threshold point
[val_a ind_a] = max(half_a);
[val_a_min ind_a_min] = min(abs(half_a(ind_v-find_exact_threshold_samples:ind_v))); 
min_ind = ind_v-(find_exact_threshold_samples+1-ind_a_min); 

%Get full width half maximum - add 2 samples since we go from acceleration to
%original data
base_trace = [min_ind+2-baselines_samps:min_ind+2]; 
fwhm_base = mean(seg_trace(base_trace));
fwhm_peak = mV_data(first_spike,min_spike_train);
half_maximum = fwhm_base+(fwhm_peak-fwhm_base)/2;

%Get amplitude
amplitude = fwhm_peak-fwhm_base;
%Get time to AP start
spike_threshold = fwhm_base;

%define samps to check for full width
spike_ind = find(spike_segs==first_spike);
first_half = spike_segs(:,1:spike_ind-1);
second_half = spike_segs(:,spike_ind+1:end);
[min_first_half_val min_first_half_ind] = mink(abs(seg_whole(first_half')-half_maximum),2);
[min_second_half_val min_second_half_ind] = mink(abs(seg_whole(second_half')-half_maximum),2);
min_first_half = sort(first_half(min_first_half_ind));
min_second_half = sort(second_half(min_second_half_ind));

%get time of maximum
%calculate with interp1, spline 
xi = [min_first_half(1)+0.001:0.001:min_first_half(2)];
yi = interp1(min_first_half, seg_whole(min_first_half),xi,'spline');
[val_fwhm ind_fwhm] = min(abs(yi-half_maximum));
secs_betw_samps_first = sec_per_samp*(ind_fwhm/length(xi));
time_first = (min_first_half(1)-1)/rate + secs_betw_samps_first;

%do the same for second sample
xi = [min_second_half(1)+0.001:0.001:min_second_half(2)];
yi = interp1(min_second_half, seg_whole(min_second_half),xi,'spline');
[val_fwhm ind_fwhm] = min(abs(yi-half_maximum));
secs_betw_samps_second = sec_per_samp*(ind_fwhm/length(xi));
time_second = (min_second_half(1)-1)/rate + secs_betw_samps_second;

%calculate duration
spike_duration = time_second - time_first; 

%% save values in struct
cellStruct.APDuration.sweepNo(cellCount) = min_spike_train;
cellStruct.APDuration.spike_ind(cellCount) = first_spike;
cellStruct.APDuration.spike_duration(cellCount) = spike_duration;
cellStruct.APDuration.spike_threshold(cellCount) = spike_threshold;
cellStruct.APDuration.spike_amplitude(cellCount) = amplitude;

%% Sanity check - make sure threshold, amplitude, and FWHM are in right
% place
if flag == 1
    figure(cellCount*10+5)
    plot(t, mV_data(:,min_spike_train))
    hold on
    plot((first_spike-1)/rate, mV_data(first_spike,min_spike_train), 'ro')
    hold on
    plot((spike_segs-1)/rate, mV_data(spike_segs,min_spike_train))
    hold on
    plot(time_first, half_maximum,'bo')
    hold on
    plot(time_second, half_maximum,'bo')
    title(strrep(['AP Duration ',cell_name],'_',' '))
    xlim([spike_segs(1)/rate (spike_segs(end)-1)/rate])
end

end



