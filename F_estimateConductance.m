function [cellStruct] = F_estimateConductance(data,sampling_interval,header_parameters,cellStruct, cellCount, I_rec_cnt, flag, cell_name)
%Note: Must run estimate resistance function before running this function
%on a cell

%% Initiate constant values
avgPeriod = 30; %in ms, for getting baseline and stimulus period average
constantSamps = 5;
rate = 1000/sampling_interval; 

%% Calculate conductance for cell

%Arrange data
mV_channel = find(strcmp(header_parameters.recChNames,'Vm_prime'));
mV_data = squeeze(data(:,mV_channel,:));
I_channel = find(strcmp(header_parameters.recChNames,'Im_sec'));
I_data = squeeze(data(:,I_channel,:)); 

%Get start and end point
[squareBeg, squareEnd] = squarePulseDetect(I_data(:,:));
trigPoints = sort([squareBeg,squareEnd]);
trigPoints = [trigPoints(1),trigPoints(2)];
onLengths = length(trigPoints(1):trigPoints(2));

%Identify sweep of interest - most negative sweep
baseline = mean(I_data(trigPoints(1)-constantSamps-avgPeriod:trigPoints(1)-constantSamps,:),1);
stimulus = mean(I_data(trigPoints(1)+constantSamps:trigPoints(1)+constantSamps+avgPeriod,:),1);
sweepNo_Comp_Neg = find(stimulus < baseline-1);
[min_val min_ind] = min(abs(stimulus(sweepNo_Comp_Neg)));
mV_data_Neg = mV_data(:,min_ind); 
I_data_Neg = I_data(:,min_ind);

%Get baseline
mV_bl_samps = [(trigPoints(1)-constantSamps-avgPeriod):(trigPoints(1)-constantSamps)];
mV_bl_Neg = mean(mV_data_Neg(mV_bl_samps),1);

%Get fit vectors - taking fitPoints from F_estimateResistance and
%normalize baseline
fitPoints = cell2mat(cellStruct.Resistance.fitSamps_Neg(cellCount, I_rec_cnt));
fitTime = ([0:fitPoints(2)-fitPoints(1)]/rate);
fitSamps_Neg = mV_data_Neg(fitPoints(1):fitPoints(2));
fitSamps_Neg_normed = fitSamps_Neg-mV_bl_Neg;

% Set up fittype and options.
[xData, yData] = prepareCurveData( fitTime, fitSamps_Neg_normed );
ft = fittype( 'a*(1-exp(-x/b))+C', 'independent', 'x', 'dependent', 'y' ); %positive values
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -100 -Inf];
opts.Upper = [Inf 0 Inf];

% Fit model to data and extract tau
[fitresult, gof] = fit( xData, yData, ft, opts );
coefficients = coeffvalues(fitresult);
tau_time_constant = coefficients(:,3);

%% Save in struct
membrane_resistance = cellStruct.Resistance.membrane(cellCount);
membrane_capacitance = tau_time_constant*(1/membrane_resistance);
cellStruct.Capacitance.membrane(cellCount,I_rec_cnt) = membrane_capacitance;
cellStruct.Capacitance.time(cellCount,I_rec_cnt) = cellStruct.Resistance.time(cellCount, I_rec_cnt);

%% Sanity check - check fit
if flag == 1 & I_rec_cnt == 1
    figure(cellCount*10+2)
    h = plot( fitresult, xData, yData );
    legend(h, 'Data', 'Fit', 'Location', 'SouthEast' )
    text(mean(xData), mean(yData)+1, ['Rsquare: ',num2str(round(gof.rsquare,3))])
    xlabel('Time in ms')
    xlim([0 xData(end)])
    ylabel('Activity in mV')
    sgtitle(strrep(['Time constant fit mV over time',cell_name],'_',' '))
end

end



