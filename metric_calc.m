%% 

% Metric calculation 
% Implement: SNR MSE t-test p-score 
% Purporse: to compare with grand average to obtain objective 'comparison
% metrics' 


% grand average = 'grand average' dataset
% filtered = filtered dataset (128*410*trialCount)


function [] = metric_calc(grand_average, filtered)
trialCount = trialCount = size(filtered, 3);     % getting total count of generated trials
%{

!!! REVAMP THIS - FUNCTIONS DO NOT WORK WITHIN NESTED FUNCTIONS


% SNR (single-trial) 
function [inputSNR, outputSNR, noiseReductionFactor] = snr_single(signal, unfiltNoise, filtSignal, filtNoise)    
    sigVar = var(signal);   % variance of GA (marked as signal) 
    noiseVar = var(unfiltNoise);  % variance of noise SEPERATE from sig+noise
    filtSigVar = var(filtSignal);    % variance of 'filtered GA' (marked as filtSignal)
    filtNoiseVar = var(filtNoise);   % variance of filtered noise SEPERATE
    
    inputSNR = sigVar/noiseVar; 
    outputSNR = filtSigVar/filtNoiseVar; 
    noiseReductionFactor = noiseVar/filtNoiseVar; 

end 

% SNR (average)
function [meanInputSNR, meanOutputSNR, meanNoiseReductionFactor, snrImprovement] = snr_avg(signal, noiseMatrix, filtSigMatrix, filtNoiseMatrix)

    meanSigVar = var(mean(signal(:)));    % variance of GA (marked as signal) -> typically this is the same value
    meanNoiseVar = var(mean(noiseMatrix(:), 3));    % variance of noise SEPERATE from sig+noise 
    meanFiltSigVar = var(mean(filtSigMatrix(:)));    % variance of filtered GA matrix
    meanFiltNoiseVar = var(mean(filtNoiseMatrix(:), 3)); % variance of filtered noise SEPERATE from sig+noise 
    
    
    meanInputSNR = meanSigVar/meanNoiseVar;
    meanOutputSNR = meanFiltSigVar/meanFiltNoiseVar;
    meanNoiseReductionFactor = meanNoiseVar/meanFiltNoiseVar;
    snrImprovement = meanOutputSNR/meanInputSNR

    fprintf('Statisitcs after averaging %g trials:\n', size(filtNoiseMatrix,3))
    fprintf('Noise Reduction factor: %g.\n', meanNoiseReductionFactor)
    fprintf('Input SNR: %g.  Output SNR: %g, SNR Improvement factor: %g \n', meanInputSNR, meanOutputSNR, meanOutputSNR/meanInputSNR)

end 

%} 




% MSE average 
% MSE tests how close the predictions are to the actual values 
% in this case predictions = filtered, actual values = ground truth 
% this will output 128, trialCount.


mse_avg = zeros(size(filtered, 1), size(filtered, 3));   % should be 128, trialCount
for j = 1:size(filtered, 3) %   for each trial
    for i = 1:size(filtered, 1) % for each electrode
    
        filt_temp = filtered(i, :, j); 
        grand_average_temp = grand_average(i, :);
        mse_temp = mean((filt_temp - grand_average_temp).^2);   % mean squared error
        mse_avg(i,j) = mse_temp; % Store the MSE
    end 
end



% RMSE 
% this is simply the square root of MSE

rmse_avg = sqrt(mse_avg)

% need to ask Justin if implementation of Bonferroni is necessary 




% unpaired samples ttest
% outputs a 128*trialNo matrix, 
% representing the significance value for each electrode (128) in each
% trial (trialNo)

T_mat = []
for j = 1:size(filtered, 3) %   for each trial

    filt_temp_ttest = filtered(:, :, j);  % 1 trial worth of ttest 
    [T df] = ttest_cell(grand_average, squeeze(filt_temp_ttest)); % EEGLAB ttest_cell
    T_mat(:,j) = squeeze(F); 
end 

% p-score 
signif = []; 

for j = 1:size(T_mat, 2)    % trial
    for i = 1:size(T_mat, 1)    % electrode
        signif_temp = tcdf(T_mat(i,j), df); 
        signif(i,j) = signif_temp; 
    end 
end 

% tally up significant p-scores
% trialNo*1
% this is a very dumb way to do this

sig_channels = []
for j = 1:size(signif, 2)   % trial
    s = 0; 
    for i = 1:size(signif, 1)   % electrode 

        if signif(i,j) <= 0.05
            s = s+1; 
        end 
    end 
    sig_channels(j) = s; 
end 
