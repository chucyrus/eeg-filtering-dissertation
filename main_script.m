% Main file - script for EEG testing goes here


% init 
%participantList = [1:Np];
 participantList = 1:4;
Ns = length(participantList);
sessionList = [1:10];
Ns = length(sessionList);

%eeg_folder = fullfile('..','data','ft10_eeg')
eeg_folder = '../ft10_eeg/';
%beh_folder = fullfile('..','data','ft10_behaviour')
beh_folder = 'ft10_behaviour/';
cedFile = fullfile('Glasgow_BioSemi_132.ced');
exptname = 'ft10';

load('noise.mat')

% load generated signal 
EEG_signal = pop_loadset('eeg_signal.set'); 

% load grand average
load("grand_average.mat"); 

% load noise 
EEG_noise = pop_loadset('eeg_noise.set');

trialCount = size(trialMatrix, 3);     % getting total count of generated trials




%%

% highpass (baseline)

EEG_filt_highpass = pop_eegfilt(EEG_signal, 0.5,0, 135);    % manually set the filter order
EEG_filt_highpass.data = reshape(EEG_filt_highpass.data, 128,410,trialCount)

%%
% OSF

% 1. Filtering signal + noise
% modified pop_optimalsubspacefilter such that it actually outputs the
% weights :) 
[EEG_filt_osf, com, filterWeights] = pop_optimalsubspacefilter(EEG_signal, 'mse', [EEG_signal.xmin 0] ,[0 EEG_signal.xmax] ,'channel');

% 2. Use filter weights to filter grand average 

grand_average_filtered = filterWeights*grand_average(:,:); 

% 3. Use filter weights to filter noise 
noise_filtered = zeros(size(EEG_noise.data)); 

for i = 1:size(noiseMatrix, 3)
    noise_filtered(:,:,i) = filterWeights*EEG_noise.data(:,:,i); 
end 



% metric evaluation 
[meanInputSNR_osf, meanOutputSNR_osf, meanNoiseReductionFactor_osf, snrImprovement_osf] = snr_avg(grand_average, EEG_noise.data, grand_average_filtered, noise_filtered); 




%% 

% EEGLAB pipeline
% from https://github.com/sccn/eeg_pipelines/blob/master/eeglab/process_eeglab_template.m

removeChans = {};conditions  = { 'oddball_with_reponse' 'standard' }; % conditions to extract
 % channels to ignore

% end of parameters ************

% load data
EEG_eeglab = EEG_signal
EEG_eeglab.data = reshape(EEG_eeglab.data, 128,410*trialCount)


% filter data
EEG_eeglab = pop_eegfiltnew(EEG_eeglab, 'locutoff',0.5);
EEG_eeglab = pop_select( EEG_eeglab, 'nochannel',removeChans); % list here channels to ignore
chanlocs = EEG_eeglab.chanlocs;


%{
% remove bad channels and bad portions of data
EEG_eeglab = pop_clean_rawdata(EEG_eeglab, 'FlatlineCriterion',4,'ChannelCriterion',0.85,'LineNoiseCriterion',4,'Highpass','off',...
    'BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
%}


% Run ICA and IC Label
EEG_eeglab = pop_runica(EEG_eeglab, 'runica'); % Use mode standard for Infomax
EEG_eeglab = pop_iclabel(EEG_eeglab);
EEG_eeglab = pop_icflag(EEG_eeglab, [NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
EEG_eeglab = pop_subcomp( EEG_eeglab, [], 0);

% Interpolate removed channels 
EEG_eeglab = pop_interp(EEG_eeglab, chanlocs);
EEG_eeglab.data = reshape(EEG_eeglab.data, 128,410, trialCount); 



%{

% epoch extraction and saving of datasets


[~,fileBase] = fileparts(EEG_eeglab.filename);
EEGcond1 = pop_epoch( EEG_eeglab, {  conditions{1} }, [-0.3 0.7], 'newname', 'Cond1', 'epochinfo', 'yes');
EEGcond1 = pop_saveset( EEGcond1, 'filename',[fileBase '_cond1_eeglab.set'],'filepath',EEG_eeglab.filepath);

EEGcond2 = pop_epoch( EEG_eeglab, {  conditions{2}  }, [-0.3 0.7], 'newname', 'Cond2', 'epochinfo', 'yes');
EEGcond2 = pop_saveset( EEGcond2, 'filename',[fileBase '_cond2_eeglab.set'],'filepath',EEG_eeglab.filepath);

%}

%%
% Data visualization 
% The x axis is seconds, -0.3 to 0.5 seconds
% creating a series such that plot exists on a -0.3 to 0.5 scale
series = linspace(-0.3, 0.5, 410);

% Data plotting 

chanToPlot = 27; 

for i = 1:trialCount
    plot(series, squeeze(EEG_signal.data(chanToPlot, :, i)), 'k')
    hold on; 
    plot(series, squeeze(EEG_filt_osf.data(chanToPlot, :, i)), 'b')

end 

% Misc figure labels
xlabel('Time (seconds)');
ylabel('Amplitude (microvolts)'); 


