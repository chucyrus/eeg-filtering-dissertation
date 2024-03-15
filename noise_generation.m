% generating noise component of simulated signal 

Np = 1; % number of participants
%participantList = [1:Np];
participantList = 1:4;
Ns = length(participantList);
sessionList = 1:10;
Ns = length(sessionList);


%eeg_folder = fullfile('..','data','ft10_eeg')
eeg_folder = 'ft10_eeg/';
%beh_folder = fullfile('..','data','ft10_behaviour')
beh_folder = 'ft10_behaviour/';
cedFile = fullfile('Glasgow_BioSemi_132.ced');
exptname = 'ft10';
ft10_eeg_trials = zeros(Np,Ns); % save number of trials
ft10_eeg_beh_match = zeros(Np,Ns);
filepathout = pwd; % where the epoch file is saved

% init variables
avg_noise = [];
trialMatrix = []; 
noiseMatrix = [];
cycle = 100; % each 'cycle' generates 10 noisy events, eg 10 'cycles' = 100 'events'
eegTimes = [];
trialNo = []; 
global times;  
times = []; 
trial = zeros(128, 410); 
uniqueTrials = []; 


% loading grand average file
% load('grand_average.mat')

% specify grand average to use 
grand_average = grand_average; 

% test file
% EEG = pop_readbdf('ft10_p1s1.bdf', [],137, 1:128);

 for t = 1:cycle 
        error = true; 
        flag = 0; 
        while error == true
            try
        
         
        P = randi([1 4]); % no participants (1-4) 
        ses = randi([1 10]); % no sessions (1-10) 

        % load random data and set reference to the average of electrodes [1:128]
        id = [exptname,'_p',num2str(P),'s',num2str(ses)];
        bdffile = fullfile(eeg_folder,[id,'.bdf']);
        %EEG = pop_readbdf(bdffile, [],137, 1:128);
        EEG_trial = pop_biosig(bdffile,'importannot','off');
        EEG_trial = eeg_checkset( EEG_trial );
        
        
        % remove electrodes [133:143]
        EEG_trial = pop_select(EEG_trial, 'channel',1:132);
        % electrodes 1 to 128 are the standard Biosemi electrodes
        % electrodes 129 to 132 are extra-ocular channels
        % see readme file for details
        
        %Average reference
        % EEG_trial = pop_reref(EEG_trial,1);
        
        
        % epoch continuous data
        % originally the epoch time was set to -1.5 1.5

        % this step occasionally throws up errors during this step, might be 
        % due to the fact that is it attempting to epoch missing trials 
        % ''Index exceeds the number of array elements. Index must not
        % exceed 1.''
        % in the following code: 
        % if there is an error, retry with a new participant and session
        % number 
        % terminates when 5 errors occur


                EEG_trial = pop_epoch( EEG_trial, {'11' '12' '21' '22' '31' '32' '41' '42' ...
            '51' '52' '61' '62' '71' '72' '81' '82' '91' '92' ...
            '101'  '102'  '111'  '112'}, [-1.101 -0.300], ...    % prestimulus 0.8s interval -> in seconds  
            'newname', 'ft10 epoched file', 'epochinfo', 'yes'); % for some unknown reason -1.1 - -0.3 extracts 409 datapoints      
                EEG_trial = eeg_checkset( EEG_trial );
                error = false; % carry on 
            catch ME 
                disp ([ME.message, id])
                flag = flag + 1; 
                if flag == 5
                    break;  
                end 

            end 

        end 
     
            temp = []; 
            temp(1) = P; 
            temp(2) = ses;  

        % choose x random trials
        while(true)
         
            trialNo = []; 
            uniqueTrials = []; 

                for x = 1:10
                    trialNo = cat(1, trialNo, randi([1 size(EEG_trial.data, 3)]));
                end 
            % internal check 
            % if no duplicates occur, break out of the while loop - do until loop  
            
            uniqueTrials = unique(trialNo'); 
            if length(uniqueTrials) == length(trialNo') 
                break; 
            end 
     
        end 


        for y = 1:length(trialNo)    
            % noise generation for 128 channels
            noise = EEG_trial.data(1:128, :, trialNo(y)); 
            scale = rand;
            if scale < 0.2 
                scale = scale + 0.2; 
            end 

            noise = noise * scale; 

            for e = 1:128 
                t = grand_average(e,:) + noise(e,:); 
             trial(e,:) = squeeze(t); 
          
            end 

         % saving simulated signal+noise into trialMatrix 
         trialMatrix = cat(3, trialMatrix, trial); 
         % saving simulated noise into noiseMatrix 
         noiseMatrix = cat(3, noiseMatrix, noise); 

        end 

         % save trialNo alongside parti/ses number 
         % number #1 = participant number 
         % number #2 = session number 
         % numbers #3 - #12 = times 1 through 10 
         times_temp = [];
         times_temp(1,1) = P; 
         times_temp(2,1) = ses;  
         times_temp(3:12,1) = trialNo;  
         times = cat(2, times, times_temp); 

 end 
        


save('noise.mat', 'noiseMatrix', 'trialMatrix', 'times');

%% 

% EEGLAB packaging 
% required for OSF 
% saves in .set format within current dir
% technically only EEG_signal is required to be in .set but 

elec_locs = readlocs('Glasgow_BioSemi_132.ced'); 
elec_locs = elec_locs(1:128);

EEG_signal = eeg_emptyset;

EEG_signal.data = trialMatrix;  
EEG_signal.nbchan = size(trialMatrix, 1);  % 128
EEG_signal.pnts = size(trialMatrix, 2);  % 410 - no. of datapoints
EEG_signal.trials = size(trialMatrix, 3);  % Number of trials
EEG_signal.srate = 512;  
EEG_signal.times = linspace(-0.3, 0.5, 410);  % timepoints
EEG_signal.xmin = EEG_signal.times(1);
EEG_signal.xmax = EEG_signal.times(end); 
EEG_signal.chanlocs = elec_locs; 


EEG_noise = eeg_emptyset;

EEG_noise.data = noiseMatrix;  
EEG_noise.nbchan = size(noiseMatrix, 1);  % 128
EEG_noise.pnts = size(noiseMatrix, 2);    % 410
EEG_noise.trials = size(noiseMatrix, 3);  % Number of trials
EEG_noise.srate = 512; % sampling rate
EEG_noise.times = linspace(-0.3, 0.5, 410);  % timepoints 
EEG_noise.xmin = EEG_noise.times(1);
EEG_noise.xmax = EEG_noise.times(end);
EEG_noise.chanlocs = elec_locs(1:128); 



EEG_grand_average = eeg_emptyset;

EEG_grand_average.data = grand_average;  
EEG_grand_average.nbchan = size(grand_average, 1);  % 128
EEG_grand_average.pnts = size(grand_average, 2);    % 410
EEG_grand_average.trials = 1;  % Number of trials
EEG_grand_average.srate = 512; % sampling rate
EEG_grand_average.times = linspace(-0.3, 0.5, 410);  % timepoints 
EEG_grand_average.xmin = EEG_grand_average.times(1);
EEG_grand_average.xmax = EEG_grand_average.times(end);
EEG_grand_average.chanlocs = elec_locs(1:128); 




EEG_signal = pop_rmbase( EEG_signal, [-0.3    0]);
EEG_noise = pop_rmbase( EEG_noise, [-0.3    0]);

EEG_grand_average = pop_rmbase( EEG_grand_average, [-0.3    0]);



pop_saveset(EEG_signal)
pop_saveset(EEG_noise)

%% 

%{
faceminusnoise = grand_avg_face - grand_avg_noise; 



h_mat = []
p_mat = []
sig = []

for t = 1:size(faceminusnoise, 2)   % time
    for c = 1:size(faceminusnoise, 1) % chans
        [h,p,ci,stats] = ttest(grand_avg_face(c,t), grand_avg_noise(c,t));  % H0: mean = 0, alpha 0.05
        p_mat(c, t) = p;

    end 
end


%}

%% 

% Data visualization 
% The x axis is seconds, -0.3 to 0.5 seconds
% creating a series such that plot exists on a -0.3 to 0.5 scale
series = linspace(-0.3, 0.5, 410);

% Data plotting 
plot(series, squeeze(trialMatrix(1, :, :)));
hold on; 
plot(series, grand_average(channelToUse,:)); 


% Misc figure labels
xlabel('Time (seconds)');
ylabel('Amplitude (microvolts)'); % its probably in microvolts but i dont want to make any assumptions

%% 


% Check trial information 

checktrial(23)

function [] = checktrial(trialNo)
    global times; 
    if trialNo > (size(times,2) * 10)
        disp('Out of range!')
        return
    end 

    if trialNo > 0 && trialNo < 10
        x = 1;
        y = trialNo +2;
    else
        x = floor(trialNo/10)+1; 
        y = mod(trialNo, 10) + 2; % first 2 numbers are participant & session

    end 

    % Syntax: times(row, column) 
    disp(['Trial',num2str(trialNo)]); 
    disp(['Participant ', num2str(times(1,x))]);    % Participant 
    disp(['session ', num2str(times(2,x))]);        % Ses number                                                
    disp(['noise snippet ', num2str(times(y,x))]);  % prestimulus 'snippet'



end 


 
 