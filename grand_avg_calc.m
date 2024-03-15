% data plotting 



Np = 1; % number of participants
%participantList = [1:Np];
 participantList = 1:4;
Ns = length(participantList);
sessionList = [1:10];
Ns = length(sessionList);

% EEG = pop_readbdf('ft10_p1s1.bdf', [],137, 1:128);


%eeg_folder = fullfile('..','data','ft10_eeg')
eeg_folder = 'ft10_eeg/';
%beh_folder = fullfile('..','data','ft10_behaviour')
beh_folder = 'ft10_behaviour/';
cedFile = fullfile('Glasgow_BioSemi_132.ced');
exptname = 'ft10';
ft10_eeg_trials = zeros(Np,Ns); % save number of trials
ft10_eeg_beh_match = zeros(Np,Ns);
filepathout = pwd; % where the epoch file is saved
grand_avg_face = [] 
grand_avg_noise = [] 
coef = []  % for regression methods
rejected = {}

elec_locs = readlocs('Glasgow_BioSemi_132.ced'); 

for P = participantList % for each participant
    
    for ses = sessionList % for each session
        

        % load data and set reference to the average of electrodes [1:128]
        id = [exptname,'_p',num2str(P),'s',num2str(ses)];
        bdffile = fullfile(eeg_folder,[id,'.bdf']);
        EEG = pop_readbdf(bdffile, [],137, 1:128);
        % EEG = pop_biosig(bdffile,'importannot','off');
        EEG = eeg_checkset( EEG );
        
        % remove electrodes [133:143]
        EEG = pop_select( EEG,'channel',1:132);
        % electrodes 1 to 128 are the standard Biosemi electrodes
        % electrodes 129 to 132 are extra-ocular channels
        % see readme file for details
        

        %}
        % filter data using ERPLAB
  %      EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Filter','bandpass','Design','butter',...
  %          'Cutoff', [1 30],'Order',4,'RemoveDC','on');
        
  
%         EEG = pop_runica(EEG);
        
        
        
        % epoch continuous data
        EEG = pop_epoch( EEG, {'11' '12' '21' '22' '31' '32' '41' '42' ...
            '51' '52' '61' '62' '71' '72' '81' '82' '91' '92' ...
            '101'  '102'  '111'  '112'}, [-1.5 1.5], ...    % -.3 to .5 -> in seconds 
            'newname', 'ft10 epoched file', 'epochinfo', 'yes');
            EEG = eeg_checkset( EEG );
        EEG.chanlocs = elec_locs; 



        % ones digit 
        % 1 = face
        % 2 = not face
        
        % use elec 1 (Cz) as reference
        % EEG = pop_reref(EEG,1);
    
        % in future, use REST or similar 
        % i cannot get it to work 
        % EEG_newnew = pop_REST_reref(EEG.data);


        %Remove DC from each epoch
        EEG.data=EEG.data-mean(EEG.data,2);
        
        % subtract mean of channel baseline from epoch
        EEG = pop_rmbase( EEG, [-300    0]);
        
        EEG = eeg_checkset( EEG );
        



        %{
        % remove eyeblink artifacts
        % linear regression approach
        % set of 1's = intercept

        EEG.data = EEG.data(1:128, :, :); 
        EOG = EEG.data(129:132, :, :); 

    
        for j = 1:size(EEG.data, 3)     % trialNo 
            for i = 1:size(EEG.data, 1) % shld be 128
                b = regress(squeeze(EEG.data(i, :, j)'), [ones(size(EOG, 2), 1), squeeze(EOG(:,:,j)')]); % cat set of 1 to EOG
                coef(:,i,j) = b(2:end);
            end 
        end 


        for j = 1:size(EEG.data, 3) 
            for i = 1:size(EEG.data, 1) 
                EEG.data(i, :, j) = EEG.data(i, :, j) - (coef(:, i, j)' * squeeze(EOG(:,:,j)));
            end
        end
        
        %}


        % load the channel location file


        % EEG = eeg_checkset( EEG );
        
        % remove events for potential extra triggers
        % happens once in a blue moon - doesn't hurt to check
        events_off = zeros(length(EEG.event),1);
        for ev = 1 : length(EEG.event)
            if EEG.event(ev).type >= 120
                events_off(ev) = 1;
            end
        end
        EEG = pop_editeventvals(EEG, 'delete', find(events_off) );
        %    EEG = eeg_checkset( EEG );
        
        ft10_eeg_trials(P,ses) = EEG.trials;
        
        % ********************************************
        % add behavioural information to EEG structure
        
        % load mat file containing behavioural results
        load(fullfile(beh_folder,id))
        correct = results.correct; % 1 for correct, 0 for incorrect response
        alternative = results.alternative; % face trial (1) or noise trial (2)
        rt = results.rt; % reaction times
        identity = results.identity; % face id [1:10]
        
        % -----------------------------------------------------------------
        % remove potential missing triggers - recordings crashed on 3 occasions
        % we left plenty of time after the crash and after resuming
        % recordings to avoid filtering distortions due to the
        % discontinuity. The indices below tell you which trials were not
        % recorded, so you could filter the segments before and after the,
        % separately if needed.
        if P == 2 && ses == 1
            gap = 978:993;
        elseif P == 3 && ses ==6
            gap = 1001:1012;
        elseif P == 3 && ses == 9
            gap = 692:699;
        else
            gap = [];
        end
        
        expected = 1:numel(results.rt);
        correct = correct(setdiff(expected,gap));
        alternative = alternative(setdiff(expected,gap));
        rt = rt(setdiff(expected,gap));
        identity = identity(setdiff(expected,gap));
        
        % check rejection
        EV = rem([EEG.event.type],10); % triggers recoded as 1-2
        if sum(alternative == EV) == numel(EV)
            fprintf('Participant %i, session %i: triggers & behaviour match...\n',P,ses)
            ft10_eeg_beh_match(P,ses) = 1;
        else
            fprintf('Participant %i, session %i: triggers & behaviour DO NOT match...\n',P,ses)
            ft10_eeg_beh_match(P,ses) = 0;
        end
   
        %EEG = pop_saveset( EEG,  'filename', filename, 'filepath', filepathout);
        


        
        %Remove Eyeblinks
        % ask if this or if the regression is better
        % modified such that it has a more lenient definition of eyeblinks
        % (average of 129-132) instead of just 129 
        % which changes depending on reference due to rank change

        disp('Removing eyeblinks...')
        [~, blinkWinIdx]=min(((EEG.times-[-500 300]').^2),[],2);
        blinkRange = blinkWinIdx(1):blinkWinIdx(2);
        % EOG_dat = squeeze(mean(EEG.data(129:132, :,:), 1)); 
        % isBlink = any(EOG_dat<-70);
        isBlink = any(squeeze(EEG.data(130,blinkRange,:))<-60);
        EEG     = pop_rejepoch(EEG,isBlink,0);
        nTrialPerSes(ses) = size(EEG.data,3);
        EEG = pop_select( EEG,'channel',1:128);



                % re-epoch continuous data
        EEG = pop_epoch( EEG, {'11' '12' '21' '22' '31' '32' '41' '42' ...
            '51' '52' '61' '62' '71' '72' '81' '82' '91' '92' ...
            '101'  '102'  '111'  '112'}, [-0.3 0.5], ...    % -.3 to .5 -> in seconds 
            'newname', 'ft10 epoched file', 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        
      
        % ****************************************************************

        try 
            % check ERPs without any clean up
            nev=rem([EEG.epoch.eventtype],10);
            X=EEG.times;
            c1Indices = find(nev==1); %Face
            c2Indices = find(nev==2); %noise
    
            nTrials2Keep = 300;  % takes 300 trials of each of face & noise
            c1Trials2Keep = randperm(length(c1Indices));
            c2Trials2Keep = randperm(length(c2Indices));
            
            c1Indices = c1Indices(c1Trials2Keep(1:nTrials2Keep));
            c2Indices = c2Indices(c2Trials2Keep(1:nTrials2Keep));
            EEG = pop_selectevent(EEG,'event',[  c1Indices c2Indices]);
           
            nev=rem([EEG.epoch.eventtype],10);
            X=EEG.times;
            c1Indices = find(nev==1); %Face
            c2Indices = find(nev==2); %noise
        catch 
            disp(['Skipping due to large number of rejected trials']); 
            rejected{end+1} = id;   % add id to 'rejected' 
            continue; 
        end 
        mferp = squeeze(mean(EEG.data(1:128,:,c1Indices),3));
        mnerp = squeeze(mean(EEG.data(1:128,:,c2Indices),3));
        % mdata = mferp - mnerp;

   
        % allMdata(:,:,P,ses) = mdata;            
        % allMFdata(:,:,P,ses) = mferp;           % guessing this is face
        % allMNdata(:,:,P,ses) = mnerp;           % guessing this is ..noise?


        grand_avg_face = cat(3, grand_avg_face, mferp);  %  seperate for face and noise 
        grand_avg_noise = cat(3, grand_avg_noise, mnerp); 

        varName = ['p', num2str(P), 's', num2str(ses)];
        disp(['Saved ', varName, ' with size: ', num2str(size(mferp))]);


    end 

end

%% 

grand_avg_face = mean(grand_avg_face, 3); 
grand_avg_noise = mean(grand_avg_noise, 3);


% output of this should be 2 dimensions, of 128*410
grand_average = mean(cat(3, grand_avg_face, grand_avg_noise), 3);
save('grand_average.mat', 'grand_average', 'grand_avg_face', 'grand_avg_noise');






%%


% Data visualization 
% The x axis is seconds, -0.3 to 0.5 seconds
n = 27
% creating a series such that plot exists on a -0.3 to 0.5 scale
series = linspace(-0.3, 0.5, 410)

% Data plotting 
plot(series, grand_avg_face(n, :), 'b');
hold on; 
plot(series, grand_avg_noise(n,:), 'r');


% Misc figure labels
xlabel('Time (seconds)');
ylabel('Amplitude (microvolts)');
legend('Face', 'Noise'); 



hold off; 






