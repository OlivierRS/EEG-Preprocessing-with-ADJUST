% Read and preprocess Raw EEG signals in EEGLAB format
% Notch filter at 60 Hz to remove power line artefact
% Band pass filtering 0.1 Hz to 50 Hz
% Channel removal
% ICA + Adjust

eeglab_folder = 'eeglab2021.1\';
addpath(eeglab_folder);

clear all;
close all;
clc;

%% Data directory
data_folder = '.\signals\';
raw_data_folder = strcat(data_folder, 'raw\');
prep_folder = strcat(data_folder, 'prep\');
try_mkdir(prep_folder);

% Data from subjects 2, 8 and 12 have been removed
subjects = [1, 3, 4, 5, 6, 7, 9, 10, 11];
sessions = 1:3;

n_subjects = length(subjects);
subjects_done = zeros(n_subjects, 1);

%% Configure preprocessing settings
param.raw_folder = raw_data_folder;
param.eeglab_folder = eeglab_folder;
param.rerun = false;
param.aarType = {'adjust'};

prep.enabled = true;
prep.in_folder = prep_folder;
prep.out_folder = prep_folder;
prep.notch = true;
prep.cleanLine = false;
prep.automaticChanReject = true;
param.prep = prep;

run_ica.enabled = true;
run_ica.in_folder = prep_folder;
run_ica.out_folder = prep_folder;
run_ica.removeEOG = true;
run_ica.icaType = 'runica'; % or 'amica'
param.ica_step = run_ica;

aar_clean.enabled = true;
aar_clean.in_folder = prep_folder;
aar_clean.out_folder = prep_folder;
param.aar_clean = aar_clean;

export_step.in_folder = prep_folder;
export_step.out_folder = prep_folder;
param.export_step = export_step;
param.should_epoch = true;
param.baseline_as_class = false; % _base ?
param.epoch_interval = [-1 15];
param.baseline_interval = [0 30];


%% Preprocess each EEG recording
for ix_sub=subjects
    for ix_ses=sessions
        sub_name = sprintf('s%02.f', ix_sub);
        ses_name = strcat(sub_name, '_', num2str(ix_ses));
        fprintf('Preprocessing subject %d session %d\n', ix_sub, ix_ses);
        try
            run_prep(param, ix_sub, ix_ses);
        catch
            fprintf('Error with %s\n', ses_name);
        end
    end
end


%% -------------------------------------------------------------------------------------
function run_prep(param, subject, session)

    sub_name = sprintf('s%02.f', subject);
    ses_name = strcat(sub_name, '_', num2str(session));
    class_names = {'ROT', 'WORD', 'SING', 'SUB', 'NAV', 'FACE', 'MI', 'trial', 'BASELINE'};
    
    %% EEGLAB analysis
    % See also:
    % http://sccn.ucsd.edu/wiki/Preprocess_pipeline
    % http://sccn.ucsd.edu/wiki/Quick_Rejection_Tutorial
    % Move to EEGLAB folder and add relevant paths
    % eeglabPath = 'C:\Users\Hubert\Documents\MATLAB\eeglab13_4_4b';
    % cd(eeglabPath);
    % addpath(genpath(eeglabPath))
    % rmpath(genpath('C:\Users\Hubert\Documents\MATLAB\eeglab13_3_2b\plugins\Biosig2.88\biosig\maybe-missing')); % This folder has functions that clash with MATLAB in-built functions.
    % rmpath(genpath('C:\Users\Hubert\Documents\MATLAB\eeglab13_3_2b\functions\octavefunc'));

    filename1 = sprintf('%s_step_1_filt.set', ses_name);
    filename2 = sprintf('%s_step_2_ica.set', ses_name);
    filename3 = sprintf('%s_step_3_aar_%s.set', ses_name, strjoin(param.aarType, '-'));
    filename4 = ses_name;
    
    EOG_channels = {'EXG1', 'EXG2', 'EXG5', 'EXG6'};

    %% PREPROCESSING
    file_exist = isfile(fullfile(param.prep.out_folder, filename1));
    
    if param.prep.enabled && (~file_exist || param.rerun)
        % Load dataset and rereference with Cz (ch48) or Mastoids (ch65 or ch66)
        % (pop_biosig())
        % Rereferencing and keeping the channel: http://sccn.ucsd.edu/pipermail/eeglablist/2008/002310.html
        % Biosemi referencing: http://www.biosemi.com/faq/cms&drl.htm
        % Adding back a reference channel in the data
        refChan = 48;
        
        raw_filename = strcat(param.raw_folder, ses_name, '.bdf');
        EEG = pop_biosig(raw_filename, 'ref', refChan);
        EEG = eeg_checkset(EEG);
        
        if subject == 5 && session == 1
            EEG.event = EEG.event(31:end);
        end
        
        if subject == 8 && ismember(session, [1, 3])
            EEG.event = EEG.event(3:end);
        end
        
        if subject == 9 && session == 2
            EEG.event = EEG.event(2:end);
        end
        
        if subject == 3 && session == 1
            EEG.event = EEG.event(7:end);
        end
        
        lat = extractfield(EEG.event, 'latency');
        lat = lat/EEG.srate;
        y = extractfield(EEG.event, 'edftype');
        
        %% Markers visualization
        %{
        figure;bar(lat, y);
        ylim([0, 10]);
        title(raw_filename);
        %}
        
        %%
        % Each recording session contains four sub-sessions
        [begin_index, end_index] = find_subsession_bounds_idx(EEG);
        
        y(begin_index) = 9;
        y(end_index) = 9;
        
        figure;bar(lat, y);
        ylim([0, 11]);
        title(raw_filename);
        
        y = num2cell(y);
        [EEG.event(:).edftype] = y{:};
        
        % Update dataset info
        % 1 - Subject number
        % 2 - Session number
        % 3 - Channel locations

        EEG.setname = ses_name;
        EEG.subject = sub_name;
        EEG.session = session;
        
        fname_ch_info = 'plugins\dipfit4.3\standard_BESA\standard-10-5-cap385.elp';
        fullname_ch_info = strcat(param.eeglab_path, fname_ch_info);
        EEG = pop_chanedit(EEG, 'lookup', fullname_ch_info);
        EEG = eeg_checkset(EEG);
        EEG.event = rmfield(EEG.event, 'type');
        
        [EEG.event(:).type] = EEG.event(:).edftype;
        
        % Make sure the total number of events is right
        events = zeros(length(EEG.event),1);
        latencies = zeros(length(EEG.event),1);
        for i = 1:length(EEG.event)
            if ischar(EEG.event(1,i).type) % For some reason sometimes the event types are char...
                 value = regexp(EEG.event(1,i).type,'[0-9]','match');
                 events(i) = str2num(value{1});
            else
                events(i) = EEG.event(1,i).type;
            end
            latencies(i) = EEG.event(1,i).latency;
        end
        
        %% Visual consistency check
        for i = 1:7
            nbTrials = sum(events == i);
            if nbTrials ~= 16 % There is supposed to be 16 trials of each task.
                figure(10)
                hist(events,8)
                title('Number of trials per task (there should be 16 of each)')
                xlabel('Task number')
                ylabel('Number of trials')
                warning('preprocess:wrongtrialnb',['There is ',num2str(nbTrials), ' trials of task ', num2str(i), ' instead of 16.']);
            end
        end

        % Rename events for the baselines at the beginning and at the end of
        % each subsession
        %{
        % Magic numbers = dangerous errors
        longBaselineInd = [1,58,59,116,117,174,175,232];
        for i = 1:length(longBaselineInd)
            try
                EEG = pop_editeventvals(EEG, 'changefield', {longBaselineInd(i) 'type' 0});
            catch
                warning(['Long baseline index ',num2str(longBaselineInd(i)),' could not be reached...']);
            end
        end
        %}
        EEG = eeg_checkset(EEG);

        % Rename the EOG channels and give them new locations
    %     for i = 1:EEG.nbchan % Find EXG3-EXG6 (EOG + nose channels)
    %         if strcmp(EEG.chanlocs(i).labels,'EXG3')
    %             eogCh = i:i+3;
    %         end
    %     end
    %     eogChLoc = [68.69, 49.71, -35; 68.72, -49.67, -35; 75, 33, -45; 89, 0, -45];
    %     for i = 1:length(eogCh)
    %         EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'changefield', {eogCh(i) 'labels' eogChNames{i}});
    %         EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'changefield', {eogCh(i) 'X' eogChLoc(i,1)});
    %         EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'changefield', {eogCh(i) 'Y' eogChLoc(i,2)});
    %         EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'changefield', {eogCh(i) 'Z' eogChLoc(i,3)});
    %     end
        % Convert from (X,Y,Z) to polar coordinates
        EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'convert', 'cart2all');

        % Remove unused channels (AF7, AF8, MastoidL, MastoidR, EXG7, EXG8, Ana1-8)
        % unusedCh = [2, 35, 65, 66, 71:87];
        
        unusedChNames = {'AF7','AF8','EXG7','EXG8','GSR1','GSR2',...
            'Erg1','Erg2','Resp','Plet','Temp','Ana1','Ana2','Ana3','Ana4',...
            'Ana5','Ana6','Ana7','Ana8','EXG3','EXG4'};
        EEG = pop_select(EEG, 'nochannel', unusedChNames);
        EEG = eeg_checkset(EEG);

        % Save the chanlocs structure
        EEG.goodChChanlocs = EEG.chanlocs;
        
        % Notch filter around 60 Hz
        if param.prep.notch
            EEG = pop_eegfiltnew(EEG, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1);
        end
        
        % Highpass at 0.1 Hz (pop_eegfiltnew())
        EEG = pop_eegfiltnew(EEG, [], 0.1, [], true, [], false);
        
        % NOTE
        % Resample to 256 Hz (pop_resample())
        %EEG = pop_resample(EEG, 256);
        

        % Remove line noise using cleanline()
        if param.prep.cleanLine
            EEG = cleanline(EEG, 'Bandwidth',2,                                               ...
                             'SignalType','Channels','ComputeSpectralPower',false,            ... %'ChanCompIndices',[1:EEG.nbchan],
                             'LineFrequencies',[60 120] ,'NormalizeSpectrum',false,           ...
                             'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,          ...
                             'ScanForLines',true,'SmoothingFactor',100,'VerboseOutput',1,     ...
                             'SlidingWinLength',4,'SlidingWinStep',1); %'SlidingWinLength',EEG.pnts/EEG.srate,'SlidingWinStep',EEG.pnts/EEG.srate);
        end
        

        % Lowpass at 100 Hz
        %EEG = pop_eegfiltnew(EEG, [], 100, [], false, [], false);

        % Use trimOutliers() to remove paroxysmal artifacts
        % 1- Reject channels based on their STD
        % 2- Reject data points based on their amplitude
        % EEG = pop_trimOutlier(EEG);

        % Ask for the visually identified bad channels using a dialog window
    %     prompt = {'Enter channel(s) name'};
    %     dlg_title = 'Manual rejection of bad channels';
    %     num_lines = [1, 50];
    %     def = {'Elec1, elec2, ...'};
    %     badChanNames = inputdlg(prompt,dlg_title,num_lines,def,'on');
    %     % Parse the output of the prompt
    %     badChanNames = strtrim(strsplit(badChanNames{1},','));
    
        % These channels are removed based on the script
        % 'eeg_preprocessing_telequebec'.
        badChanNames = {'P8','O1'};
        EEG = pop_select(EEG, 'nochannel', badChanNames);
        EEG = eeg_checkset(EEG);
        
    %     [EEG] = findAlphaContamination(EEG, 0.9, 98, 9439, 9825);
    %     EEG = pop_saveset( EEG, 'filename', [sessionNbStr,'_resampled_hp_lp_epoched_alpha.set'], 'filepath', savepath);

        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG.rejmanual = EEG.reject.rejmanual;
        EEG = pop_rejepoch( EEG, EEG.reject.rejmanual ,0);

        EEG.setname=[ses_name,' step 1: filt + resampled'];
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename', filename1, 'filepath', param.prep.out_folder);
    end
    %}
    

    %% COMPUTING ICA
    file_exist = isfile(fullfile(param.ica_step.out_folder, filename2));
    
    if param.ica_step.enabled && (~file_exist || param.rerun)
        if ~exist('EEG','var')
            % Load the preprocessed dataset for ICA decomposition...
            EEG = pop_loadset('filename', filename1, 'filepath', param.ica_step.in_folder); 
        end
        
        aarType = param.aarType;
        
        % Systematically remove Cz as this channel is used as reference and ICA
        % does not account for it in the decomposition
        % If the number of channels if not the same as the number of ICs
        % (the case when Cz was kept although it's used as reference),
        % remove Cz
        EEG = pop_select(EEG, 'nochannel', {'Cz'});
        EEG = eeg_checkset(EEG);

        % Run ICA on the good channels
        tic
        icaType = param.ica_step.icaType;
        if strcmp(icaType,'runica')
            EEG = pop_runica(EEG,'icatype', 'fastica', 'extended',1,'interupt','on');
        elseif strcmp(icaType,'amica')
            [EEG.icaweights, EEG.icasphere, ~] = runamica12(EEG.data(:,:));
        end
        elapsedTime = toc;
        disp(['The ICA decomposition took ', num2str(elapsedTime/60), ' min.']);

        EEG = eeg_checkset(EEG);
        EEG.setname=[ses_name,' step 2: ICA'];
        EEG = pop_saveset(EEG, 'filename', filename2, 'filepath', param.ica_step.out_folder);
    end
    

    %% ARTEFACT CORRECTION
    file_exist = isfile(fullfile(param.aar_clean.out_folder, filename3));
    
    if param.aar_clean.enabled && (~file_exist || param.rerun)
        if ~exist('EEG','var')
            % Load the preprocessed dataset with the computed ICA for artefact correction
            EEG = pop_loadset('filename', filename2, 'filepath', param.aar_clean.in_folder);
            % Make sure the ICA was computed
            assert(isfield(EEG,'icaweights'), 'The ICA was not found in the EEG dataset.');
        end
        
        suffix = '';
        aarType = param.aarType;
        savepath_report = param.aar_clean.out_folder;
        
        if any(ismember(aarType,'adjust'))
            event = EEG.event;
            % Run ADJUST and remove the unwanted ICs
            %EEG_ADJUST = interface_ADJ(EEG, [savepath, filesep, 'ADJreport_',ses_name]);
            %EEG_ADJUST = pop_subcomp(EEG_ADJUST);
            [EEG_ADJUST, art] = interface_ADJ_mod(EEG, [savepath_report, filesep, 'ADJreport_',ses_name]);
             EEG_ADJUST = pop_subcomp(EEG_ADJUST, art);
    %         [EEG, art] = interface_ADJ_HJB(EEG, [savepath, filesep, ...
    %         'ADJreport_',ses_name]); % Warning: function not found
    %         EEG = pop_subcomp(EEG, art);

            EEG_ADJUST = eeg_checkset(EEG_ADJUST);
            EEG_ADJUST.setname=[ses_name,' resampled hp lp epoched ica ADJUST_corrected'];
            suffix = strcat(suffix, '_ADJUST');
            EEG_ADJUST = pop_saveset( EEG_ADJUST, 'filename', [ses_name, suffix, '.set'], 'filepath', param.aar_clean.out_folder);
            EEG_ADJUST.event = event;
            vis_artifacts(EEG_ADJUST, EEG);
            EEG = EEG_ADJUST;
        end

        
        % Reject artifacted epochs (pop_autoref())
        disp('You should reject artifacted epochs here... As much as 10% according to http://sccn.ucsd.edu/wiki/Chapter_09:_Decomposing_Data_Using_ICA');

        % Re-interpolate bad channels from neighbouring channels
        %   Interpolate the channels that were deleted (pop_interp())
        EEG = pop_interp(EEG, EEG.goodChChanlocs, 'spherical');
        EEG.setname=[ses_name, suffix];
        %TODO uncomment
        %EEG = pop_saveset( EEG, 'filename', [ses_name,'_aar_corrected.set'], 'filepath', param.aar_clean.out_folder);

        % Remove EOG channels before re-referencing
        EEG = pop_select(EEG, 'nochannel', EOG_channels);
        EEG = eeg_checkset(EEG);
        % Re-reference to common average (pop_reref()) and keep the previous
        % reference
        EEG = pop_reref(EEG, [], 'keepref', 'on');

        % SPATIAL FILTERING instead of CAR?
        % Current Source Density interpolation toolbox
        % (http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/index.html)
        % OR using this function http://sccn.ucsd.edu/wiki/Flt_laplace
        %% -- NEW -- Apply Laplacian filter
        if any(ismember(aarType,{'lp', 'CSD', 'Laplacian'}))
            EEG_CSD = EEGLAB_CSD(EEG, EEG.chanlocs);
            label = 1;
            str_title = sprintf('EEG - %s', class_names{label});
            topo_class_average(EEG, {label}, param.epoch_interval, str_title);
            str_title = sprintf('EEG CSD - %s', class_names{label});
            topo_class_average(EEG_CSD, {label}, param.epoch_interval, str_title);
            vis_artifacts(EEG_CSD, EEG);
            
            suffix = strcat(suffix, '_CSD');
            EEG_CSD.setname=[ses_name, suffix];
            EEG_CSD = pop_saveset(EEG_CSD, 'filepath', param.aar_clean.out_folder, 'filename', [ses_name, suffix, '.set']);
            EEG = EEG_CSD;
        end
        
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset( EEG, 'filename', filename3, 'filepath', param.aar_clean.out_folder);
    end
    
    %% --- Final preprocessing step ---
    if ~exist('EEG','var')
        % Load the corrected dataset
        EEG = pop_loadset('filename', filename3, 'filepath', param.export_step.in_folder);
        % Make sure the ICA was computed
        assert(isfield(EEG,'icaweights'), 'The ICA was not found in the EEG dataset.');
    end
    
    %% Remove unused channels
    unusedChNames = {'EXG3','EXG4','EXG5','EXG6'};
    unusedChNames = cat(2, unusedChNames, EOG_channels);
    EEG = pop_select(EEG, 'nochannel', unusedChNames);
    
    EEG = pop_resample(EEG, 100);
    plot_event(EEG);
    
    % Slice into epochs without removing baseline values (pop_epoch())
    % I take all the mental tasks, and use the interval [-5 20] or [-0.3 15]
    % If using more data before the trigger, chances are blinks and
    % movements will be picked up more often
    if param.baseline_as_class
        eventsOfInterest = {'9'};
        interval = param.baseline_interval;
    else
        eventsOfInterest = {'1','2','3','4','5','6','7'};
        interval = param.epoch_interval;
    end
    
    EEG.setname = ses_name;
    
    %% Keep trials and baseline signal segments exclusively
    
    % Marging to add before and after subsessions to prevent filtering
    % artefacts and to keep end of subsessions baselines
    pre_subses_margin = 10;
    post_subses_margin = 30+10;
    
    latency = extractfield(EEG.event, 'latency');
    latency = latency/EEG.srate; % convert latency in seconds
    [begin_index, end_index] = find_subsession_bounds_idx(EEG);
    begin_latency = latency(begin_index)-pre_subses_margin;
    end_latency = latency(end_index)+post_subses_margin;
    
    subses_intervals = cat(1, begin_latency, end_latency);
    subses_intervals = subses_intervals';
    EEG = pop_select(EEG, 'time', subses_intervals);

    % Select valid event only
    y = extractfield(EEG.event, 'type');
    EEG.event = EEG.event(ismember(y, eventsOfInterest));
    
    % Add class names
    y = extractfield(EEG.event, 'type');
    y = str2double(y);
    y_names = class_names(y);
    [EEG.event(:).task] = deal(y_names{:});
    [EEG.event(:).type] = deal(y_names{:});
    
    % Epoching?
    if param.should_epoch
        EEG = pop_epoch(EEG, [], interval);
    end
    
    %% Save EEGLAB dataset
    pop_saveset( EEG, 'filename', filename4, 'filepath', param.export_step.out_folder);
        
end

function [begin_subses_index, end_subses_index] = find_subsession_bounds_idx(EEG)
	gap_min_duration = 60; % in seconds
    
    latency = extractfield(EEG.event, 'latency');
    latency = latency/EEG.srate;
    
    % First event is necessarely the beginning of the first subsession
    lat_diff_begins = cat(2, -gap_min_duration, latency);
    lat_diff_ends = cat(2, latency, EEG.xmax+gap_min_duration);
    % subsessions are separated with more than 30s from each others
    lat_diff_begins = abs(diff(lat_diff_begins));
    lat_diff_ends = abs(diff(lat_diff_ends));
    % Find gaps that longer than 1 minute
    begin_subses_mask = lat_diff_begins > gap_min_duration;
    end_subses_mask = lat_diff_ends > gap_min_duration;
    
    % Event index arrays of beginning and end of each subsessions
    begin_subses_index = find(begin_subses_mask);
    end_subses_index = find(end_subses_mask);
    
    
    %% Verify validity of subsession begin
    % End of gap = Begin of subsession
    % Long baseline are last approx. 30s after each subsessions beginning
    begin_baseline = latency(begin_subses_index);
    end_baseline = begin_baseline+30;
    
    % End of baseline should be align with its corresponding end of
    % baseline event
    mask2 = zeros(1, length(begin_subses_mask));
    for ix=1:length(end_baseline)
        % Check if closest event from end is in acceptable variability
        % window (empirically 5s)
        mask2(begin_subses_index(ix)) = min(abs(end_baseline(ix)-latency)) < 5;
    end
    begin_subses_mask = begin_subses_mask & mask2;
    
    
    %% Verify validity of subsession end
    % Beginning of gap = End of subsession
    % Each end of subsession should be at least 'gap_min_duration' min away from a beginning
    begin_latency = latency(begin_subses_mask);
    mask2 = zeros(1, length(end_subses_mask));
    for ix=1:length(end_subses_index)
        lat = latency(end_subses_index(ix));
        mask2(end_subses_index(ix)) = min(abs(lat - begin_latency)) > gap_min_duration;
    end
    end_subses_mask = end_subses_mask & mask2;

    begin_subses_index = find(begin_subses_mask);
    end_subses_index = find(end_subses_mask);
end


function plot_event(EEG)
    labels = extractfield(EEG.event, 'type');
    mask = ismember(labels, 'boundary');
    
    lat = extractfield(EEG.event, 'latency');
    lat = lat/EEG.srate;
    lat = lat(~mask);
    y = extractfield(EEG.event, 'edftype');
    
    figure;bar(lat, y);
    ylim([0, 11]);
end

function topo_class_average(EEG, target_events, interval, str_title)
    EEG_epoched = pop_epoch(EEG, target_events, interval);
    
    begin_epoch = interval(1);
    end_epoch = interval(2)-1;
    timestamp = linspace(begin_epoch, end_epoch, end_epoch-begin_epoch+1)';
    timestamp = timestamp * 1000;
    pop_topoplot(EEG_epoched, 1, timestamp);
    title(str_title);
end



