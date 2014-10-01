function WEEG_MakeEEGLABDataset(Participant, SaveDir)
% function WEEG_MakeEEGLABDataset(Participant, SaveDir)
% Creates EEGLAB datasets for all series of the supplied participant 
% These datasets are notated with event times, in the EEGLAB format
% such as time of liftoff, LEDon, etc.
% This file also generates and adds "special" events
% (denoted special since they do not occur in every lift)
% regarding unexpected and expected weight and surface 
% on first touch events (for surface) and liftoff events (for weight)
%
% To run this script, EEGLAB must be on your MATLAB path
%
% INPUTS
%   Participant - Participant number
%   SaveDir     - Saving directory for EEGLAB dataset (default=current dir)
%  
%   for placing electrodes on the scalp map,
%   chanlabels_32channel.xyz must be in this directory
%   (you can also change the directory location below)
%
%   how to load into EEGLAB: this code creates a .set file. 
%   Start EEGLAB, then File->Load Existing Dataset... and select the .set file
%
% This MATLAB script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava



if nargin < 2
    SaveDir = '.';
end

ChanLocDir = '';    %directory containing channel position files %
% Get the channel location names
    try
        Chanlocs = readlocs([ChanLocDir 'chanlabels_32channel.xyz'], 'filetype', 'xyz');  %using .xyz format so it fits in the scalp map
    catch
        error('ERROR: EEGLab does not seem to be on the path...');
    end
   
%"Special" event names to add
WeightEventNames = [{'Unexp. heavy weight lift'}, ...
                    {'Unexp. light weight lift'}, ...
                    {'Expected weight lift'}];

SurfEventNames =  [{'Unexp. higher friction touch'}, ...
                   {'Unexp. lower friction touch'}, ...
                   {'Expected friction touch'}];

               
%% Load Events 
fprintf('\nLoading Data and Events ...\n');
AllSeriesInfo = WEEG_GetEventsInHS(Participant, 4);

if ((strcmp(computer,'PCWIN')) || (strcmp(computer, 'PCWIN64')))
    dirLoc = ['..\P', int2str(Participant), '\'];
else
    dirLoc = ['../P', int2str(Participant), '/'];
end
   
try
    load([dirLoc, 'P', int2str(Participant), '_AllLifts.mat'])
catch
    error('ERROR: PX_AllLifts should be in directory ../PX');
end

for Series = 1:9
    Events = []; 
    EEG = [];
    
    %% Load the EEG data
    fprintf(['\nLoading Series ' num2str(Series) '\n']);
    load([dirLoc, 'HS_P', num2str(Participant), '_S', num2str(Series), '.mat']);
        
    %% Create the dataset and lookup channel locations
        EEG = pop_importdata(...
        'setname', ['WAY-GAL Participant ' num2str(Participant) ' Series ' num2str(Series)], ...
        'data', 0.1*hs.eeg.sig', ...
        'dataformat', 'array', ...
        'subject', hs.name, ...
        'nbchan', size(hs.eeg.names, 2), ...
        'chanlocs', Chanlocs, ...
        'srate', hs.eeg.samplingrate);
    
    %% Add events
    % Get all events and names for this series
    SeriesEvents = AllSeriesInfo.Events(Series);

    EventNames = fieldnames(SeriesEvents);
    
    for EventNumber = 1:size(EventNames, 1)
        % For each non-special event name(type), make the Events cell array and import
        EventName = char(EventNames(EventNumber));
        EventTimes = SeriesEvents.(EventName);
        Events = [repmat({EventName}, size(EventTimes, 1), 1) num2cell(EventTimes)];
        EEG = pop_importevent(EEG, 'event', Events, 'fields', {'type', 'latency'}, 'append', 'yes');
    end
    
    %special events (expected/unexpected weight/friction)
    %AllSeriesInfo(Series)
    SeriesType = AllSeriesInfo.SeriesType(Series);  %type of series   
    
    %get rows in P.AllLifts for this series
    Lifts = find(P.AllLifts(:,2)==Series); 
    
    if (SeriesType == 1)  %it is a weight series, let's look at weight
        CurrWeights = P.AllLifts(Lifts,4);
        PrevWeights = P.AllLifts(Lifts,6);
    
        UxHeavy = (find([CurrWeights > PrevWeights]==1));  %unexpected heavy
        UxLight = (find([CurrWeights < PrevWeights]==1));  %unexpected light
        Expected = (find([CurrWeights == PrevWeights]==1));  %same as before 
    
        LiftoffEvents = SeriesEvents.tLiftOff;  %times of liftoff 
        SpecialEvents(1).times = LiftoffEvents(UxHeavy);
        SpecialEvents(2).times = LiftoffEvents(UxLight);
        SpecialEvents(3).times = LiftoffEvents(Expected);
    
        %import
        for EventNumber = 1:3
            % For each event name(type), make the Events cell array and import
            EventTimes = SpecialEvents(EventNumber).times; 
            Events = [repmat({WeightEventNames{EventNumber}}, ...
                size(EventTimes, 1), 1) num2cell(EventTimes)];
            EEG = pop_importevent(EEG, 'event', Events, ...
                'fields', {'type', 'latency'}, 'append', 'yes');
        end
    
    else   %it is a surface/mixed series, let's look at friction
        CurrSurf = P.AllLifts(Lifts,5);
        PrevSurf = P.AllLifts(Lifts,7);
    
        UxHigher = (find([CurrSurf > PrevSurf]==1));   %unexpected higher friction
        UxLower =  (find([CurrSurf < PrevSurf]==1));  %unexpected lower friction
        Expected = (find([CurrSurf == PrevSurf]==1));  %same as before 
    
        TouchEvents = SeriesEvents.tTouch;  %times of touch 

        SpecialEvents(1).times = TouchEvents(UxHigher);
        SpecialEvents(2).times = TouchEvents(UxLower);
        SpecialEvents(3).times = TouchEvents(Expected);
        
        %import 
        for EventNumber = 1:3
            % For each event name(type), make the Events cell array and import
            EventTimes = SpecialEvents(EventNumber).times; 
            Events = [repmat({SurfEventNames{EventNumber}}, ...
                     size(EventTimes, 1), 1) num2cell(EventTimes)];
            EEG = pop_importevent(EEG, 'event', Events, 'fields', {'type', 'latency'}, 'append', 'yes');
        end
    
    end
    
      
    %% Save and exit
    FileName = [[SaveDir '/HS_P'], num2str(Participant), '_S', num2str(Series), '.set'];
    pop_saveset(EEG, ...
        'filename', FileName, ...
        'check', 'off', ...
        'version', '7.3');

    fprintf('\nSaved EEGLAB dataset %s\n', FileName);
end
return