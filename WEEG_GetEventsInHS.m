function [Runs] = WEEG_GetEventsInHS(participant, seriesType)
% function Runs = WEEG_GetEventsInHS(participant, seriesType)
% Identifies the time of events in HS*.mat files of a specific type, e.g., 
% 'Weight series'
%
% INPUTS
%   participant    - Participant number (default = 1)
%   seriesType - 1=Weight series (default), 2=Friction series, 
%                3=Mixed series, 4=All series            
%
% OUTPUT
%   Runs with the following fields: 
%     HSs()     - runs (from HS*.mat files)
%     Events().   tLEDOn         - time of LED on
%                 tLEDOff        - time of LED off
%                 tHandStart     - time of HandStart
%                 tObjectContact - time of first digit contact
%                 tLoadPhaseOn   - time when load phase began on both digits
%                 tLiftOff       - time of liftoff of the object
%                 tReplace       - time of replace of the object
%                 tRelease       - time when both digits have released object
%     States().   WeightExpect   - comparing current weight with previous: 
%                                  0 = equal; +1 = larger; -1 = smaller
%                 SurfaceExpect  - comparing current friction with previous: 
%                                  0 = equal; +1 = larger; -1 = smaller
%     SeriesType  - 1: Weight, 2 Friction, 3: Mixed
%
% (Two code examples are provided at the end of this m file)
%
% This MATLAB / Octave script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava


if nargin < 1, participant = 1; disp('WARNING: default Participant=1'); end
if nargin < 2, seriesType = 1; disp('WARNING: default seriesType=1 (weight series)'); end; 

if ((strcmp(computer,'PCWIN')) || (strcmp(computer, 'PCWIN64')))
    dirLoc = ['..\P', int2str(participant), '\'];
else
    dirLoc = ['../P', int2str(participant), '/'];
end
    
try
    load([dirLoc, 'P', int2str(participant), '_AllLifts.mat'])
catch
    error('ERROR: PX_AllLifts should be in directory ../PX');
end

cFile    = LOC_GetVarColumns('Run', P.ColNames); 
cCurW    = LOC_GetVarColumns('CurW', P.ColNames); 
cCurS    = LOC_GetVarColumns('CurS', P.ColNames); 
cPrevW   = LOC_GetVarColumns('PrevW', P.ColNames); 
cPrevS   = LOC_GetVarColumns('PrevS', P.ColNames); 
cSTs     = LOC_GetVarColumns('StartTime', P.ColNames); 

c_seriesType = LOC_GetVarColumns('BlockType', P.ColNames); 
if seriesType < 4,  q_seriesType = P.AllLifts(:,c_seriesType) == seriesType;
else q_seriesType = P.AllLifts(:,c_seriesType) < seriesType;
end

Files = unique(P.AllLifts(q_seriesType,cFile));

cLEDOn = LOC_GetVarColumns('LEDOn', P.ColNames); 
cLEDOff = LOC_GetVarColumns('LEDOff', P.ColNames); 
ctFirstDigitTouch = LOC_GetVarColumns('tFirstDigitTouch', P.ColNames); 
ctBothStartLoadPhase = LOC_GetVarColumns('tBothStartLoadPhase', P.ColNames); 
ctLiftOff = LOC_GetVarColumns('tLiftOff', P.ColNames); 
ctReplace = LOC_GetVarColumns('tReplace', P.ColNames); 
ctBothReleased = LOC_GetVarColumns('tBothReleased', P.ColNames); 
ctHandStart = LOC_GetVarColumns('tHandStart', P.ColNames); 
cCurWeight = LOC_GetVarColumns('CurW', P.ColNames); 
cPrevWeight = LOC_GetVarColumns('PrevW', P.ColNames); 


for File = 1 : numel(Files)   
    fname = [dirLoc, 'HS_P', int2str(participant), '_S', int2str(Files(File)), '.mat'];
    fprintf(1, '%s\n', fname);
    load(fname);
    Runs.HSs(File) = hs;
    
    %note the type of series
    temp = find(P.AllLifts(:,cFile) == Files(File));
    Runs.SeriesType(File) = P.AllLifts(temp(1),c_seriesType);
    
    % Find the rows in P.AllLifts
    if (seriesType < 4)
        q = find( (P.AllLifts(:,c_seriesType) == seriesType) & ...
              (P.AllLifts(:,cFile) == Files(File)));
    else
        q = find( (P.AllLifts(:,c_seriesType) < seriesType) & ...
              (P.AllLifts(:,cFile) == Files(File))); 
    end
    
    %% Events
    STs = P.AllLifts(q,cSTs); % Offset time of each single lift trial     
    Runs.Events(File).LEDOn = STs + P.AllLifts(q,cLEDOn);
    Runs.Events(File).LEDOff = STs + P.AllLifts(q,cLEDOff);
    Runs.Events(File).tTouch = STs + P.AllLifts(q,ctFirstDigitTouch);
    Runs.Events(File).tStartLoadPhase = STs + P.AllLifts(q,ctBothStartLoadPhase);
    Runs.Events(File).tLiftOff = STs + P.AllLifts(q,ctLiftOff);  %-P.AllLifts(q,cLEDOn);
    Runs.Events(File).tReplace = STs + P.AllLifts(q,ctReplace);
    Runs.Events(File).tRelease = STs + P.AllLifts(q,ctBothReleased);
    Runs.Events(File).tHandStart = STs + P.AllLifts(q,ctHandStart);
    %% Lift off with an object heavier then expected
    %q = find(P.AllLifts(:,cCurWeight)>P.AllLifts(:,cPrevWeight) & ...
    %         P.AllLifts(:,cFile)==Files(File));
    %if ~isempty(q) 
    %  STs = P.AllLifts(q,cSTs);
    %  Runs.Events(File).tUnExpectedHeavy = ...
    %        STs + P.AllLifts(q, ctBothStartLoadPhase);
    %end
    %q = find(P.AllLifts(:,cCurWeight) == P.AllLifts(:,cPrevWeight) & ...
    %         P.AllLifts(:,cFile)==Files(File));
    %if ~isempty(q) 
    %  STs = P.AllLifts(q,cSTs);
    %  Runs.Events(File).tExpectedWeight = ...
    %        STs + P.AllLifts(q, ctBothStartLoadPhase);
    %end

    %% States    
    % CurW : 1=165g, 2=330g, 3=660g
    % CurS : 1=silk, 2=suede, 3=sandpaper
    Runs.States(File).WeightExpect = sign(P.AllLifts(q,cCurW) - ...
                                  P.AllLifts(q,cPrevW));
    Runs.States(File).SurfaceExpect = - sign(P.AllLifts(q,cCurS) - ...
                                   P.AllLifts(q,cPrevS));
end

return

%% LOCAL utility functions
function [Col] = LOC_GetVarColumns(P, ColNames)
% function [Col] = LOC_GetVarColumns(P, ColNames)
%   Returns the index of ColNames in P

Col = -1;
if numel(P) == 1 % Mistaken input; structure as second input?
    Temp = P; P = ColNames; ColNames = Temp;
end
for C = 1:numel(ColNames)
    if strcmp(P, ColNames(C)), Col = C; break, end
end

return

%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Example usage

%% 1.
% Get all weight series for Participant #1 
% and mark moment of liftoff on the vertical object position

Runs = WEEG_GetEventsInHS(1,1);  % Subj 1 and weight block
close all, figure
for R = 1:numel(Runs.HSs); 
    subplot(6, 1, R), 
    t = 0.002*(1:size(Runs.HSs(R).kin.sig, 1));
    plot(t, Runs.HSs(R).kin.sig(:,27)); hold on; % 27 = Vertical object pos
    tP = round(500*Runs.Events(1,R).tLiftOff);
    plot(t(tP), Runs.HSs(R).kin.sig(round(tP),27), '.r', 'MarkerSize', 10);
end

%% 2. 
% Get all series for Participant #1 
% - Mark lift force with a red dot if heavier weight than previous and with
%   a green dot if lighter than previous 
% - Mark grip force with a red circe if slippier than previous and with a 
%   green circle if less slippery than previous
% (The GF and LF changes between the lifts are due to the experimenter
% manipulating the surface plates)

Runs = WEEG_GetEventsInHS(1,4);  % Subj 1 and all series
close all, figure
for R = 1: min([6 numel(Runs.HSs)]); 
    subplot(6, 1, R), 
    t = 0.002*(1:size(Runs.HSs(R).kin.sig, 1));
    GF = -0.5*(Runs.HSs(R).kin.sig(:,17) + Runs.HSs(R).kin.sig(:,18));
    LF = Runs.HSs(R).kin.sig(:,13) + Runs.HSs(R).kin.sig(:,14);
    plot(t, LF, t, GF); hold on; 
    legend('LF', 'GF');
    tP = Runs.Events(1,R).tLiftOff;
    q = find(Runs.States(1,R).WeightExpect == 1);
    plot(tP(q), LF(round(500*tP(q))), '.r', 'MarkerSize', 10);
    q = find(Runs.States(1,R).WeightExpect == -1);
    plot(tP(q), LF(round(500*tP(q))), '.g', 'MarkerSize', 10);
    q = find(Runs.States(1,R).SurfaceExpect == 1);
    plot(tP(q), LF(round(500*tP(q))), 'or', 'MarkerSize', 10);
    q = find(Runs.States(1,R).SurfaceExpect == -1);
    plot(tP(q), LF(round(500*tP(q))), 'og', 'MarkerSize', 10);
end
