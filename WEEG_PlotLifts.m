function WEEG_PlotLifts(participants, series)
%% function WEEG_PlotLifts(participants, series)
% Creates participant specific figures illustrating all single trials performed
%
% INPUTS
%   participants: participant number such as 1, 
%   a range, e.g., [1:12], or a set [1 3 4 5]
%
%   series: series numbers e.g., 1, [1:10], [1 4 5]
%
% This MATLAB / Octave script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava


global ColNames
global ws

global SGF_WIN
SGF_WIN = 31; % Window size of sgolayfilt(x, 3, SGF_WIN)

if nargin < 2, series = 1:9; disp('Warning: Using default series (all) '); end
if nargin < 1, participants = 1:13; disp('Warning: Using default Participants (all) '); end


%% For each participant
for i = participants

    if ((strcmp(computer,'PCWIN')) || (strcmp(computer, 'PCWIN64')))
        dirLoc = ['..\P', int2str(i), '\'];
    else
        dirLoc = ['../P', int2str(i), '/'];
    end

    %load participant AllLifts structure
    AllLiftsFName = [dirLoc, 'P', int2str(i), '_AllLifts.mat'];
    try
        load(AllLiftsFName);
    catch
        error('ERROR: WEEG_PlotLifts.m should be in the same directory as PX_AllLifts');
    end
    
    ColNames = P.ColNames;  %get column names
    
    fprintf('\nLoaded file %s', AllLiftsFName); 
        
    
    %% For each file/series
    for CFile = series
        
            
        f = figure; set(f, 'Name', sprintf('Participant: %d -- Series %d', i, CFile))
        
        %load WS file
        WSFName = [dirLoc, 'WS_P', int2str(i), '_S', int2str(CFile),'.mat'];
        load(WSFName);
        fprintf('\nLoaded file %s', WSFName)
        
        NumLifts = numel(ws.win);
               
        %row in P.AllLifts of the first lift in this series
        SeriesLiftRows = find(P.AllLifts(:,2)==CFile);
        FirstLiftRow = SeriesLiftRows(1);
        
        %column indices
        cGF = LOC_GetVarColumns('GF', ws.names.kin);
        cLF = LOC_GetVarColumns('LF', ws.names.kin);
        cWristPosx = LOC_GetVarColumns('Px4 - position x sensor 4', ws.names.kin);
        cWristPosy = LOC_GetVarColumns('Py4 - position y sensor 4', ws.names.kin);
        cWristPosz = LOC_GetVarColumns('Pz4 - position z sensor 4', ws.names.kin);
        cLEDOn = LOC_GetVarColumns('LEDOn', P.ColNames);
        cLEDOff = LOC_GetVarColumns('LEDOff', P.ColNames);
        cOnsetHand = LOC_GetVarColumns('tHandStart', P.ColNames);
        cBothDigitTouched = LOC_GetVarColumns('tBothDigitTouch', P.ColNames);
        cBothReleased = LOC_GetVarColumns('tBothReleased', P.ColNames);
        cLiftOff = LOC_GetVarColumns('tLiftOff', P.ColNames);
        cReplace = LOC_GetVarColumns('tReplace', P.ColNames);
        cBothStartLoadPhase = LOC_GetVarColumns('tBothStartLoadPhase', P.ColNames);
        cCurW = LOC_GetVarColumns('CurW', P.ColNames);
        cCurS = LOC_GetVarColumns('CurS', P.ColNames);
        
        %% For each lift/trial
        for L = 1:NumLifts
            subplot(6,6, L); disp(sprintf('Participant %d Series %d Lift/Trial %d', i, CFile, L));
            
            %grip and load force
            GF = -ws.win(L).kin(:,cGF); 
            LF = ws.win(L).kin(:,cLF); 
            
            %hand velocity
            HandPosx = ws.win(L).kin(:, cWristPosx);
            HandPosy = ws.win(L).kin(:, cWristPosy);
            HandPosz = ws.win(L).kin(:, cWristPosz);
            HandVel = LOC_GetHandVel(HandPosx, HandPosy, HandPosz);
            
            %selected events
            t_LEDOn = P.AllLifts(FirstLiftRow + L - 1, cLEDOn);
            t_LEDOff = P.AllLifts(FirstLiftRow + L - 1, cLEDOff);   
            t_OnsetHand = P.AllLifts(FirstLiftRow + L - 1, cOnsetHand);
            t_BothDigitTouched = P.AllLifts(FirstLiftRow + L - 1, cBothDigitTouched);
            t_BothReleased = P.AllLifts(FirstLiftRow + L - 1, cBothReleased);
            t_LiftOff = P.AllLifts(FirstLiftRow + L - 1, cLiftOff);
            t_Replace = P.AllLifts(FirstLiftRow + L - 1, cReplace);
            %t_BothStartLoadPhase = P.AllLifts(FirstLiftRow + L - 1, cBothStartLoadPhase);
 
            %current weight and surface
            CurW = P.AllLifts(FirstLiftRow + L - 1, cCurW);
            CurS = P.AllLifts(FirstLiftRow + L - 1, cCurS);
            
            %% plot data
            %grip force
            t = 0.002*(1:numel(GF));
            plot(t, GF, 'b', 'LineWidth', 1.5); 
            hold on
            
            %load force
            plot(t, LF, 'g', 'LineWidth', 1.5);
            
            %hand velocity (scaled for visibility)
            plot(t, HandVel/10, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.3);
            
            %% plot events
            plot([t_LEDOn t_LEDOn],  [0 10], '--', 'Color',[0.8 0.5 0.2], 'LineWidth', 1.5)
            plot([t_OnsetHand t_OnsetHand],[0 10], 'k--', 'LineWidth', 1.5);
            plot([t_BothDigitTouched t_BothDigitTouched],[0 10], 'b--', 'LineWidth', 1.5);
            %plot([t_BothStartLoadPhase t_BothStartLoadPhase],[0 10], 'r');
            plot([t_LiftOff t_LiftOff],[0 10], '--', 'Color', [0.9 0.2 0.3], 'LineWidth', 1.5);
            plot([t_LEDOff t_LEDOff],[0 10], '--', 'Color', [0.8 0.3 0.2], 'LineWidth', 1.5);
            plot([t_Replace t_Replace],[0 10], '--', 'Color', [0.8 0.4 0.8], 'LineWidth', 1.5);         
            plot([t_BothReleased t_BothReleased],[0 10], '--', 'Color', [0.3 0.4 0.8], 'LineWidth', 1.5);
            
            %axis scaling
            axis([0 length(GF)/500 0 10]);
            
            %plot title
            title([num2str(L), ' - W ', num2str(CurW)', ' - S ', num2str(CurS)]);
        end
        legend('Grip Force', 'Load Force', 'Hand Velocity (scaled)',...
            'Event - LED On', 'Event - Hand Starts', 'Event - Touch', 'Event - LiftOff', ...
            'Event - LED Off', 'Event - Placed Down', 'Event - Release', ...
            'Location', 'BestOutside')
    end
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

function [X, dX, ddX] = LOC_GetDerivatives(X, dt)
%LOC_GetDerivatives

global SGF_WIN
    if nargin < 2 
        dt = 0.002; % 500 Hz sampling rate
    end
    dX = sgolayfilt([0; diff(X)], 3, SGF_WIN)/dt;
    ddX = sgolayfilt([0; diff(dX)], 3, SGF_WIN)/dt;
return

function [HandVel] = LOC_GetHandVel(X, Y, Z)
%LOC_GetHandVel
   
    [~, dX] = LOC_GetDerivatives(X);
    [~, dY] = LOC_GetDerivatives(Y);
    [~, dZ] = LOC_GetDerivatives(Z);
    HandVel = sqrt(dX.^2 + dY.^2 + dZ.^2);
return
