function WEEG_FindEvents(participants, savefiles, debug)
%% function WEEG_FindEvents(participants, save, debug)
% Creates a participant specific matrix describing all single trials performed
%
% INPUTS
%   InPartNums participant number (1-12) or a range, e.g., 1:12; default = 1
%   savefiles if true, PX_AllLifts.mat file(s) are created
%   debug  if true, debug info is displayed
%
% OUTPUTS
%   Generates a file in this folder with the name
%   Px_AllLifts.mat where x is Participant number, 
%   e.g, P5_AllLifts.mat. This file
%   contain a structure with two fields: 
%       AllLifts a matrix (number of trials, 39 columns)
%       ColNames a vector (39 column names)
%
% This MATLAB / Octave script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava

global ColNames
global ws
global SGF_WIN
SGF_WIN = 31; % Window size of sgolayfilt(x, 3, SGF_WIN)

%% Algorithm
% 1. Find t_LiftOff and t_Replace
%     - unequivocal events that can ground both initial touch, load
%       phase and relase 
% 2. Find touch
% 3. Find release
% 4. Identify preload phase, load phase and release phase
% 5. Hand related variables

%% Set default parameters primarily for debugging purposes  
if nargin < 1, participants = 1; disp('Warning: Using default Participant (1) '); end
if nargin < 2, savefiles = true; end
if nargin < 3, debug = false; end
    
ColNames = {'Part', 'Run', 'Lift', ...
            'CurW', 'CurS', 'PrevW', 'PrevS', ...
            'StartTime', 'LEDOn', 'LEDOff', 'BlockType', ...
            'tIndTouch', 'tThumbTouch', 'tFirstDigitTouch', 'tBothDigitTouch', ...
            'tIndStartLoadPhase', 'tThuStartLoadPhase', 'tBothStartLoadPhase', ...
            'tLiftOff', 'tReplace', ...
            'tIndRelease', 'tThuRelease', 'tBothReleased', ...
            'GF_Max', 'LF_Max', 'dGF_Max', 'dLF_Max', ...
            'tGF_Max', 'tLF_Max', 'tdGF_Max', 'tdLF_Max', ...
            'GF_Hold', 'LF_Hold', ...
            'tHandStart', 'tHandStop', ...
            'tPeakVelHandReach', 'tPeakVelHandRetract', ...
            'GripAparture_Max', 'tGripAparture_Max', ...
            'Dur_Reach', 'Dur_Preload', 'Dur_LoadPhase', 'Dur_Release'}';

if debug, close all, end    
    
%for each participant
for i = participants
    %% Load WS file to be processed
    if ((strcmp(computer,'PCWIN')) || (strcmp(computer, 'PCWIN64')))
        DirLoc = ['..\P', int2str(i), '\'];
    else
        DirLoc = ['../P', int2str(i), '/']; 
    end
    
    Files = dir([DirLoc, 'WS_P', int2str(i), '_*.mat']);

    
    AllLifts = [];
    
    if (numel(Files) == 0)
        error('ERROR: The WS_PX_SY.mat files are not in the directory ../PX');
    end
    
    %for each file
    for CFile = 1:numel(Files)
        load([DirLoc, Files(CFile).name]);
        nLifts = numel(ws.win);
        Temp_AllLifts = zeros(nLifts, numel(ColNames));
        Temp_AllLifts(:,LOC_GetVarColumns('Part', ColNames)) = i;
        Temp_AllLifts(:,LOC_GetVarColumns('Run', ColNames)) = CFile;
        disp(sprintf('Participant %2.d: Processing file %s (%d/%d)', ...
            i, Files(CFile).name, CFile, numel(Files)));
        
        if debug
            f = figure; set(f, 'Name', sprintf('Participant: %d -- File %d', i, CFile))
        end
        
        %for each lift
        for L = 1:nLifts
            if debug
                subplot(6,6, L); disp(sprintf('File %d Trial %d', CFile, L));
            end;
            c = LOC_GetVarColumns('Lift', ColNames); 
                Temp_AllLifts(L,c) = L;
            % Weights and surface
            c = LOC_GetVarColumns('CurW', ColNames); 
                Temp_AllLifts(L,c) = ws.win(L).weight;
                CurW = Temp_AllLifts(L,c);
            c = LOC_GetVarColumns('CurS', ColNames);
                Temp_AllLifts(L,c) = ws.win(L).surf;
                CurS = Temp_AllLifts(L,c);
            c = LOC_GetVarColumns('PrevW', ColNames);
                if ws.win(L).weight_prev < 0, Temp_AllLifts(L,c) = NaN;
                else Temp_AllLifts(L,c) = ws.win(L).weight_prev; end
            c = LOC_GetVarColumns('PrevS', ColNames);
                if ws.win(L).surf_prev < 0, Temp_AllLifts(L,c) = NaN;
                else Temp_AllLifts(L,c) = ws.win(L).surf_prev; end
            c = LOC_GetVarColumns('StartTime', ColNames);
                Temp_AllLifts(L,c)=ws.win(L).trial_start_time;
            c = LOC_GetVarColumns('LEDOn', ColNames);
                Temp_AllLifts(L,c)=ws.win(L).LEDon;
            c = LOC_GetVarColumns('LEDOff', ColNames);
                Temp_AllLifts(L,c)=ws.win(L).LEDoff;
            
            % Find data column in ws.names.kin
            cIndGF = LOC_GetVarColumns('FZ1 - force z plate 1', ws.names.kin);
            cThuGF = LOC_GetVarColumns('FZ2 - force z plate 2', ws.names.kin);
            cIndLF = LOC_GetVarColumns('FX1 - force x plate 1', ws.names.kin);
            cThuLF = LOC_GetVarColumns('FX2 - force x plate 2', ws.names.kin);
            cObjVPos = LOC_GetVarColumns('Pz1 - position z sensor 1', ws.names.kin);
            cWristPosx = LOC_GetVarColumns('Px4 - position x sensor 4', ws.names.kin);
            cWristPosy = LOC_GetVarColumns('Py4 - position y sensor 4', ws.names.kin);
            cWristPosz = LOC_GetVarColumns('Pz4 - position z sensor 4', ws.names.kin);
            cIndPosx = LOC_GetVarColumns('Px2 - position x sensor 2', ws.names.kin);
            cIndPosy = LOC_GetVarColumns('Py2 - position y sensor 2', ws.names.kin);
            cIndPosz = LOC_GetVarColumns('Pz2 - position z sensor 2', ws.names.kin);
            cThuPosx = LOC_GetVarColumns('Px3 - position x sensor 3', ws.names.kin);
            cThuPosy = LOC_GetVarColumns('Py3 - position y sensor 3', ws.names.kin);
            cThuPosz = LOC_GetVarColumns('Pz3 - position z sensor 3', ws.names.kin);
      
             
            % Get the data streams to be used for trial data
            GF_Index = - ws.win(L).kin(:, cIndGF);
            GF_Thumb = - ws.win(L).kin(:, cThuGF);
            LF_Index = ws.win(L).kin(:, cIndLF);
            LF_Thumb = ws.win(L).kin(:, cThuLF);
            VObjPos = ws.win(L).kin(:, cObjVPos);
            HandPosx = ws.win(L).kin(:, cWristPosx);
            HandPosy = ws.win(L).kin(:, cWristPosy);
            HandPosz = ws.win(L).kin(:, cWristPosz);
            IndPosx = ws.win(L).kin(:, cIndPosx);
            IndPosy = ws.win(L).kin(:, cIndPosy);
            IndPosz = ws.win(L).kin(:, cIndPosz);
            ThuPosx = ws.win(L).kin(:, cThuPosx);
            ThuPosy = ws.win(L).kin(:, cThuPosy);
            ThuPosz = ws.win(L).kin(:, cThuPosz);
            WristPosy = ws.win(L).kin(:,cWristPosy);
            t = 0.002*(1:numel(GF_Index));
            
            
            %% STEP 1: Find t_LifOff and t_Replace
            [t_LiftOff, t_Replace] = LOC_GetLiftOffReplace(VObjPos);
            Temp_AllLifts(L, LOC_GetVarColumns('tLiftOff', ColNames)) = t_LiftOff;
            Temp_AllLifts(L, LOC_GetVarColumns('tReplace', ColNames)) = t_Replace;
            
            %% STEP 2: Digit touches
            [t_IndTouch] = LOC_GetTouches(GF_Index);
            [t_ThuTouch] = LOC_GetTouches(GF_Thumb);
            t_FirstDigitTouch = min(t_IndTouch, t_ThuTouch);
            t_BothDigitTouch = max(t_IndTouch, t_ThuTouch);
            Temp_AllLifts(L, LOC_GetVarColumns('tIndTouch', ColNames)) = t_IndTouch;
            Temp_AllLifts(L, LOC_GetVarColumns('tThumbTouch', ColNames)) = t_ThuTouch;
            Temp_AllLifts(L, LOC_GetVarColumns('tFirstDigitTouch', ColNames)) = t_FirstDigitTouch;
            Temp_AllLifts(L, LOC_GetVarColumns('tBothDigitTouch', ColNames)) = t_BothDigitTouch;
            
            %% STEP 3: Timing of digit releases
            [t_IndRelease] = LOC_GetRelease(GF_Index, t_Replace);
            [t_ThuRelease] = LOC_GetRelease(GF_Thumb, t_Replace);
            t_BothReleased = max(t_IndRelease, t_ThuRelease);
            Temp_AllLifts(L, LOC_GetVarColumns('tIndRelease', ColNames)) = t_IndRelease;
            Temp_AllLifts(L, LOC_GetVarColumns('tThuRelease', ColNames)) = t_ThuRelease;
            Temp_AllLifts(L, LOC_GetVarColumns('tBothReleased', ColNames)) = t_BothReleased;
            
            %% STEP 4: Identify preload phase, load phase and release phase
            [t_IndStartLoadPhase] = LOC_GetPreload(LF_Index, t_IndTouch);
            [t_ThuStartLoadPhase] = LOC_GetPreload(LF_Thumb, t_ThuTouch);
            t_BothStartLoadPhase = max(t_IndStartLoadPhase, t_ThuStartLoadPhase);
            Dur_Preload = min(t_IndStartLoadPhase, t_ThuStartLoadPhase) - t_FirstDigitTouch;
            Dur_LoadPhase = t_LiftOff - t_BothStartLoadPhase;
            Dur_Release = t_BothReleased - t_Replace;
            Temp_AllLifts(L, LOC_GetVarColumns('tIndStartLoadPhase', ColNames)) = t_IndStartLoadPhase;
            Temp_AllLifts(L, LOC_GetVarColumns('tThuStartLoadPhase', ColNames)) = t_ThuStartLoadPhase;
            Temp_AllLifts(L, LOC_GetVarColumns('tBothStartLoadPhase', ColNames)) = t_BothStartLoadPhase;
            Temp_AllLifts(L, LOC_GetVarColumns('Dur_Preload', ColNames)) = Dur_Preload;
            Temp_AllLifts(L, LOC_GetVarColumns('Dur_LoadPhase', ColNames)) = Dur_LoadPhase;
            Temp_AllLifts(L, LOC_GetVarColumns('Dur_Release', ColNames)) = Dur_Release;
            
            [MaxGF, tMaxGF, MaxdGF, tMaxdGF] = LOC_GetPeak(GF_Index, GF_Thumb);
            [MaxLF, tMaxLF, MaxdLF, tMaxdLF] = LOC_GetPeak(LF_Index, LF_Thumb); MaxLF = 2*MaxLF;
            Temp_AllLifts(L, LOC_GetVarColumns('GF_Max', ColNames)) = MaxGF;
            Temp_AllLifts(L, LOC_GetVarColumns('tGF_Max', ColNames)) = tMaxGF;
            Temp_AllLifts(L, LOC_GetVarColumns('LF_Max', ColNames)) = MaxLF;
            Temp_AllLifts(L, LOC_GetVarColumns('tLF_Max', ColNames)) = tMaxLF;
            Temp_AllLifts(L, LOC_GetVarColumns('dGF_Max', ColNames)) = MaxdGF;
            Temp_AllLifts(L, LOC_GetVarColumns('tdGF_Max', ColNames)) = tMaxdGF;
            Temp_AllLifts(L, LOC_GetVarColumns('dLF_Max', ColNames)) = MaxdLF;
            Temp_AllLifts(L, LOC_GetVarColumns('tdLF_Max', ColNames)) = tMaxdLF;
            
            %% STEP 5: Hand related variables
            [t_OnsetHand, HandVel, pHandVelReach, pHandVelRetract] = ...
                LOC_GetHand(HandPosx, HandPosy, HandPosz, ...
                t_FirstDigitTouch, t_Replace);
            [t_HandStop] = LOC_HandStop(IndPosy, ThuPosy, WristPosy, ws.win(L).LEDoff);
            [MGA, tMGA] = LOC_GetMGA(t_BothDigitTouch, ...
                IndPosx, IndPosy, IndPosz, ...
                ThuPosx, ThuPosy, ThuPosz);
            DurReach = t_FirstDigitTouch - t_OnsetHand;
            Temp_AllLifts(L, LOC_GetVarColumns('tHandStart', ColNames)) = t_OnsetHand;
            Temp_AllLifts(L, LOC_GetVarColumns('tHandStop', ColNames)) = t_HandStop;
            Temp_AllLifts(L, LOC_GetVarColumns('tPeakVelHandReach', ColNames)) = pHandVelReach;
            Temp_AllLifts(L, LOC_GetVarColumns('tPeakVelHandRetract', ColNames)) = pHandVelRetract;
            Temp_AllLifts(L, LOC_GetVarColumns('Dur_Reach', ColNames)) = DurReach;
            Temp_AllLifts(L, LOC_GetVarColumns('GripAparture_Max', ColNames)) = MGA;
            Temp_AllLifts(L, LOC_GetVarColumns('tGripAparture_Max', ColNames)) = tMGA;
            
            %% STEP 6: Hold phase related variables
            GF_Index = - ws.win(L).kin(:, cIndGF);
            GF_Thumb = - ws.win(L).kin(:, cThuGF);
                GF = (GF_Index + GF_Thumb)/2;
            LF_Index = ws.win(L).kin(:, cIndLF);
            LF_Thumb = ws.win(L).kin(:, cThuLF);
                LF = LF_Index + LF_Thumb;
            t_LEDOff = Temp_AllLifts(L,LOC_GetVarColumns('LEDOff', ColNames));
            t_LEDOn = Temp_AllLifts(L,LOC_GetVarColumns('LEDOn', ColNames));
            CutOut = [(round(t_LEDOff*500)-300):(round(t_LEDOff*500)-100)];
            Temp_AllLifts(L, LOC_GetVarColumns('GF_Hold', ColNames)) = ...
                mean(GF(CutOut));
            Temp_AllLifts(L, LOC_GetVarColumns('LF_Hold', ColNames)) = ...
                mean(LF(CutOut));
            if debug
                if L == 6
                    [GF, dGF] = LOC_GetDerivatives(GF);
                    [LF, dLF] = LOC_GetDerivatives(LF);
                    QQQQ = 1; % Suitable debug stop before manually running 
                              % the code at the bottom of the file  
                end
                
                %%plotting 
                %grip force
                t = 0.002*(1:numel(GF));
                plot(t, GF, 'b', 'LineWidth', 1.5); 
                hold on
            
                %load force
                plot(t, LF, 'g', 'LineWidth', 1.5);
            
                %hand velocity (scaled for visibility)
                plot(t, HandVel/10, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.3);
            
                %plot (some) events
                plot([t_LEDOn t_LEDOn],  [0 10], '--', 'Color',[0.8 0.5 0.2], 'LineWidth', 1.5)
                plot([t_OnsetHand t_OnsetHand],[0 10], 'k--', 'LineWidth', 1.5);
                plot([t_BothDigitTouch t_BothDigitTouch],[0 10], 'b--', 'LineWidth', 1.5);
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
            
        end % for L = 1:nLifts
        
        if (debug)
            legend('Grip Force', 'Load Force', 'Hand Velocity (scaled)',...
                'Event - LED On', 'Event - Hand Starts', 'Event - Touch', 'Event - LiftOff', ...
                'Event - LED Off', 'Event - Placed Down', 'Event - Release', ...
                'Location', 'BestOutside')
        end
        
        %% STEP 6: Set block type (1=Weight, 2=Surface and 3=Mixed)
        c = LOC_GetVarColumns('CurW', ColNames)';
            Ws = unique(Temp_AllLifts(:,c));
        c = LOC_GetVarColumns('CurS', ColNames)';
            Ss = unique(Temp_AllLifts(:,c));
        c = LOC_GetVarColumns('BlockType', ColNames);
            Temp_AllLifts(:,c) = (numel(Ws)>1) + 2*(numel(Ss)>1);
        
        %% Add Temp_AllLifts to AllLifts
        AllLifts = [AllLifts; Temp_AllLifts];
        
    end % for CFile = 1:numel(Files)
    
    %% Write output file
    clear P
    P.AllLifts = AllLifts;
    P.ColNames = ColNames;
    FileName = ['P', int2str(i), '_AllLifts.mat'];
    
    if (savefiles)
        save(FileName, 'P');
        disp(sprintf('Saved data in %s', FileName))
        disp(' ')
    end
    
end % for i = participants
return

%% :::::::::::::::::::::::::::::::::
%  LOCAL utility functions
%  :::::::::::::::::::::::::::::::::

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

function [t_LiftOff, t_Replace] = LOC_GetLiftOffReplace(VP)
    [VP, dVP] = LOC_GetDerivatives(VP);
    t = 0.002*(1:numel(VP));
    ZeroVP = mean(VP(1:500));
    MaxVP = max(VP); 
    % There may be a trailing replace from previous trial, so start at index 500
    t_MaxVP = find(VP(500:end) == max(VP(500:end))); t_MaxVP = 499 + t_MaxVP(1);
    t_50 = find(VP(500:t_MaxVP)-ZeroVP < 0.5*(MaxVP-ZeroVP));
    t_50 = 499 + t_50(end);
    P = find(dVP(t_50:-1:1)>1); % Check 'backward in time'
    dP = diff(P); q = find(dP > 1);
    if ~isempty(q), P = P(q(1));
    else P = P(end); end
    t_LiftOff = t_50 - P;  
    
    t_50 = find(VP(t_MaxVP:end)-ZeroVP < 0.5*(MaxVP-ZeroVP));
    t_50 = t_MaxVP + t_50(1) - 1;
    P = find(dVP(t_50:end)<-2); 
    dP = diff(P); q = find(dP > 1);
    if ~isempty(q), P = P(q(1));
    else P = P(end); end
    t_Replace = t_50 + P;  
    t_LiftOff = t(t_LiftOff); t_Replace = t(t_Replace);
return

function [MGA, tMGA] = LOC_GetMGA(tTouch, IndPosx, IndPosy, IndPosz, ...
                                 ThuPosx, ThuPosy, ThuPosz)

    GA = sqrt( (IndPosx-ThuPosx).^2 + ...
               (IndPosy-ThuPosy).^2 + ...
               (IndPosz-ThuPosz).^2 );
    t = 0.002*(1:numel(IndPosx));
    MGA = max(GA(500:round(500*tTouch)));
    q = find(GA == MGA);
    tMGA = t(q(1));
return                             

function [MaxF, tMaxF, MaxdF, tMaxdF] = LOC_GetPeak(F1, F2)
    t = 0.002*(1:numel(F1));
    F = (F1+F2)/2;
    [F, dF] = LOC_GetDerivatives(F);
    MaxF = max(F(500:end));
    q = find(MaxF == F);
    tMaxF = t(q(1));
    MaxdF = max(dF(500:end));
    q = find(MaxdF == dF);
    tMaxdF = t(q(1));
return    

function [t_Touch, MaxF] =  LOC_GetTouches(F)

    [~, dF] = LOC_GetDerivatives(F);
    MaxF = max(F(800:round(end*3/4)));
    t_MaxF = find(F(800:end) == MaxF); t_MaxF = 799 + t_MaxF(1);
    t_50 = find(F(800:t_MaxF) < 0.25*max(F(800:end)));
    
    LO = 799 + t_50(end);
    t = 0.002*(1:numel(F));
    m_dF = mean(dF(LO+(-500:-400))); sd_dF = std(dF(LO+(-500:-400)));
    q = find( dF(LO:-1:1) > m_dF + 4*sd_dF);
    qq =find(diff(q)>1);
    if ~isempty(qq)
        q = q(1:qq(1));
    end
    t_Touch = t(LO-q(end) + 1);
return

function [t_Release] = LOC_GetRelease(F, t_Replace)
    [F, dF] = LOC_GetDerivatives(F);
    RE = round(t_Replace*500);
    t = 0.002*(1:numel(F));
    m_dF = mean(dF(700:1000)); sd_dF = std(dF(700:1000));
    q = find( dF(RE:end) < m_dF - 2*sd_dF);
    qq =find(diff(q)>1);
    if ~isempty(qq)
        q = q(1:qq(1));
    end
    if isempty(q)
        a = 2;
    end;
    t_Release = t(RE + q(end) - 1);
return
    
function [T, MaxLF] = LOC_GetPreload(LF, t_Touch)
    t = 0.002*(1:numel(LF));
    MaxLF = max(LF(500:end));
    START = round(t_Touch*500);
    [LF, ~, ddLF] = LOC_GetDerivatives(LF);
    % 1. Find the peak LF
    P = find(LF(START:end) == max(LF(START:end))); P = START + P(1) - 1 ;
    % 2. Find LF > 0.2
    qP = find(LF(P:-1:START) > 0.2);  
    Q = find(diff(qP) > 1); 
    if ~isempty(Q), 
        qP = qP(1:(Q(1))); 
    end   
    P = START + length(P:-1:START) - qP(end);
    % 3. Find zero crossing of ddLF 
    qq = find( ddLF(P:-1:START) > 0);
    Q = find(diff(qq) > 1);
    if ~isempty(Q)
        qq = qq(1:(Q(1)));
    end
    if isempty(qq)
        T = P-length(START:P)+1;
    else
        T = P-qq(end)+1;
    end
    if LF(T)<0
        if LF(T) > min(LF(START:P))
            q = find(LF(START:P) == min(LF(START:P)));
            T = START + q(1);
        end
    end
    T = t(T);
return

function [T, HandVel, pHandVelReach, pHandVelRetract] = ...
                               LOC_GetHand(X, Y, Z, t_Touch, t_Replace)
    t = 0.002*(1:numel(X));
    P = round(500*t_Touch) - 200;    
    [~, dX] = LOC_GetDerivatives(X);
    [~, dY] = LOC_GetDerivatives(Y);
    [~, dZ] = LOC_GetDerivatives(Z);
    HandVel = sqrt(dX.^2 + dY.^2 + dZ.^2);
    
    % First find an approximate start, i.e., Vel > 5 cm/s
    q = find(HandVel(P:-1:1) > 5); 
    qq =find(diff(q)>1);
    if ~isempty(qq)
        q = q(1:qq(1));
    end
    T = P-q(end) + 1; 
    
    % Then search for for mean + 2·SD
    mWin = (T-600):(T-100);
    m_HV = mean(HandVel(mWin)); SD_HV = std(HandVel(mWin));
    q = find(HandVel(P:-1:1) > m_HV + 2*SD_HV); 
    qq =find(diff(q)>1);
    if ~isempty(qq)
        q = q(1:qq(1));
    end
    pHandVelReach = max(HandVel(T:round(500*t_Touch)));
    pHandVelRetract = max(HandVel(round(500*t_Replace):end));
    if pHandVelRetract == HandVel(end)
        pHandVelRetract = NaN;
    end
    
    T = t(T);
return

function [t_HandStop] = LOC_HandStop(cIndPosy, cThuPosy, cWristPosy, LEDoff)

    ForwardArmPos = mean([cIndPosy cThuPosy cWristPosy]');
    WinLength = length(cIndPosy);
    meanFAP = mean(ForwardArmPos(ceil(WinLength/2):end));
    maxFAP = max(ForwardArmPos(ceil(WinLength/2):end));
    FAPThresh = maxFAP * 0.90 + meanFAP * 0.1;

    EventFound = false;
    j=round(LEDoff * 500); 
    while (EventFound == false)
      if (ForwardArmPos(j) > FAPThresh)
        EventFound = true;
      end
      j = j + 1;
      if (j > numel(ForwardArmPos))
          QQQQ = 1   
      end
    end
    t_HandStop = (j/500);   
            
return

function [X, dX, ddX] = LOC_GetDerivatives(X, dt)
global SGF_WIN
    if nargin < 2 
        dt = 0.002; % 500 Hz sampling rate
    end
    dX = sgolayfilt([0; diff(X)], 3, SGF_WIN)/dt;
    ddX = sgolayfilt([0; diff(dX)], 3, SGF_WIN)/dt;
return


