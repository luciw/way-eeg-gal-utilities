function WEEG_PlotStats(participants, inOUTLIER_SD)
% function WEEG_PlotStats(participants, inOUTLIER_SD)
% Generates histograms for the following 
%  - Time difference between digit contacts & Preload phase duration
%  - Load phase duration (3x3)
%  
% INPUTS
%   participants        - Participant(s) to be displayed (default = 1)
% 	inOUTLIER_SD  - Defines outliers eliminated from histograms = x SD
%                   (default = 4)
%
% This MATLAB / Octave script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava


global OUTLIER_SD 

if nargin < 1, participants = 1; end
if nargin < 2, inOUTLIER_SD = 4; end
close all
OUTLIER_SD = inOUTLIER_SD;
WEIGHT_LIGHTEST_OBJECT = 165; %g

for i = participants
    f = figure;
    set(f, 'Name', sprintf('Participant %d', i)); 
    
    
    if ((strcmp(computer,'PCWIN')) || (strcmp(computer, 'PCWIN64')))
        dirLoc = ['..\P', int2str(i), '\'];
    else
        dirLoc = ['../P', int2str(i), '/'];
    end

    try
        load([dirLoc, 'P', int2str(i), '_AllLifts.mat'])
    catch
        error('ERROR: PX_AllLifts should be in directory ../PX');
    end
    
    % Digit time differences
    subplot(1,2,1);
    x = LOC_GetVarColumns('tIndTouch', P.ColNames);
    y = LOC_GetVarColumns('tThumbTouch', P.ColNames);
    X = 1000*(P.AllLifts(:,x) - P.AllLifts(:,y));
    Range = -250:5:250;
    Ti = ['Index contact relative thumb contact ({\mu}=' ...
        sprintf('%0.0f ms)', mean(X))];
    YL = 'Count'; XL = 'ms';
    LOC_MakeHisto(X, Ti, YL, XL, Range)
    
    % Preload phase duration
    subplot(1,2,2);
    x = LOC_GetVarColumns('Dur_Preload', P.ColNames);
    X = 1000*(P.AllLifts(:,x));
    Ti = sprintf('Preload phase duration)',mean(X));
    YL = 'Count'; XL = 'ms';
    Range =0:5:400;
    LOC_MakeHisto(X, Ti, YL, XL, Range)
        
    % 'Dur_LoadPhase'
    f = figure;
    set(f, 'Name', sprintf('Participant %d - LoadPhase', i));
    x = LOC_GetVarColumns('Dur_LoadPhase', P.ColNames);
    clear X;
    Weis = [1 2 4];
    for WP = 1:3
        for WC = 1:3
            q = find((P.AllLifts(:, 6) == Weis(WP)) & ...
                (P.AllLifts(:, 4) == Weis(WC)));
            X = 1000*P.AllLifts(q, x);
            subplot(3,3,WP + (WC-1)*3);
            Ti = '';
            YL = ''; XL = '';
            Range =0:10:1000;
            LOC_MakeHisto(X, Ti, YL, XL, Range)
        end
    end
    for WP = 1:3
        subplot(3,3,WP);
        title(sprintf('Expected weight %dg', ...
            WEIGHT_LIGHTEST_OBJECT*Weis(WP)));
    end
    
    for WC = 1:3
        subplot(3,3,1 + (WC-1)*3);
        ylabel(sprintf('Current %dg', WEIGHT_LIGHTEST_OBJECT*Weis(WC)));
    end
        
end

return

%% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Local functions
%  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

%LOC_EliminateOutliers
function [X, R] = LOC_EliminateOutliers(X)
global OUTLIER_SD 
    SD = std(X);
    R = [mean(X)-OUTLIER_SD*SD mean(X)+OUTLIER_SD*SD];
    q = find(abs(X-mean(X)) < OUTLIER_SD*SD);
    X = X(q);
return

%LOC_MakeHisto
function LOC_MakeHisto(X, Ti, YL, XL, Range)
    [X, R] = LOC_EliminateOutliers(X);
    hist(X, Range)
    m = mean(X);
    if numel(Ti) > 0, title(Ti); end
    if numel(YL) > 0, ylabel(YL); end
    if numel(XL) > 0, xlabel(XL); end
    hold on; 
    ax = axis;
    axis([[Range(1) Range(end)] ax(3:4)]);
    plot([m m], [ax(3:4)], 'r')
    plot([R(1) R(1)], [ax(3:4)], ':r')
    plot([R(2) R(2)], [ax(3:4)], ':r')
return
