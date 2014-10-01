function WEEG_MakeAllEEGLABDatasets(SaveDir)
% function WEEG_MakeAllEEGLABDatasets(SaveDir)
% Creates EEGLAB datasets for all participants and series.
% To run this, EEGLAB must be on your MATLAB path
%
% INPUTS
%   SaveDir - Saving directory for EEGLAB datasets (default=current dir)
%
% This MATLAB script was created for use with the WAY-EEG-GAL dataset. 
% Copyright (c) 2014, Benoni Edin, Matthew Luciw, and Rupesh Srivastava


if nargin < 1
    SaveDir = '.';
end

Participants = 1:12;

for Participant = Participants
    fprintf(['\nLoading data for Participant ' num2str(Participant) '\n']);
    WEEG_MakeEEGLABDataset(Participant, SaveDir);
end

return