way-eeg-gal-utilities
=====================

A set of Matlab scripts for working with the WAY-EEG-GAL dataset.  

WEEG_GetEventsInHS.m returns the times of various events within the series, instead of within the windowed trials. One can select the participant and a particular series type – Weight, Friction, Mixed, or All. 
WEEG_PlotLifts.m enables a per-window visualization of a participant’s activities. One can select the participants and series to plot. Each subplot shows a different trial and displays three signals: the grip force, the load forces and the hand velocity. Seven events are indicated by dotted vertical lines. The events shown are the time of: LEDon, when the hand starts moving, first contact, liftoff, LEDoff, the object is placed down, and object release. Above each subplot, the weight and surface type are indicated.

WEEG_PlotStats.m displays histograms indicating, per-participant, the time of index finger contact relative to thumb contact, the duration of the preload phase, and the duration of the load phase. The load phase duration is broken into 9 subplots, shown in the 3 x 3 grid. The current weight is shown on the y-axis, while the weight in the previous trials is shown on the x-axis.

WEEG_FindEvents.m is the script used to determine event timings and lift characterizations, and was used to generate the P.AllLifts structure.

The MATLAB data files and scripts described above can be loaded and run with Octave 3.8 or higher (http://www.gnu.org/software/octave/).

The open-source MATLAB software EEGLab can be used to assist in processing the EEG signals. We provide two scripts for importing the data to EEGLab. 

WEEG_MakeEEGLABDataset.m imports the EEG and event data for all series for one participant, and saves it as a 'set'; 

WEEG_MakeAllEEGLABDatasets.m does the same but for all participants. The file chanlabels_32channel.xyz can be used to localize the electrode positions in EEGLab, for topographic plotting.

We also provide brief instructions for using EEGLAB for ERSP plotting with the WAY-EEG-GAL data.

Also see Usage.txt in the Utilities folder.