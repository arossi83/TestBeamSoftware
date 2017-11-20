## Beam Test Data Analysis Software
## Installation

Set up ROOT(works with ROOT 5, requires ROOT to be compiled with libminut2)
##On lxplus, you can setup gcc and root as follows
source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

git clone https://github.com/sroychow/TestBeamSoftware.git

cd TestBeamSoftware

git checkout may16dev

#Set to NOV_15 or MAY_16 or OCT_16
setenv  BT_ERA OCT_16

make cint
make

## Usage
##Running on tuple

#Step1: Alignment Analysis

./alignmentReco \<jobcardName\>

Apart from a root file with all the histograms, the application also produces a text file with all the alignment parameters.*** 
This application can be run in "production" or "new" mode(can be set in the job card). If you run in production mode, the 
alignment output parameter file will be opened in append mode. It can be used for analysing many runs and dunping the 
output in the same file. In the "new" mode, a new alignment output file will be created and if a file exists with the same name,
 it will be overwritten.

#JobCard Description

Only the parameters required for this application are described.

Run=\<RUN-NUMBER\>

inputFile=AnalysisTree_\<RUN-NUMBER\>.root

outputFile=\<output-file name\>

fei4Z=\<FEI4-distance\>

doTelescopeMatching=1 #set to true; will look for telescope data; if =0, this step will not work

doChannelMasking=0 #=0 if no channel masked; =1 if channel masking required

channelMaskFile=\<filename\> #set name of file with channel mask; required if doChannelMasking=1 

nStrips=254 # No. of strips

pitchDUT=0.09 # Pitch

alignmentOutputFile=\<filename\> #Filename where the alignment parameters will be written

isProductionmode=0 #if =0 new mode; =1 append mode 
#Step2: Baseline Analysis to study detector performance

Only the parameters required for this application are described.

Run=\<RUN-NUMBER\>

inputFile=AnalysisTree_\<RUN-NUMBER\>.root

outputFile=\<output-file name\>

fei4Z=\<FEI4-distance\>

readAlignmentFromfile=1 #if=1, alignment parameters will be read from alignmentOutputFile; if=0, telescope matching will not work(or you have to hardcode values in your code! Not recommended!!)

doTelescopeMatching=1 #set to true; will look for telescope data; if =0, detector only performance

doChannelMasking=0 #=0 if no channel masked; =1 if channel masking required

channelMaskFile=\<filename\> #set name of file with channel mask; required if doChannelMasking=1 

nStrips=254 # No. of strips

pitchDUT=0.09 # Pitch

residualSigmaDUT=0.026 #residual(4times this value) in mm used for matching DUT hit/cluster/stub to Track

alignmentOutputFile=\<filename\> #Filename from where the alignment parameters will be read

#Alignment Paremter file format

Each line in the file corresponds to a Run and is a ":" separated list of all alignment parameters required by our analysis written out in the following order.

Run:offsetFEI4X:offsetFEI4Y:residualSigmaFEI4X:residualSigmaFEI4Y:zD0:offsetD0:deltaZ:angle

Run=\<Run Number\>

offsetFEI4X/Y=Offset of track with respect to FeI4 hits in X and Y

residualSigmaFEI4X/Y=Pitch at FeI4 plane obtained from telescope-Fei4 residual fits

zD0=Position of DUT0 w.r.t FeI4

offsetD0=Offset of track-DUT hits obtained from residual fit at DUT planr

deltaZ=difference in Z pos between DUT0 and DUT1

angle=DUT angle w.r.t beam

**If alignment parameters are read from file, the code searches for the Run Number and takes the alignment parameters from that line.
***For TDR, use the alignment parameters under data directory
