=======
## Beam Test Data Analysis 
##Merging Code: Takes Ntuple from EDM, Telescope and merges into one Analysis Tree for user
The code fill the Track vector inside the EDM ntuple
## Installation

Set up ROOT(works with ROOT 6)

After checking out the code:-
make cint

make 

## Usage
##
./ntuplemerger \<EDMTupleName\> \<TelescopeTupleName\> \<RunNumber\>

It produces 1 file AnalysisTree_RunNumber.root
