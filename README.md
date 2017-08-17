# Zebranalysis

Analysis pipeline for animal swimming parameters obtained from the MBD program.
Started on 2017-08-11.

Two sample bin files are provided. Bin file 8_9_2017_S3 contains only frame#, xy coords and angle; 
8_11_2017_S1 contains cameraTime, deltaAngle, and StimuliTrigger in addition to the ones available in 8_9_2017_3.

myBinReader can be used to read all these parameters from the bin file and save them to a csv file. 
The current version is adapted for the 8_11_2017_S1.bin file. But it can be easily downgraded to work with 8_9_2017_S3.bin.

freeswim_pipeline will be the first analysis project looking at the swimming kinematics of the fish.

%this is a test! %test123