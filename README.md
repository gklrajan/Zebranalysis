# Zebranalysis

Analysis pipeline for animal swimming parameters obtained from the MBD program.
Started on 2017-08-11.

A sample binary file is provided: 2017-10-17_20.42.24_swim18_trial10.bin contains frame#, cameraTime(us), systemTime(ms), xcoords, ycoords, fishAngle(rads) (in that precise order). This acquisition was done carried out at 750 Hz. (It is a 13 minutes acquisition)

binreader_pix.m will be the first analysis project looking at the swimming kinematics of the fish.
