# Zebranalysis

Analysis pipeline for animal swimming parameters obtained from the MBD program.
Started on 2017-08-11.

A sample binary file is provided: 2017-10-31_15.20.14_swim_1_Sample15mins_1_.bin contains frame#, cameraTime(us), systemTime(ms), xcoords, ycoords, fishAngle(rads) (in that precise order). This acquisition was done carried out at 750 Hz. (It is a 15 minutes acquisition)

binreader_pix.m will be the first analysis project looking at the swimming kinematics of the fish.
