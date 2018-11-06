# Open automation.log file
# Needs to be run from an *experiment* window in TopSpin
# (i.e. will not work if called from edte/stdisp/etc)
# @Yar 2016

# import routines to run shell commands
from subprocess import call
import os.path

# get info about current dataset
curdat = CURDATA()

# create path to file, using dataset name and TopSpin data folder
target_file = curdat[3]+"/"+curdat[0]+"/automation.log"

# open text editor with the file
call(["gedit", target_file])

#print("test message to TopSpin console")
#MSG(datasetNotes)