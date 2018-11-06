## Open/create file with sample list for automation runs.
## Needs to be run from an *experiment* window in TopSpin
## (i.e. will not work if called from edte/stdisp/etc)
## @Yar 2016

# import routines to run shell commands
from subprocess import call
import os.path

# get info about current dataset
NAME, expno, procno, CURDIR = CURDATA()

# create path to file, using dataset name and TopSpin data folder
sampleList = CURDIR+"/"+NAME+"/sample_list.py"

# If the file doesn't exist yet - create one
if not(os.path.isfile(sampleList)):
	
	if CONFIRM("Warning",
		'<html>About to create sample list in this folder: <br><br>'+
		'<font color="red">'+NAME+'</font>'+
		'<br><br> Are you sure its the TARGET, not the template dataset?! <br> Continue?!</html>') == 0:
		EXIT()
	
	configPlaceholders = (
	'templateSet = "XXX"\n\n'
	'solvent = "H2O+D2O"\n\n'
	'## Current scripts (e.g. nc_v008.yn.py) assume first expt is sample title, second - calibrations (zg.ethz)\n'
	'sourceExpts = [10,11]\n'
	'\n'
	'## If samples are incremented with 10 (e.g. 10, 20, 30) - use 10.\n
	'## If difference between samples is >=100 - set ZERO!\n'
	'subtract_to_set_expt_numbers = 10\n\n'
	'## If want to use only a subset of sourceExpts for certain samples, add such sublist after comma after sample name, like:\n'
	'## \t(100, "ABC"),\n'
	'## \t(200, "XYZ", [10,11,17]),\n'
	
	'expset = [\n\t(00, "XYZ"),\n]'
	)
	
	f = open(sampleList, 'w')
	f.write( configPlaceholders )
	f.close()

# open text editor with the file
call(["gedit", sampleList])

#print("test message to TopSpin console")
#MSG(datasetNotes)
