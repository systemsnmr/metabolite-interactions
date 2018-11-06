# nmrcase.py - v.20181106
# Automation routines for operation of NMR Case (potentially other sample changers which switch samples by ej/ij command).

# Based on nc_v008c.yn.py
# Development version!! May not function properly. Needs testing on individual machines before execution.
# Written by Yar Nikoaev at ETH Zurich. 

# We assume no responsibility, whatsoever, and shall not be liable for any problems, damages or losses that may result from running this code on your computers and/or spectrometers.

# Known issues:
# - To stop execution of this script in TopSpin - need to restart TopSpin itself (did not find way to terminate the script w/o restarting TopSpin)
# - Sometimes spectrometer generates errors during topshim, ATMA/wobb or during lock.

# --------------------------------------------------------------------------------

## Helper scripts
## Create/edit sample list - slist.py
## View automation log - autolog.py
# In ANY script - need to have CALIBRATION/ref expt first!! (to avoid over-writing!)

from __future__ import division # to avoid errors in float division
from TopCmds import *
from subprocess import call
import os.path
from os import remove
import math
import logging as log

### "TOP-LEVEL" variables (technically NOT GLOBAL)
DEBUG = 0 # For code debugging - set to zero for normal runs.

NAME, expno, procno, CURDIR = CURDATA()

### Config of logging

log_file = CURDIR+'/'+NAME+'/'+'automation.log'

log.basicConfig(
	filename=log_file,
	filemode='a', # overwrite file contents instead of appending
	level=log.DEBUG, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
	format='%(asctime)s %(levelname)-10s %(message)s',
	#datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
	)

## Add log display to the console
log.getLogger().addHandler(log.StreamHandler())

#####################################

### Use this header to debug small things
#if DEBUG:

#EXIT()

#####################################
def copyTemplate(sExpts,tExpts,sSet,tSet):
	""" Copies (sources expts) to (target expts) from (template set) to (target set).
	"""
	# TODO - check that target set exists!
	for i in range(0,len(sExpts)): 
		#print(sourceExpts[i])
		source = [sSet, str( sExpts[i] ), "1", CURDIR]
		target = [tSet, str( tExpts[i] ), "1", CURDIR]

		RE(source, show="n") # read the dataset w/o showing it. Syntax/defaults: RE(dataset = None, show = "y")

		# KEEP IN MIND - WR copies also the data!! (not only the parameters)
		# Unlike other functions, WR() gives error if override parameter in not just <"n"> but <override="n">
		# TODO - check if can add WRAPARAM in TopCmds?
		#	WR(target, "n") # write/copy the current dataset. Ask to override if exists. Syntax/defaults: WR(Dataset = None, override = "y") 
		WR(target) # overwrites by default
		# Example of doing ternary if in one line:
		# 'Yes' if fruit == 'Apple' else 'No'
		log.debug('Copied experiment #'+str(i+1)+"\n"+'  '.join(target))

	log.info('Finished copying expts\n')


def readExpt(n,e,p,dir):
	fullpath = [n, str(e), str(p), dir]
	RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")
	log.debug('RE experiment: \n'+'  '.join(fullpath)+'\n')
	
	
def calibrateP1():
	log.debug('Starting P1 calibration...')
	
	exptDB = float( GETPAR('PLdB 1') )
	log.debug("exptDB = %.2f" % exptDB)
	exptW = 10**(exptDB/-10.0) # refDB = -10*math.log10(refW)
	log.debug("exptW = %.4f" % exptW)
	refDB = -9.0
	refW = 10**(refDB/-10.0) # refDB = -10*math.log10(refW)
	PUTPAR('PLW 1', str(refW))
	PUTPAR('P 1', str(5)) # set p1 to low-us, otherwise pulsecal may fail!
	log.debug('refW for pulsecal = %s' % GETPAR("PLW 1"))
	log.debug('Starting P1 for pulsecal = %s' % GETPAR("P 1"))
	XCPR("xau pulsecal quiet same", WAIT_TILL_DONE)
	P90H = float( GETPAR("P 1") )
	log.debug('P1 @ %.4f W = %.2f' % (refW, P90H))
	P90H = math.sqrt(refW / exptW) * P90H
	log.debug('P1 @ %.4f W = %.2f' % (exptW, P90H))
	PUTPAR('PLW 1', str(exptW))
	PUTPAR('P 1', str(P90H))
	
	log.info('Finished calibration, P1 = %.2f' % P90H)
	return P90H

def calibrateO1():
	log.debug('Starting O1 calibration...')
#	XCPR("xau o1calib", WAIT_TILL_DONE)
	XCPR("o1calib_silent.yn.au", WAIT_TILL_DONE)
	H_Hz = GETPAR("O1")
	log.info('Finished calibration, H_Hz = %s' % str(H_Hz))
	return H_Hz


def lockTuneShim(expno, solvent):
	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around
	
	# TODO: add options for EXACT and tunebxyz

	## LOCKING
	lock_string = "lock "+solvent
	lock_result = -1
	
	lock_attempt = 0
	while lock_result == -1 and lock_attempt < 20:
		log.info('Locking '+lock_string+' ...')
		result = XCPR(lock_string, WAIT_TILL_DONE)
		lock_result = result.getResult()
		log.debug('lock_attempt / lock_result = %s / %s' % (str(lock_attempt), str(lock_result)))
		lock_attempt += 1

	## ATMA
	RE(target, show="y") # in case user has switched elsewhere
	log.info('ATMA...')
	result = XCPR("atma", WAIT_TILL_DONE)
#	result = XCPR("atma exact", WAIT_TILL_DONE)


	log.debug('autophase')
	XCPR("autophase", WAIT_TILL_DONE)
	log.debug('autogain')
	XCPR("autogain", WAIT_TILL_DONE)

	log.info('Shimming...')
	RE(target, show="y") # in case user has switched elsewhere
#	XCPR("topshim ls tunebxyz", WAIT_TILL_DONE)
#	XCPR("topshim ls tuneb", WAIT_TILL_DONE)
	XCPR("topshim ls", WAIT_TILL_DONE)
	log.debug('autogain')
	XCPR("autogain", WAIT_TILL_DONE)
### If want to use loopadj - make a SILENT version of it!
### But better d n use it - It optimises for best long-term stability, 
# but not for best lineshape, resolution or homogeneity.
#	log.debug('loopadj')
#	XCPR("loopadj", WAIT_TILL_DONE)
	log.info('lockTuneShim() finished')
	# TODO - set AUTOSHIM ON somehow

																		
def recO1andH2Oshim(expno,P90H,H_Hz):
	"""
	Records H2O lineshape - after tap and 360 pulses
	(for shim and o1 checks respectively).
	"""

	target = [NAME,str(expno),"1",CURDIR] # save so that can re-read if ppl are switching around

	RE(target, show="y") # in case user has switched elsewhere
		
	PUTPAR('O1', str(H_Hz))
	PUTPAR('P 1', '0.02')
	log.debug('--- Set P1 0.02 on:')
	log.debug(CURDATA())
	try:
		os.remove(os.path.join(CURDIR,NAME,str(expno),'fid'))
	except OSError:
		pass

	# TODO JUST USE ZG(), FT(), PK() here?!
	XCPR('zg', WAIT_TILL_DONE)

	RE(target, show="y") # REREAD - in case user has switched elsewhere
	XCPR('fp', WAIT_TILL_DONE) # combined 'qu zgfp' doesn't work!
	XCPR('apk0', WAIT_TILL_DONE)
	WR([NAME,str(expno),"2",CURDIR]) # overwrites procno 2 (complained if was using 
	
	XCPR('rep 1', WAIT_TILL_DONE) # save for quality checks of shimming
	log.info('Recorded and phased water with tip pulse. Starting 360..')

	P360H = float(P90H)*4
	PUTPAR('P 1', str(P360H)) 
	log.debug('--- Set P1 (360deg)='+str(P360H)+' on dset:')
	log.debug(CURDATA())
	log.debug('Deleted the fid (otherwise may ask user for confirmation):')
	log.debug(os.path.join(CURDIR,NAME,str(expno),'fid'))
	try:
		os.remove(os.path.join(CURDIR,NAME,str(expno),'fid'))
	except OSError:
		pass

	XCPR('zg', WAIT_TILL_DONE)
	RE(target, show="y") # REREAD - in case user has switched elsewhere
	XCPR('fp', WAIT_TILL_DONE)
	log.info('recO1andH2Oshim() finished')
	
	
def recExp(expno,P90H,H_Hz):

	# TODO:
	# RPAR ...
	# If TITLE - just set title and return
	# Else:
  #  - set o1, p1
	#  - run additional commands (if any)
	#  - run zg(), etc

	log.debug('Setting up expno = '+str(expno))
	target = [NAME,str(expno),"1",CURDIR]
	RE(target, show="y") # read the dataset w/o showing it. Syntax/defaults: RE(dataset = None, show = "y")
	PUTPAR('P 1', str(P90H))
	log.debug('--- Set P1 (90deg)='+str(P90H)+' on dset:')
	log.debug(CURDATA())
	PUTPAR('O1', str(H_Hz))
	log.info('Recording experiment # '+str(expno))
	XCPR('zg', WAIT_TILL_DONE)
	
	### TODO - process only if DIM = 1D
	if 0:
		XCPR('efp', WAIT_TILL_DONE)
		XCPR('apk0', WAIT_TILL_DONE)
		
	RE(target, show="y") # read dataset, just in case if users are switching around during run
	log.info('Finished recExp() ')


def nextSample():
	#### based on nc.yn.py
	log.info('Ejecting sample ..')
	XCPR("ej", WAIT_TILL_DONE)
	SLEEP(25)   # @M.O.Ebert used 20s.
	log.info('Injecting new sample ..')
	XCPR("ij", WAIT_TILL_DONE)
	log.info('Temperature equilibration (30s) ..')
	SLEEP(30)
	log.info('Finished nextSample()')


def set_title(dataset,text):
	outpath = dataset[3]+'/'+dataset[0]+'/'+dataset[1]+'/pdata/'+dataset[2]
	f = open(outpath+'/title','w')
	f.write(text)
	f.close()


def get_timetag():
	timevar = time.localtime()
	timetag = '%d%02d%02d_%02d%02d'%(timevar[0],timevar[1],timevar[2],timevar[3],timevar[4])
	return timetag


def read_sample_list():
	dirpath = CURDIR+'/'+NAME
	sys.path.append(dirpath)
	#print('dirpath='+dirpath)
	try:
		from sample_list import solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset
	except:
		MSG("Sample list (sample_list.py) not found or has errors")
		EXIT()
	log.debug('Solvent = '+solvent)
	log.debug('TemplateSet = '+templateSet)
	log.debug('sourceExpts = '+str(sourceExpts))
	log.debug('expset = '+str(expset))
	
	return solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset

####################################
# Actual run commands
def main():
	targetSet = NAME
	
	log.info('== Logging (Re)start for new run / expset')
	
	if len(sys.argv) > 1 and sys.argv[1] == 'noconfirm':
		print("noconfirm option specified - recording w/o asking for user confirmation")
	else:
		if CONFIRM("Warning",
			'<html>The automation script will execute in: <br><br>'+
			'<font color="red">'+targetSet+'</font>'+
			'<br><br> This should be the TARGET dataset (any data will be overwritten)! <br> Continue?!</html>') == 0:
			EXIT()
	
	solvent, templateSet, sourceExpts, subtract_to_set_expt_numbers, expset = read_sample_list()
	
	## If multiple samples - check that increment size will not cause conflict
	## This assumes (when calc s_expno_deltas) that experiments are IN ORDER)
	## TODO: SORT numbers before calculating s_expno_deltas
	if len(expset) > 1:
		s_start_nums = [int(item[0]) for item in expset]
		log.debug("s_start_nums=" + str(s_start_nums))
		s_expno_deltas = [abs(s_start_nums[i+1]-s_start_nums[i]) for i in range(0,len(s_start_nums)-1)]
		log.debug("s_expno_deltas=" + str(s_expno_deltas))
		log.debug("len(sourceExpts)=" + str(len(sourceExpts)))
		log.debug("max(sourceExpts)-min(sourceExpts)=" + str(max(sourceExpts)-min(sourceExpts)))
		log.debug("min(s_expno_deltas)=" + str(min(s_expno_deltas)))
		if (len(sourceExpts) or max(sourceExpts)-min(sourceExpts)) > min(s_expno_deltas):
			MSG("ERROR: Min increment between starting sample EXPNOS needs to be \n"+
			" more than number of expts per sample")
			EXIT()
	
  # TODO(@YN): 
  # - Check experiment times (assuming first and second entries are title and calibr)

	# for our "standard" parameter namings - check: edau hsqc15N.all
	BF1 = GETPAR("BF1")
	P90H = 5
	H_Hz = 4.7*float(BF1);
	
	### Go over each SAMPLE, do calibrations (once), record all expts
	for i in range(len(expset)):

		### Inject new sample
		nextSample()
		
		### Copy expts and set title
		log.info('=== Sample # %s - %s:' % (str(i+1),expset[i][1]))
		tIncrement = expset[i][0]-subtract_to_set_expt_numbers
		
		### Enabling the use of additional settings:
		### 1. Use individual experiment list - based on 3rd element of sample definition
		sourceExptsFinal = sourceExpts
		### Loop if sample has additional elements (experiments, ...)
		if len(expset[i]) > 2:
			log.info(' Sample definition has extra elements (len(expset[i]) > 2)')
			### Check 3rd element
			if (type(expset[i][2]) is list) and expset[i][2] and all(isinstance(item, int) for item in expset[i][2]):
				log.info(' 3rd element is a non-empty list of integers - updating sourceExpts list')
				sourceExptsFinal = expset[i][2]
			else:
				log.info(' 3rd element is NOT a non-empty list of integers - using default sourceExpts')
			log.info('sourceExptsFinal = '+str(sourceExptsFinal))
			
		targetExpts = [x+tIncrement for x in sourceExptsFinal]

		copyTemplate(sourceExptsFinal,targetExpts,templateSet,targetSet)
		title = '*********  %s  *********' % expset[i][1]
		titleExpt = sourceExptsFinal[0]+tIncrement
		set_title([NAME,str(titleExpt),"1",CURDIR], title)
	
		### Do all calibrations in the first expt after title
		log.info('= Optimization and calibrations')
		calibrExpt = targetExpts[1]
		readExpt(targetSet,calibrExpt,"1",CURDIR)
		lockTuneShim(calibrExpt, solvent)

		## Read again, just in case user was switching around TopSpin
		readExpt(targetSet,calibrExpt,"1",CURDIR)
		## Run o1 calibration, catching and suppressing errors if such occur.
		H_Hz_calibrated = None
		while H_Hz_calibrated is None:
			try:
				H_Hz_calibrated = calibrateO1()
			except: # don't care whats the error here
				log.debug('Got exception during o1calib')
				pass
		
		H_Hz = H_Hz_calibrated
		log.debug('P90H(manual set)='+str(P90H))
	
		dBH = GETPAR("PLdB 1")
		# Read again, just in case user was switching around TopSpin
		readExpt(targetSet,calibrExpt,"1",CURDIR)
		P90H = calibrateP1()
		log.debug('P90H(returned)='+str(P90H))
		P90H = GETPAR("P 1")
		log.debug('P90H(after function)='+str(P90H))

		log.info('= Record o1/shim quality checks')
		recO1andH2Oshim(calibrExpt, P90H, H_Hz)
	
		### Run actual experiments
		for e in range(2,len(targetExpts)):
			recExp(targetExpts[e], P90H, H_Hz)
	
	log.info('Automation script finished.\n\n')

#EXIT()

### Own checks for getResult
### getResult returns DIFFERENT THINGS - depending on what 
### function is called. I.e. watchout if implementing this!
#	while result.getResult() < 0:
#		print('result in while loop = '+str(result.getResult()))
#		SLEEP(5)

#####################
# Fore debugging use constructions like:
#MSG(VARNAME) # displays in a pop-up window
#print("VARNAME = " + str(VAR))
#EXIT() # make a break-point

#####################
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
		main()
