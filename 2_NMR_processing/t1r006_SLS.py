## Script for processing T1rho "Short-Long-Short" (10ms-200ms-10ms) data.
## Input:
## - 3*n expno arguments as input - corresponding to T1r_LONG spectra
##   (for Target, Target+Ligand, free Ligand).

### MAY NEED TO SET REFERENCE PEAK SEARCH WIDTH to 0.5 ppm - in EDLOCK table!
### Turned it down to 0.1 - for IVTNMR datasets

# test

## Requirements:
## - each T1r_LONG spectrum is expected to be sandwiched by two T1r_SHORT spectra.
## - an empty (unused) expno after the last T1r_SHORT - where the SUM of two T1r_SHORT
##   spectra (== FSI reference) will be stored. If expno already exists - it will be overwritten!

## Output:
## - The final T1r_LONG PM-P-M difference is saved as procno #4 in the Target+Ligand T1r_LONG expno.
## - The final T1r_SHORT PM-P difference is saved as procno #2 in the Target+Ligand T1r_SHORTsum expno. 
## - The two resulting spectra (main and FSI reference) are saved into ./spectra subfolder
##   (file names are T1r_prot_mix.txt and T1r10PMP_prot_mix.txt) respectively.
## - Also a T1rLONG spectrum of free metabolite is saved as secondary reference (T1r200freeM_PROT_MIX.txt) -
##   this is useful in cases when compound peaks are significantly shifted in presence of protein
##   (this is often manifested by FSI value > 1).


## Options of the script:
#=========================
## -n, --nosref     Do not use SREF during processing (if used w/o -k, the SR is reset to 0)
## -k, --keepsr     Keep SR as in the data, and do not use SREF during processing (overrides -n)
## -q, --quiet      Drop any user-confirmation requests
## -s, --separate   Export T1r_PMP and T1r_M spectra as separate files
## -t, --sep2       Export all 10ms/200ms and short1,short2 as separate files (adds T1r200PMP_, T1r10freeM, T1rS1/S2_)

## Some notes / observations:
# - Once got weird errors while processing - potentially because switched to another program while processing was going on. Some quirks in Java memory management when the app is out of focus? Better keep focus on TopSpin window while processing is going on. 


### TODO
#=========================
# - Make option -a --auto: automatically figure out the list of expts to process
#  (like in matlab - check which expts were indeed recorded - in specified range)
# - on below two points - see more in "NMR Automation for screens" - Task 06
# - Change the naming of export: PMexpno_T1r_PMPM, PMexpno_T1r_short_PMP, PMexpno_T1r_long_M
# -------
# - Check that all expts from the triplet have data in them - and pass if not!
# - decide whether to fix ZOOM?!

## Versions
#=========================
# v006_SLS - 2016-09-06. Notes are in "NMR Automation for Screens".
#    - Branched off for 0-200-0 (short-long-short) T1rho High-Throughput analysis
#    - Uses APK instead of APKS, APKM
#    - Added parsing of OPTIONS - nosref, quiet, ..
#    - Solved issue with setting of SR (need to do both SF = BF1 and XCMD("sr "+str(0)))
#    - Saves the input OPTIONS as string, and passes them during recursive calls
#    - Added --keepsr option
#    - Added ph180 (negative multiply) correction for LONG-diff
#    - Added SSL "triple-diff" generation.
#    - Automatic extraction of mix and protein names
#    - Save T1r_long_freeM as txt - for referencing of FSI_freeM
# 2017-03-27
#    - Turned logging back on.
#    - Option -s --separate: to export T1r_PMP and T1r_M spectra separately
#    - Moved T1rLONG txt export to AFTER it is processed! (on the first pass won't export otherwise!)
#    - 2017-04-05: turned export of T1r_LONG on again
#    - 2017-05-02: added export of 10ms and 200ms as separate files
#    - 2017-07-04: Adds short1/short2 export for M and PM
#    - 2017-08-28: Sets short1/short2 export under -sep2 flag (was not under option flag before)

# v006 - 2016-09-05. Notes are in "NMR Automation for Screens".
#    - Functionalized the code. When functions modify shared variables (e.g. PHC0) - explicitly assigned those as GLOBAL.

# v005 - 2016-08-15. Notes are in "NMR Automation for Screens".
#    - Taking ANY expnos as argument (P,PM,M)
#    - Saves to spectra/*.txt automatically (T1r-short PMP as the reference for FSA).
#    - Can take in 3*n expnos (processing multiple datasets)
# v004 - 2016-08-11. Notes are in "NMR Automation for Screens" - debug and decisions
#    on which processing params to use. (Edits: overwrite data, 
#    use XCPR (so that the script actually waits fro exec!), Overwrite in duadd, Do automatic phase
#    Use ADDFID for short-long, Use ONE phase for all expts!
# v003 - 2016-02-09. added 0th-order baseline correction (for many spectra baseline was slightly negative...
#    which should not be a problem when comparing RELATIVE intensities with other spectra,...
#    but for ABSOLUTE ones - it better be zero.
# v002 - added option for 240, 440, 540 (was only 1xx and 3xx allowed).
#====================================================

# @Yar.Nikolaev, 2014-2018

from TopCmds import *
import math
import logging as log
import ntpath
import os
import time
import re # regex to get sample name from title

from optparse import OptionParser
## OPTIONS - parsed options
## INPUT_OPTIONS_STRING - options before parsing
global OPTIONS, INPUT_OPTIONS_STRING

start_time = time.time()

######################################
## Initial checks & global variables
######################################
if CURDATA() == None:
    MSG("Please open a 1D NMR dataset!")
    EXIT()

[NAME, expno, procno, CURDIR] = CURDATA()
user = ''

PHC0_default = 270
PHC1_default = 0

LB = 2
SI = '128k'
BCMOD = 'qfil' # qpol / qfil (BCFW=0.2 normally - set for qfil is below)
APK_setting = 'apk' # apk apks apkm
NEG_PROCNO = '10' # procno where negative version of the spectrum will be stored (for subtractions)

## Procnos to store t1r subtractions
PM_P_PROCNO = '2'
PM_M_PROCNO = '3'
PM_P_M_PROCNO = '4'


##############################################
## Functions
##############################################
def start_option_parser():
    parser = OptionParser()
    parser.add_option("-n", "--nosref", dest="sref", action="store_false", default=True,
    	help="Set SR=0 and do not use SREF during processing")
    parser.add_option("-k", "--keepsr", dest="keepsr", action="store_true", default=False,
    	help="Keep existing SR, and do not use SREF")
    parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False,
    	help="Drop any user-confirmation requests")
    parser.add_option("-s", "--separate", dest="separate", action="store_true", default=False,
        help="Export T1r_PMP and T1r_M spectra as separate files")
    parser.add_option("-t", "--sep2", dest="sep2", action="store_true", default=False,
        help="Export 10ms and 200ms spectra as separate files (both for M and PM)")
    	
    return parser

def start_log():
	#### - uncomment such lines below if want to log into a file
	[NAME, expno, procno, CURDIR] = CURDATA()
	log_file = CURDIR+'/'+NAME+'/'+'processing.log'
	log.basicConfig(
	    filename=log_file,
	    #filemode='w', # overwrite file contents instead of appending (default mode is 'a')
        level=log.DEBUG, # minimal level of severity to keep in log: DEBUG INFO WARNING ERROR CRITICAL
        ## if want to debug something specific - use this level + log.warning('...') in the code
        # level=log.WARNING,
		format='%(asctime)s %(levelname)-10s %(message)s',
		#datefmt='%Y%m%d_%H%M%S' # format w/o dashes, but cant easily add millisec
		)
    
		## Add log display to the console
	log.getLogger().addHandler(log.StreamHandler())

		# log example:
		# log.debug('Copied experiment #'+str(i+1)+"\n"+'  '.join(target))

def directory_check(NAME):
    if CONFIRM("Warning",
        '<html>The automation script will execute in: <br><br>'+
        '<font color="red">'+NAME+'</font>'+
        '<br><br>Processed data will be overwritten! <br> Continue?!</html>') == 0:
        EXIT()


def check_n_split_args(args):
    n_required_args = 3
    
    log.debug("len(args)="+str(len(args)))
    log.debug(args)

    if len(args) == n_required_args:
        t1r = [int(i) for i in args] # defined in order: target, target+ligand, ligand
        log.debug("t1r:")
        log.debug(t1r)
        
        return t1r

    elif not len(args) % n_required_args:
        for i in range(0, len(args), n_required_args):
            log.debug("xpy %s %s %s %s %s" % (ntpath.basename(sys.argv[0]), str(args[i]), str(args[i+1]), str(args[i+2]), INPUT_OPTIONS_STRING))
            XCMD("xpy %s %s %s %s %s" % (ntpath.basename(sys.argv[0]), str(args[i]), str(args[i+1]), str(args[i+2]), INPUT_OPTIONS_STRING), wait = WAIT_TILL_DONE) # XCPR would not understand this
        EXIT()
    
    else:
        ERRMSG("Error: Script expects 3 (or 3*n) arguments - expnos for 1.Target; 2.Ligand+Target; 3.Ligand")
        EXIT()


def read_expt(e, p=1, n=NAME, dir=CURDIR):
    ''' Wrapper to read expts - assumes default values for p, n, dir
    '''
    fullpath = [n, str(e), str(p), dir]
    RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")
    #log.debug('RE experiment: \n'+'  '.join(fullpath)+'\n')


def zoom1D(left,right):
    ''' To zoom specific region for export. Currently not used.
    '''
    log.info('Zoom called')
    fullrange = putil.DataChecks.getNMRDataOfSelectedFramePrintMsg().getFullPhysicalRange()
    print(fullrange)
    newRange = fullrange
    newRange[0].setStart(float( left ))
    newRange[0].setEnd(float( right ))
    newRange[0].setUnit("ppm")
    return newRange


def copy_expt(exp_source,exp_target):
    CLOSEWIN()
    read_expt(exp_source)
    log.debug('Copy expts on:')
    log.debug(CURDATA())

    ## Make a copy of expt - to store new FID
    RE([NAME, str(exp_source), '1', CURDIR])
    WR([NAME, str(exp_target), '1', CURDIR], "y") # overwrite
    RE([NAME, str(exp_target), '1', CURDIR])


def add_fids(exp1,exp2,DC):
    CLOSEWIN()
    read_expt(exp1)

    ## ADDFID - adds RAW FIDs, writes into CURRENT one:
    ## EDC1(current) = EDC2 + EDC3*DC
    PUTPAR("DC", str(DC))
    SET_CURDATA2(CURDIR, user, NAME, str(exp1), '1') # EDC2
    SET_CURDATA3(CURDIR, user, NAME, str(exp2), '1') # EDC3
    XCPR('addfid y',wait = WAIT_TILL_DONE) # overwrites the data    


def duadd(e1p1,e2p2,e3p3,*add):
    ''' Adds spectra ppm/Hz-wise.
        EDC3 = EDC1 + EDC2
        # DUADD seems to work like EDC3 = EDC1+EDC2 (result stored in EDC3)
        # DUADD seems to NOT TAKE into account the DC parameter.
        # "Work-around" is to create negative versions of spectra using NM()).
    '''
    
    # For SET_CURDATA2/SET_CURDATA3 to work properly - need to use
    # CLOSEWIN() function before reading datasets!

    # If adding m ADD (using cmd line) - writes into current dataset (EDC1 = EDC2+EDC3*DC).
    # Not sure if DC is taken into account by the ADD command when run from cmd line.

    # Syntax for DATASET2/3 definition:
    # see /opt/topspin3.2/classes/lib/topspin_py/py/pycmd/TopCmds.py:
    # SET_CURDATA2(dir2, user2, name2, expno2, procno2)
    # SET_CURDATA3(dir3, user3, name3, expno3, procno3)
    # SET_CURDATA3 did not exist in the default TopCmds.py
    # > had to manually add it by duplicating SET_CURDATA2.
        
    ## Unpack arguments
    e1,p1 = e1p1
    e2,p2 = e2p2
    e3,p3 = e3p3
    
    CLOSEWIN()
    log.info('Adding spectra ppm-wise (duadd): \
        e1/p1=%s/%s, e2/p2=%s/%s, e3/p3=%s/%s' % \
        (str(e1),str(p1),str(e2),str(p2),str(e3),str(p3)))
    
    read_expt(e1,p1) # EDC1

    PUTPAR("DC", str(1)) # just in case if tests were wrong and DC param IS USED by duadd
    
    log.debug('Have read dataset:')
    log.debug(CURDATA())
    SET_CURDATA2(CURDIR, user, NAME, str(e2), str(p2)) # WHAT to add
    SET_CURDATA3(CURDIR, user, NAME, str(e3), str(p3)) # WHERE to store
    ## 'duadd y' - adds Hz/ppm-wise, 'add y' - adds point-wise (spectra need be same size)
    if not add:
        XCPR('duadd y',wait = WAIT_TILL_DONE) # see header for notes on duadd
    else:
        XCPR('add y',wait = WAIT_TILL_DONE) # see header for notes on duadd

def process(expno,*phases):
    '''
    Processes spectrum using specified phases.
    If no phase specified - phases automatically & returns the phase values.
    '''
    global OPTIONS    
    
    if phases:
        if len(phases)==2:
            PHC0, PHC1 = phases
        else:
            ERRMSG("process() requires 1 or 3 params: expno, or expno + two phases (PHC0, PHC1)")
            EXIT()
    else:
        PHC0, PHC1 = PHC0_default, PHC1_default

    read_expt(expno)
    log.warning('expno= %s, procno = ..' % str(expno))
    # log.warning('SF / SR (before PUTPAR)= %s / %s' % (GETPAR('SF'),GETPAR('SR')))

    if not OPTIONS.keepsr:
        PUTPAR("SF", GETPAR("BF1")) # resets SR to zero
        XCMD('sr '+str(0), wait = WAIT_TILL_DONE) # resets another SR to zero    
    log.warning('SF / SR (after PUTPAR)= %s / %s' % (GETPAR('SF'),GETPAR('SR')))

    # log.warning('SF / SR (before SI)= %s / %s' % (GETPAR('SF'),GETPAR('SR')))
    PUTPAR("SI", str(SI))
    # log.warning('SF / SR (after SI)= %s / %s' % (GETPAR('SF'),GETPAR('SR')))
    PUTPAR("PHC0", str(PHC0))
    PUTPAR("PHC1", str(PHC1))

    PUTPAR("BCFW", str(0.2)) # in case if using BC_mod qfil
    PUTPAR("BC_mod", BCMOD)

    PUTPAR("WDW", "EM")
    PUTPAR("LB", str(LB))

    EFP()
    
    ## Determine and return phases if these were not provided by user
    if not phases:        
        # In initial T1r - APK() and APKS() give similar results!
        # wL and 0-200-0 - APK() looked better!
        if APK_setting == 'apk':
            APK()
        elif APK_setting == 'apks':
            APKS()
        elif APK_setting == 'apkm':
            APKM()
                                    
        PHC0 = GETPAR("PHC0")
        PHC1 = GETPAR("PHC1")
        return PHC0, PHC1

    ## Alternatively - actual processing: incl baseline, sref and saving a negative.
    ## TODO: this is perhaps not optimal to put all in same function, but for
    ## the moment keeping as it was - to avoid having to re-read expnos for each action.
    else:
        # v003 - add 0th-order baseline corr
        PUTPAR("ABSG", str(0))        
        XCMD('abs n', wait = WAIT_TILL_DONE) # if want direct ABSN() - add it to TopCmds.py
        
        if OPTIONS.sref and not OPTIONS.keepsr: SREF()
        log.warning('SF / SR (after SREF)= %s / %s' % (GETPAR('SF'),GETPAR('SR')))
        
        if OPTIONS.sref and not OPTIONS.keepsr: log.warning('--- SREF WAS DONE')
        else: log.warning('--- SREF NOT DONE')
        
        # Create a negative version of the spectra for subtraction operations:
        WR([NAME, str(expno), NEG_PROCNO, CURDIR], "y")
        RE([NAME, str(expno), NEG_PROCNO, CURDIR])
        XCMD('nm', wait = WAIT_TILL_DONE) # if want direct NM() - add it to TopCmds.py

def save_as_txt(expno, procno, prefix):
    """ Save results as TXT in /spectra subfolder """
    ## Create folder for spectra if its not there yet:
    spectra_path = os.path.join(CURDIR,NAME,'spectra')
    if not os.path.exists(spectra_path):
    	os.makedirs(spectra_path)
    	
    read_expt(expno, procno)
    # save_path = os.path.join(spectra_path, prefix + str(expno) + '.txt')
    save_path = os.path.join(spectra_path, prefix + '.txt')
    XCMD('totxt ' + save_path, wait = WAIT_TILL_DONE) # seems XCPR d n understand this command!
    
def get_sample_name(expno):
    filepath = os.path.join(CURDIR,NAME,str(expno),'pdata/1/title') 
    textfile = open(filepath, 'r')
    # filetext = textfile.read() # read full file
    first_line = textfile.readline() # read first line
    textfile.close()

    name = re.sub('[ *]', '', first_line) # remove stars and spaces
    return name
    
        
###############################
## Main execution
###############################
def main():
    start_log()

    parser = start_option_parser()
    global OPTIONS, INPUT_OPTIONS_STRING
    
    ## For recursive calls:
    ## Filters sys.argv[1:] through an anonymous function (lambda),
    ## leaving only the items which have dash as their first symbol.
    ## The filtered items are joined with spaces again - so can pass in recursive calls.
    INPUT_OPTIONS_STRING = ' '.join( filter(lambda x: x[0] == '-', sys.argv[1:]) )
    
    (OPTIONS, args) = parser.parse_args()
    log.info('OPTIONS:')
    log.info(OPTIONS)
    log.info('args:')
    log.info(args)
    
    ## Get experiment numbers
    t1r = check_n_split_args(args)
    t1rSHORT1 = [x-1 for x in t1r]
    t1rSHORT2 = [x+1 for x in t1r]
    t1rSHORT_sum = [x+2 for x in t1r]
    t1rAll = t1r + t1rSHORT_sum
    
    t1rSSL = [x+1 for x in t1rSHORT_sum]
    
    if OPTIONS.quiet:
        pass
    else:
        directory_check(NAME)    

    ## Get protein and mix names
    #=================================
    # 0,1,2 - P, PM, M
    name_p = get_sample_name(str(t1r[0])[:-1]+'0')
    name_m = get_sample_name(str(t1r[2])[:-1]+'0')
    log.info('Names extracted: Protein = %s, Mix = %s' % (name_p, name_m))
    
    log.info('Creating sum of two T1r_short expts ...')
    for e in t1rSHORT1:
        copy_expt(e,e+3) # copies to the expno position of t1rSHORT_sum
        add_fids(e+3,e+2,1) # current(is read) + DC*second


    ## Processing and exporting PM-P (SHORT) - for FSI
    #=================================
    log.info('Determining phase for free ligand in T1r_short ...')
    PHC0, PHC1 = process(t1rSHORT_sum[2])

    log.info('Processing T1r_SHORT ...')
    for exp in t1rSHORT_sum:
        process(exp, PHC0, PHC1)

    ## Calc the difference PM-P
    log.info('PM-P reference (max Ligand intensity) for FSI calc ...')
    ## EDC3 = EDC1 + EDC2
    duadd([t1rSHORT_sum[1],'1'], [t1rSHORT_sum[0],NEG_PROCNO], [t1rSHORT_sum[1],PM_P_PROCNO])

    ## Save results        
    # save_as_txt(t1rSHORT_sum[1], PM_P_PROCNO, 'T1rSHORT_PMP_') ## SHORT - PM-P reference

    ## Save PMP reference
    save_as_txt(t1rSHORT_sum[1], PM_P_PROCNO, ('T1r10PMP_%s_%s' % (name_p, name_m)))


    ## Processing and exporting SHORT1/2 PM and M spectra - for compound stability checks
    #=================================
    for exp in t1rSHORT1:
        log.info('Processing Short1 and Short2 with same phase as Short_sum ...')
        process(exp, PHC0, PHC1)
        process(exp+2, PHC0, PHC1)
    
    if OPTIONS.sep2:
        # Save PM (index [1])
        save_as_txt(t1rSHORT1[1], '1', ('T1rS1_PM_%s_%s' % (name_p, name_m)))
        save_as_txt(t1rSHORT2[1], '1', ('T1rS2_PM_%s_%s' % (name_p, name_m)))
        # Save M (index [2])
        save_as_txt(t1rSHORT1[2], '1', ('T1rS1_M_%s_%s' % (name_p, name_m)))
        save_as_txt(t1rSHORT2[2], '1', ('T1rS2_M_%s_%s' % (name_p, name_m)))

    ## Processing and exporting PM-P-M (LONG)
    #=================================
    log.info('Determining phase for free ligand in T1r LONG ...')
    PHC0, PHC1 = process(t1r[2])

    log.info('Processing T1r LONG ...')
    for exp in t1r:
        process(exp, PHC0, PHC1)

    ## Make difference PM-P-M
    log.info('PM-P-M differences for T1r LONG ...')
    ## EDC3 = EDC1 + EDC2
    duadd([t1r[1],'1'], [t1r[0],NEG_PROCNO], [t1r[1],PM_P_PROCNO]) # PM-P
    duadd([t1r[1],'1'], [t1r[2],NEG_PROCNO], [t1r[1],PM_M_PROCNO]) # PM-M
    duadd([t1r[1],PM_P_PROCNO], [t1r[2],NEG_PROCNO], [t1r[1],PM_P_M_PROCNO]) # PM-P-M

    ## Invert the phase for LONG-diff
    RE([NAME, str(t1r[1]), PM_P_M_PROCNO, CURDIR])
    XCMD('nm', wait = WAIT_TILL_DONE) # if want direct NM() - add it to TopCmds.py
    
    ## Was not saving after checking compound stability bug - 
    ## decided using pure LONG-difference is not optimal!
    ## 2017-04-05 - returned saving again, to do QC control of compound stability
    ## > unstable compound should show as negative peaks in L-diff
    ## Save LONG PM-P-M
    save_as_txt(t1r[1], PM_P_M_PROCNO, ('T1r_LONG_%s_%s' % (name_p, name_m)))

    #### # Save T1r_long_freeM
    ### save_as_txt(t1r[2], '1', 'T1r_')

    ## Save freeM - T1rLONG (for FSI freeM)
    save_as_txt(t1r[2], '1', ('T1r200freeM_%s_%s' % (name_p, name_m)))


    ## Processing and exporting S+S-L
    #=================================
    log.info('Copy S+S, and make S+S-L diff ...')
    for e in t1rSSL:
        copy_expt(e-1,e)
        ## EDC3 = EDC1 + EDC2.
        duadd([e,'1'], [e-3,NEG_PROCNO], [e,'1']) # overwrite initial spectrum
        # duadd([e,'1'], [e-3,NEG_PROCNO], [e,'1'], 'true') # USE ADD option 'true'
        ## Create a negative version of the spectra for subtraction operations:
        RE([NAME, str(e), '1', CURDIR])
        WR([NAME, str(e), NEG_PROCNO, CURDIR], "y")
        RE([NAME, str(e), NEG_PROCNO, CURDIR])
        XCMD('nm', wait = WAIT_TILL_DONE) # if want direct NM() - add it to TopCmds.py

    ## Make difference PM-P-M
    log.info('Make PM-P-M differences for SSL ...')
    ## EDC3 = EDC1 + EDC2
    duadd([t1rSSL[1],'1'], [t1rSSL[0],NEG_PROCNO], [t1rSSL[1],PM_P_PROCNO]) # PM-P
    duadd([t1rSSL[1],'1'], [t1rSSL[2],NEG_PROCNO], [t1rSSL[1],PM_M_PROCNO]) # PM-M
    duadd([t1rSSL[1],PM_P_PROCNO], [t1rSSL[2],NEG_PROCNO], [t1rSSL[1],PM_P_M_PROCNO]) # PM-P-M

    ## Save main expt
    save_as_txt(t1rSSL[1], PM_P_M_PROCNO, ('T1r_%s_%s' % (name_p, name_m)))

    ## Save T1r_PM-P and T1r_M separately (for separate peak/maxima picking)
    if OPTIONS.separate:
        save_as_txt(t1rSSL[1], PM_P_PROCNO, ('T1r_%s_%s_PMP' % (name_p, name_m))) # PM-P
        save_as_txt(t1rSSL[2], '1', ('T1r_%s_%s_M' % (name_p, name_m))) # just M spectrum


    # ## Add export of remaining 10ms and 200ms as separate files (2017-05-02)
    if OPTIONS.sep2:
        # M_short (M_10)
        # save_as_txt(t1r[2], '1', ('T1r200freeM_%s_%s' % (name_p, name_m)))
        save_as_txt(t1rSHORT_sum[2], '1', ('T1r10freeM_%s_%s' % (name_p, name_m)))
        
        # PMP_long (PMP_200)
        # save_as_txt(t1rSHORT_sum[1], PM_P_PROCNO, ('T1r10PMP_%s_%s' % (name_p, name_m)))
        save_as_txt(t1r[1], PM_P_PROCNO, ('T1r200PMP_%s_%s' % (name_p, name_m)))


    ## Script complete
    log.info('\n T1r processing of %s completed in %.4f sec.' % (str(NAME), time.time()-start_time))
    
#####################
## This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
		main()