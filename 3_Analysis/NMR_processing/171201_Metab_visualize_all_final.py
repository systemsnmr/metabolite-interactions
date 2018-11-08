## Just overlay all Metab spectra.
## @Yar.Nikolaev. 2017-..

## TODO:
## - Make to take parameters for the MIX

## Made for Evn "Manual checks" of 2016 CCM hits - A02

# Try-catch to allow testing this script as normal Python too (outside of TopSPin)
try:
    from TopCmds import *
    # get info about current dataset
    NAME, expno, procno, CURDIR = CURDATA()
    #MSG(CURDIR)
    #EXIT()
    RUNNING_IN_TOPSPIN = 1
except:
    RUNNING_IN_TOPSPIN = 0
    pass

import sys

####################################
## Actual run commands
def main():

    ### Set here which mix you want to see:
    # ccm1-ccm2-ccm3-ccm4 - 10-20-30-40
    mix_decimals = [20]

    if RUNNING_IN_TOPSPIN:
        name, expno, procno, curdir = CURDATA()
        XCMD('.ret', wait = WAIT_TILL_DONE) # return from .md mode (if still in it)

    # Order of expts:
    # procnos:
    PM_P_PROCNO = '2'
    PM_M_PROCNO = '3'
    PM_P_M_PROCNO = '4'
    # expnos:
    # 12 - spinecho
    # 13 - S
    # 14 - L // procno 1,2,3,4 == raw, -P, -M, -PM
    # 15 - S
    # 16 - S+S
    # 17 - SLS difference
    procno_to_show = PM_P_M_PROCNO
    last_num_of_expno = 7
    
    dsets = [
    # '160915_Metab_TalB_GltA_ccm1-4_600-2',
    '161117_Metab_eno_BSA_ccm1-4_600-2',
    '161123_Metab_TktA_RpiA_ccm1-4_600-2',
    '161214_Metab_acs_cra_etc_ccm1-4_600-2',
    '170321_Metab_acea_fbp_gpma_etc_ccm1-4_600-2',
    '170516_Metab_mdh_pckA_pta_ccm1-4_600-2',
    '170628_Metab_pgi_pfka_etc_ccm1-4_600-2',
    '170823_Metab_acnb_fuma_ccm1-4_600-2',
    '170824_Metab_fis_lrp_ccm1-4_400ms_600-2',
    '170830_Metab_pgk_pgl_rpe_tpia_ccm1-4_400ms_600-2',
    ]
    
    prot_expnos = [
    [200],
    [200,300],
    [200,300,400,500,600,700,800],
    [200,300,400,500,600],
    [200,300,400],
    [200,300,400,500,600,700],
    [200,300],
    [200,300],
    [200,300,400,500],
    ]
        
    n_sets = len(dsets)
    # print n_sets
    
    # exit()
    flag_first_set = 1
    
    for mix in mix_decimals:
        print '----'
        for iSet in range(0,n_sets):
            # print dsets[iSet]
            for prot_expno in prot_expnos[iSet]:
                expno_final = prot_expno+mix+last_num_of_expno
                
                if not RUNNING_IN_TOPSPIN:
                    # print dsets[iSet], prot_expno, mix
                    print dsets[iSet], expno_final
                # read_expt(p_set, sls_expno, PM_P_M_PROCNO, curdir) # read main SLS diff
                else:    
                    read_expt(dsets[iSet], expno_final, procno_to_show, curdir) # read main SLS diff
                    if flag_first_set:
                        XCMD('.md no_load', wait = WAIT_TILL_DONE) # multiple display
                        flag_first_set = 0
        
    # if RUNNING_IN_TOPSPIN:
    #     ### XCMD('re '+str(expno), wait = WAIT_TILL_DONE) # read
    #     read_expt(p_set, sls_expno, PM_P_M_PROCNO, curdir) # read main SLS diff
    #     XCMD('.md no_load', wait = WAIT_TILL_DONE) # multiple display
    #     read_expt(p_set, p_expno+7, '1', curdir) # overlay P sls
    #     read_expt(p_set, 100+mix_dec+7, '1', curdir) # overlay M
    #     read_expt(p_set, sls_expno, PM_P_PROCNO, curdir) # overlay PMP check before subtracting compound
    #     # overlay PM_short ? (ref for FSA)        
    #     XCMD(read_compound_command, wait = WAIT_TILL_DONE) # overlay free compound
    #     # XCMD('zaro.py', wait = WAIT_TILL_DONE) # zoom
    #     # XCMD('zmet.py', wait = WAIT_TILL_DONE) # zoom

#####################################
## Helper functions

def read_expt(n,e,p,dir):
    fullpath = [n, str(e), str(p), dir]
    if RUNNING_IN_TOPSPIN:
        RE(fullpath, show="y") # Syntax/defaults: RE(dataset = None, show = "y")
        
    # Example
    #name = '_sref_test'
    #dir = '/Volumes/Data/yar/_eth2/data_NMR/spectra_tests'
    #readExpt(name,1,1,dir)
                        
#####################
# Fore debugging use constructions like:
#print("VARNAME = " + str(VAR))
#print "VARNAME = ", str(VAR)
#print("VARNAME = %.4f" % VAR)
#EXIT() # make a break-point

#####################
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
		main()