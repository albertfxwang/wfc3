#!/usr/bin/env python
"""
Parse APT file

- G. Brammer (July 2015)

"""
import os

import xmltodict

def go(apt_file='test.apt', program=None):
    """
    Parse an APT file to print a report like
    
==  (77) Visit: 4R  (ORIENT = 51.5  - 51.5 ) ==

   1 FIGS-GN2  GRISM1024  G102   SPARS100/15      (+1.762, +1.483) / (  0.00,   0.00)
   2 FIGS-GN2  GRISM1024  G102   SPARS100/14      (+1.558, +1.786) / ( -1.51,   2.50)
   3 FIGS-GN2  GRISM1024  F105W  SPARS25/11       (+1.558, +1.786) / ( -1.51,   2.50)
   4 FIGS-GN2  GRISM1024  G102   SPARS100/15      (-1.558, -1.786) / (-24.50, -26.99)
   5 FIGS-GN2  GRISM1024  G102   SPARS100/14      (-1.762, -1.483) / (-26.01, -24.49)
   6 FIGS-GN2  GRISM1024  F105W  SPARS25/13       (-1.762, -1.483) / (-26.01, -24.49)
    
    To make a pdf, pipe to an output file and run through a2ps:
    a2ps 13779.report -1 -f 7 --borders=no -s 2 -E -o 13779.ps
    
    """
    #### Plate scale
    ps = {'WFC3/UVIS':{'@X':0.039750, '@Y':0.039591}, 
          'WFC3/IR':{'@X':0.135603, '@Y':0.121308}} 
    
    if program:
        apt_file = '%d.apt' %(program)
        if not os.path.exists(apt_file):
            os.system('wget http://www.stsci.edu/hst/phase2-public/%d.apt' %(program))
            
    xml = ''.join(open(apt_file).readlines())
    d = xmltodict.parse(xml)
    obs = d['HSTProposal']['Observations']
    visits = d['HSTProposal']['Visits']['Visit']
    
    print ''
    
    pattern_dict = {}
    patterns = d['HSTProposal']['Patterns']
    if patterns is not None:
        patterns = patterns['Pattern']
        if not isinstance(patterns, list):
            patterns = [patterns]
        #
        for pattern in patterns:
            pattern_dict[pattern['Number']] = pattern 
            print 'Pattern |%2s|: %26s  N-points=%s' %(pattern['Number'], pattern['PatternType'], pattern['NumberOfPoints'])           
    else:
        patterns = []
        
    if not isinstance(visits, list):
        visits = [visits]
        
    #visit = visits[0]
    for i, visit in enumerate(visits):
        orient_str = ''
        if 'Orient' in visit.keys():
            orients = visit['Orient']
            if not isinstance(visit['Orient'], list):
                orients = [orients]
            #
            orient_str = '(ORIENT = %s' %(', '.join(['%11s' %(('%5s-%5s' %(orient['From']['@Degrees'], orient['To']['@Degrees'])).replace(' ','')) for orient in orients]))
            if 'FromVisit' in orients[0].keys():
                orient_str += ', from %s)' %(orients[0]['FromVisit']['@Visit'])
            else:
                orient_str += ')'
                
        print '\n==  (%2d) Visit: %s  %s ==\n' %(i+1, visit['@Number'], orient_str)
        if 'ExposureGroup' not in visit.keys():
            grps = []
        else:
            grps = visit['ExposureGroup']
        
        if not isinstance(grps, list):
            grps = [grps]
            
        i_exp = 0
        has_postarg = 0
        
        all_exps = []
        grp_pattern = []
        if 'Exposure' in visit.keys():
            for exp in visit['Exposure']:
                if 'WFC3' in exp['@Config']:
                    all_exps.append(exp)
                    grp_pattern.append('')
            
        for grp in grps:
            if '@Pattern' in grp.keys():
                pattern_str = ' |%2s|' %(grp['@Pattern'])
            else:
                pattern_str = ''
                
            if 'ExposureGroup' in grp.keys():
                sub_grps = grp['ExposureGroup']
            else:
                sub_grps = [grp]
            
            if not isinstance(sub_grps, list):
                sub_grps = [sub_grps]
            
            for sub_grp in sub_grps:
                exps = sub_grp['Exposure']
                if not isinstance(exps, list):
                    exps = [exps]
                #
                for exp in exps:
                    if 'WFC3' in exp['@Config']:
                        all_exps.append(exp)
                        grp_pattern.append(pattern_str)
                        
        for exp, exp_pattern in zip(all_exps, grp_pattern):
            if 'IR' in exp['@Config']:
                #sampseq = '.'.join(exp['InstrumentParameters'].values()[::-1])
                sampseq = '%8s %-2s' %(exp['InstrumentParameters']['@SAMP-SEQ'], exp['InstrumentParameters']['@NSAMP']) 
            else:
                sampseq = exp['ExposureTime']['@Secs']+'s'
                flash = ''
                if isinstance(exp['InstrumentParameters'], dict):
                    if '@FLASH' in exp['InstrumentParameters'].keys():
                        flash = ' FLASH-%-2s' %(exp['InstrumentParameters']['@FLASH'])
                sampseq += flash
                
            sampseq = '%-11s %s' %(sampseq, exp_pattern)
            
            #### POS-TARGs
            offset_str = ''
            if 'PosTarg' in exp.keys():
                pt_xy = ', '.join(exp['PosTarg'].values())
                pt_asec = []
                for key in ['@X', '@Y']:
                    pt_asec.append(float(exp['PosTarg'][key]) / ps[exp['@Config']][key])
                #
                if has_postarg == 0:
                    pt_ref = [pt_asec[0], pt_asec[1]]
                    has_postarg = 1
                #
                pt_asec = ', '.join(['%6.2f' %(pt_asec[i]-pt_ref[i]) for i in range(2)])
                offset_str = '(%14s) / (%s)' %(pt_xy, pt_asec)
            #
            if '@SamePosAs' in exp.keys():
                offset_str += ' SamePosAs-%s' %(exp['@SamePosAs'])
            #
            details = '  %2s %s  %9s %-16s  %-5s  %-20s  %s' %(exp['@Number'], exp['@TargetName'], exp['@Config'], exp['@Aperture'], exp['@SpElement'], sampseq, offset_str)
            print details #
            i_exp += 1

if __name__ == '__main__':
    import sys
    if 'apt' not in sys.argv[1]:    
        go(program=int(sys.argv[1]))
    else:
        go(apt_file=sys.argv[1])
        
