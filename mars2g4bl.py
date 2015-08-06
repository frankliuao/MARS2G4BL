"""Convert the MARS output to G4Beamline beam format

Usage:
python mars2g4bl.py marsfile

The MARS file should be in the following format:
Every row contains:
i x y z px py pz w (t)
i - particle type, 1=proton 3=pi+
x,y,z - particle coordinate
pz,py,pz - particle momentum projections
w - particle weight
t - time of record (In some cases, the t column is not included.)

"""
from __future__ import division
import os
import sys

import numpy as np

from beam import beam

__version__ = "1.02a"
__history__ ={"02/09/2015": "antiproton's PDGID should have been -2212; was "
                            "2212. changed.",
              "07/07/2015": "Changed the script so that it can be used as a "
                            "module",
              }
__help__ = "Usage: python mars2g4bl.py marsfile"
__extrahelp__ = """Explanation of each key:

'p_cut':[p_min, p_max], a list/ndarray with 2 floats numbers, i.e. minimum and
maximum cut on the desired momentum range. If no momentum cut is desired,
use [0, 0];

'mult': float, a factor to be multiplied on the particle weights. Notice
that this procedure is done before duplicating the particles;

'duplicate': logical, a logical value to determine whether each particle
will be duplicated by its weight. e.g. A particle with weight 2.3 will be
duplicated by 2 times with a 30% probability, and 3 times with a 70%
probability;

'frac': float, the fraction of the beam to be outputted. The choice of
the particles is completely random. Notice this step
is done AFTER the beam is duplicated (if set). e.g. If this value = 0.1,
and there are 1000 particles in the beam, only 100 particles will be
outputted.

"""
__all__ = ['MARS2G4BL']

MARS2G4BL_PDG = {1:2212,2:-2112,3:211,4:-211,5:321,6:-321,7:-13,8:13,9:22,
                10:11,11:-11}


def main(marsfile,inputdict):


    g4blbeam = beam()

    try:
        print 'Loading MARS beam file'
        marsbeam_mat = g4blbeam.loadtxt_fast(marsfile,0)
        print 'MARS beam file loaded!'
    except Exception as ee:
        print ee
        print '---Error: Can not load file '+marsfile + '---'
        sys.exit()

    # See if the time column is included. If yes, the # of columns should
    # be 9.
    if marsbeam_mat.shape[1]==9:
        time_included = True
    elif MARSbeamMat.shape[1]==8:
        time_included = False
    else:
        print __help__
        print '---Error: Not enough number of columns provided to get the ' \
          'particle information! ---'

    # Read the MARS pathname and filename without extension
    marsfile_dir = os.path.dirname(os.path.abspath(marsfile))
    marsfile_base = os.path.basename(marsfile)
    marsfile_noext = '.'.join(marsfile_base.split('.')[:-1])
    g4blfile = marsfile_dir+'/'+marsfile_noext+'.beam'
    print "The name of the output will be "+g4blfile

    # The total momentum of the particles:
    total_p = np.sqrt(marsbeam_mat[:,4]**2 +
                      marsbeam_mat[:,5]**2 +
                      marsbeam_mat[:,6]**2)

    # Momentum cut
    if inputdict['p_cut'] != [0,0]:
        marsbeam_mat = marsbeam_mat[np.all([total_p<=inputdict['p_cut'][1]/1e3,
                                            total_p>=inputdict['p_cut'][0]/1e3],
                                           axis=0),:]
    g4blbeam_mat = np.zeros([marsbeam_mat.shape[0],12])
    # MARS uses cm, GeV/c, and second. Do unit conversion.
    g4blbeam_mat[:,0:3] = marsbeam_mat[:,1:4]*10
    g4blbeam_mat[:,3:6] = marsbeam_mat[:,4:7]*1000
    if time_included: g4blbeam_mat[:,6] = marsbeam_mat[:,8]*1e9
    # Convert the PDGid in MARS to G4BL:
    parse_pdg = np.vectorize(lambda x: MARS2G4BL_PDG[int(x)])
    g4blbeam_mat[:,7] = parse_pdg(marsbeam_mat[:,0])
    # Weight multiplication
    marsbeam_mat[:,7] *= inputdict['mult']

    # Duplication
    if inputdict['duplicate']:
        particle_weight_int = np.floor(marsbeam_mat[:,7])
        particle_weight_float = marsbeam_mat[:,7] - particle_weight_int
        add_weight_logic = np.random.random(marsbeam_mat.shape[0]) >= \
                           particle_weight_float
        marsbeam_mat[:,7] = particle_weight_int + add_weight_logic
        g4blbeam_mat = np.repeat(g4blbeam_mat,
                                 np.array(marsbeam_mat[:,7], dtype=int),
                                 axis=0)
        g4blbeam_mat[:,11] = 1
    else:
        g4blbeam_mat[:,11] = marsbeam_mat[:,7]

    # Fraction of the beam.
    g4blbeam_index = np.arange(g4blbeam_mat.shape[0])
    np.random.shuffle(g4blbeam_index)
    index_keep = g4blbeam_index[0:np.floor(inputdict[
                                'fraction']*g4blbeam_mat.shape[0])]
    g4blbeam_mat = g4blbeam_mat[index_keep,:]

    # Re-assign the eventID.
    g4blbeam_mat[:,8] = np.arange(g4blbeam_mat.shape[0])+1

    g4blbeam.loadBeam(g4blbeam_mat)
    g4blbeam.writeBeam(g4blfile)
    print 'The conversion is done!'
    sys.exit(0)


def mars2g4bl(marsfile,inputdict):
    """Call the main function.

    Will be imported through "from MARS2G4BL import *"
    Equivalent to the main function.
    """
    if inputdict.has_key('p_cut') and inputdict.has_key('fraction') and \
            inputdict.has_key('mult') and inputdict.has_key('duplicate'):
        main(marsfile,inputdict)
    else:
        print "The inputdict dictionary must have all the following keys:\n" \
              "'p_cut','fraction','mult','duplicate'."
        print __extrahelp__
        exit(1)


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print "Provide the MARS file as the script argument!"
        print __help__
        sys.exit()
    else:
        if sys.argv[1].lstrip('-') == 'help':
            print __help__
            sys.exit(0)
    # User input: Momentum cut
    p_cut = raw_input('\nEnter the Momentum Cut boundaries in MeV/c '
                         'with the form [p_low, p_high]\nEnter 0 if you do '
                         'not need a cut\n')
    if p_cut != '0':
        try:
            p_cut = eval(p_cut)
            if CutLower >= CutUpper:
                print "---Error: The lower cut is greater than the upper cut!"
                sys.exit(1)
        except Exception as e:
            print e
            print '---Error: Wrong format: Enter the Momentum Cut ' \
                  'boundaries in MeV/c with the form [p_low, p_high] ---'
            sys.exit(1)
    else:
        p_cut = [0, 0]
    # User input: the multiplication factor on the particle weights.
    mult = raw_input('\nDo you want to multiply the particle weights '+
                        'by a factor?\n'+
                        'Enter factor: (Enter 1 if no multiplication)\n' +
                        'If duplication is needed later, the particle will '+
                        'be duplicated based on the new weights\n')
    mult = float(mult)
    if mult <= 0:
        print '---Error: You can only multiply the weight by a positive ' \
              'number! ---'
        sys.exit(1)
    # User input: whether to duplicate the particles based on the weights.
    duplicate = raw_input('\nDo you want to duplicate the particles based '
                           'on their current weight?\ny or n? (n)\n')
    duplicate = True if duplicate.lower() == 'y' else False
    # User input: the fraction of the beam to be outputted.
    fraction = raw_input('\nDo you want to take only a fraction of the beam? '
          'The pick will be random.\n'+'Enter the fraction n: (1 for all the '
                                     'beam)\n')
    if float(fraction)>0:
        fraction = float(fraction) if fraction<=1 else 1
    else:
        print '---Error: Could not take a negative ' \
              'fraction of the beam! ---'
        sys.exit(1)

    main(sys.argv[1], {'p_cut':p_cut, 'duplicate':duplicate, 'mult':mult,
                      'fraction':fraction})
