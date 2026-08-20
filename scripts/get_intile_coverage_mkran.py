'''
Given some input file (assumed to be fits) and input tiling file (in ecsv format), determine the number of times each
ra,dec in the input will be observed by the grism with the assumed wavelength coverage
requires that the optical model be installed and the environment variable ROMAN_GDPS_OPTICAL_MODEL_CONFIG be set to the path of the yaml config file for the optical model
Example run, TBC, some details are still NERSC specific,:
module load conda
conda activate /global/common/software/m4943/grizli1
export github_dir=/global/common/software/m4943/grizli0/
export PYTHONPATH=$PYTHONPATH:$github_dir/observing-program/py/:$github_dir/optical_model_tools/py/:$github_dir/GDPS_optical_model/
export ROMAN_GDPS_OPTICAL_MODEL_CONFIG=$github_dir/GDPS_optical_model/roman_gdps_optical_model/config/Roman_grism_OpticalModel_v0.8.yaml
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account m4943 python $github_dir/observing-program/scripts/get_intile_coverage_inDESIran.py --nran 1 --ramin 49 --ramax 51 --decmin -11 --decmax -9
to run the full survey 
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account m4943 python $github_dir/observing-program/scripts/get_intile_coverage_inDESIran.py --fullsurvey
exposure ID and detector number should be added to the output for a file containing each observation of each point, but this is not implemented yet
A bigger piece would be to get all of the pixel information along each trace. This would be a bit of work and hard to make fast.
'''
from rstgrs_footprint import sky_coords, trace_on_det
import logging
import argparse
from astropy.table import Table, unique
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from time import time


# create logger
logname = 'Roman_coverage'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


parser = argparse.ArgumentParser()
parser.add_argument(
    "--wficen", help="if y, positions are detector center", default='y')
parser.add_argument(
    "--wavmin", help="set minimum wavelength, if not None", default=None)
parser.add_argument(
    "--wavmax", help="set maximum wavelength, if not None", default=None)
parser.add_argument(
    "--chunksize", help="objects to process per chunk", default=10000000, type=int)
parser.add_argument("--tilefile", help="full path to tile file",
                    default=os.environ['github_dir']+'Roman_GRS_footprint/data/994-hlwas-Feb26_fixed_orients.sim.ecsv')
parser.add_argument(
    "--nran", help="number of randoms to use", default=10000, type=int)
parser.add_argument(
    "--outroot", help="root directory for output", default=os.getcwd())
parser.add_argument(
    "--outname", help="additional name for output files", default='test')
parser.add_argument("--ramin", help="minimum ra", default=49, type=float)
parser.add_argument("--ramax", help="maximum ra", default=51, type=float)
parser.add_argument("--decmin", help="minimum dec", default=-11, type=float)
parser.add_argument("--decmax", help="maximum dec", default=-9, type=float)
parser.add_argument("--target", help="which targets to include based on TARGET_NAME in tile file (ignored if not in tile file), options are medium (default, includes all except XMM and COSMOS), xmm, cosmos",
                    default='medium', choices=['medium', 'xmm', 'cosmos'])
parser.add_argument(
    "--fullsurvey", help="If set, override any ra,dec bounds and simulate the whole survey", action='store_true')
parser.add_argument(
    "--palist", help="list of PA values for repeated RA,DEC", default=None, type=str, nargs='*')
parser.add_argument(
    "--radiff", help="diff in RA for the repeated values", default=0, type=float)
parser.add_argument(
    "--decdiff", help="diff in DEC for the repeated values", default=0, type=float)
parser.add_argument(
    "--decpa", help="add diff in DEC for the flipped roll angles", default=0, type=float)
parser.add_argument(
    "--par", help="whether to use the parallel processing", action='store_true')
args = parser.parse_args()

if args.fullsurvey:
    logger.info(
        'fullsurvey flag is set, will ignore any ra,dec bounds and simulate the whole survey')
    # no Roman footprint seems to go out of these bounds
    decm = -60
    decx = 10
    ram = -10
    rax = 180

else:
    decm = args.decmin
    decx = args.decmax
    ram = args.ramin
    rax = args.ramax
    logger.info('will apply ra,dec bounds of ra > '+str(ram)+' and ra < '+str(rax) +
                ' and dec > '+str(decm)+' and dec < '+str(decx))
data = []

ran4degfn = args.outroot+'/DESIran' + \
    str(args.nran)+str(ram)+str(rax) + \
    str(decm)+str(decx)+'.fits'

data = Table()
ra, dec = sky_coords.generate_randoms(
    nran=args.nran, ra_bounds=(ram, rax), dec_bounds=(decm, decx))
data['RA'] = ra
data['DEC'] = dec
data['ID'] = np.arange(len(data)).astype(int)
logger.info('number of randoms used is '+str(len(data)))
# logger.info('apply ra,dec bounds, data has been cut from '+str(inputsize)+' to '+str(len(data)))

minwav = 1
maxwav = 1.9
wavstr = ''
wfistr = ''
if args.wficen != 'y':
    wfistr = 'notwficen'

if args.fullsurvey:
    outdir = args.outroot+'/test' + \
        '/fullsurvey'+wfistr+'/'
else:
    outdir = args.outroot+'/test' + \
        '/ramin'+str(ram)+'decmin'+str(decm)+wfistr+'/'
logger.info('results will be written to '+outdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)

fstr = 'Genran'+str(args.nran)+'_lam'+str(minwav)+str(maxwav)
outf = outdir+'nobs'+fstr+'grid_'+args.outname+'.ecsv'
logger.info('will save results to '+outf)


tiles = Table.read(args.tilefile)
tiles['EXPID'] = np.arange(len(tiles)).astype(int)
gtiles = tiles['BANDPASS'] == 'GRISM'
tls = tiles[gtiles]
if 'TARGET_NAME' in tiles.dtype.names:
    selxmm = tls['TARGET_NAME'] == 'XMM-LSS'
    selcos = tls['TARGET_NAME'] == 'COSMOS'
    if args.target == 'medium':
        sel = ~selxmm & ~selcos
        tls = tls[sel]
    if args.target == 'xmm':
        tls = tls[selxmm]
    if args.target == 'cosmos':
        tls = tls[selcos]
pad = 0.2
if args.wficen != 'y':
    pad = 0.5

logger.info('number of grism tiles in tiling file is '+str(len(tls)))

if args.palist is not None:  # this will only work for medium targets
    palist = np.array(args.palist).astype(float)
    logger.info(
        'palist provided, will set PA of repeated RA,DEC values to '+str(palist))
    tlsu = unique(tls, keys=['RA', 'DEC'])
    tral = []
    tdecl = []
    tpal = []
    for i in range(0, len(tlsu)):
        # for pa in palist:
        for j in range(0, len(palist)):
            pa = palist[j]
            tpal.append(pa)
            tral.append(tlsu['RA'][i]+j*args.radiff)
            deci = tlsu['DEC'][i]+j*args.decdiff
            if pa-palist[0] > 100:
                deci += args.decpa
            tdecl.append(deci)
    tls = Table()
    tls['RA'] = tral
    tls['DEC'] = tdecl
    tls['PA'] = tpal

if args.fullsurvey is False:
    selreg = tls['RA'] > ram-2*pad/np.cos(args.decmin*np.pi/180)
    selreg &= tls['RA'] < rax+2*pad/np.cos(args.decmin*np.pi/180)
    selreg &= tls['DEC'] > decm-pad
    selreg &= tls['DEC'] < decx+pad
    # print(len(tls[selreg]))
    tls = tls[selreg]
logger.info(
    'number of grism tiles to simulate after ra/dec cuts is '+str(len(tls)))
# nobs = np.zeros(len(ral_tot))


dets = np.arange(1, 19)
nr = 0
ra_all = []
dec_all = []

cnts_all = []
indx_all = []  # will contain unique IDs for each point, so will be shorter than the re other arrays
indx_re = []  # will contain repeated IDs for each observation of each point, so will be longer than the all other arrays
det_bits_re = []
expid_re = []
tottl = len(tls)

Nchunk = len(data)//args.chunksize + 1
Nchunk = int(Nchunk)
logger.info('will go through '+str(Nchunk)+' chunks')
for chunk in range(0, Nchunk):
    rand_indx = []
    det_bits_chunk = []
    expid_chunk = []
    min_indx = int(chunk*args.chunksize)
    max_indx = int((chunk+1)*args.chunksize)
    if max_indx > len(data):
        max_indx = len(data)

    t0 = time()

    ral_tot = data[min_indx:max_indx]['RA']
    decl_tot = data[min_indx:max_indx]['DEC']
    ran_indices = data[min_indx:max_indx]['ID']
    logger.info('cut data to chunk '+str(chunk))

    def get_idx_tl(tl):
        ra0 = tls['RA'][tl]
        if ra0 > 180:
            ra0 -= 360
        dec0 = tls['DEC'][tl]
        pa = tls['PA'][tl]
        eid = tls['EXPID'][tl]

        idx = []
        det_bits = []  # for encoding which detectors see each point
        ddec = decl_tot-dec0
        dra = ral_tot-ra0
        sel1deg = ddec > -1
        sel1deg &= ddec < 1
        dfac = np.cos(dec0*np.pi/180)
        sel1deg &= dra < 1/dfac
        sel1deg &= dra > -1/dfac
        ral_tl = ral_tot[sel1deg]
        decl_tl = decl_tot[sel1deg]
        ran_indices_tl = ran_indices[sel1deg]

        if len(ran_indices_tl) > 0:

            xfpa, yfpa = sky_coords.tangent_plane(
                np.array(ral_tl), np.array(decl_tl), ra0, dec0, pa)

            for det in dets:
                xpixl, ypixl = test_det.optmod.coords.convert_fpa_to_sca(
                    xfpa, yfpa, sca=det)
                selp = xpixl > -1000
                selp &= xpixl < 5088
                selp &= ypixl > -1000
                selp &= ypixl < 5088
                for i in range(0, len(xpixl[selp])):
                    xfpai = xfpa[selp][i]
                    yfpai = yfpa[selp][i]
                    test = 0

                    test = test_det.test_foot(
                        xfpai, yfpai, det=det, min_lam_4foot=minwav, max_lam_4foot=maxwav)
                    if test == 1:
                        idx_det = ran_indices_tl[selp][i]
                        idx.append(idx_det)
                        # append the detector number to the list of bits for this point
                        det_bits.append(int(det))
        return idx, det_bits, eid

    if args.par:
        from concurrent.futures import ProcessPoolExecutor
        tl_idx = list(np.arange(len(tls)).astype(int))
        with ProcessPoolExecutor() as executor:
            for idx, det_bits, eid in executor.map(get_idx_tl, tl_idx):
                rand_indx.append(idx)
                det_bits_chunk.append(det_bits)
                expid_chunk.append(np.ones(len(idx), dtype=int) * eid)
    else:
        for tl in range(0, len(tls)):
            idx, det_bits, eid = get_idx_tl(tl)
            rand_indx.append(idx)
            det_bits_chunk.append(det_bits)
            expid_chunk.append(np.ones(len(idx), dtype=int) * eid)
            print(str(tl)+' completed')

    logger.info('completed chunk '+str(chunk))
    # logger.info('length of list of ids '+str(len(rand_indx)))
    rand_indx = np.concatenate(rand_indx)
    det_bits_chunk = np.concatenate(det_bits_chunk)
    expid_chunk = np.concatenate(expid_chunk)
    logger.info('length of concatenated array of ids '+str(len(rand_indx)))
    rans, indices, cnts = np.unique(
        rand_indx, return_counts=True, return_inverse=True)
    selobs = np.isin(ran_indices, rans)
    if np.array_equal(ran_indices[selobs], rans):
        logger.info('input/final ids are matched in order')
    else:
        sys.exit('ids are not matched')
    # this doesn't work as a bitmask because some points are observed by the same detector in multiple tiles, so cannot reduce to a mapping of the set of detectors that see each point, would need to keep track of the number of times each detector sees each point to be able to encode as a bitmask, which is more complicated and may not be worth it for this application
    # det_bits_u = np.zeros(rans.size)
    # np.add.at(det_bits_u, indices, det_bits_chunk)
    racut = ral_tot[selobs]
    deccut = decl_tot[selobs]
    ra_all.append(racut)
    dec_all.append(deccut)
    cnts_all.append(cnts)
    indx_all.append(rans)
    indx_re.append(rand_indx)
    det_bits_re.append(det_bits_chunk)
    expid_re.append(expid_chunk)

    tf = time()
    logger.info('finished chunk '+str(chunk)+' in ' +
                str(tf-t0)+' out of '+str(Nchunk))


ra_all = np.concatenate(ra_all)
dec_all = np.concatenate(dec_all)
cnts_all = np.concatenate(cnts_all)
indx_all = np.concatenate(indx_all)
det_bits_re = np.concatenate(det_bits_re)
expid_re = np.concatenate(expid_re)
indx_re = np.concatenate(indx_re)
logger.info('concatenated data '+str(len(ra_all)) +
            ' data points with at least one observation')

tout = Table()
tout['RA'] = ra_all
tout['DEC'] = dec_all
tout['ID'] = indx_all
tout['NOBS'] = np.array(cnts_all, dtype=int)
# tout['DET_BITS'] = np.array(det_bits_all, dtype=int)
logger.info('about to write output to '+outf)
tout.write(outf, overwrite=True)
logger.info('wrote output for unique points to '+outf)

tout_re = Table()
tout_re['ID'] = indx_re
tout_re['SCA'] = det_bits_re
tout_re['EXPID'] = expid_re
outfre = outdir+fstr+'repeated_'+args.outname+'.ecsv'
logger.info('length of repeated table is '+str(len(tout_re)))
logger.info('about to write output for repeated points to '+outfre)
tout_re.write(outfre, overwrite=True)
logger.info('finished successfully!')
