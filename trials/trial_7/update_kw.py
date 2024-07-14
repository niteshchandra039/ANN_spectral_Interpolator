import numpy as np
import pandas as pd
import glob
from astropy.io import fits

TGM = np.loadtxt('/home/nitesh/nitesh/PhD/GSL_Interpolator/data/grid_params.txt')


file_name = glob.glob('*.fits')[0]


dict_ = {}
dict_['name'] = ['Teff', 'Logg', 'Fe/H']
dict_['unit'] = ['K', 'dex', 'dex']
dict_['lower bound'] = [ min(TGM[:,0]), min(TGM[:,1]), min(TGM[:,2])]
dict_['higher bound'] = [ max(TGM[:,0]), max(TGM[:,1]), max(TGM[:,2])]
dict_['step for numerical derivation'] = np.array([100, 0.5, 1.0])/100

dict_['mean'] = [TGM[:,0].mean(), TGM[:,1].mean(), TGM[:,2].mean()]
dict_['std'] = [TGM[:,0].std(), TGM[:,1].std(), TGM[:,2].std()]
_temp_df = pd.DataFrame(dict_)


new_hdul = fits.open(file_name)

for i in range(len(new_hdul)):
    del new_hdul[i].header['I_PR_*']
    for j in range(len(_temp_df)):
        new_hdul[i].header['I_P_NM{}'.format(j+1)]='{}'.format(_temp_df["name"][j])
        new_hdul[i].header.comments['I_P_NM{}'.format(j+1)]='{}'.format(_temp_df.columns[0])

        new_hdul[i].header['I_P_UN{}'.format(j+1)]='{}'.format(_temp_df['unit'][j])
        new_hdul[i].header.comments['I_P_UN{}'.format(j+1)]='{}'.format(_temp_df.columns[1])

        new_hdul[i].header['I_P_LB{}'.format(j+1)]= _temp_df['lower bound'][j]
        new_hdul[i].header.comments['I_P_LB{}'.format(j+1)]='{}'.format(_temp_df.columns[2])

        new_hdul[i].header['I_P_HB{}'.format(j+1)]= _temp_df['higher bound'][j]
        new_hdul[i].header.comments['I_P_HB{}'.format(j+1)]='{}'.format(_temp_df.columns[3])

        new_hdul[i].header['I_P_ST{}'.format(j+1)]= _temp_df['step for numerical derivation'][j]
        new_hdul[i].header.comments['I_P_ST{}'.format(j+1)]='{}'.format(_temp_df.columns[4])

        new_hdul[i].header['I_PR_B{}'.format(j+1)] = _temp_df['mean'][j]
        new_hdul[i].header.comments['I_PR_B{}'.format(j+1)] = 'Mean of the {}'.format(_temp_df["name"][j])

        new_hdul[i].header['I_PR_S{}'.format(j+1)]= _temp_df['std'][j]
        new_hdul[i].header.comments['I_PR_S{}'.format(j+1)] = 'Std of the {}'.format(_temp_df["name"][j])

    new_hdul[i].header['I_VERSIO'] = '0.0.3'
    new_hdul[i].header.comments['I_VERSIO'] = 'Version of the perceptron file format'

    new_hdul.writeto(file_name[0:22]+'3'+file_name[23:], overwrite=True)
print("DONE")
