import sys
import scipy.optimize
import decoratio
import numpy as np

bedfile = sys.argv[1]
outfile = sys.argv[2]

prefix = bedfile.replace('.bed', '')

print(f"Sample:\t{prefix}")

with open(outfile) as f:
    lines = f.read().split('\n')
for line in lines[-15:-2]:
    if line.startswith('Number'):
        print(line[:-1].replace(': ', ':\t'))

# do decoratio (copied from its __main__ function)

LOG_CLASSES_BASE = np.e ** 0.02
ratio_min = 0.02


with open(bedfile) as f:
    dups = [int(x.split()[-1]) for x in f]
u, n = np.unique(dups, return_counts=True)
n = n/u
cz_d = dict(zip(u, n))
n_clones = sum(cz_d.values())
cz_d = [
        cz_d[z] if z in cz_d.keys() else 0
        for z in range(max(cz_d.keys())+1) ]
obs_cz_d = [ n / sum(cz_d) for n in cz_d ]
obs_duprate = decoratio.calc_duprate(obs_cz_d)
ratio_max = 1/(1-obs_duprate)
model = decoratio.AmplModel.factory('logskewnormal')
obs_log_cz_d = decoratio.radinitio_log_convert_clone_dist(obs_cz_d, LOG_CLASSES_BASE)
logratio_bounds = (np.log(ratio_min), np.log(ratio_max))
model_bounds = model.get_optimizable_parameters_bounds()
def f(x):
    ratio = np.exp(x[0])
    model.set_optimizable_parameters(x[1:])
    ampl_cz_d = model.get_ampl_cz_d(LOG_CLASSES_BASE)
    seq_log_cz_d = decoratio.radinitio_log_convert_clone_dist(
    decoratio.generate_cz_distrib(ratio, ampl_cz_d, LOG_CLASSES_BASE),
            LOG_CLASSES_BASE)
    return decoratio.score_log_cz_distrib(obs_log_cz_d, seq_log_cz_d, 'readfrac_sumsquares')
    
optim = scipy.optimize.shgo(f, [logratio_bounds] + model_bounds,
        sampling_method='sobol')

model.set_optimizable_parameters(optim.x[1:])
ratio, model = np.exp(optim.x[0]), model

ampl_cz_d = model.get_ampl_cz_d(LOG_CLASSES_BASE)
pred_cz_d = decoratio.generate_cz_distrib(ratio, ampl_cz_d, LOG_CLASSES_BASE)
model_fit = decoratio.cz_d_overlap(pred_cz_d, obs_cz_d)
pred_duprate = decoratio.calc_duprate(pred_cz_d)
saturation = ratio * (1.0 - decoratio.calc_duprate(pred_cz_d))
noisiness = decoratio.AmplModel.calc_noisiness(ampl_cz_d, LOG_CLASSES_BASE)
skeweness = decoratio.AmplModel.calc_skewness(ampl_cz_d, LOG_CLASSES_BASE)

print(f"D/C Ratio:\t{ratio}")
print(f"Model:\t{model}")
print(f"Observed duplicates:\t{obs_duprate}")
print(f"Predicted duplicates:\t{pred_duprate}")
print(f"Saturation:\t{saturation}")
print(f"Noisiness:\t{noisiness}")
print(f"Skeweness:\t{skeweness}")

