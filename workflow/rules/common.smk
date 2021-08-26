import pandas as pd
import numpy as np

sample_info = pd.read_csv(config['samples'], sep='\t')
recipe = pd.read_csv(
    config['twoclus_recipe'],
    sep='\t',
    dtype={'prop1': float, 'prop2': float, 'seed1': int, 'seed2': int})

# input checks
assert len(sample_info) == 2
assert(recipe.sample1.isin(sample_info.sample_name).all())
assert(recipe.sample2.isin(sample_info.sample_name).all())

s1_samples = recipe['sample1'].values
s2_samples = recipe['sample2'].values

s1_props = [str(p).split('.')[1] for p in recipe['prop1'].values]
s2_props = [str(p).split('.')[1] for p in recipe['prop2'].values]

s1_seeds = recipe['seed1'].values
s2_seeds = recipe['seed2'].values

samples = np.concatenate((s1_samples, s2_samples))
props = np.concatenate((s1_props, s2_props))
seeds = np.concatenate((s1_seeds, s2_seeds))

# calculate preprocessing subsampling
# required to normalise both bams for
# consistent tumour coverage
nrpcc = (sample_info['coverage'] / sample_info['ploidy']) * sample_info['purity']
sub1 = round(nrpcc[0] / nrpcc[1], 4)
sub2 = round(nrpcc[1] / nrpcc[0], 4)

pseeds = sample_info['seed'].values
norms = [str(s).split('.')[1] if s < 1 else '0' for s in [sub1, sub2]]
input_samples = sample_info.sample_name.values
sample_info['norms'] = norms

# add these values back into the
# recipe dataframe
recipe = recipe.merge(
    sample_info,
    how='left',
    left_on='sample1',
    right_on='sample_name'
)
recipe = recipe.merge(
    sample_info,
    how='left',
    left_on='sample2',
    right_on='sample_name',
    suffixes=('_s1', '_s2')
)

print(
'''
#------------------------ RECIPE -------------------------- 
'''
)
for idx, row in recipe.iterrows():
    print(
        'sample 1: %s, sample 2: %s, prop1: %s, prop2: %s' % \
        (row['sample1'], row['sample2'], row['prop1'], row['prop2'])
        )
print(
'''
#----------------------------------------------------------
'''
)
