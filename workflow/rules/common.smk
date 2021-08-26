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

# calculate preprocessing subsampling
# required to normalise both bams for
# consistent tumour coverage
cov = sample_info.coverage.values
ploidy = sample_info.ploidy.values
purity = sample_info.purity.values
nrpcc = (cov / ploidy) * purity

# we must normalise based on the ratio
# or effective tumour coverage between
# the samples. only one sample will need
# normalisation, the other will the the
# constant
norm1 = round(nrpcc[0] / nrpcc[1], 4)
norm2 = round(nrpcc[1] / nrpcc[0], 4)
sample_info['norms'] = [norm1, norm2]

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

recipe['norms_s1'] = [str(n).split('.')[1] if n < 1 else '0' for n in recipe.norms_s1.values]
recipe['norms_s2'] = [str(n).split('.')[1] if n < 1 else '0' for n in recipe.norms_s2.values]
recipe['prop1'] = [str(p).split('.')[1] for p in recipe.prop1.values]
recipe['prop2'] = [str(p).split('.')[1] for p in recipe.prop2.values]

samples = np.concatenate((recipe.sample1.values, recipe.sample2.values))
norm_seeds = np.concatenate((recipe.seed_s1.values, recipe.seed_s2.values))
norms = np.concatenate((recipe.norms_s1.values, recipe.norms_s2.values))
seeds = np.concatenate((recipe.seed1.values, recipe.seed2.values))
props = np.concatenate((recipe.prop1.values, recipe.prop2.values))

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
