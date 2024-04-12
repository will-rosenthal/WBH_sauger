import numpy
from numpy import array

import dadi
import demographic_models
import matplotlib as plt
import nlopt

dadi.cuda_enabled(True)

data = dadi.Spectrum.from_file("/project/ysctrout/sauger_walleye/dadi/output/dadi/Bighorn_River-Wind_River.sfs")

ns = data.sample_sizes
pts = [390, 400, 410]
func = demographic_models.split_migration_mod


upper_bound = [1, 1, 10, 5, 5,]
lower_bound = [1e-10, 0, 1e-10, 1e-20, 1e-20]

p0 = [0.5, 0.5, 0.01, 0.5, 0.5]

print('Starting optimization')

func_ex = dadi.Numerics.make_extrap_func(func)

p2 = dadi.Misc.perturb_params(p0, fold=1, upper_bound = upper_bound, lower_bound = lower_bound)
popt, ll_model = dadi.Inference.opt(p2, data, func_ex, pts, lower_bound=lower_bound, upper_bound=upper_bound,
                                      verbose=len(p0), fixed_params=[None, None, None, 0, 0],
                                      algorithm=nlopt.LN_COBYLA, log_opt=True)
print('Finished optimization, divergence w/out gene flow')
print(popt)

p0 = popt
#calculate best-fit model AFS
model = func_ex(popt, ns, pts)

# Likelihood of the data given the model AFS.
#ll_model = dadi.Inference.ll_multinom(model, data)
simp_ll = ll_model
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)

print('Optimal value of theta: {0}'.format(theta))


#ref = 2 * (((theta/4)/3.5e-9)/560118) #corrected for sauger
#ref = 2 * (((theta/4)/3.5e-9)/1043696) #corrected for sauger with MAF 0.01
ref = 2 * (((theta/4)/3.5e-9)/25442985) #corrected for sauger with no MAF filter
ref_simp_human = 2 * (((theta/4)/2.5e-8)/25442985)
ref_simp_herring = 2 * (((theta/4)/2.0e-9)/25442985)

print('ref = ')
print(ref)
ref_simp = ref

import pylab
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                   pop_ids =('Bighorn River','Wind River'), show=False)
# Save the figure
pylab.savefig('WBH_diverge_0.3_nomig_nomaf_mod.png', dpi=300)


print('#------------------------------NEW MODEL ---------------------------#')
p2 = dadi.Misc.perturb_params(p0, fold=1, upper_bound = upper_bound, lower_bound = lower_bound)
popt, ll_model = dadi.Inference.opt(p2, data, func_ex, pts, lower_bound=lower_bound, upper_bound=upper_bound,
                                      verbose=len(p0), fixed_params=[None, None, None, None, None],
                                      algorithm=nlopt.LN_BOBYQA, log_opt=True)
print('Finished optimization, divergence with migration')
print(popt)

#calculate best-fit model AFS
model = func_ex(popt, ns, pts)

# Likelihood of the data given the model AFS.
#ll_model = dadi.Inference.ll_multinom(model, data)
comp_ll = ll_model
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.

mig_theta = dadi.Inference.optimal_sfs_scaling(model, data)

print('Optimal value of theta: {0}'.format(mig_theta))

mig_ref = 2 * (((mig_theta/4)/3.5e-9)/25442985) #corrected for sauger
mig_ref_human = 2 * (((mig_theta/4)/2.5e-8)/25442985)
mig_ref_herring = 2 * (((mig_theta/4)/2.0e-9)/25442985)

print('ref = ')
print(ref)

import pylab
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('Bighorn River','Wind River'), show=False)
# Save the figure
pylab.savefig('WBH_diverge_0.3_mig_nomaf_mod.png', dpi=300)



print('#------------------------------NEW MODEL ---------------------------#')
p4 = dadi.Misc.perturb_params(p0, fold=1, upper_bound = upper_bound, lower_bound = lower_bound)
upstream_popt, ll_model = dadi.Inference.opt(p4, data, func_ex, pts, lower_bound=lower_bound, upper_bound=upper_bound,
                                      verbose=len(p0), fixed_params=[None, None, None, None, 0],
                                      algorithm=nlopt.LN_BOBYQA, log_opt=True)
print('Finished optimization, divergence with upstream migration')
print(upstream_popt)


#calculate best-fit model AFS
model = func_ex(upstream_popt, ns, pts)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
upstream_ll = ll_model
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
upstream_theta = dadi.Inference.optimal_sfs_scaling(model, data)

print('Optimal value of theta: {0}'.format(upstream_theta))


upstream_ref = 2 * (((upstream_theta/4)/3.5e-9)/25442985) #corrected for sauger
upstream_ref_human = 2 * (((upstream_theta/4)/2.5e-8)/25442985)
upstream_ref_herring = 2 * (((upstream_theta/4)/2.0e-9)/25442985)

print('ref = ')
print(upstream_ref)

pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('Bighorn River','Wind River'), show=False)
# Save the figure
pylab.savefig('WBH_diverge_0.3_upstreammig_nomaf_mod.png', dpi=300)



print('#----------------------------------------------------------#')
print('LRT and GIM uncert')

import dadi.Godambe
dict = dadi.Misc.make_data_dict_vcf("./sauger_miss0.3_filt.vcf", "./sauger_popinfo.txt")
chunk_size = 1e6
nboot = 100
chunks =  dadi.Misc.fragment_data_dict(dict, chunk_size)

all_boot = dadi.Misc.bootstraps_from_dd_chunks(chunks, nboot, ['Bighorn_River', 'Wind_River'], ns)


lrt = dadi.Godambe.LRT_adjust(func_ex, pts, all_boot, p0, data, nested_indices=[3,4], multinom= True)								   
adj = lrt*2*(comp_ll - simp_ll)
pval = dadi.Godambe.sum_chi2_ppf(adj, weights=(0.5,0.5))
print('Wind-Bighorn, no mig vs mig divergence p value = ')
print(pval)

if (pval > 0.05):
    print('simple model favored, no migration at all')
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts, all_boot, p0, data, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    div_time_gen =  p0[2]*ref_simp
    div_time_yrs = div_time_gen/0.33
    div_time_gen_human = p0[2]*ref_simp_human
    div_time_yrs_human = div_time_gen_human/0.33
    div_time_gen_herring = p0[2]*ref_simp_herring
    div_time_yrs_herring = div_time_gen_herring/0.33
    print('Divergence time in yrs [cichlid, human, herring] = ')
    print(div_time_yrs, div_time_yrs_human, div_time_yrs_herring)
    print('Divergence time in gen [cichlid, human, herring] = ')
    print(div_time_gen, div_time_gen_human, div_time_gen_herring)
else:
    print('complex model (with mig) favored')
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts, all_boot, popt, data, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    div_time_gen =  popt[2]*mig_ref
    div_time_yrs = div_time_gen/0.33
    div_time_gen_human = popt[2]*ref_simp_human
    div_time_yrs_human = div_time_gen_human/0.33
    div_time_gen_herring = popt[2]*ref_simp_herring
    div_time_yrs_herring = div_time_gen_herring/0.33
    print('Divergence time in yrs [cichlid, human, herring] = ')
    print(div_time_yrs, div_time_yrs_human, div_time_yrs_herring)
    print('Divergence time in gen [cichlid, human, herring] = ')
    print(div_time_gen, div_time_gen_human, div_time_gen_herring)


lrt = dadi.Godambe.LRT_adjust(func_ex, pts, all_boot, p0, data, nested_indices=[3], multinom= True)								   
adj = lrt*2*(upstream_ll - simp_ll)
pval = dadi.Godambe.sum_chi2_ppf(adj, weights=(0.5,0.5))
print('Wind-Bighorn, no mig vs upstream mig divergence p value = ')
print(pval)

if (pval > 0.05):
    print('simple model favored, no migration at all')
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts, all_boot, p0, data, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    div_time_gen =  p0[2]*ref_simp
    div_time_yrs = div_time_gen/0.33
    div_time_gen_human = p0[2]*ref_simp_human
    div_time_yrs_human = div_time_gen_human/0.33
    div_time_gen_herring = p0[2]*ref_simp_herring
    div_time_yrs_herring = div_time_gen_herring/0.33
    print('Divergence time in yrs [cichlid, human, herring] = ')
    print(div_time_yrs, div_time_yrs_human, div_time_yrs_herring)
    print('Divergence time in gen [cichlid, human, herring] = ')
    print(div_time_gen, div_time_gen_human, div_time_gen_herring)
else:
    print('complex model (with upstream mig) favored')
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts, all_boot, upstream_popt, data, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    div_time_gen =  upstream_popt[2]*upstream_ref
    div_time_yrs = div_time_gen/0.33
    div_time_gen_human = upstream_popt[2]*upstream_ref_human
    div_time_yrs_human = div_time_gen_human/0.33
    div_time_gen_herring = upstream_popt[2]*upstream_ref_herring
    div_time_yrs_herring = div_time_gen_herring/0.33
    print('Divergence time in yrs [cichlid, human, herring] = ')
    print(div_time_yrs, div_time_yrs_human, div_time_yrs_herring)
    print('Divergence time in gen [cichlid, human, herring] = ')
    print(div_time_gen, div_time_gen_human, div_time_gen_herring)