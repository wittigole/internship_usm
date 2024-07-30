# Project Overview on my internship at USM (LMU Munich, 2024)

Note: As this also shows the progress I made during the internship, some plots in the beginning might be outdated. For now, I decided to keep them so I can see what improved over time. I might delete those figures at the end of the internship. Values like R<sub>500</sub> always refer to the critical density, not the mean density, i.e. R<sub>500</sub> = R<sub>500,crit</sub>, M<sub>200</sub> = M<sub>200,crit</sub> etc.

## Getting familiar with the topic and the data

To get an impression for the kind of data I will be working with, my first task was to update estimates for the stellar masses in several clusters. [Chiu et al. (2018)](https://academic.oup.com/mnras/article/478/3/3072/4996803) first computed stellar masses for ~90 clusters from SED fitting. However, this relies on correct halo masses (M<sub>500</sub>). These were updated since the paper release. That is why I'm computing new stellar masses with the updated halo masses.

For this, I'm using stellar mass profiles from [Hennig et al. (2017)](https://academic.oup.com/mnras/article/467/4/4015/2939812). Specifically,

$$ \cfrac{m_\ast^{new}}{m_\ast^{old}} = \cfrac{\int_0^{R_{500}^{new}} \rho_\ast r^2\rm{d}r}{\int_0^{R_{500}^{old}} \rho_\ast r^2\rm{d}r} $$

This means, the ratio of new and old stellar masses depend on the NFW stellar mass profiles. We use the parameters (concentration parameter c) for those from [Hennig et al. (2017)](https://academic.oup.com/mnras/article/467/4/4015/2939812).

$$ \rho = \cfrac{\rho_0}{\left(\cfrac{r}{R_S}\right)\left(1 + \cfrac{r}{R_S}\right)^2} $$

Since they used R<sub>200</sub>, we have to convert c to the c for R<sub>500</sub>. This is done for every halo separately: c<sub>500</sub> = c<sub>200</sub> * R<sub>500</sub> / R<sub>200</sub>.

Since we don't know R<sub>200</sub> but only R<sub>500</sub>, we have to compute the former. We use `scipy.optimize.root` on the following expression.

$$ \left(\cfrac{R_{200}}{R_{500}}\right)^3 \cdot \cfrac{\log\left(\cfrac{R_{200} + cR_{500}}{R_{200}} - \cfrac{cR_{500}}{R_{200} + cR_{500}}\right)}{\log(1+c) - \cfrac{c}{1+c}} - 2.5 = 0$$

Here we used the following <a href="https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile#Density_distribution" target="_blank">integral</a> for the enclosed mass within R<sub>max</sub>.

$$ M = 4 \pi\int_0^{R_{max}}\rho(r)r^2\rm{d}r = 4\pi\rho_0 R_s^3 \left[\log\left(\cfrac{R_s + R_{max}}{R_s}\right) - \cfrac{R_{max}}{R_s + R_{max}}\right] = 4\pi\rho_0 \left(\cfrac{R_{vir}}{c}\right)^3 \left[\log\left(1 + c\right) - \cfrac{c}{1+c}\right]
 $$

for R<sub>max</sub> = R<sub>vir</sub> = cR<sub>s</sub>.

This gives us an estimate for R<sub>200</sub> only from the (fixed) concentration parameter and R<sub>500</sub>. We can now compute c<sub>500</sub> and thus find the desired ratio.

[//]: # (Comment test)

![First results](./plots/stellar_mass_ratio_distribution.jpg)
![newwer plot](./plots/ratio_comparison.jpg)

## Comparison to Simulations [^1]

In the following figure, I'm showing the total stellar mass as a function of their total cluster masses (yellow data points). For comparison, the best fit from the old stellar mass-halo mass relation from [Chiu et al. (2018)](https://academic.oup.com/mnras/article/478/3/3072/4996803) is shown. As the masses did not change too much, the fit still describes the data quite well (by eye). Additionally, the stellar mass-halo mass relation from the IllustrisTNG300-1 simulation is shown. There is a significant offset to the data. One aspect, although not important enough to explain the deviation, is the intra-cluster light (ICL) which is included for the TNG300 curve, but not for the observational data (as there is currently no way to determine its contribution).

![Comparison of data with best fit from Chiu et al. (2018) and TNG300-1.](./plots/stellar_vs_halo_mass.jpg)

The next step is to fit the scaling relation

$$ M_\ast = A_\ast  \left(\cfrac{M_{500}}{M_{\rm piv}}\right)^{B_\ast}\left(\cfrac{1+z}{1+z_{\rm piv}}\right)^{C_\ast} $$

to my data. I used an MCMC python module named <a href = "https://johannesbuchner.github.io/PyMultiNest/index.html#" target="_blank">PyMultiNest</a> for this purpose. It evaluates different sets of parameters (A<sub>* </sub>,B<sub>* </sub>,C<sub>* </sub>) based on some priors using the Log-Likelihood.

The first step to finding the best parameters is to choose a set of parameters (A<sub>* </sub>,B<sub>* </sub>,C<sub>* </sub>). Since we have some expectations on what the parameters will be, we can limit the range from which parameters are chosen. With these parameters and the M<sub>500</sub> and z data, we can use the scaling relation to compute a prediction for M<sub>* </sub>. This prediction together with a scatter quantity are then used as mean µ and standard deviation of a lognormal distribution:

$$ P(x|\mu,\sigma) = \cfrac{1}{x\sigma \sqrt{2\pi}}\exp\left(\cfrac{(\ln{x} - \mu)^2}{2\sigma^2}\right) $$

We need a lognormal distribution here because M<sub>* </sub> is normally distributed. However, since we're working with log(M<sub>* </sub>), this quantity is lognormally distributed. We can now compute the value of the lognormal distribution for every data point and add them up. The result is the Log-Likelihood for this set of parameters. The goal is to maximize this quantity by varying the parameters within the range in a smart way (so it does not take too long). The results from my fit to the new stellar masses can be found (not yet) [here](./files/chain_no_measurement_error_1/_1_stats.dat). The second block of values gives the results of the best fit for the parameters in order and the standard deviation. The latter is computed by evaluating different sets of parameters that have similar Log-Likelihoods. The results are:

+ normalization: A<sub>* </sub> = (3.85 +- 0.30) * 10^12 M<sub>sun</sub>
+ halo mass trend: B<sub>* </sub> = 0.703 +- 0.207
+ redshift trend: C<sub>* </sub> = 0.185 +- 0.516
+ scatter: D<sub>* </sub> = 0.057 +- 0.044

[//]: # (Clarify this matter with Aditya!)

From this one receives the following plot where I set z to z<sub>piv</sub>=0.6 and compare the data to several simulation cases.

![my_fit](./plots/stellar_vs_halo_mass_my_fit.jpg)

The next task is now to show the trend with halo mass and redshift separately. To exclude one or the other from the data, one simply divides by the respective term from the scaling relation with the best fit parameter:

$$ \left(\cfrac{M_{500}}{M_{\rm piv}}\right)^{B_\ast} \hspace{0.5cm}\text{or}\hspace{0.5cm}\left(\cfrac{1+z}{1+z_{\rm piv}}\right)^{C_\ast} $$

The fits are then altered by using M<sub>piv</sub>=4.8* 10^14 M<sub>sun</sub> and z<sub>piv</sub>=0.6 for M<sub>500</sub> and z, respectively. This procedure yields the following plots.

![mass trend](./plots/stellar_vs_halo_mass_wo_z.jpg)
![redshift trend](./plots/stellar_vs_halo_mass_wo_m500.jpg)

For all the simulations below except TNG300 and TNG-Clusters, I only have stellar masses within R<sub>500</sub> or stellar mass fractions. Consequently, it is not clear whether the ICL is include or exclude in these values. However, the latter is more likely, as one simply takes the sum of the masses of all stellar particles in the halo (or R<sub>500</sub>). One has to make corrections to account for that and a missing lower mass limit as discussed in the [TNG300 section](https://github.com/wittigole/internship_usm?tab=readme-ov-file#tng-clusters).


Another correction has to be taken into account, namely the different cosmologies applied in the simulations. In our data and [Chiu et al. (2018)](https://academic.oup.com/mnras/article/478/3/3072/4996803?login=false), a flat $\Lambda\text{CDM}$ cosmology is used with the following cosmological parameters (only the ones relevant for this analysis):

$$ \Omega_{\rm m} = 0.309, \hspace{0.5cm} \Omega_{\rm bar} = 0.048$$

Since the total amount of matter in the universe impacts the amount of stellar matter formed, we correct for the different baryonic matter densities by multiplying the stellar masses in each simulation with

$$f_{\rm corr} = \frac{f_{\rm bar}^{\rm data}}{f_{\rm bar}^{\rm sim}} \hspace{0.5cm}\text{where}\hspace{0.5cm} f_{\rm bar} = \frac{\Omega_{\rm bar}}{\Omega_{\rm m}}$$

Thereby, we receive the following results for the best simulation variation runs.

![best_sims](./plots/stellar_vs_halo_mass_wo_z_best_sims_rescaled.jpg)

Another correction one can do is to subtract the simulated BCG stellar masses from the total stellar masses of the clusters. This step is necessary as e.g. in TNG, the BCG stellar masses are too high when compared to observations. For this reason, we fit the scaling relation to the reduced stellar masses and add the scaling relation fit of the observed BCG masses on top of it. We use the BCG masses summarized by [Stott et al. (2010)](https://iopscience.iop.org/article/10.1088/0004-637X/718/1/23#fnref-apj340771r41). Their halo masses are given within R<sub>200</sub>, so they are converted to M<sub>500</sub> using the abovementioned methods. One obtains e.g. the dotted line in the TNG plot. The `scipy.optimize.odr` fitting routine yields the following for the normalization and power law index.

$$ M_{\rm norm} = (8.3\pm 1.0)\times 10^{11}\, \rm M_\odot \hspace{0.5cm}\text{and}\hspace{0.5cm} A_\ast = 0.50\pm 0.15 $$

In the following, we will discuss the performance of each of the simulations individually and compare (if applicable) different variation runs.

### TNG300 [^2]

![tng_data_comparison](./plots/stellar_vs_halo_mass_wo_z_tng.jpg)

We can see in the mass trend plot that there is a deviation between all TNG300 lines. They were computed as the fit of the scaling relation to the set of individual halo data points. However, especially for larger halo masses (essentially above 4* 10<sup>14</sup> solar masses) this fit might not be accurate since there are no clusters that massive in TNG300. The dashed line captures all stellar mass in satellites (within R<sub>500</sub>) as well as twice the stellar half-mass radius of the BCG to exclude the extended ICL. Obviously, this is not enough as there is a large deviation to our data points. The next correction was to only add the stellar masses of galaxies above a certain threshold. I chose 10<sup>10</sup> M<sub>sun</sub> in stellar mass, as this is comparable to what Chiu et al. (2018) did in their study. They integrated their stellar mass function with this value as the lower limit. Unfortunately, this brings the line only marginally closer to our data points. For comparison - although not meaningful - I also show the stellar mass halo mass relation for TNG300 without considering the contribution of BCGs (dash-dotted). This line is in the regime where we want it to be as it roughly captures all data points. Since [Chiu et al. (2018)](https://academic.oup.com/mnras/article/478/3/3072/4996803?login=false) did include BCGs in their analysis, this is just to evaluate the amount of stellar mass in the BCG.

As there is no widely accepted method for determining the contribution of the ICL to the total cluster stellar mass, both lines for TNG300 rather give an upper and lower limit. As described above, taking the total stellar mass of a cluster includes the ICL whereas only considering stars contained within 2 R<sub>0.5,* </sub> of each galaxy's center might miss outer parts of some galaxies. The true value is thus between the two lines. Nevertheless, even the lower boundary deviates strongly from our data.

When replacing the BCG stellar masses of TNG300 with those from [Stott et al. (2010)](https://iopscience.iop.org/article/10.1088/0004-637X/718/1/23#fnref-apj340771r41), we obtain significantly better results for the total scaling relation (compare dotted red line). TNG300 still has to large stellar masses, but all the corrections pay off.

IllustrisTNG is also based upon the Planck 2016 cosmological parameters: 

$$ \Omega_{\rm m} = 0.3089, \hspace{0.5cm} \Omega_{\rm bar} = 0.0486, \hspace{0.5cm} h = 0.6774 $$

### Magneticum [^3]

![magneticum_data_comparison](./plots/stellar_vs_halo_mass_wo_z_magneticum.jpg)

The only data available to me is the one listed on their website. It contains the M<sub>500</sub> and M<sub>star</sub> values for cluster in the Box 2 with high resolution, i.e. a box of length 352 MPc/h and $1584^3$ resolution elements (for both gas and DM respectively). There is also the larger Box 2b with side length of 640 MPc/h and $2\cdot 2880^3$ resolution elements ([Hirschmann et al. (2014)](https://academic.oup.com/mnras/article/442/3/2304/1039443), [Ragagnin et al. (2017)](https://academic.oup.com/mnras/article/486/3/4001/5475127)). In the plot above one cannot really distinguish both lines. The reason for this is that, despite one box being roughly eight times larger in volume than the other, the same physical models are used. As the simulation parameters are calibrated to yield certain relations, one would expect both lines to be very similar.

Nevertheless, both lines deviate strongly from our data points. Since there is no easy access to other variations of Magneticum, one cannot give a reason for this. Additionally, I only have stellar mass fractions given for each cluster. Presumably, they correspond to all stars within R<sub>500</sub>. Thus, we again have the problem with excluding the ICL.

Magneticum is based upon the WMAP7 cosmology:

$$ \Omega_{\rm m} = 0.272, \hspace{0.5cm} \Omega_{\rm bar} = 0.0456, \hspace{0.5cm} h = 0.704 $$

### FLAMINGO [^4]

![flamingo_data_comparison](./plots/stellar_vs_halo_mass_wo_z_flamingo.jpg)

The FLAMINGO ([Schaye et al. (2023)](https://academic.oup.com/mnras/article/526/4/4978/7246074)) clusters have stellar masses closer to our data than Magneticum. As for the other simulations, I don't know exactly what the stellar mass includes, only that it is limited to R<sub>500</sub>. Thus, the shown stellar masses might be too large which would lower the deviations to our data. There are 8 variations of FLAMINGO differing in supernovae (SN) and active galactic nuclei (AGN) feedback. They are displayed in two different colors where the color is meaningless. The most prominent line (dashed red) is the one closest to the data points whereas the others are plotted with $\alpha = 0.3$. It is from the variation with stronger SN feedback. Also quite good are the one with increased SN _and_ AGN feedback as well as the one with strong jets. Stronger and strongest AGN feedback actually does not perform as well as the beforementioned variation runs. Thus there is a limit to how strong AGN should be or SN are far more important.

What is good about the FLAMINGO simulations is that they possess clusters with masses above $10^{15} \rm{M}_\odot$. This is important as it matches our data better than TNG300 for example. Additionally, the fit of the scaling relation is more accurate at higher masses.

All of the eight FLAMINGO variations use the D3A cosmology with the following parameters:

$$ \Omega_{\rm m} = 0.306, \hspace{0.5cm} \Omega_{\rm bar} = 0.0486, \hspace{0.5cm} h = 0.681 $$

### BAHAMAS

![bahamas_data_comparison](./plots/stellar_vs_halo_mass_wo_z_bahamas.jpg)

All of the BAHAMAS ([McCarthy et al. (2016)](https://academic.oup.com/mnras/article/465/3/2936/2417021)) variations match our data better than any of the FLAMINGO runs does. However, there is still a significant deviation. The best agreement is reached when increasing the AGN feedback strength by raising the AGN heating temperature. In the two variations of the fiducial model of BAHAMAS, the heating temperatur of AGN is lowered (raised) by 0.2 dex. This means, that neighboring gas cells receive less (more) thermal energy, which decreases (increases) the AGN feedback strength. As visible in the plot, this has no strong impact on the stellar masses of galaxy clusters though.

The BAHAMAS simulations are based upon the WMAP9 cosmological parameters, namely:

$$ \Omega_{\rm m} = 0.279, \hspace{0.5cm} \Omega_{\rm bar} = 0.0463, \hspace{0.5cm} h = 0.700 $$

### C-OWLs

![cowls_data_comparison](./plots/stellar_vs_halo_mass_wo_z_cowls.jpg)

We have the reference run without AGN feedback (but _with_ SN feedback) and three models where the energy output of BHs is varied. In C-OWLs ([Le Brun et al. (2014)](https://academic.oup.com/mnras/article/441/2/1270/1064000)), the BH stores energy it gains from accreting mass until it reaches the threshold of heating up n<sub>heat</sub>=1 neighbouring gas cells by $\Delta T = 10^8 \rm K$. This is done in the AGN 8.0 model, in AGN 8.5 $\Delta T = 3\cdot 10^8 \rm K$ and in AGN 8.7 $\Delta T = 5\cdot 10^8 \rm K$. Since it takes longer to acquire the necessary amount of energy, the feedback is generally more bursty and energetic in case of higher $\Delta T$.

Just as Magneticum, C-OWLs is based upon the WMAP7 cosmological parameters (see there).

### TNG-Clusters [^5]



### mTNG [^6]


## Discussion

Now that we have seen the final (as of right now) [plot](https://github.com/wittigole/internship_usm/blob/main/plots/stellar_vs_halo_mass_wo_z_best_sims_rescaled.jpg), some remarks have to be adressed and discussed.

1. Some cluster catalogs from simulations are taken at $z=0$ whereas others are taken at redshifts close to 0.6. Ideally, one would prefer $z\approx 0.6$ as it represents z<sub>piv</sub> which is roughly the median redshift of our cluster sample. However, I couldn't access other snapshots for some simulations and had to settle for redshift zero. Though, this might not be as bad since the redshift evolution of our data (see respective plot) is almost negligible. So it wouldn't matter, at which redshift we consider the data. We can compare with the simulations basically at all redshifts.

2. For TNG300, we used mass-to-light ratios in the z-band to confirm that there are no huge deviations between the simulation and our clusters. We plotted the zband magnitude of a galaxy against its stellar mass and found agreement. Unfortunately, I had no access to z-band magnitudes of cluster galaxies for the other simulations either. The stellar masses are determined from the luminosities of mainly elliptical, red galaxies as they constitute most of the stellar mass of a cluster. If the simulations overestimate the luminosities in red bands, the stellar masses will alikely not agree with pur observational cluster sample.

3. contents of the stellar mass: ICL, BCG, satellites (mass cut??)



[^1]: Note: From here on, instead of using the abovementioned procedure to compute R<sub>200</sub>, we directly take the M<sub>200</sub> measurements and compute R<sub>200</sub> directly from this. By doing so, we avoid numerical inaccuracies that could arise as M<sub>500</sub> itself was not measured but calculated from the M<sub>200</sub> measurements.

[^2]: https://www.tng-project.org/

[^3]: http://www.magneticum.org/

[^4]: https://flamingo.strw.leidenuniv.nl/

[^5]: https://www.tng-project.org/cluster/

[^6]: https://www.mtng-project.org/