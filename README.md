# Project Overview on my internship at USM (LMU Munich, 2024)

## Getting familiar with the topic and the data

To get an impression for the kind of data I will be working with, my first task was to update estimates for the stellar masses in several clusters. Chiu et al (2018) first computed stellar masses for ~90 clusters from SED fitting. However, this relies on correct halo masses (M<sub>500</sub>). These were updated since the paper release. That is why I'm computing new stellar masses with the updated halo masses.

For this, I'm using stellar mass profiles from Hennig et al (2017). Specifically,

![equation](https://latex.codecogs.com/gif.latex?%5Cdpi%7B110%7D%5Cbg%7Bwhite%7D%5Ccfrac%7Bm_%5Cast%5E%7Bnew%7D%7D%7Bm_%5Cast%5E%7Bold%7D%7D=%5Ccfrac%7B%5Cint_0%5E%7BR_%7B500%7D%5E%7Bnew%7D%7D%5Crho_%5Cast%20r%5E2%5Crm%7Bd%7Dr%7D%7B%5Cint_0%5E%7BR_%7B500%7D%5E%7Bold%7D%7D%5Crho_%5Cast%20r%5E2%5Crm%7Bd%7Dr%7D)

This means, the ratio of new and old stellar masses depend on the NFW stellar mass profiles. We use the parameters (concentration parameter c) for those from Hennig et al (2017).

![equation](https://latex.codecogs.com/gif.image?%5Cinline%20%5Cdpi%7B110%7D%5Cbg%7Bwhite%7D%5Crho(r)=%5Ccfrac%7B%5Crho_0%7D%7B%5Cleft(%5Ccfrac%7Br%7D%7BR_s%7D%5Cright)%5Cleft(1&plus;%5Ccfrac%7Br%7D%7BR_s%7D%5Cright)%5E2%7D)

Since they used R<sub>200</sub>, we have to convert c to the c for R<sub>500</sub>. This is done for every halo separately: c<sub>500</sub> = c<sub>200</sub> * R<sub>500</sub> / R<sub>200</sub>.

Since we don't know R<sub>200</sub> but only R<sub>500</sub>, we have to compute the former. We use `scipy.optimize.root` on the following expression.

![insert equation](https://latex.codecogs.com/gif.image?%5Cdpi%7B110%7D%5Cbg%7Bwhite%7D%5Cleft(%5Ccfrac%7BR_%7B200%7D%7D%7BR_%7B500%7D%7D%5Cright)%5E3%5Ccfrac%7B%5Clog%5Cleft(%5Ccfrac%7BR_%7B200%7D&plus;cR_%7B500%7D%7D%7BR_%7B200%7D%7D%5Cright)-%5Ccfrac%7BcR_%7B500%7D%7D%7BR_%7B200%7D&plus;cR_%7B500%7D%7D%7D%7B%5Clog(1&plus;c)-%5Ccfrac%7Bc%7D%7B1&plus;c%7D%7D-2.5)

Here we used the following [integral](https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile#Density_distribution) for the enclosed mass within R<sub>max</sub>.

![equation](https://latex.codecogs.com/gif.image?%5Cinline%20%5Cdpi%7B110%7D%5Cbg%7Bwhite%7DM=4%5Cpi%5Cint_0%5E%7BR_%7Bmax%7D%7D%5Crho(r)r%5E2%5Crm%7Bd%7Dr=4%5Cpi%5Crho_0%20R_s%5E3%5Cleft%5B%5Clog%5Cleft(%5Ccfrac%7BR_s&plus;R_%7Bmax%7D%7D%7BR_s%7D%5Cright)-%5Ccfrac%7BR_%7Bmax%7D%7D%7BR_s&plus;R_%7Bmax%7D%7D%5Cright%5D=4%5Cpi%5Crho_0%5Cleft(%5Ccfrac%7BR_%7Bvir%7D%7D%7Bc%7D%5Cright)%5E3%5Cleft%5B%5Clog%5Cleft(1&plus;c%5Cright)-%5Ccfrac%7Bc%7D%7B1&plus;c%7D%5Cright%5D)

for R<sub>max</sub> = R<sub>vir</sub> = cR<sub>s</sub>.

This gives us an estimate for R<sub>200</sub> only from the (fixed) concentration parameter and R<sub>500</sub>. We can now compute c<sub>500</sub> and thus find the desired ratio.

[//]: # (Comment test)

![First results](./plots/stellar_mass_ratio_distribution.jpg)

## Comparison to Simulations

In the following figure, I'm showing the total stellar mass as a function of their total cluster masses (yellow data points). For comparison, the best fit from the old stellar mass-halo mass relation from Chiu et al. (2018) is shown. As the masses did not change too much, the fit still describes the data quite well (by eye). Additionally, the stellar mass-halo mass relation from the IllustrisTNG300-1 simulation is shown. There is a significant offset to the data. One aspect, although not important enough to explain the deviation, is the intra-cluster light (ICL) which is included for the TNG300 curve, but not for the observational data (as there is currently no way to determine its contribution).

![Comparison of data with best fit from Chiu et al. (2018) and TNG300-1.](./plots/stellar_vs_halo_mass.jpg)



