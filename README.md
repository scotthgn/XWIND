# XWIND

Analytic model for calculating line-profiles from outflowing winds in AGN. Specifically, this is designed for slow(ish) winds at BLR scales. It is originally designed with the Fe-K $\alpha$ complex in mind, but also allows for the convolution of the same transfer functions for any line species (under the strong assumption that their emissivity traces density in the same way as netral Fe-K $\alpha$). This code is specifically designed for use in XSPEC. For a non-XSPEC version, see [PyXWIND](https://github.com/scotthgn/PyXWIND).

If this model is used in your work, please cite: **Hagen et al. (2026, submitted to A&amp;A)**

**Note 1**: Depending on your version of XSPEC you will need to use one of two available branches:
1. For XSPEC version $\geq$ v.15.1, use **main** branch
2. For XSPEC version < v.15.1, use **old_xspec** branch

This comes from recent changes in how XSPEC handles and compiles fortran code

**Note 2**: Within XWIND there are three submodels: `xwindline`, `xwindconv`, and `xindfe`. After installation these are called within XSPEC directly (i.e `model xwindline`). See below for a parameter description of each.



## Requirements

A working [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/) installation. Must be compiled from Source code (otherwise local models do not work).
If HEASOFT installs properly, then you will also fulfull the requierements for fortran compilors

## Installation
1. Clone the repository
2. Option 1 (manual installation): </br>
&emsp; &emsp; 2.1 Open XSPEC, and `cd` into the `/src` directory </br>
&emsp; &emsp; 2.2 Within XSPEC type: `initpackage xwind lmod_xwind.dat .` This will compile the code, and assumes you are currectly within the directory containing the source code files. </br>
&emsp; &emsp; 2.3 Still within XSPEC type: `lmod xwind .` </br>
3. Option 2 (quick automatic installation) </br>
&emsp; &emsp; 3.1 Open XSPEC in the top level directory </br>
&emsp; &emsp; 3.2 Run the compile script. Type: `@compile.xcm` (this does steps 2.2 and 2.3 for you) </br>
4. (OPTIONAL EXTRA): By default you will have to run the `lmod xwind .` step each time you start a new XSPEC session. To avoid this, locate your `xspec.rc` file (typically located in the `.xspec/` directory within the home directory). Edit the `xspec.rc` file to contain the line: `lmod xwind /path/to/XWIND/src` where /path/to/XWIND/src is the full path to the source code (i.e what you get by typing `pwd` within the `XWIND/src` directory).
    

## Model Descrption
Within XWIND there are three submodels. These are: `xwindline`, `xwindconv`, and `xwindfe`. Note, while a breif overview of parameters is given here, for a full model desciption including a more physical meaning of each parameter, see **Hagen et al. (2026, submitted to A&amp;A)**

### xwindline
The base model. This simply focuses on the line-shape as given by a wind. As such the res-frame energy is treated as a free-parameter, and the rest-frame emission is considered a simple delta-function. 

**Par 1. &ensp; $\log_{10} \dot{m}_{w}$** </br>
&emsp; &emsp; &#9656; **Units:** $\dot{M}_{w} / \dot{M}\_{\rm{Edd}}$ </br>
&emsp; &emsp; &#9656; **Description:** Wind mass-outflow rate, scaled by the Eddington accretion rate 

**Par 2. &ensp; $r_{\rm{in}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Inner launch radius

**Par 3. &ensp; $r_{\rm{out}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Outer launch radius

**Par 4. &ensp; $d_{f}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Distance from origin to wind focus (see Fig. 1 in Hagen et al. 2026)

**Par 5. &ensp; $f_{\rm{cov}}$** </br>
&emsp; &emsp; &#9656; **Units:** $\frac{\Omega}{4 \pi}$ </br>
&emsp; &emsp; &#9656; **Description:** Covering fraction fo the wind as seen from the central (illuminating) source

**Par 6. &ensp; $\log_{10} v_{\infty}$** </br>
&emsp; &emsp; &#9656; **Units:** $c$ </br>
&emsp; &emsp; &#9656; **Description:** Outflow velocity at infinity

**Par 7. &ensp; $r_{v}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Velocity scale length. i.e. The distance along the streamline where the wind reaches half $v_{\infty}$

**Par 7. &ensp; $\beta$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description**: Wind velocity exponent. Determines the acceleration along a streamline

**Par 8. &ensp; $v_{\rm{turb}}$** </br>
&emsp; &emsp; &#9656; **Units:** km/s </br>
&emsp; &emsp; &#9656; **Description:** Turbulents velocity. Assumed constant throughout the wind. Sets the width of the Gaussian smoothing kernel used to emulate tubulence

**Par 9. &ensp; $\kappa$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description:** Sets the weighting for wind launching efficiency as function of radius

**Par 10. &ensp; Inc.** </br>
&emsp; &emsp; &#9656; **Units:** Degrees </br>
&emsp; &emsp; &#9656; **Description:** Observer inclination (measured from the z-axis)

**Par 11. &ensp; $E_{0}$** </br>
&emsp; &emsp; &#9656; **Units:** keV </br>
&emsp; &emsp; &#9656; **Description:** Rest frame line energy

**Par 12. &ensp; Norm** </br>
&emsp; &emsp; &#9656; **Units:** photons/s/cm$^2$ </br>
&emsp; &emsp; &#9656; **Description:** Normalisation. Sets the total number of photons within the line


### xwindconv
Convolution model. This takes the line profiles from `xwindline` and uses them as a convolution kernel that is then applied to an input spectrum. This can then be used on a series of lines or continuum. Note, the internal normalisation is set to conserve **photon number**. i.e. If the input spectrum being convolved contains 100 photons, then the output specturm will also contain 100 photons.

**Par 1. &ensp; $\log_{10} \dot{m}_{w}$** </br>
&emsp; &emsp; &#9656; **Units:** $\dot{M}_{w} / \dot{M}\_{\rm{Edd}}$ </br>
&emsp; &emsp; &#9656; **Description:** Wind mass-outflow rate, scaled by the Eddington accretion rate 

**Par 2. &ensp; $r_{\rm{in}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Inner launch radius

**Par 3. &ensp; $r_{\rm{out}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Outer launch radius

**Par 4. &ensp; $d_{f}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Distance from origin to wind focus (see Fig. 1 in Hagen et al. 2026)

**Par 5. &ensp; $f_{\rm{cov}}$** </br>
&emsp; &emsp; &#9656; **Units:** $\frac{\Omega}{4 \pi}$ </br>
&emsp; &emsp; &#9656; **Description:** Covering fraction fo the wind as seen from the central (illuminating) source

**Par 6. &ensp; $\log_{10} v_{\infty}$** </br>
&emsp; &emsp; &#9656; **Units:** $c$ </br>
&emsp; &emsp; &#9656; **Description:** Outflow velocity at infinity

**Par 7. &ensp; $r_{v}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Velocity scale length. i.e. The distance along the streamline where the wind reaches half $v_{\infty}$

**Par 7. &ensp; $\beta$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description**: Wind velocity exponent. Determines the acceleration along a streamline

**Par 8. &ensp; $v_{\rm{turb}}$** </br>
&emsp; &emsp; &#9656; **Units:** km/s </br>
&emsp; &emsp; &#9656; **Description:** Turbulents velocity. Assumed constant throughout the wind. Sets the width of the Gaussian smoothing kernel used to emulate tubulence

**Par 9. &ensp; $\kappa$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description:** Sets the weighting for wind launching efficiency as function of radius

**Par 10. &ensp; Inc.** </br>
&emsp; &emsp; &#9656; **Units:** Degrees </br>
&emsp; &emsp; &#9656; **Description:** Observer inclination (measured from the z-axis)


### xwindfe
Additive model, specific for the Fe-K $\alpha$ complex. This calculates self-consistently the equivalent wdith of the line, such that the normalisation within XSPEC **should always be fixed to unity**. Here the absolute line flux is calculated from the number of photons absorbed (assuming an input spectrum) and then re-emitted by the wind, which is fundamentally governed by the wind density profile. The input spectrum is a power-law, which should always be tied to a correspinding fit of the broad-band continuum. Additionally, `xwindfe` uses the 7-Lorentzian Holzer (1997) profile for the rest frame emission. This naturally gives Lorentzian wings, as well as the spin doublet (Fe-K $\alpha_1$ and Fe-K $\alpha_2$).

**Par 1. &ensp; $\log_{10} \dot{m}_{w}$** </br>
&emsp; &emsp; &#9656; **Units:** $\dot{M}_{w} / \dot{M}\_{\rm{Edd}}$ </br>
&emsp; &emsp; &#9656; **Description:** Wind mass-outflow rate, scaled by the Eddington accretion rate 

**Par 2. &ensp; $r_{\rm{in}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Inner launch radius

**Par 3. &ensp; $r_{\rm{out}} $** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Outer launch radius

**Par 4. &ensp; $d_{f}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Distance from origin to wind focus (see Fig. 1 in Hagen et al. 2026)

**Par 5. &ensp; $f_{\rm{cov}}$** </br>
&emsp; &emsp; &#9656; **Units:** $\frac{\Omega}{4 \pi}$ </br>
&emsp; &emsp; &#9656; **Description:** Covering fraction fo the wind as seen from the central (illuminating) source

**Par 6. &ensp; $\log_{10} v_{\infty}$** </br>
&emsp; &emsp; &#9656; **Units:** $c$ </br>
&emsp; &emsp; &#9656; **Description:** Outflow velocity at infinity

**Par 7. &ensp; $r_{v}$** </br>
&emsp; &emsp; &#9656; **Units:** $R_{G}$ </br>
&emsp; &emsp; &#9656; **Description:** Velocity scale length. i.e. The distance along the streamline where the wind reaches half $v_{\infty}$

**Par 7. &ensp; $\beta$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description**: Wind velocity exponent. Determines the acceleration along a streamline

**Par 8. &ensp; $v_{\rm{turb}}$** </br>
&emsp; &emsp; &#9656; **Units:** km/s </br>
&emsp; &emsp; &#9656; **Description:** Turbulents velocity. Assumed constant throughout the wind. Sets the width of the Gaussian smoothing kernel used to emulate tubulence

**Par 9. &ensp; $\kappa$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description:** Sets the weighting for wind launching efficiency as function of radius

**Par 10. &ensp; Inc.** </br>
&emsp; &emsp; &#9656; **Units:** Degrees </br>
&emsp; &emsp; &#9656; **Description:** Observer inclination (measured from the z-axis)

**Par 11. $ensp; $A_{\rm{Fe}}$** </br>
&emsp; &emsp; &#9656; **Units:** $[\rm{Fe}]/[\rm{Fe}_{\odot}]$ </br>
&emsp; &emsp; &#9656; **Description:** Iron abundance relative to solar. Uses the abundance values from Anders & Grevesse 1989

**Par 12. &ensp; $N_{0}$** </br>
&emsp; &emsp; &#9656; **Units:** photons/s/cm^2 at 1keV </br>
&emsp; &emsp; &#9656; **Description:** Normalisation of the incident X-ray power0-law emission

**Par 13. &ensp; $\Gamma$** </br>
&emsp; &emsp; &#9656; **Units:** Dimensionless </br>
&emsp; &emsp; &#9656; **Description:** Photon index of the incident power-law spectrum
