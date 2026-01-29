c      Contains the individual model subroutines that XSPEC sees
c
c      Included model routines (note they do not calculate errors)
c      -----------------------------------------------------------------
c      xwindline:
c          Additive model
c          Single line profile from xwind, with rest energy as an input.
c          Always normalised such that the total line flux is 1 photon/s/cm^2.
c          Can in theory be used for any emission line (though note that 
c          this assumes the same fluoresence profile as calculated for
c          neutral Fe-Kalpha as this is needed to give the emissivity
c          through the wind)
c
c      xwindfe:
c          Additive model
c          Line profile specific to Fe-Kalpha.
c          This includes the intrinsic Holzer profile for the rest frame
c          emission (giving both Kalpha1 and 2) as well as a physical 
c          normalisation calculated from the illuminating spectrum and
c          wind density profile.
c          For this model XSPEC norm parameter should always be fixed to 
c          1, since normalisation calcualted internally.
c          Note the illuminating spectrum is assumed a power-law between
c          7 and 40keV (beyond 40keV and below 7keV the interaction cross
c          section is 0 or negligible). Convention used is identical to 
c          XSPEC powerlaw model, so parameters can be tied if you wish
c          to jointly fit both the line profile and continuum
c
c      xwindconv:
c          Convolution model
c          Convolves the xwindline profile with any input spectrum.
c          Model preserves photon count (i.e convolution kernel is 
c          normalised to units)
c      -----------------------------------------------------------------
c
c      If you use XWIND please cite: Hagen et al. (2026, submitted to A&A)
c
c      -----------------------------------------------------------------
c      -----------------------------------------------------------------
      

       subroutine xwindline(ear,ne,param,ifl,photar,photer)
c      Gives line profile centred on some energy 0
c      Re-normsalised to unity
       implicit none

       integer npars, cpars
       parameter(npars=12)
       parameter(cpars=14)

       integer ne, ifl
       real ear(0:ne),photar(ne),param(npars),photer(ne),cparam(cpars)
       real oldpars(npars)

       integer Nnew
       parameter(Nnew=21000)    !sets internal resolution to 1.2eV at 6keV
       real e_enew(0:Nnew), ph(Nnew), Enew(0:Nnew)

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)

       logical parchange
       save oldpars

       integer i, n

c      param(1):    log10 mdot_w, wind mass-outflow rate, Units: Mdot/Mdot_edd
c      param(2):    r_in, inner launch radius, Units: Rg
c      param(3):    r_out, outer launch radius, Units: Rg
c      param(4):    d_foci, distance to wind focus (from BH), Units: Rg
c      param(5):    fcov, wind covering fraction (as seen from BH), Units: Omega/4pi
c      param(6):    log10 vinf, outflow velocity at infinity, Units: c   
c      param(7):    rv, wind velocity scale length, Units: Rg
c      param(8):    beta, wind velocity exponent
c      param(9):    vturb, intrinsic turbulence, Units: km/s
c      param(10):    kappa, radial density law exponent
c      param(11):   inc, observer inclination (measured from z-axis), Units: deg
c      param(12):   E0, rest frame line energy, Units: keV
c
c      START ROUTINE
c      -----------------------------------------------------------------

c      First checking if parameters have changed
       parchange = .false.
       do i=1,npars,1
          if (param(i).ne.oldpars(i)) then
             parchange=.true.
          end if
          oldpars(i) = param(i)
       end do

       
c      generating internal energy bins (see xw_utils.f)
       call init_egrid(e_enew, Enew, Nnew, param(12))

c      performing main calculation IF parameters have changed
       if (parchange) then
c         Filling parameter array for line-calculation
c         Actual line calculation done by calc_windline in xw_core.f
          do i=1, 8, 1          !stopping loop before vturb, since not done by calc_windline
             cparam(i) = param(i) 
          end do
          cparam(9) = param(10) !kappa
          cparam(10) = param(11) !inc
          cparam(11) = 1.0      !Afe, doesn't matter since re-normalise later
          cparam(12) = param(12) !E0
          cparam(13) = 1.0      !inspec norm - arbitrary since re-norm later
          cparam(14) = 2.0      !inspec gamma - arbitrary here

c         calculating line
          call calc_windline(e_enew, Nnew, cparam, ph) !local, in xw_core.f

c         Adding tubulence
          if (param(9).gt.0.0) then
             call do_turbulence(Enew,Nnew,ph,e_enew,param(9)) !local, in xw_utils.f
          end if
       end if
       
c      re-binning and re-normalising
       call inibin(Nnew, Enew, ne, ear, istart, iend, fstart, fend, 0) !xspec routine
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar) !xspec routine
       call renorm_line(ear, ne, photar) !local, in xw_utils.f

       return
       end



       subroutine xwindfe(ear,ne,param,ifl,photar,photer)
c      Caclulates windline profile specifically for Fe-Kalpha.
c      Uses the physical normalisation based on wind profile + incident spectrum.
c      Also uses Holzer et al. (1997) intrinsic line profile to account
c      for intrinsic line width and Kalpha1/2 doublet
       implicit none

       integer npars, cpars
       parameter(npars=14)
       parameter(cpars=14)

       integer ne, ifl
       real ear(0:ne),photar(ne),photer(ne),param(npars),cparam(npars)
       real oldpars(npars)

       integer Nnew
       parameter(Nnew=21000)
       real e_enew(0:Nnew),Enew(0:Nnew),ph(Nnew),ph_wnd(Nnew)

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)

       logical parchange
       save oldpars

       integer i

c      param(1):    log10 mdot_w, wind mass-outflow rate, Units: Mdot/Mdot_edd
c      param(2):    r_in, inner launch radius, Units: Rg
c      param(3):    r_out, outer launch radiu, Units: Rg
c      param(4):    d_foci, distance to wind focus, Units: Rg
c      param(5):    fcov, wind covering fraction (as seen from BH), Units: Omega/4pi
c      param(6):    log10 vinf, outflow velocity at infinity, Units: c
c      param(7):    rv, wind velcity scale length, Rg
c      param(8):    beta, wind velocity exponent,
c      param(9):    vturb, intrinsic turbulence, Units: km/s
c      param(10):   kappa, radial density law exponent
c      param(11):   inc, observer inclination, Units: deg
c      param(12):   Afe, iron abundance, Units: [Fe]/[Fe solar]
c      param(13):   N0, incident power law normalisation, Units ph/s/cm^2/keV at 1 keV
c      param(14):   Gamma, incident power law photon index
c
c      START ROUTINE
c      -----------------------------------------------------------------

c      checking if model parameters have changed
       parchange=.false.
       do i=1,npars,1
          if (param(i).ne.oldpars(i)) then
             parchange=.true.
          end if
          oldpars(i) = param(i)
       end do
       
c      initiating internal energy bins
       call init_egrid(e_enew, Enew, Nnew, 6.4) !local, in xw_utils.f

c      doing calculation IF parameters have changed
       if (parchange) then
c          generating holzer profile for rest frame emission
          call holzer_profile(Enew, Nnew, ph)

c         filling parameter array for calculation of line profile
          do i=1, 8, 1
             cparam(i) = param(i)
          end do
          cparam(9) = param(10)
          cparam(10) = param(11)
          cparam(11) = param(12)
          cparam(12) = 6.4
          cparam(13) = param(13)
          cparam(14) = param(14)
       
c         generating xwindline profile, including correct norm
          call calc_windline(e_enew, Nnew, cparam, ph_wnd) !local, in xw_core.f

c         convolving wind profile with holzer profile to get observed frame Fe-Kalpha
          call do_logconv(Enew,Nnew,ph,e_enew,ph_wnd) !local, in xw_utils.f

c         Inclduing turbulence if non 0
          if (param(9).gt.0.0) then
             call do_turbulence(Enew,Nnew,ph,e_enew,param(9)) !local, in xw_utils.f
          end if
       end if
          
c      re-binning back to xspec energy grid
       call inibin(Nnew,Enew,ne,ear,istart,iend,fstart,fend,0)
       call erebin(Nnew,ph,ne,istart,iend,fstart,fend,photar)

       return
       end



       subroutine xwindconv(ear,ne,param,ifl,photar,photer)
c      Convolves an input spectrum (ph) with the xwindline profile
c      Normalised to conserve photon flux of the input specrum
       implicit none

       integer npars, cpars
       parameter(npars=11)
       parameter(cpars=12)

       integer ne, ifl
       real ear(0:ne), photar(ne), param(npars), photer(ne)

       real photar_wnd(ne), cparam(cpars)
       real emid

       integer i

c      param(1):    log10 mdot_w, wind mass-outflow rate, Units: Mdot/Mdot_edd   
c      param(2):    r_in, inner launch radius, Units: Rg
c      param(3):    r_out, outer launch radius, Units: Rg
c      param(4):    d_foci, distance to wind focus, Units: Rg
c      param(5):    fcov, wind covering fraction, Units: Omega/4pi
c      param(6):    log10 vinf, outflow velocity at infinity, Units: c
c      param(7):    rv, wind velocity scale length, Units: Rg
c      param(8):    beta, wind velocity exponent
c      param(9):    vturb, wind trubulence, Units: km/s
c      param(10):   kappa, radial density law
c      param(11):   inc, observer inclination
c
c      START ROUTINE
c      -----------------------------------------------------------------
       
c      filling parameter array for xwindline
       do i=1, npars, 1
          cparam(i) = param(i)
       end do
       cparam(12) = 0.5*(ear(ne) + ear(0)) !E0, arbitrary so set to centre of energy array

c      generating line profile and then convolving with input spec
       call xwindline(ear,ne,cparam,ifl,photar_wnd,photer)
       call do_logconv(ear,ne,photar,ear/cparam(12),photar_wnd)

       return
       end
