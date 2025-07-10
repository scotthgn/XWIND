c      Model for calculating the line-emission from a outflowing wind.     
c      Assumes a bi-conical structure, launched between r_in and r_out
c      reacing some velcotiy v_inf at infinity. The azmiuthal velocity
c      is assumed Keplerian at the base, and calculated by conserving
c      angular momentum along a trajectory
c
c      To calculate the line profile, the code will first calculate the
c      density within a cell in a wind, which is used to estimate the
c      fluorescent emission (i.e by integrating the input spectrum 
c      through the cell). It then calculates the expected energy shift
c      in the cell, applies relativstic boosting effects, and then 
c      adds photon flux to relevant energy bin
c     
c      For an input spectru, a simple power-law, calculated between
c      7 and 40 keV, is used.
c
c      Note: internal resolution set to give 3eV at 6keV
c      (oversample XRISM resolution)



       subroutine windline(ear,ne,param,ifl,photar,photer)
c      Returns line profile centered on input energy E0
c      Re-normalised to unity
       implicit none

       integer npars, cpars
       parameter(npars=12)
       parameter(cpars=15)
       
       integer ne, ifl
       real ear(0:ne), photar(ne),param(npars),photer(ne),cparam(cpars)
       real oldpars(npars)

       integer Nnew
       parameter(Nnew=8400)
       real e_enew(0:Nnew), ph(Nnew), Enew(0:Nnew) !e_enew defined as Eobs_Eem

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)
       
       logical parchange, echange
       save oldpars

       integer i, n

c      param(1):    Mbh, Black hole mass, Msol
c      param(2):    mdot_w, wind mass outflow rate, Mdot/Mdot_edd 
c      param(3):    r_in, Inner launch radius, Rg
c      param(4):    r_out, Outer launch radius, Rg
c      param(5):    d_foci, Distance to wind foci, Rg
c      param(6):    fcov, Covering fraction of wind 
c      param(7):    vinf, outflow velocity at infinty, c
c      param(8):    rv, Wind velocity scale-length, Rg
c      param(9):    vexp, Wind velocity exponent
c      param(10):   kappa, Radial density law exponent
c      param(11):   inc, Observer inclination, deg
c      param(12):   E0, Rest frame line energy, keV

c      Defining internal energy grid
c      Linearly spaced in Eobs_Eem, to give 3eV at 6keV
c      Implies delta Eobs_Eem = 5e-4
c      Range of 0.2 -> 4.4 gives max bulk velocity 0.9c
       e_enew(0) = 0.2
       Enew(0) = e_enew(0) * param(12)
       do n=1, Nnew, 1
          e_enew(n) = e_enew(0) + 5e-4*float(n)
          Enew(n) = e_enew(n) * param(12)
       end do

c      Filling parameter array for line-calculation
       do i=1, npars-1, 1
          cparam(i) = param(i)
       end do
       cparam(12) = 1.0 !Afe
       cparam(13) = param(12) !E0
       cparam(14) = 1.0 !incident spec norm
       cparam(15) = 2.0 !incident spec gamma
       
c      calculating line
       call calc_windline(e_enew,Nnew,cparam,ifl,ph)

c      re-binning to xspec grid
       call inibin(Nnew, Enew, ne, ear, istart, iend, fstart, fend, 0)
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar)

c      re-normalising
       call re_norm(ear,ne,photar)

       return
       end
       
       
       subroutine windlinefk(ear,ne,param,ifl,photar,photer)
       implicit none

       integer npars, cpars
       parameter(npars=14)
       parameter(cpars=15)

       integer ne, ifl
       real ear(0:ne),photar(ne),param(npars),photer(ne),cparam(cpars)
       real oldpars(npars)

       integer Nnew
       parameter(Nnew=8400)
       real e_enew(0:Nnew), ph(Nnew), Enew(0:Nnew) !e_enew defined as Eobs_Eem

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)
       
       logical parchange, echange
       save oldpars

       integer i, n

c      param(1):    Mbh, Black hole mass, Msol
c      param(2):    mdot_w, wind mass outflow rate, Mdot/Mdot_edd 
c      param(3):    r_in, Inner launch radius, Rg
c      param(4):    r_out, Outer launch radius, Rg
c      param(5):    d_foci, Distance to wind foci, Rg
c      param(6):    fcov, Covering fraction of wind 
c      param(7):    vinf, outflow velocity at infinty, c
c      param(8):    rv, Wind velocity scale-length, Rg
c      param(9):    vexp, Wind velocity exponent
c      param(10):   kappa, Radial density law exponent
c      param(11):   inc, Observer inclination, deg
c      param(12):   Afe, iron abundance, relative to solar
c      param(13):   N0, incident power-law normalisation at 1 keV
c                       has units: ph/s/cm^2/keV
c      param(14):   Gamma, incident power-law photon index


c      Defining internal energy grid
c      Linearly spaced in Eobs_Eem, to give 3eV at 6keV
c      Implies delta Eobs_Eem = 5e-4
c      Range of 0.2 -> 4.4 gives max bulk velocity 0.9c
       e_enew(0) = 0.2
       Enew(0) = e_enew(0) * 6.4
       do n=1, Nnew, 1
          e_enew(n) = e_enew(0) + 5e-4*float(n)
          Enew(n) = e_enew(n) * 6.4
       end do

c      Filling param array
       do i=1, 12, 1
          cparam(i) = param(i)
       end do
       cparam(13) = 6.4 !central energy
       cparam(14) = param(13)
       cparam(15) = param(14)
       
c      calling model
       call calc_windline(e_enew, Nnew, cparam, ifl, ph)

c      Re-binning to xspec grid
       call inibin(Nnew, Enew, ne, ear, istart, iend, fstart, fend, 0)
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar)

       return
       end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Main subroutine for calculating wind line profile and emission
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       subroutine calc_windline(ees, nn, param, ifl, ph)
       implicit none

       integer nn, ifl
       real ees(0:nn), ph(nn), param(*)
       integer ebin_idx

c      Input Parameters
       double precision Mbh, mdot_w, r_in, r_out, d_foci
       double precision fcov, vinf, rv, vexp, kappa, inc, Afe
       double precision E0, N0, Gamma

c      System (internal) parameters
       double precision rmax, lw, r_base, dr_lin, cosB
       double precision v_l, v_phi, v_r, v_z
       double precision ndens, rbase, dm_dA
       double precision cos_th_grd, r_grd, phi_grd, rmax_th, rmid

c      Incidient emission (power-law)
       integer nePow
       parameter(nePow=300)
       double precision esPow(0:nePow), ph_pow(nePow) 
       double precision dlogE_pow, emid

c      Absorption and transmittance
       double precision absCoeff(nePow), crossSec(nePow)
       double precision phIn(nePow)
       double precision NfeK, Nfe_tot, NfeK_obs, eshift
       double precision tpars(3)
       
c      Function values
       double precision r_intercept, cos_beta, l_stream, energy_shift
       
c      Constants
       double precision pi, sigma0, Ek, alpha, c, G, mp
       double precision Mdot_edd, Rg
       
c      Resolution controls
       double precision dcos_theta, dphi, dlog_r

c      iteration indices
       integer i, j, k
       
c      Strings for writing to terminal
       character(30) fc_str

       
c      Setting physical constants
       pi = 4.0*atan(1.0)
       alpha = 2.67          !Absorption cross-section exponent
       sigma0 = 3.37d-20     !Thompson cross-section, cm^-2 
       Ek = 7.1              !Fe-K absorptin edge, keV
       c = 3.0d+10           !light speed, cm/s
       G = 6.67d-8 * 1.99d33 !Grav const, cm^-3 s^-1 Msol^-1
       mp = 1.67d-24         !proton mass, g

c      Setting wind resolution
       dlog_r = 0.01
       dcos_theta = 0.002
       dphi = 0.001

c      Reading parameters
       Mbh = dble(param(1))
       mdot_w = dble(param(2))
       r_in = dble(param(3))
       r_out = dble(param(4))
       d_foci = dble(param(5))
       fcov = dble(param(6))
       vinf = dble(param(7))
       rv = dble(param(8))
       vexp = dble(param(9))
       kappa = -1.0*dble(param(10))
       inc = (pi * dble(param(11)))/180.0
       Afe = 4.68d-5 * dble(param(12)) !uses abund from Anderson & Grevesse 1989
       E0 = dble(param(13))
       N0 = dble(param(14))
       Gamma = dble(param(15))


c      Defining system pars
       Rg = (G*Mbh)/c**2
       Mdot_edd = (1.39d38*Mbh)/(0.057*c**2)


c      Checking geometric limits and setting boundary
       if (acos(fcov).le.atan(r_in/d_foci)) then
          fcov = cos(atan(r_in/d_foci) + 0.1)

          write(fc_str, '(f10.2)') fcov
          call xwrite('Input fcov never reached for inpu wind!!', 15)
          call xwrite('Re-setting to: '//fc_str//'', 15)
       end if
       rmax = r_intercept(acos(fcov), r_in, d_foci)

c      Defining incident power-law over range 7.1-40 keV
c      Also generating asborption cross-sec
       dlogE_pow = log10(100.0/7.0)/float(nePow)
       esPow(0) = 7.0
       do i=1, nePow, 1
          esPow(i) = 10**(log10(7.0) + dlogE_pow * float(i))
          emid = 10**(log10(esPow(i)) - 0.5*dlogE_pow) !geometric mid
          ph_pow(i) = N0*emid**(-Gamma) !ph/s/cm^2/keV
          if (emid.lt.Ek) then
             crossSec(i) = 0.0
          else
             crossSec(i) = sigma0 * (emid/Ek)**(-alpha) !cm^2
          end if
       end do

c      Initiating photon array (setting 0 since additive model)
       do i=1, nn, 1
          ph(i) = 0.0
       end do
       
c      -----------------------------------------------------------------
c      ACTUAL CALCULATION
c      -----------------------------------------------------------------
       cos_th_grd = dcos_theta*0.5 !evaluating at midpoint
       Nfe_tot = 0.0
       do while(cos_th_grd.lt.fcov)
          !re-setting incident spectrum
          do i=1, nePow, 1
             phIn(i) = ph_pow(i)
          end do

          !Calcalating along a sight-line
          r_grd = r_intercept(acos(cos_th_grd), r_in, d_foci)
          rmax_th = min(rmax,r_intercept(acos(cos_th_grd),r_out,d_foci))
          do while(r_grd.lt.rmax_th)
             dr_lin = 10**(log10(r_grd)+dlog_r) - r_grd
             rmid = r_grd + 0.5*dr_lin !evaluate at centre of bin

             !streamline pars
             cosB = cos_beta(cos_th_grd, rmid, d_foci)
             r_base = d_foci * tan(acos(cosB))
             lw = l_stream(acos(cos_th_grd), rmid, r_base)
             
             !Calculating wind velocities
             v_phi = (r_base/rmid) * sqrt(1/r_base)
             v_l = vinf * (1 - (rv/(rv+lw)))**vexp
             v_r = v_l * sin(acos(cosB))
             v_z = v_l * cosB

             !Wind density (at grid)
             dm_dA = mdot_w*Mdot_edd*(kappa + 2) * r_base**(kappa)
             dm_dA = dm_dA/(4*pi * (r_out**(kappa+2) - r_in**(kappa+2)))
             dm_dA = dm_dA * Rg**(-2)

             ndens = (1/(1.23*mp*v_l*c)) * dm_dA !cm^-3

             !Updating transmitted spectrum and calculating fluoresence
             tpars(1) = Afe
             tpars(2) = ndens
             tpars(3) = dr_lin*Rg
             call do_tran(esPow,nePow,phIn,crossSec,NfeK,tpars)
             NfeK = NfeK*dcos_theta*dphi*0.3 !K-alpha yield in one cell, ph/s/cm^2
             !Nfe_tot = Nfe_tot + NfeK

             !Looping over phi, and calculating energy shift at each azimuth
             phi_grd = dble(0.0) + 0.5*dphi !evaluate at centre of bin
             do while(phi_grd.lt.2.0*pi)
                eshift = energy_shift(inc,rmid,phi_grd,v_r,v_z,v_phi)
                NfeK_obs = eshift**(3.0) * NfeK/(4.0*pi)

                !Placing in bin
                !Since grid exactly known, bin idx is analytic
                ebin_idx = ceiling((eshift - 0.2)/(5.0d-4))
                ph(ebin_idx) = ph(ebin_idx) + NfeK_obs !ph/s/cm^2/bin

                phi_grd = phi_grd + dphi
             end do
             r_grd = r_grd + dr_lin
          end do
          cos_th_grd = cos_th_grd + dcos_theta
       end do

       !write(*,*) Nfe_tot
       
       return
       end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Helper routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      

       subroutine do_tran(es,ne,ph_in,cross_sec,N_feK,pars)
c      Does the transfer of a input spectrum through a wind cell
c      Updates ph_in to be the transmitted spectrum 
c      Also updates the fluoresence (i.e difference between in and out spec)
c
c      pars(1): Afe
c      pars(2): ndens
c      pars(3): dR
       implicit none

       integer ne
       double precision es(0:ne), ph_in(ne), pars(*)
       double precision abs_coeff(ne), cross_sec(ne)
       double precision N_feK
       double precision Afe, ndens, dR

       integer i

       Afe = pars(1)
       ndens = pars(2)
       dR = pars(3)

       N_feK = 0.0
       do i=1, ne, 1
          abs_coeff(i) = exp(-cross_sec(i)*Afe*ndens*dR)
          N_feK = N_feK + (es(i)-es(i-1))*ph_in(i)*(1 - abs_coeff(i))
          ph_in(i) = ph_in(i)*abs_coeff(i)
       end do
       
       return
       end


       subroutine re_norm(ear,ne,ph)
c      Re-normalises such that the integrated line-profile is unity
       implicit none

       integer ne
       real ear(0:ne), ph(ne)
       double precision ph_int
       
       integer i

       ph_int = 0.0
       do i=1, ne, 1
          ph_int = ph_int + ph(i)
       end do
       ph = ph*(1.0/ph_int)
       
       return
       end


      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Model functions
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      

       function energy_shift(inc, r, phi, vr, vz, vphi)
c      Calculates the fractional energy shift of a photon between the 
c      emitted frame and observed frame (i.e Eobs/Eem)
c
c      Derived in the Schwarzchild metric, but still assuming weak field
c      s.t photons traevl in straight (ish) lines
       implicit none
       double precision inc, phi, r, vr, vz, vphi
       double precision gamma, eshift_inv, energy_shift

       gamma = 1/sqrt(1 - (vr**2 + vz**2 + vphi**2)) !Lorentz factor

       eshift_inv = sin(inc) * (vr*cos(phi) - vphi*sin(phi))
       eshift_inv = eshift_inv + vz*cos(phi)
       eshift_inv = eshift_inv*(-1.0) + (1.0-(1.0/(4.0*r**2.0)))**(2.0)
       eshift_inv = eshift_inv*gamma*(1.0+(1.0/(2.0*r)))**(4.0)

       energy_shift = 1/eshift_inv
       return
       end
       

       function r_intercept(theta, r_l, d_foci)
c      Function to calculate the intercept between a line of sight 
c      and a streamline launched at r_l
       implicit none
       double precision theta, r_l, d_foci
       double precision r_intercept

       r_intercept = r_l/(sin(theta) - (r_l/d_foci)*cos(theta))
       return
       end


       function cos_beta(cos_theta, r_grd, d_foci)
c      Function to calculate the laucnh radius of a streamline
c      passing through the point (r_grd, cos_theta) in the wind
       implicit none
       double precision cos_theta, r_grd, d_foci
       double precision cos_beta, root_term

       root_term = sqrt(r_grd**2 + 2*r_grd*d_foci*cos_theta + d_foci**2)
       cos_beta = (d_foci + r_grd*cos_theta)/root_term

       return
       end


       function l_stream(theta, r_grd, r_l)
c      Calculates length along streamline for the point (r_grd,theta)
       implicit none
       double precision theta, r_grd, r_l
       double precision l_stream

       l_stream = sqrt(r_grd**2 + r_l**2 - 2*r_grd*r_l*sin(theta))
       return
       end
