c      Convolves line profile (including strength) with Holtzer profile
c      for full Fe-Kalpha line shape.
c
c      Generates Holtzer profile here, uses windline for shape and
c      strength, then convolution routine in windconv to calculate
c      final line shape


       subroutine windfe(ear,ne,param,ifl,photar,photer)
       use xsfortran
c      Calculates line profile, convolves with Holtzer profile, returnes
c      Fe Kalpha line
       implicit none
       
       integer npars, cpars
       parameter(npars=14) !pars for this subflavour of model
       parameter(cpars=14) !pars needed for the calc_windline

       integer ne, ifl
       real ear(0:ne),photar(ne),param(npars),photer(ne),cparam(cpars)
       real oldpars(npars)

       integer Nnew
c       parameter(Nnew=8400)     !num new Ebins
       parameter(Nnew=21000)
       real e_enew(0:Nnew),ph(Nnew),ph_wnd(Nnew),Enew(0:Nnew)
       real ph_err(Nnew)

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)

       logical parchange, echange
       save oldpars

       integer i,n

       integer num_bins
       
c      param(1):    mdot_w, Wind mass outflow rate, Mdot/Mdot_edd
c      param(2):    r_in, inner luanch radius, Rg
c      param(3):    r_out, outer launch radius, Rg
c      param(4):    d_foci, distance to wind focus, Rg
c      param(5):    fcov, wind covering fraction
c      param(6):    vinf, outflow velocity at infinity
c      param(7):    rv, wind velocity scale length, Rg
c      param(8):    vexp, wind velocity exponent
c      param(9):    vturb, turbulent velcoity, km/s
c      param(10):   kappa, radial density law exponent
c      param(11):   inc, observer inclination, deg
c      param(12):   Afe, iron abundance, relative to solar
c      param(13):   N0, incident power-law normalisation at 1keV
c                       units: ph/s/cm^2/keV
c      param(14):   Gamma, incident power law photon index


c      Defining internal energy grid
c      same as in windline.f

       e_enew(0) = 0.2
       Enew(0) = e_enew(0) * 6.4 !centre for now
       do n=1, Nnew, 1
*          e_enew(n) = e_enew(0) + 5.0e-4*float(n)
          e_enew(n) = e_enew(0) + 2e-4*float(n)
          Enew(n) = e_enew(n) * 6.4
       end do


c      generating holtzer profile
       call holtzer(Enew, Nnew, ph)

c      calculating windline profile (inc. normalisation)
c      First fillinf parameter array
       do i=1, 8, 1
          cparam(i) = param(i)
       end do
       cparam(9) = param(10)
       cparam(10) = param(11)
       cparam(11) = param(12)
       cparam(12) = 6.4         !central energy
       cparam(13) = param(13)
       cparam(14) = param(14)

       call calc_windline(e_enew, Nnew, cparam, ifl, ph_wnd)
       
c      convolving windline with holtzer profile
       call logconv_line(Enew,Nnew,ph,e_enew,ph_wnd)

c      Adding turbulence if greater than 0
c      log convolves a gaussian (so some constant velocity width)
c      Writing own as gmsooth in c++ (argh!!) 
       if (param(9).gt.0.0) then
          call do_turbulence(Enew,Nnew,ph,e_enew,param(9))
       end if
          
       
       call inibin(Nnew, Enew, ne, ear, istart, iend, fstart, fend, 0)
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar)


       return
       end


       subroutine holtzer(ear, nn, ph)

c      Caclulates intrinsic line profile using seven Lorentzian fit from
c      Holtzer et al. 1997
       implicit none

       integer nn
       real ear(0:nn), ph(nn)

       real Ei(7), Wi(7), Ii(7)
       
       integer i, n !iteration indices
       real norm

c      Filling Lorentzian pars with values from Holtzer et al 1997
       Ei(1) = 6.404148    !E11, keV
       Ei(2) = 6.403295    !E12
       Ei(3) = 6.400653    !E13
       Ei(4) = 6.402077    !E14
       Ei(5) = 6.391190    !E21
       Ei(6) = 6.389106    !E22
       Ei(7) = 6.390275    !E23

       Wi(1) = 1.613e-3    !W11, keV
       Wi(2) = 1.965e-3    !W12
       Wi(3) = 4.833e-3    !W13
       Wi(4) = 2.803e-3    !W14
       Wi(5) = 2.487e-3    !W21
       Wi(6) = 2.339e-3    !W22
       Wi(7) = 4.433e-3    !W23

       Ii(1) = 0.697       !I11, Intensity
       Ii(2) = 0.376       !I12
       Ii(3) = 0.088       !I13
       Ii(4) = 0.136       !I14
       Ii(5) = 0.339       !I21
       Ii(6) = 0.060       !I22
       Ii(7) = 0.102       !I23


c      Initialising (i.e zeroing) output array
       do n=1,nn,1
          ph(n) = 0.0
       end do

       
c      Generating seven Lorentzian line component (in intensity/energy units!)
       do i=1,7,1
          do n=1,nn,1
             ph(n) = ph(n) + Ii(i)/(1+((ear(n)-Ei(i))/Wi(i))**2.0)
          end do
       end do

       
c      Changin to units of ph/s/keV (since profile in energy units)
c      The normalising to unity (since windline calculation handles total flux)
c      Doing as sum, since want units of ph/keV (so this gives the dE)
       
       norm = 0.0
       do n=1,nn,1
          norm = norm + ph(n)
       end do
       ph = ph*(1.0/norm)
          
       return
       end



       subroutine do_turbulence(ear, ne, ph, e_ear, vturb)
c      Generates a gaussian with some fixed energy independent velocity
c      width (with norm unity). Then log convolves with windline
       implicit none

       integer ne
       real ear(0:ne), e_ear(0:ne), ph(ne), gau(ne)
       real ee_mid
       real vturb, sig, pi, norm

       integer n

       pi = 4.0*atan(1.0)
       
       sig = vturb/3.0e5        !v/c - since units of E/E0 for gaussian kernel
       do n=1,ne,1
          ee_mid = e_ear(n-1) + 0.5*(e_ear(n) - e_ear(n-1))
          gau(n)=(1/sig*sqrt(2.0*pi))*exp(-0.5*(((ee_mid-1)**2)/sig**2))
          norm = norm + gau(n)
       end do
       gau = gau*(1.0/norm)

       call logconv_line(ear,ne,ph,e_ear,gau)

       return
       end
