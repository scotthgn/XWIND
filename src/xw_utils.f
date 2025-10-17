c      Contains generic utility routines for XWIND models
c      
c      Included routines:
c      -----------------------------------------------------------------
c      init_egrid:
c         Generates internal energy bins used for calculations.
c         Sets the model energy resolution to Delta E/E0 = 2e-4.
c         Bins are linearly spaced in E/E0.
c         This gives roughly 1.2eV resolution at 6keV.
c
c      do_logconv:
c         Does a direct (i.e no FFTs) log convolution between an input
c         spectrum (ph) and emission line kernel (kern)
c         (or any other kernel for that matter).
c         This treats the photon count in each bin (in the input spectrum)
c         as a delta function centred on E. The output spectrum in bin i
c         is then sum(kern(Ei/Ej) * ph(Ej))
c
c      direct_bin_search:
c         Helper routine for do_logconv. If the energy bins in the input
c         spectrum are not evenly spaced (in either lin or log) then
c         cannot calculate relevant bin indices analytically. So this
c         routine will instead do a direct search through the energy
c         grid during the convolution routine.
c         NOTE: This is ONLY called IF input energy bins are not evenly
c         spaced, as it quickly becomes slow for large energy arrays
c      
c      renorm_line
c         Re-normalised a line profile to unity (i.e such that the total
c         photon flux in a line is 1 photon/s/cm^2
c      
c      do_turbulence
c         Convolves the windline profile with a Gaussian of some velocity
c         width (in km/s) to emulate the effects of trubulence
c      
c      holzer_profile
c         Generates the 7 Lorentzian Holzer profile for the intrinsic 
c         Fe-Kalpha line
c      -----------------------------------------------------------------
c
c      For details on XWIND, see: Hagen, Done & Matzeu (in prep)
c
c      -----------------------------------------------------------------
c      -----------------------------------------------------------------


       subroutine init_egrid(e_ear, ear, ne, E0)
c      Initiates the energy bins used unternally for calculations
       implicit none

       integer ne
       real e_ear(0:ne), ear(0:ne), E0

       real d_ee !Delta E/E0

       integer n !iteration index
       
c      inargs
c      ------
c      e_ear : array
c         Empty bins to be filled with E/E0
c      ear : array
c         Empty bins to be filled with E
c      ne : int
c         Number of energy bins
c      E0 : float
c         Rest frame energy
c       
c      outargs
c      -------
c      e_ear : array
c         Bin edges in E/E0
c      ear : array
c         Bin edges in E (keV)
c
c      START ROUTINE
c      -----------------------------------------------------------------

c      Filling bins
c      Linearly spaced in range E/E0: 0.2 -> 4.4
c      Gives max bulk velocity of 0.9c 
       e_ear(0) = 0.2
       ear(0) = e_ear(0) * E0
       d_ee = 4.2/float(ne)
       do n=1,ne,1
          e_ear(n) = e_ear(0) + d_ee*float(n)
          ear(n) = e_ear(n) * E0
       end do

       return
       end

      

       subroutine do_logconv(ear, ne, ph, e_kern, ph_kern)
c      Does a direct log-convolution between an input spectrum and the
c      emission line kernel
c      Treats the photon count in each bin in the input spectrum as a 
c      delta function centred on E. The output photon count in bin i
c      is then sum(kern(Ei/Ej) * ph(Ej))
c
c      In order to conserve photon number, the windline normalisation 
c      is set such that its sum is equal to 1
       implicit none

       integer ne
       real ear(0:ne), e_kern(0:ne)
       real ph(ne), ph_in(ne), ph_kern(ne)

       real emid_i, emid_j, ee_ij
       real ee_min, ee_max
       integer idx_tran

       logical is_log, is_lin
       real precision
       real erat, crat, dini, di, fchange, dchange
       real ldE, deps

       integer i, j, k

c      inargs
c      ------
c      ear : array
c         Energy bins of input spectrum
c      ne : int
c         Number of energy bins
c      ph : array
c         Input spectrum (i.e photomn flux)
c      e_kern : array
c         Bins in E/E0 for convolution kernel
c         Must have same number (ne) as ear
c      ph_kern : array
c         Convolution kernel (i.e line profile)
c
c      outargs
c      -------
c      ph : array
c         Spectrum after convolution
c
c      START ROUTINE
c      -----------------------------------------------------------------
       
       !Copying input spec, and zeroing output array
       do i=1, ne, 1
          ph_in(i) = ph(i)
          ph(i) = 0.0
       end do

       !checking if input energy grid is log 
       is_log = .true.
       precision = 1e-6 !since using reals
       erat = e_kern(1)/e_kern(0)
       aitr: do i=1, ne, 1
          crat = e_kern(i)/e_kern(i-1)
          fchange = abs(erat-crat)/erat
          if (fchange.gt.precision) then
             is_log = .false.
             exit aitr
          else
             continue
          end if
       end do aitr

       !if not log, then checkin if lin
       if (is_log) then
          ldE = log10(e_kern(1)) - log10(e_kern(0))
       else
          is_lin = .true.
          dini = e_kern(1) - e_kern(0)
          bitr: do i=1, ne, 1
             di = e_kern(i) - e_kern(i-1)
             dchange = abs(dini - di)
             if (dchange.gt.precision) then
                is_lin = .false.
                exit bitr
             else
                continue
             end if
          end do bitr

          if (is_lin) then
             ldE = e_kern(1) - e_kern(0)
          end if
       end if

       !Note - if energy grid is neither log nor lin, then will simply
       !search directly through energy grid for re-binning
       !WARNGIN! This will be VERY slow for large energy arrays!!!
       
       !performing convolution
       !Looping over first energy array
       ee_min = e_kern(0)
       ee_max = e_kern(ne)
       do i=1, ne, 1
          emid_i = 0.5*(ear(i) + ear(i-1))
          do j=1, ne, 1
             emid_j = 0.5*(ear(j) + ear(j-1))
             ee_ij = emid_i/emid_j

             if ((ee_ij.ge.ee_max).or.(ee_ij.le.ee_min)) then
                continue
             else
                if (is_log) then
                   idx_tran = ceiling((log10(ee_ij)-log10(ee_min))/ldE)
                   ph(i) = ph(i) + ph_in(j) * ph_kern(idx_tran)
                
                else if (is_lin) then
                   idx_tran = ceiling((ee_ij-ee_min)/ldE)
                   ph(i) = ph(i) + ph_in(j) * ph_kern(idx_tran)
                
                else
                   call direct_bin_search(e_kern,ne,ph,ph_in,ph_kern,
     $                  ee_ij,i,j)
                end if
             end if
          end do
       end do

       
       return
       end      



       subroutine direct_bin_search(ee_ar,ne,ph,ph_in,ph_kern,ee_ij,i,j)
c      Subroutine for doing direct search through bins and then
c      applying to correct bin
c      This is slow for large energy arrays, so will ONLY be used when
c      input energy grid is not equally spaced in either lin or log
       implicit none

       integer ne
       real ee_ar(0:ne), ph(ne), ph_in(ne), ph_kern(ne)
       real ee_ij
       integer i,j,k

c      inargs
c      ------
c      ee_ar : array
c         Bins on E/E0 for kernel (i.e line profile)
c      ne : int
c         Number of energy bins
c      ph : array
c         Contains convolving spectrum
c      ph_in : array
c         Raw input spectrum (i.e no operations have been acted on it)
c      ph_kern : array
c         Convolution kernel
c      ee_ij : float
c         Current energy bin (in Ei/Ej)
c      i : int
c         Index of first energy bin 
c      j : int
c         Index of second energy bin
c
c     outargs
c     -------
c     ph : array
c         Updated output spectrum
c
c     START ROUTINE
c     ------------------------------------------------------------------
       
       citr: do k=1, ne, 1
          if ((ee_ij.lt.ee_ar(k)).and.(ee_ij.gt.ee_ar(k-1))) then
             ph(i) = ph(i) + ph_in(j)*ph_kern(k)
             exit citr
          else
             continue
          end if
       end do citr

       return
       end



       subroutine renorm_line(ear, ne, ph)
c      Re-normalises such that the integrated line-profile is unity
       implicit none

       integer ne
       real ear(0:ne), ph(ne)
       double precision ph_int
       
       integer i

c      inargs
c      ------
c      ear : array
c         Energy bins 
c      ne : int
c         Number of energy bins
c      ph : array
c         Line profile in photons/s/cm^2/bin
c
c      outargs
c      -------
c      ph : array
c         Normalised line profile in photons/s/cm^2/bin
c
c      START ROUTINE
c      -----------------------------------------------------------------
       
       ph_int = 0.0
       do i=1, ne, 1
          ph_int = ph_int + ph(i)
       end do
       ph = ph*(1.0/ph_int)
       
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

c      inargs
c      ------
c      ear : array
c          Energy bins in keV
c      ne : int
c          Number of bins
c      ph : array
c          Line profile
c      e_ear : array
c          Energy bins in E/E0
c      vturb : float
c          Turbulent velocity in km/s
c
c      outargs
c      -------
c      ph : array
c          Line profile including turbulence
c
c      START ROUTINE
c      -----------------------------------------------------------------

       
       pi = 4.0*atan(1.0)
       
       sig = vturb/3.0e5        !v/c - since units of E/E0 for gaussian kernel
       do n=1,ne,1
          ee_mid = e_ear(n-1) + 0.5*(e_ear(n) - e_ear(n-1))
          gau(n)=(1/sig*sqrt(2.0*pi))*exp(-0.5*(((ee_mid-1)**2)/sig**2))
          norm = norm + gau(n)
       end do
       gau = gau*(1.0/norm)

       call do_logconv(ear,ne,ph,e_ear,gau)

       return
       end
      

       subroutine holzer_profile(ear, nn, ph)
c      Caclulates intrinsic line profile using seven Lorentzian fit from
c      Holtzer et al. 1997
       implicit none

       integer nn
       real ear(0:nn), ph(nn)

       real Ei(7), Wi(7), Ii(7)
       
       integer i, n !iteration indices
       real norm

c      inargs
c      ------
c      ear : array
c         Energy bins (keV)
c      nn : int
c         Number of energy bins
c      ph : array
c         Empty photon array
c
c      outargs
c      -------
c      ph : array
c          Filled photon array, norm unity
c          Holzer profile
c
c      START ROUTINE
c      ----------------------------------------------------------------- 

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

       
c      Renormalising line to unity, since actual wind normalisation 
c      handled by calc_windline (in xw_core.f)
       
       call renorm_line(ear, nn, ph)
          
       return
       end
