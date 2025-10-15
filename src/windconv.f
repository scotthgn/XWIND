c      Routine to smooth an input spectrum
c      Uses windline to calculate a line-profile at some arbitrary
c      rest energy. Then convolves in log energy, since
c      \Delta log E = \Delta log epsilon, where epsilon = E/E0
c
c      CAVEAT! The line-profile is re-normalised to conserve photon
c      number. Hence Mbh and mdot are ONLY used to set the density
c      profile, required for calculating the emissivity


       subroutine windconv(ear,ne,param,ifl,photar,photer)
       use xsfortran
       implicit none
       
       integer npars
       parameter(npars=10)

       integer ne, ifl
       real ear(0:ne), photar(ne), param(npars), photer(ne)
       real oldpars(npars)

       real photar_wnd(ne), wnd_pars(14) !w is wnd_pars. Short since fixed form
       real emid !midpoint in energy range - used for initial windline call
       
       logical parchange, echange
       save oldpars

       integer i
       
c      param(1):    Mbh, Black hole mass, Msol
c      param(2):    mdot_w, wind mass outflow rate, Mdot/Mdot_edd
c      param(3):    r_in, Inner Launch radius, Rg
c      param(4):    r_out, Outer launch radius, Rg
c      param(5):    d_foci, Distance to wind foci, Rg
c      param(6):    fcov, wind covering fraction, Omega/4pi
c      param(7):    vinf, outflow velcoity at infinity, c
c      param(8):    rv, wind velocity scale length, Rg
c      param(9):    vexp, Wind velocity exponent
c      param(10):   kappa, radial density law exponent
c      param(11):   inc, observer inclination, deg

c      Uses windline code to generate line profile for convolution
c      Start by filling internal parameter array
       do i=1, npars, 1
          wnd_pars(i) = param(i)
       end do
       wnd_pars(11) = 0.5*(ear(ne) + ear(0)) !E0 - arbitrary, set to centre of E range

       call windline(ear,ne,wnd_pars,ifl,photar_wnd,photer)
       call logconv_line(ear,ne,photar,ear/wnd_pars(11),photar_wnd)
       
       return
       end


       subroutine logconv_line(ear, ne, ph, e_kern, ph_kern)
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
       real ph(ne), ph_in(ne), ph_kern(ne), ph_kern_tp(ne)

       real emid_i, emid_j, ee_ij
       real ee_min, ee_max
       integer idx_tran

       logical is_log, is_lin
       real precision
       real erat, crat, dini, di, fchange, dchange
       real ldE, deps

       integer i, j, k
       
       !Copying input arrays, and zeroing output array
       do i=1, ne, 1
          ph_in(i) = ph(i)
          ph_kern_tp(i) = ph_kern(i)
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
                   ph(i) = ph(i) + ph_in(j) * ph_kern_tp(idx_tran)
                
                else if (is_lin) then
                   idx_tran = ceiling((ee_ij-ee_min)/ldE)
                   ph(i) = ph(i) + ph_in(j) * ph_kern_tp(idx_tran)
                
                else
                   call direct_bin_search(e_kern,ne,ph,ph_in,ph_kern_tp,
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
       

       


