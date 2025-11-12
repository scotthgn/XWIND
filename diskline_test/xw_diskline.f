c      Routine for doing line profile for a Keplerian disc
c      Exists as a validation routine for the general line-profile
c      calculations
c      -----------------------------------------------------------------


       subroutine xwdisk(ear,ne,param,ifl,photar,photer)
       use xsfortran

       implicit none

       integer npars
       parameter(npars=4)

       integer ne, ifl
       real ear(0:ne),photar(ne),param(npars),photer(ne)

       integer Nnew
       parameter(Nnew=21000)
       real e_enew(0:Nnew),ph(Nnew),Enew(0:Nnew)

       real fstart(ne), fend(ne)
       real istart(ne), iend(ne)

       real ph_int
       integer n

c      param(1):   r_in
c      param(2):   r_out
c      param(3):   inc
c      param(4):   E0
c

c      Generating new bins
       e_enew(0) = 0.2
       Enew(0) = e_enew(0) * param(4)
       do n=1,Nnew,1
          e_enew(n) = e_enew(0) + (4.2/float(Nnew)) * float(n)
          Enew(n) = e_enew(n) * param(4)
       end do


       call calc_xwdiskline(e_enew, Nnew, param, ph)

c      renorm
       ph_int = 0.0
       do n=1,Nnew,1
          ph_int = ph_int + ph(n)
       end do
       ph = ph*(1.0/ph_int)


       call inibin(Nnew,Enew,ne,ear,istart,iend,fstart,fend,0)
       call erebin(Nnew,ph,ne,istart,iend,fstart,fend,photar)

       return
       end



       subroutine calc_xwdiskline(ees, nn, param, ph)
       use xsfortran
       implicit none

       integer nn
       real ees(0:nn), ph(nn), param(*)
       integer ebin_idx

c      input pars
       real r_in, r_out, inc, E0

c      system pars
       real rmid, dlog_r_c, Nrbins, dr
       real phi_mid, dphi_c, Nphbins, vphi_r
       real Nphot_em, Nphot_obs, eshift

       real energy_shift
      
c      constants
       real pi, c, G

c      resolution constrols
       real dlog_r, dphi

c      iteration
       integer i, j


c      Setting physical constants]
       pi = 4.0*atan(1.0)
       c = 3.0e10               !light speed, cm/s
       G = 6.67e-8 * 1.99e30    !grav const

c      setting resolution
       dlog_r = 0.001
       dphi = 0.0001

c      reading pars
       r_in = param(1)
       r_out = param(2)
       inc = (pi * param(3))/180.0
       E0 = param(4)


c      inititaing phootn array
       do i=1,nn,1
          ph(i) = 0.0
       end do


c      finding number of evenly spaced bins
       Nrbins = ceiling((log10(r_out) - log10(r_in))/dlog_r)
       dlog_r_c = (log10(r_out) - log10(r_in))/Nrbins

       Nphbins = ceiling((2.0*pi)/dphi)
       dphi_c = (2.0*pi)/Nphbins

       
c      CALCULATION
       do i=1, int(Nrbins), 1
          rmid = 10**(log10(r_in)+float(i-1)*dlog_r_c+dlog_r_c/2.0)
          dr = 10**(log10(rmid)+dlog_r_c/2.0)
          dr = dr - 10**(log10(rmid)-dlog_r_c/2.0)

          vphi_r = sqrt(1/rmid)
          Nphot_em = 1.0 * rmid**(-3.0)
          Nphot_em = Nphot_em * rmid * dr * dphi
          do j=1, int(Nphbins), 1
             phi_mid = 0.0 + float(j-1)*dphi_c+dphi_c/2.0
             
             eshift = energy_shift(inc, rmid, phi_mid, 0.0, 0.0, vphi_r)
             Nphot_obs = eshift**(3.0) * Nphot_em

             ebin_idx = ceiling((eshift-0.2)/(2.0e-4))
             ph(ebin_idx) = ph(ebin_idx) + Nphot_obs
          end do
       end do


       return
       end   


       function energy_shift(inc, r, phi, vr, vz, vphi)
c      Calculates the fractional energy shift of a photon between the 
c      emitted frame and observed frame (i.e Eobs/Eem)
c
c      Derived in the Schwarzchild metric, but still assuming weak field
c      s.t photons traevl in straight (ish) lines
       implicit none
       real inc, phi, r, vr, vz, vphi
       real gamma, eshift_inv, energy_shift

       gamma = 1/sqrt(1 - (vr**2 + vz**2 + vphi**2)) !Lorentz factor

       eshift_inv = sin(inc) * (vr*cos(phi) - vphi*sin(phi))
       eshift_inv = eshift_inv + vz*cos(phi)
       eshift_inv = eshift_inv*(-1.0) + 1.0
       eshift_inv = eshift_inv*gamma*(1.0 - (2.0/r))**(-0.5)
       
       energy_shift = 1/eshift_inv
       return
       end
