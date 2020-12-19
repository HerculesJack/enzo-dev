#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!TODO  INCREASE THE TEMPERATURE OF THE CENTRAL GAS to 1d4
c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker11(nx, ny, nz,
     &                        d, dm, temp, u, v, w, cooltime,
     &                        dt, r, metal, zfield1, zfield2,
     &                        dx, t, z, procnum,
     &                        d1, x1, v1, t1,
     &                        nmax, xstart, ystart, zstart, ibuff,
     &                        imetal, imethod, mintdyn,
     &                        odthresh, masseff, smthresh, level, np,
     &                        xp, yp, zp, up, vp, wp,
     &                        mp, tdp, tcp, metalf,
     &                        imetalSNIa, metalSNIa, metalfSNIa, mp0,
     &                        verbose_mk)

c
c  CREATES GALAXY PARTICLES
c
c  written by: Brian O'Shea
c  date:       13 November 2002
c	This file was originally a copy of star_maker2.src,
c	which was originally written by Chris Loken.  As of
c	today, this is intended to be the unigrid version of
c	star_maker2, so the jeans mass and stochastic star
c	formation have been completely removed.  See
c	star_maker2.src for changes previous to 13 Nov. 2002.
c
c  modified1:  20 Dec 2002 by BWO
c       Stochastic star formation is added again.
c       The particle masses are averaged over several cells to avoid
c       the "runaway star particle" phenomenon
c  modified2: 27 Sept 2011 by JHW
c    added metalSNIa & metalfSNIa; included feedback from SN Ia/PN
c    (original changes by M. Ryan Joung)
c  modified2: 27 Sept 2013 by Taysun Kimm
c    added creation mass for more accurate feedback
c  modified3: 7 Nov 2020 by He Jia
c    removed Jeans criterion, added verbose mode
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold* avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
c    procnum - processor number (for output)
c    imetalSNIa - SN Ia metallicity flag (0 or 1 - none, 2 - yes)
c    verbose_mk - output some details
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    mp0      - creation mass of new particles in Msun (not density)
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    nmax     - particle array size specified by calling routine
c    metalfSNIa - metallicity fraction of particle (from SN Ia) ! MKRJ
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC procnum, imetalSNIa, verbose_mk
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz)
      R_PREC    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax), mp0(nmax)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    odthresh, masseff, smthresh, mintdyn
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      real*8 msolar
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass)
c
      ii = np
c
c  for each zone, : "star" particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c    is the gas mass greater than the Jeans mass?
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff

c
c              1) is this finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) is density greater than threshold?
c
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) is divergence negative?
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
               if (imethod .eq. 2) then
                  div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
               else
                  div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
               endif
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(d1))) *
     &                    dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

c
c  THIS IS COMMENTED OUT - NO JEANS MASS CRITERION IN THIS ALGORITHM!!!
c  BWO, 13 NOV 02 (fix 3 dec 02)
c              if (bmass .lt. jeanmass) goto 10
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c
               starfraction = min(masseff*dt/tdyn, 0.9_RKIND)
               tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)

c
c  STOCHASTIC STAR FORMATION HAS BEEN ADDED AGAIN - BWO 20 Dec 2002
c
#define STOCHASTIC_STAR_FORMATION
c
#ifdef STOCHASTIC_STAR_FORMATION
c
c                 Keep global count of "unfullfilled" star formation
c                 and when total is larger than threshold, then create
c                 a star particle with the threshold mass or 1/2 the
c                 gas in the cell, whichever is smaller.
c
               if (starfraction*bmass .lt. smthresh) then
                  sformsum = sformsum + starfraction*bmass
                  if (sformsum .lt. smthresh) goto 10
                  starfraction = smthresh/bmass
                  if (starfraction>0.5) goto 10
                  sformsum = sformsum - starfraction*bmass
               endif
#else
c
c              is star mass greater than threshold, then make it.
c              if it's less than threshold, go to the next cell.
c
               if (starfraction*bmass .lt. smthresh) goto 10
#endif

c
c              Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction * d(i,j,k)
               mp0(ii) = starfraction * d(i,j,k) *(d1*(dx*x1)**3/msolar)
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM

               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
!                 write(*,'("Setting metal fraction")')
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
!                 write(*,'("Zero metal fraction")')
                  metalf(ii) = 0._RKIND
               endif
c
c              MKRJ 2/20/08 Do the same for particle metal fraction from SN Ia
c
               if (imetalSNIa .gt. 0) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)    ! in here metal is a fraction
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1._RKIND - starfraction)*d(i,j,k)
c
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)
c
1000           format(i5,1x,6(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20

10          continue

            enddo
         enddo
      enddo
20    continue

      if (ii .ge. nmax) then
         write(6,*) 'star_maker11: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii

      if (verbose_mk .eq. 1 .and. np .ne. 0) then
         write(6,*) 'Stars created: number, time, level: ', np, t, level
      endif

      return
      end

c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback11(nx, ny, nz,
     &                   d, te, ge, u, v, w, metal,
     &                   idual, imetal, imethod, dt, r, dx, t, z,
     &                   d1, x1, v1, t1, t2, sn_param, m_eject, 
     &                   yield, imode,
     &                   npart, xstart, ystart, zstart, ibuff, level,
     &                   xp, yp, zp, up, vp, wp, mp, tdp, tcp, id0, 
     &                   metalf, imetalSNIa, metalSNIa, metalfSNIa,
     &                   type, justburn, mp0, FLEFTSN,
     &                   consv_p, consv_e, tmaxshock, 
     &                   verbose_fb, iPN, loadtype, vlimit)
c
c
c  RELEASES "STAR" PARTICLE MOMENTUM, ENERGY, MASS AND METALS
c
c  written by: Taysun Kimm
c  date:       8 Jan 2015
c              Supernova explosions are modeled as a mechanical input.
c              The amount of momentum input is determined by the stage 
c              in which the expansion is. In practice, this is determined
c              by the amount of mass swept up (See the Appendix of Kimm & Cen 2014)
c              On top of that, the main sequence life time of massive stars is modeled 
c              as realistically as possible (i.e. multiple SNe), based on the 
c              Padova track. This requires a speical field for initial mass of a 
c              star particle ('initial_mass' in my case after 'typeia_fraction')
c 
c  CONS:       1. This routine requires 4x4x4 grids,
c              and if the number of ActiveDimension in one direction is just 2,
c              momentum that is supposed to go for the missing grids will be omitted.
c              2. If SN is too close to the boundary, the center is shifted.
c              3. Notice also that Zeus solver leads to a cross-like expansion,
c              instead of a perfectly spherical one.
c
c  INPUTS:
c
c    d     - density field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metal fraction 
c    r     - refinement field (0 if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t,T to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imetalSNIa - metallicity flag for SNIa
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    level - current level of refinement
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles at t=Time
c    mp0      - creation mass of new particles in Msun (Msol should be 2d33)
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    id0      - particle index (used to generate the random number for sn2)
c    metalf   - star particle metal fraction
c    metalfSNIa - star particle SNIa metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    FLEFTSN - the fraction of *ejecta* mass left in the cell where SNe go off
c    consv_p - (NOT USED) correct the momentum so that it can be conserved 
c                  before/after the SN explosions (default:off)
c    consv_e - conserve the total energy (thermal+kinetic) before/after
c                  the SN explosions. This is to mimic shock-heating
c                  when two flows collide at the time of injection (def:on)
c    tmaxshock - maximum shock temperature
c    verbose_fb - output some details
c    iPN        - Planetary Nebulae is actually done in star_feedback_pn_ia.F
c                 This is only necessary to estimate the mass loss
c    loadtype   - parameter determining the loading from the host cell of SN
c    vlimit     - limit velocity
c    imode      - injection mode (0: immediate 6 faces, 2: 2nd layer,
c                 1: hybrid (eta<eta_crit: 0, eta>eta_crit: 1)
c                 StarFeedbackDistRadius is used for this
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod, level
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz), fl(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, t2, justburn
      P_PREC    xstart, ystart, zstart, t
      P_PREC    xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC imetalSNIa, type(npart), id0(npart), imode
      INTG_PREC consv_p, consv_e, verbose_fb, loadtype, iPN
      R_PREC    metalfSNIa(npart), mp0(npart)
      R_PREC    m_eject, sn_param, yield, FLEFTSN, tmaxshock, vlimit
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i,j,k,idx,n,displace,ic,jc,kc,nbub,cnt
      INTG_PREC nsn2(nx,ny,nz),nsnIa(nx,ny,nz),nsn_tot
      INTG_PREC nSN_sofar,nSN_new
      R_PREC dold(nx,ny,nz), Zsn(nx,ny,nz), ZsnIa(nx,ny,nz)
      R_PREC dpsn_x(nx,ny,nz), dpsn_y(nx,ny,nz), dpsn_z(nx,ny,nz)
      R_PREC eold(nx,ny,nz)
      R_PREC corr_x, corr_y, corr_z
      R_PREC psn(1:3),psn_solid_x,psn_solid_y,psn_solid_z,psn_solid
      R_PREC dpsn_solid_x,dpsn_solid_y,dpsn_solid_z
      R_PREC mform, tfactor, clight, energy, msolar_e51
      R_PREC minitial, xv1, xv2, age_yr,dt_yr
      R_PREC yr2s,one,sfreqSN2,mnow_msun,unit_msun,unit_nH
      R_PREC sn2_mass,snIa_mass,msn2_tot,msnIa_tot
      R_PREC msn2_cgs,snIa_msun,dratio,FLOADSN
      R_PREC vsn27(1:27,1:3), xdum(1:4), dq(npart,1:3)
      R_PREC wei27(1:27),nsight,rdum,maxwei
      R_PREC mejecta_tot, mejecta_solid, mswept_tot_cen
      R_PREC mswept_solid, chi_solid, chi_crit
      R_PREC m_h,Z_solar,X_h,nh_sn,Z_sn, dnew1, ebubb,Tbubb
      R_PREC t_rad,r_rad,r_shell,f_load_n,f_left_n
      R_PREC eint_old,eint_new,ekin_old,ekin_new,eint_host
      R_PREC tshock, tshock_old, fsmoo, fsolid
      R_PREC vsn_solid, mcen_solid, ESN
      R_PREC fZ_SNII,Z_Ia,dswept,dsolid0
      R_PREC uold,vold,wold,dall,zold,uej,vej,wej,dsolid
      R_PREC unew,vnew,wnew,dcen,vrad
      INTG_PREC indsn27(1:27,1:3),indsn,indv6(1:6,1:3)
      INTG_PREC nSNIa_new,nSNIa_sofar
      INTG_PREC iup,jup,kup,iad,jad,kad,ifa,jfa,kfa
      R_PREC uup,vup,wup,uad,vad,wad,met27(1:27)
      R_PREC d0,d11,d2,d3,d4,d5,d6,dtot,dpsn,dsn,d0l
      R_PREC u0,u11,u2,u3,u4,u5,u6
      R_PREC v0,v11,v2,v3,v4,v5,v6
      R_PREC w0,w11,w2,w3,w4,w5,w6
      R_PREC uavg,vavg,wavg,sigma_v,c_sound,Mach
      R_PREC vol,dd1,dd2,zd1,zd2,ektot,et
      R_PREC uv1(1:3),uv2(1:3),uv3(1:3)
      R_PREC eth0,eth1,ekk0,ekk1,eunit
      R_PREC pi,G,msolar,two,kms,pdum
      R_PREC tnow,rmf,MdotSNII,MdotPN
      R_PREC u0l,u0r,v0l,v0r,w0l,w0r
      R_PREC correction, chi_shift
      logical oksn,immediate_neighbor
      parameter (clight = c_light, msolar_e51 = 1800._RKIND)
      parameter (yr2s=31536000._RKIND, one=1.000000_RKIND)
      parameter (msolar=SolarMass)
      parameter (m_h=mass_h)
      parameter (Z_solar=0.02)
      parameter (X_h=0.76)
      parameter (two=2._RKIND)
      parameter (pi=pi_val, G=GravConst)
      if(npart.eq.0) return

      ! sanity check
      if(imethod.ne.2.and.idual.ne.1)then
         print *,'>>>ERR: Turn on DualEnergyFormalism '
         print *,'>>>to use Feedback=32768 for non-Zeus solvers'
         stop
      endif

c-----------------------------------------------------------------------
c    0. intialize and set some useful numbers
c-----------------------------------------------------------------------
      !if(imethod.ne.0)
      displace=1    ! safety device not to touch the boundary
      if(imethod.eq.2) displace=2

      Z_Ia    = 1d0  ! assume that 100% of the snIa ejecta is metal 
      ESN     = 1d51 
      FLOADSN = 1d0-FLEFTSN ! the fraction of ejecta that goes into the neighbors
      eth0=0d0
      eth1=0d0
      ekk0=0d0
      ekk1=0d0
      kms =v1/1d5
 
      do i=1,nx
      do j=1,ny
      do k=1,nz
        nsn2(i,j,k)  = 0  ! total number of sn2 in a cell
        nsnIa(i,j,k) = 0  ! total number of snIa in a cell
        dold   (i,j,k) = d(i,j,k)
        if(idual.eq.1)then
           eold(i,j,k) = ge(i,j,k)
        else
           eold(i,j,k) = te(i,j,k)
        endif
        dpsn_x (i,j,k) = 0._RKIND
        dpsn_y (i,j,k) = 0._RKIND
        dpsn_z (i,j,k) = 0._RKIND
        fl     (i,j,k) = FLOADSN         ! gas entrainment from the SN cell
        if(imetal.eq.1) Zsn(i,j,k)=0._RKIND       ! to re-distribute matter
        if(imetalSNIa.eq.2) ZsnIa(i,j,k)=0._RKIND ! to re-distribute matter
      end do
      end do
      end do
 
      do n=1,npart
        dq(n,1)=0._RKIND ; dq(n,2)=0._RKIND; dq(n,3)=0._RKIND
      end do

c     input mode for Zeus solver
      if (imode.eq.0) then
         immediate_neighbor = .true. ! correction term will be invoked
      else if (imode.eq.2)then
         immediate_neighbor = .false.
      end if 
c
c     - sfreqSN2  :  specific frequency of SN per Msun
c     - unit_msun:  code unit -> Msun
c     - sn2_mass :  individual SN2 mass (actually density, code unit)
c     - msn2_cgs : typical mass of SNe in [g] 

      sfreqSN2  = msolar_e51*sn_param 
      unit_msun = d1*(x1*dx)**3/msolar
      sn2_mass  = one/sfreqSN2*m_eject/unit_msun
      msn2_cgs  = sn2_mass*unit_msun*msolar
      snIa_msun = 1.4 
      snIa_mass = snIa_msun/unit_msun
      dt_yr     = dt*t1/yr2s
      unit_nH   = d1/m_h*X_h

      eunit=unit_msun*1.989d33*v1*v1
c-----------------------------------------------------------------------
c  1. calculate justburn
c-----------------------------------------------------------------------
      do n=1, npart
        if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c
c          Determine how much of a given star particle would have been 
c          turned into stars during this timestep.  Then calculate the mass
c          which should have formed during this timestel dt using the integral
c          form of the Cen & Ostriker formula.
           xv1 = (t      - tcp(n))/tdp(n)
           if (xv1 .gt. 12.0) goto 10     ! t-tcp >> tdp so ignore
           xv2 = (t + dt - tcp(n))/tdp(n)

c          First calculate the initial mass of the star particle 
c          in question.
           minitial = mp0(n)/unit_msun 
c
c          Then, calculate the amount of mass that would have formed in
c          this timestep.
c
           mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &                          (1._RKIND + xv2)*exp(-xv2))
           mform = max(mform, 0._RKIND)
c
c          Compute index of the cell that the star particle
c          resides in.
c 
           i = int((xp(n) - xstart)/dx,IKIND) + 1
           j = int((yp(n) - ystart)/dx,IKIND) + 1
           k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c          check bounds - if star particle is outside of this grid
c          then exit and give a warning.
c
           if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &        .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'star particle out of grid; i,j,k,level,',
     &                     i,j,k,level
           else
               justburn = justburn + mform * dt * dx**3
           endif
10         continue
        end if
      end do 


c-----------------------------------------------------------------------
c    2. feedback part
c-----------------------------------------------------------------------
c 
c     first gather the total number of SNII (not energy)
c     this is to do the SN explosion on a cell-by-cell basis
c
      do n=1, npart
         if (tcp(n) .le.0) goto 11
         if (mp(n)  .le.0) goto 11
         if (type(n).ne.2) goto 11


         minitial  = mp0(n)           ! initial mass in Msun
         mnow_msun = mp(n)*unit_msun  ! current mass
         age_yr   = (t+dt-tcp(n))*t1/yr2s

         ! Type Ia SNe
         if(imetalSNIa.eq.2)then
            call get_number_snIa(minitial,age_yr,dt_yr,
     &                           nSNIa_new,nSNIa_sofar)
            mnow_msun = mnow_msun + snIa_msun*nSNIa_sofar
            ! the right term is to get an accurate number of nSN2
         else
            nSNIa_new=0
         endif

         ! Take into account the mass loss due to Planetary Nebulare
         if(iPN.eq.1)then
            tnow = (t-tcp(n))*t1/yr2s
            call cal_rmf(tnow,rmf,MdotSNII,MdotPN) !Mdot means dM not dM/dt
            mnow_msun = mnow_msun + MdotPN*minitial
         endif

         ! fraction of mass formed until t+dt
         xv2 = (t + dt - tcp(n)) / tdp(n) 
         nsn_tot=floor((one - (one + xv2)*exp(-xv2))*minitial*sfreqSN2)
         nSN_sofar = nint((minitial - mnow_msun)/(sn2_mass*unit_msun))

         call get_number_sn2(age_yr,id0(n),metalf(n),nsn_tot,
     &                       nSN_sofar,nSN_new)
         if(nSN_new+nSNIa_new.le.0) goto 11

         mp(n)  = mp(n) - nSN_new*sn2_mass - nSNIa_new*snIa_mass 
         
         i = int((xp(n) - xstart)/dx,IKIND) + 1
         j = int((yp(n) - ystart)/dx,IKIND) + 1
         k = int((zp(n) - zstart)/dx,IKIND) + 1

c        if a particle is too close to the boundary,
c        shift the centre; this may be justifiable if you agree   
c        that the location of a particle is only accurate
c        at two cell size; zeus requires that SN is two cells away from the boundary

         !left boundaries
         if(i.lt.1+ibuff+displace)  i=1+ibuff+displace
         if(j.lt.1+ibuff+displace)  j=1+ibuff+displace
         if(k.lt.1+ibuff+displace)  k=1+ibuff+displace
         !right boundaries
         if(i.gt.nx-ibuff-displace) i=nx-ibuff-displace
         if(j.gt.ny-ibuff-displace) j=ny-ibuff-displace
         if(k.gt.nz-ibuff-displace) k=nz-ibuff-displace

         nsn2(i,j,k) =nsn2(i,j,k) +nSN_new
         nsnIa(i,j,k)=nsnIa(i,j,k)+nSNIa_new

c        mass-weighted metals of the original SN2 prog stars
         if (imetal.eq.1) then
            Zsn(i,j,k)=Zsn(i,j,k)
     &                +nSN_new*sn2_mass*metalf(n)
         endif

c        gather original momentum of the SN2 particles
         dpsn_x(i,j,k)=dpsn_x(i,j,k)+nSN_new*sn2_mass*up(n)
         dpsn_y(i,j,k)=dpsn_y(i,j,k)+nSN_new*sn2_mass*vp(n)
         dpsn_z(i,j,k)=dpsn_z(i,j,k)+nSN_new*sn2_mass*wp(n)

         if (imetalSNIa.eq.2)then
c        gather original momentum of the SNIa particles
            dpsn_x(i,j,k)=dpsn_x(i,j,k)+nSNIa_new*snIa_mass*up(n)
            dpsn_y(i,j,k)=dpsn_y(i,j,k)+nSNIa_new*snIa_mass*vp(n)
            dpsn_z(i,j,k)=dpsn_z(i,j,k)+nSNIa_new*snIa_mass*wp(n)
         endif

11       continue  ! end this star particle 
      end do



c     convert (sort of) metal density to metallicity of the SN progenitors
      do k=1+ibuff,nz-ibuff
      do j=1+ibuff,ny-ibuff 
      do i=1+ibuff,nx-ibuff
         if(nsn2(i,j,k).gt.0)then
            if(imetal.eq.1)then
               dall=nsn2(i,j,k)*sn2_mass
               Zsn(i,j,k)=Zsn(i,j,k)/dall ! metallicity
            endif
         endif
      end do
      end do
      end do


c     no supernova? return
      if(sum(nsnIa+nsn2).eq.0) goto 300

      if(verbose_fb.eq.1) print *,'NSN II, Ia =', sum(nsn2), sum(nsnIa)
      !if(z.lt.5.5) print *, 'HJ DEBUG 0', t, xstart, ystart, zstart

c     check zone size
      if((nx-ibuff*2.lt.1+displace*2).or.
     &   (ny-ibuff*2.lt.1+displace*2).or.
     &   (nz-ibuff*2.lt.1+displace*2)) then
         write(6,*) 'WARNING: you have nx,ny,nz=',
     &   int(nx-ibuff*2, kind=2),
     &   int(ny-ibuff*2, kind=2),
     &   int(nz-ibuff*2, kind=2)
         write(6,*) 'WARNING: but we need each side to be at least ',
     &   int(1+displace*2, kind=2)
         write(6,*) 'WARNING: you may consider using a larger RefineBy'
      endif

c     set useful arrays (it would've been better not to do this
c                        every time this routine is calld, but
c                        TK left this for readibility)

      nbub = 27
      idx  = 1
      do k=-1,1
      do j=-1,1
      do i=-1,1
         indsn27(idx,1)=i
         indsn27(idx,2)=j
         indsn27(idx,3)=k
         rdum = sqrt(float(i*i+j*j+k*k))
         if(i.ne.0.or.j.ne.0.or.k.ne.0)then
            vsn27(idx,1)=i/rdum
            vsn27(idx,2)=j/rdum
            vsn27(idx,3)=k/rdum
         endif
         cnt=0
         if (i==0) cnt=cnt+1
         if (j==0) cnt=cnt+1
         if (k==0) cnt=cnt+1
         if(imethod.eq.2)then
             select case (cnt)
               case (0)
                  wei27(idx)=0d0
               case (1)
                  wei27(idx)=0d0
               case (2)
                  wei27(idx)=1d0
               case (3) 
                  wei27(idx)=0d0
            end select
         else
            select case (cnt)
               case (0)
                  wei27(idx)=0.3d0
               case (1)
                  wei27(idx)=3.4d0
               case (2)
                  wei27(idx)=8d0
               case (3) 
                  wei27(idx)=0d0
            end select
         endif
         idx=idx+1
      end do
      end do
      end do
      nsight = sum(wei27)
      maxwei = maxval(wei27)

c-------------------------------------------------------
c   2.1 update the SN cells first 
c-------------------------------------------------------
c   add some fraction of SN ejecta and subtract mass to be entrained.
c   these cells will be updated later if they are near to other SNe

      !if(z.lt.5.5) print *, 'HJ DEBUG 1', t, xstart, ystart, zstart
      do k=1+ibuff,nz-ibuff
      do j=1+ibuff,ny-ibuff
      do i=1+ibuff,nx-ibuff

         nsn_tot=nsn2(i,j,k)+nsnIa(i,j,k)

        if(nsn_tot.gt.0)then

          !print *,nx,ny,nz,ibuff,i,j,k
          if(imethod.eq.2)then
             u0 = (u(i  ,j  ,k  )+u(i+1,j  ,k  ))/2
             v0 = (v(i  ,j  ,k  )+v(i  ,j+1,k  ))/2
             w0 = (w(i  ,j  ,k  )+w(i  ,j  ,k+1))/2
          else
             u0 = u(i  ,j  ,k  )
             v0 = v(i  ,j  ,k  )
             w0 = w(i  ,j  ,k  )
          endif
          msn2_tot   = sn2_mass*nsn2(i,j,k) ! [code density]
          msnIa_tot  = snIa_mass*nsnIa(i,j,k)
          mejecta_tot= msn2_tot+msnIa_tot 

          nh_sn    = dold(i,j,k)*unit_nH
          dratio   = dold(i,j,k)/mejecta_tot
          ! estimate the mass loading from the host cell of SN

          eint_old   = eold(i,j,k)
          tshock_old = eint_old*t2
          c_sound    = 11729.170*sqrt(tshock_old) !cm/s 

          d0 = dold(i  ,j  ,k  )
          d11= dold(i-1,j  ,k  ); d2 = dold(i+1,j  ,k  )
          d3 = dold(i  ,j-1,k  ); d4 = dold(i  ,j+1,k  )
          d5 = dold(i  ,j  ,k-1); d6 = dold(i  ,j  ,k+1)

          ! determine the fraction of mass to be entrained from the host cell of SN
          if (loadtype.eq.1) then ! experimental
   
             if(imethod.eq.2)then !Zeus solver
               u11= (u(i-1,j  ,k  )+u(i  ,j  ,k  ))/2
               u2 = (u(i+1,j  ,k  )+u(i+2,j  ,k  ))/2
               u3 = (u(i  ,j-1,k  )+u(i+1,j-1,k  ))/2
               u4 = (u(i  ,j+1,k  )+u(i+1,j+1,k  ))/2
               u5 = (u(i  ,j  ,k-1)+u(i+1,j  ,k-1))/2
               u6 = (u(i  ,j  ,k+1)+u(i+1,j  ,k+1))/2
   
               v11= (v(i-1,j  ,k  )+v(i-1,j+1,k  ))/2
               v2 = (v(i+1,j  ,k  )+v(i+1,j+1,k  ))/2
               v3 = (v(i  ,j-1,k  )+v(i  ,j  ,k  ))/2
               v4 = (v(i  ,j+1,k  )+v(i  ,j+2,k  ))/2
               v5 = (v(i  ,j  ,k-1)+v(i  ,j+1,k-1))/2
               v6 = (v(i  ,j  ,k+1)+v(i  ,j+1,k+1))/2
   
               w11= (w(i-1,j  ,k  )+w(i-1,j  ,k+1))/2
               w2 = (w(i+1,j  ,k  )+w(i+1,j  ,k+1))/2
               w3 = (w(i  ,j-1,k  )+w(i  ,j-1,k+1))/2
               w4 = (w(i  ,j+1,k  )+w(i  ,j+1,k+1))/2
               w5 = (w(i  ,j  ,k-1)+w(i  ,j  ,k  ))/2
               w6 = (w(i  ,j  ,k+1)+w(i  ,j  ,k+2))/2
             else
               u11= u(i-1,j  ,k  ); u2 = u(i+1,j  ,k  )
               u3 = u(i  ,j-1,k  ); u4 = u(i  ,j+1,k  )
               u5 = u(i  ,j  ,k-1); u6 = u(i  ,j  ,k+1)
   
               v11= v(i-1,j  ,k  ); v2 = v(i+1,j  ,k  )
               v3 = v(i  ,j-1,k  ); v4 = v(i  ,j+1,k  )
               v5 = v(i  ,j  ,k-1); v6 = v(i  ,j  ,k+1)
   
               w11= w(i-1,j  ,k  ); w2 = w(i+1,j  ,k  )
               w3 = w(i  ,j-1,k  ); w4 = w(i  ,j+1,k  )
               w5 = w(i  ,j  ,k-1); w6 = w(i  ,j  ,k+1)
             endif

             ! use the actual local properties to measure turbulence
             dtot = d0+d11+d2+d3+d4+d5+d6
             uavg = (d0*u0+d11*u11+d2*u2+d3*u3+d4*u4+d5*u5+d6*u6)/dtot
             vavg = (d0*v0+d11*v11+d2*v2+d3*v3+d4*v4+d5*v5+d6*v6)/dtot
             wavg = (d0*w0+d11*w11+d2*w2+d3*w3+d4*w4+d5*w5+d6*w6)/dtot
             
             sigma_v = 0.0
             sigma_v = d0*((u0-uavg)**2+(v0-vavg)**2+(w0-wavg)**2)
     &             + d11*((u11-uavg)**2+(v11-vavg)**2+(w11-wavg)**2)
     &             + d2*((u2-uavg)**2+(v2-vavg)**2+(w2-wavg)**2)
     &             + d3*((u3-uavg)**2+(v3-vavg)**2+(w3-wavg)**2)
     &             + d4*((u4-uavg)**2+(v4-vavg)**2+(w4-wavg)**2)
     &             + d5*((u5-uavg)**2+(v5-vavg)**2+(w5-wavg)**2)
     &             + d6*((u6-uavg)**2+(v6-vavg)**2+(w6-wavg)**2)
             sigma_v = sqrt(sigma_v/dtot)*v1 ! [cm/s]

             Mach    = sigma_v/c_sound
             ! Mixture of forcing modes (solenoidal+compressive): b~0.4 (Federrath+10)
             f_load_n = min(exp(-0.5*log(1+(0.4*Mach)**2)),FLOADSN)
             if(verbose_fb.eq.1) print *,'#### fload', sngl(f_load_n),
     &            sngl(sigma_v/1e5),sngl(c_sound/1e5)
          else 
             ! use a universal value (StarFeedbackFLeftover)
             f_load_n = FLOADSN
          endif

          f_left_n  = one - f_load_n
          fl(i,j,k) = f_left_n

          dsn      = mejecta_tot*FLEFTSN
          dnew1    = d(i,j,k)*f_left_n + dsn

          if(imethod.eq.2)then !Zeus solver
             !u11 = (p*+u11(d11+d0))/(d11+d0+d*)
             d0l  = d0*f_left_n

             dpsn = dpsn_x(i,j,k)*FLEFTSN
             u(i  ,j,k) = (dpsn+u(i  ,j,k)*(d11+d0l))/(d11+d0l+dsn) 
             u(i+1,j,k) = (dpsn+u(i+1,j,k)*(d2 +d0l))/(d2 +d0l+dsn) 
             dpsn = dpsn_y(i,j,k)*FLEFTSN
             v(i,j  ,k) = (dpsn+v(i,j  ,k)*(d3+d0l))/(d3+d0l+dsn) 
             v(i,j+1,k) = (dpsn+v(i,j+1,k)*(d4+d0l))/(d4+d0l+dsn) 
             dpsn = dpsn_z(i,j,k)*FLEFTSN
             w(i,j,k  ) = (dpsn+w(i,j,k  )*(d5+d0l))/(d5+d0l+dsn) 
             w(i,j,k+1) = (dpsn+w(i,j,k+1)*(d6+d0l))/(d6+d0l+dsn)

             ! bookkeep the momentum of gas entrained from the host cell
             dpsn_x(i,j,k)=dpsn_x(i,j,k)*FLOADSN+u0*d0*f_load_n
             dpsn_y(i,j,k)=dpsn_y(i,j,k)*FLOADSN+v0*d0*f_load_n
             dpsn_z(i,j,k)=dpsn_z(i,j,k)*FLOADSN+w0*d0*f_load_n
 
             u0 = (u(i,j,k)+u(i+1,j,k))/two
             v0 = (v(i,j,k)+v(i,j+1,k))/two
             w0 = (w(i,j,k)+w(i,j,k+1))/two

          else

             ! update velcity of the SN host cell
             u(i,j,k) = u0*d0*f_left_n+dpsn_x(i,j,k)*FLEFTSN
             v(i,j,k) = v0*d0*f_left_n+dpsn_y(i,j,k)*FLEFTSN 
             w(i,j,k) = w0*d0*f_left_n+dpsn_z(i,j,k)*FLEFTSN
             u(i,j,k) = u(i,j,k)/dnew1
             v(i,j,k) = v(i,j,k)/dnew1
             w(i,j,k) = w(i,j,k)/dnew1

             ! bookkeep the momentum of gas entrained from the host cell
             dpsn_x(i,j,k)=dpsn_x(i,j,k)*FLOADSN+u0*d0*f_load_n
             dpsn_y(i,j,k)=dpsn_y(i,j,k)*FLOADSN+v0*d0*f_load_n
             dpsn_z(i,j,k)=dpsn_z(i,j,k)*FLOADSN+w0*d0*f_load_n

             u0 = u(i,j,k)
             v0 = v(i,j,k)
             w0 = w(i,j,k)

          endif

c         temperature of the central cell
          eint_new = dold(i,j,k)*eint_old/(dold(i,j,k)+mejecta_tot)
          ekin_new = (u0*u0+v0*v0+w0*w0)/two

          if(idual.eq.1) then
             ge(i,j,k)=eint_new
             te(i,j,k)=eint_new + ekin_new
          else
             te(i,j,k)=eint_new
          endif
   

c         distribute metal
          !-----------------------------------------------------------------
          !here metal is a fraction
          ! m_Znew = yield*m_star 
          ! yield  = Z_new*m_ejecta
          ! original fluid metal = metal(i,j,k)  ! metal fraction
          ! original particle metal = Zsn(i,j,k)*msn2_tot
          ! SN2 processed   = msn2_tot*(one-Zsn(i,j,k))*(yield/m_eject)
          !-----------------------------------------------------------------
          if(imetal.eq.1)then
            fZ_SNII = Zsn(i,j,k)+
     &                (one-Zsn(i,j,k))*(yield/m_eject)

          ! Note:Zsn will be used as the metallicity of the loaded gas from now on
            dall = d(i,j,k)*f_load_n+mejecta_tot*FLOADSN
            Zsn(i,j,k) = d(i,j,k)*metal(i,j,k)*f_load_n
     &                  + fZ_SNII*msn2_tot*FLOADSN
     &                  + Z_Ia*msnIa_tot*FLOADSN 
            Zsn(i,j,k) = Zsn(i,j,k)/dall

          ! metallicity
            dall = d(i,j,k)*f_left_n+mejecta_tot*FLEFTSN
            zold = metal(i,j,k)
            metal(i,j,k) = d(i,j,k)*metal(i,j,k)*f_left_n
     &                   + fZ_SNII*msn2_tot*FLEFTSN
     &                   + Z_Ia*msnIa_tot*FLEFTSN
            metal(i,j,k) = metal(i,j,k)/dall

          endif

          if(imetalSNIa.eq.2)then
          ! Note:ZsnIa will be used as the Ia metallicity of the loaded gas in the below 
            dall = d(i,j,k)*f_load_n+mejecta_tot*FLOADSN
            ZsnIa(i,j,k)=d(i,j,k)*metalSNIa(i,j,k)*f_load_n
     &                       + Z_Ia*msnIa_tot*FLOADSN 
            ZsnIa(i,j,k)=ZsnIa(i,j,k)/dall

            dall = d(i,j,k)*f_left_n+mejecta_tot*FLOADSN
            metalSNIa(i,j,k)=d(i,j,k)*metalSNIa(i,j,k)*f_left_n
     &                      +Z_Ia*msnIa_tot*FLEFTSN
            metalSNIa(i,j,k)=metalSNIa(i,j,k)/dall
          endif

          d(i,j,k) = dnew1

        endif

      end do
      end do
      end do


c-----------------------------------------
c   2.2 inject mass, momentum, and energy
c-----------------------------------------
      !if(z.lt.5.5) print *, 'HJ DEBUG 2', t, xstart, ystart, zstart
      do k=1+ibuff,nz-ibuff
      do j=1+ibuff,ny-ibuff 
      do i=1+ibuff,nx-ibuff
   
        nsn_tot = nsn2(i,j,k)+nsnIa(i,j,k) 
        if(nsn_tot.eq.0) goto 230

c       total ejecta mass from SNe in code unit
        msn2_tot = sn2_mass*nsn2(i,j,k)
        msnIa_tot = snIa_mass*nsnIa(i,j,k)
        mejecta_tot = msn2_tot+msnIa_tot

        f_left_n  = fl(i,j,k)
        f_load_n  = 1._RKIND - f_left_n

        eint_host = eold(i,j,k)*dold(i,j,k)*fl(i,j,k)

c       total mass entrained from the SNcell + SN ejecta.
c       this excludes the fraction left in the SN cell. 
        mswept_tot_cen  = dold(i,j,k)*f_load_n 
     &                  + mejecta_tot*FLOADSN

c       momentum conservation if you wish
        corr_x=0._RKIND
        corr_y=0._RKIND
        corr_z=0._RKIND

c       keep the original metallicity before the injection for accuracy
        do indsn=1,nbub
c          mass swept up along a child cell (8._RKIND)
           ic = i+indsn27(indsn,1)
           jc = j+indsn27(indsn,2)
           kc = k+indsn27(indsn,3)
       
c          safety device if ndim is too small
           if((ic.lt.1+ibuff).or.(ic.gt.nx-ibuff)) goto 119
           if((jc.lt.1+ibuff).or.(jc.gt.ny-ibuff)) goto 119
           if((kc.lt.1+ibuff).or.(kc.gt.nz-ibuff)) goto 119

           met27(indsn)=metal(ic,jc,kc)
119        continue
        end do


c       keep the original velocity of 6 faces for Zeus
        u0l = u(i,j,k) ; u0r = u(i+1,j,k)
        v0l = v(i,j,k) ; v0r = v(i,j+1,k)
        w0l = w(i,j,k) ; w0r = w(i,j,k+1)
         
        ebubb = 0.0

c       estimate the momentum along some solid angle
        do indsn=1,nbub

c         mass swept up along a child cell (8._RKIND)
          ic = i+indsn27(indsn,1)
          jc = j+indsn27(indsn,2)
          kc = k+indsn27(indsn,3)

          if(i==ic.and.j==jc.and.k==kc) goto 229 ! skip the host cell

c         safety device if ndim is too small
          if((ic.lt.1+ibuff).or.(ic.gt.nx-ibuff)) goto 229
          if((jc.lt.1+ibuff).or.(jc.gt.ny-ibuff)) goto 229
          if((kc.lt.1+ibuff).or.(kc.gt.nz-ibuff)) goto 229

          fsolid = wei27(indsn)/nsight

          if(fsolid.lt.1d-10) goto 229

c         ejecta mass along some solid angle
          mejecta_solid = mejecta_tot*FLOADSN*fsolid

c         total mass swept up along some solid angle
          mswept_solid = dold(ic,jc,kc)*wei27(indsn)/maxwei
     &                 + mswept_tot_cen*fsolid

c         mass loading factor
          chi_solid = mswept_solid / mejecta_solid

c         transition ratio between swept-up and ejecta mass
          dd1  = dold(i,j,k)*fsolid*fl(i,j,k) ! volume-filling density in the host cell of SN 
          dd2  = dold(ic,jc,kc)*wei27(indsn)/maxwei! density in the immediate neighbor
          vol  = wei27(indsn)/maxwei+fsolid*fl(i,j,k)
          nh_sn = (dd1+dd2)/vol*unit_nH

          zd1 = dd1*metal(i,j,k) 
          zd2 = dd2*met27(indsn)
          Z_sn = (zd1+zd2)/(dd1+dd2)
          Z_sn = max(Z_sn/Z_solar,0.01_RKIND)

          chi_crit = (900._RKIND/0.676/(msn2_cgs/msolar))
     &             * nsn_tot**(-2._RKIND/17._RKIND)
     &             * nh_sn**(-4._RKIND/17._RKIND)
     &             * Z_sn**(-0.28_RKIND)

c         compute to-be-injected momentum
          if (chi_solid.ge.chi_crit)then
c           snowplow phase
            fsmoo=0.676_RKIND
            vsn_solid = sqrt(fsmoo*2*ESN*chi_crit/msn2_cgs)
     &                  /chi_solid/FLOADSN ! cm/s
            if(imode.eq.1) immediate_neighbor=.false.
          else
c           ~adabatic phase
            fsmoo=1._RKIND-(1._RKIND-0.676_RKIND)
     &             *(chi_solid-1._RKIND)/(chi_crit-1._RKIND) ! for smooth transition
            vsn_solid = sqrt(fsmoo*2*ESN/chi_solid/msn2_cgs)/FLOADSN
            if(imode.eq.1) immediate_neighbor=.true.
          endif


          if(imode.eq.0)then !experimental
             chi_shift  = -one/5.*log(one/0.05-one)*chi_crit
             correction = -5.*(chi_solid-chi_crit+chi_shift)/chi_crit
             correction = one/(one+exp(correction))
             vsn_solid  = vsn_solid*(one+correction)
          endif

c         expansion speed
          vsn_solid = vsn_solid/v1  ! [code unit]

          ! momentum due to explosion
          psn_solid_x=chi_solid*mejecta_solid*vsn_solid*vsn27(indsn,1)
          psn_solid_y=chi_solid*mejecta_solid*vsn_solid*vsn27(indsn,2)
          psn_solid_z=chi_solid*mejecta_solid*vsn_solid*vsn27(indsn,3)
          psn_solid  = psn_solid_x**2+psn_solid_y**2+psn_solid_z**2
          psn_solid  = sqrt(psn_solid)

          ekin_old = psn_solid*vsn_solid/two
          if(imode.eq.0)then !experimental
             ekin_old=ekin_old/(one+correction)**2
          endif

          ! initial momentum of a star particle
          dpsn_solid_x = dpsn_x(i,j,k)*fsolid
          dpsn_solid_y = dpsn_y(i,j,k)*fsolid
          dpsn_solid_z = dpsn_z(i,j,k)*fsolid

          ! initial energy of the gas entrained + original motion of a star particle
          dswept = mswept_tot_cen*fsolid
          ekin_old = (dpsn_solid_x**2+dpsn_solid_y**2+dpsn_solid_z**2)
     &               /two/dswept + ekin_old
          
!          if(consv_p.eq.1)then ! experimental
!             psn_solid=chi_solid*mejecta_solid*vsn_solid ! re-defined!!!
!             corr_x=corr_x+psn_solid*vsn27(indsn,1)
!             corr_y=corr_y+psn_solid*vsn27(indsn,2)
!             corr_z=corr_z+psn_solid*vsn27(indsn,3)
!          endif

          dsolid = d(ic,jc,kc) + dswept  ! this will be the updated density of d(ic,jc,kc)

          ! update TypeII metallicity
          if(imetal.eq.1)then
            zold = metal(ic,jc,kc)
            metal(ic,jc,kc) = (d(ic,jc,kc)*metal(ic,jc,kc)
     &             +Zsn(i,j,k)*dswept)/dsolid
          endif

          ! update TypeIa metallicity
          if(imetalSNIa.eq.2)then 
            metalSNIa(ic,jc,kc) = (d(ic,jc,kc)*metalSNIa(ic,jc,kc)
     &         + ZsnIa(i,j,k)*dswept)/dsolid
          endif

c         update mass (density), momentum, and energy !!!!
          if(imethod.eq.2)then

             ! ic,jc,kc   : index for density update (immediate neighbor)
             ! iup,jup,kup: index for velocity update

             iup=ic;jup=jc;kup=kc
             if(immediate_neighbor)then
                if(indsn27(indsn,1)>0)then
                   iup=ic;iad=ic;ifa=ic+1
                endif
                if(indsn27(indsn,1)<0)then
                   iup=ic+1;iad=ic+1;ifa=ic-1
                endif
                if(indsn27(indsn,2)>0)then
                   jup=jc;jad=jc;jfa=jc+1
                endif
                if(indsn27(indsn,2)<0)then
                   jup=jc+1;jad=jc+1;jfa=jc-1
                endif
                if(indsn27(indsn,3)>0)then
                   kup=kc;kad=kc;kfa=kc+1
                endif
                if(indsn27(indsn,3)<0)then
                   kup=kc+1;kad=kc+1;kfa=kc-1
                endif
             else
                if(indsn27(indsn,1)>0)then
                   iup=ic+1;iad=ic;ifa=ic+1
                endif
                if(indsn27(indsn,1)<0)then
                   iup=ic;iad=ic+1;ifa=ic-1
                endif
                if(indsn27(indsn,2)>0)then
                   jup=jc+1;jad=jc;jfa=jc+1
                endif
                if(indsn27(indsn,2)<0)then
                   jup=jc;jad=jc+1;jfa=jc-1
                endif
                if(indsn27(indsn,3)>0)then
                   kup=kc+1;kad=kc;kfa=kc+1
                endif
                if(indsn27(indsn,3)<0)then
                   kup=kc;kad=kc+1;kfa=kc-1
                endif
             endif

             u0=0;v0=0;w0=0

             ! kinetic energy before injection - ignore the central cell term
             if(immediate_neighbor)then
                if(indsn27(indsn,1)>0)then
                   u2  = u(iup  ,jc,kc)
                   u3  = u(iup+1,jc,kc)
                   d0 =  d(iup,jc,kc)
                   u0 = (u2+u3)/two
                   ektot = d0*u0*u0/two
                else if(indsn27(indsn,1)<0)then
                   u11 = u(iup-1,jc,kc)
                   u2  = u(iup  ,jc,kc)
                   d0 =  d(iup-1,jc,kc)
                   u0 = (u11+u2)/two
                   ektot = d0*u0*u0/two
                else if (indsn27(indsn,2)>0)then
                   v2  = v(ic,jup  ,kc)
                   v3  = v(ic,jup+1,kc)
                   d0 =  d(ic,jup,kc)
                   v0 = (v2+v3)/two
                   ektot = d0*v0*v0/two
                else if (indsn27(indsn,2)<0)then
                   v11 = v(ic,jup-1,kc)
                   v2  = v(ic,jup  ,kc)
                   d0 =  d(ic,jup-1,kc)
                   v0 = (v11+v2)/two
                   ektot = d0*v0*v0/two
                else if (indsn27(indsn,3)>0)then
                   w2  = w(ic,jc,kup)
                   w3  = w(ic,jc,kup+1)
                   d0 =  d(ic,jc,kup)
                   w0 = (w2+w3)/two
                   ektot = d0*w0*w0/two
                else if (indsn27(indsn,3)<0)then
                   w11 = w(ic,jc,kup-1)
                   w2  = w(ic,jc,kup)
                   d0 =  d(ic,jc,kup-1)
                   w0 = (w11+w2)/two
                   ektot = d0*w0*w0/two
                end if
             else
                ! kinetic energy before injection
                if(indsn27(indsn,1).ne.0)then
                   d0 =  d(iup,jc,kc)
                   u0 = (u(iup,jc,kc)+u(iup+1,jc,kc))/two
                   ektot = d0*u0*u0/two
                   d0 =  d(iup-1,jc,kc)
                   u0 = (u(iup-1,jc,kc)+u(iup,jc,kc))/two
                   ektot = ektot + d0*u0*u0/two
                else if (indsn27(indsn,2).ne.0)then
                   d0 =  d(ic,jup,kc)
                   v0 = (v(ic,jup,kc)+v(ic,jup+1,kc))/two
                   ektot = d0*v0*v0/two
                   d0 =  d(ic,jup-1,kc)
                   v0 = (v(ic,jup-1,kc)+v(ic,jup,kc))/two
                   ektot = ektot + d0*v0*v0/two
                else if (indsn27(indsn,3).ne.0)then
                   d0 =  d(ic,jc,kup)
                   w0 = (w(ic,jc,kup)+w(ic,jc,kup+1))/two
                   ektot = d0*w0*w0/two
                   d0 =  d(ic,jc,kup-1)
                   w0 = (w(ic,jc,kup-1)+w(ic,jc,kup))/two
                   ektot = ektot + d0*w0*w0/two
                end if
             endif

             ekin_old = ekin_old+ektot
             eint_old = te(ic,jc,kc)*d(ic,jc,kc) ! ZEUS only works with idual=0


             ! update velocities - need to be careful when injecting momentum
             ! because velocities are face-centered!
             ! if we change one value, it will always affect two cells
             ! unew = uold + (2psnx-(uup+uold)*dej)/(d1+d2+dej) 

             unew=0;vnew=0;wnew=0
 
             ! update mass
             if(.not.immediate_neighbor) d(ic,jc,kc) = dsolid

             !x-direction
             if(indsn27(indsn,1).ne.0)then
                dall = d(ic,jc,kc)+d(ifa,jc,kc) !includes dswept
                uad  = u(iad,jc,kc) ! velocity of the adjacent face to the host cell of SN
                uup  = u(iup,jc,kc) ! velocity of the face we are updating
                 
                ! Here 2 is a fudge factor that is necessary because momentum is cancelled for the gas in the host cell
                if(immediate_neighbor)then
                   unew = (u(iup,jc,kc)*d(ic,jc,kc)+psn_solid_x
     &                  + dpsn_solid_x)/dsolid
                else
                   unew = uup + (2*(psn_solid_x+dpsn_solid_x)
     &                           -(uup+uad)*dswept)/dall
 
                endif
                if(unew*kms> vlimit) unew=vlimit/kms  ! velocity limiter..
                u(iup,jc,kc)=unew
             endif

             !y-direction
             if(indsn27(indsn,2).ne.0)then
                dall = d(ic,jc,kc)+d(ic,jfa,kc) !includes dswept
                vad  = v(ic,jad,kc) ! velocity of the adjacent face to the host cell of SN
                vup  = v(ic,jup,kc) ! velocity of the face we are updating
                if(immediate_neighbor)then
                   vnew = (v(ic,jup,kc)*d(ic,jc,kc)+psn_solid_y
     &                  + dpsn_solid_y)/dsolid
                else
                   vnew = vup+(2*(psn_solid_y+dpsn_solid_y)
     &                         -(vup+vad)*dswept)/dall
                endif
                if(vnew*kms> vlimit) vnew=vlimit/kms  ! velocity limiter..
                v(ic,jup,kc)=vnew
             endif

             !z-direction
             if(indsn27(indsn,3).ne.0)then
                dall = d(ic,jc,kc)+d(ic,jc,kfa) !includes dswept
                wad  = w(ic,jc,kad) ! velocity of the adjacent face to the host cell of SN
                wup  = w(ic,jc,kup) ! velocity of the face we are updating
                if(immediate_neighbor)then
                   wnew = (w(ic,jc,kup)*d(ic,jc,kc)+psn_solid_z
     &                  + dpsn_solid_z)/dsolid
                else
                   wnew = wup+(2*(psn_solid_z+dpsn_solid_z)
     &                          -(wup+wad)*dswept)/dall
                endif
                if(wnew*kms> vlimit) wnew=vlimit/kms  ! velocity limiter..
                w(ic,jc,kup)=wnew
             endif


             ! update mass
             if(immediate_neighbor) d(ic,jc,kc) = dsolid


             ! To update thermal energy
             ! kinetic energy after the injection
             if(immediate_neighbor)then
                if(indsn27(indsn,1)>0)then
                   u2  = u(iup  ,jc,kc)
                   u3  = u(iup+1,jc,kc)
                   d0 =  d(iup,jc,kc)
                   u0 = (u2+u3)/two
                   ektot = d0*u0*u0/two
                else if(indsn27(indsn,1)<0)then
                   u11 = u(iup-1,jc,kc)
                   u2  = u(iup  ,jc,kc)
                   d0 =  d(iup-1,jc,kc)
                   u0 = (u11+u2)/two
                   ektot = d0*u0*u0/two
                else if (indsn27(indsn,2)>0)then
                   v2  = v(ic,jup  ,kc)
                   v3  = v(ic,jup+1,kc)
                   d0 =  d(ic,jup,kc)
                   v0 = (v2+v3)/two
                   ektot = d0*v0*v0/two
                else if (indsn27(indsn,2)<0)then
                   v11 = v(ic,jup-1,kc)
                   v2  = v(ic,jup  ,kc)
                   d0 =  d(ic,jup-1,kc)
                   v0 = (v11+v2)/two
                   ektot = d0*v0*v0/two
                else if (indsn27(indsn,3)>0)then
                   w2  = w(ic,jc,kup)
                   w3  = w(ic,jc,kup+1)
                   d0 =  d(ic,jc,kup)
                   w0 = (w2+w3)/two
                   ektot = d0*w0*w0/two
                else if (indsn27(indsn,3)<0)then
                   w11 = w(ic,jc,kup-1)
                   w2  = w(ic,jc,kup)
                   d0 =  d(ic,jc,kup-1)
                   w0 = (w11+w2)/two
                   ektot = d0*w0*w0/two
                end if
             else
                ! kinetic energy before injection
                if(indsn27(indsn,1).ne.0)then
                   d0 =  d(iup,jc,kc)
                   u0 = (u(iup,jc,kc)+u(iup+1,jc,kc))/two
                   ektot = d0*u0*u0/two
                   d0 =  d(iup-1,jc,kc)
                   u0 = (u(iup-1,jc,kc)+u(iup,jc,kc))/two
                   ektot = ektot + d0*u0*u0/two
                else if (indsn27(indsn,2).ne.0)then
                   d0 =  d(ic,jup,kc)
                   v0 = (v(ic,jup,kc)+v(ic,jup+1,kc))/two
                   ektot = d0*v0*v0/two
                   d0 =  d(ic,jup-1,kc)
                   v0 = (v(ic,jup-1,kc)+v(ic,jup,kc))/two
                   ektot = ektot + d0*v0*v0/two
                else if (indsn27(indsn,3).ne.0)then
                   d0 =  d(ic,jc,kup)
                   w0 = (w(ic,jc,kup)+w(ic,jc,kup+1))/two
                   ektot = d0*w0*w0/two
                   d0 =  d(ic,jc,kup-1)
                   w0 = (w(ic,jc,kup-1)+w(ic,jc,kup))/two
                   ektot = ektot + d0*w0*w0/two
                end if
             endif

             ekin_new = ektot
             eint_new = max(ekin_old-ekin_new,tiny)
  
             ! consider shock progation through the central cell
             ! i.e. shock-heated gas is likely to exist in the central cell as well if eint_new is larger
             ebubb = ebubb + eint_new*wei27(indsn)

             eint_new = eint_new+eint_old+eint_host*fsolid
             tshock=eint_new/d(ic,jc,kc)*t2
             !print *,'TKSN',log10(tshock),eint_new,v(ic,jup,kc),jup,ny
             if(tshock>tmaxshock)then
                if(verbose_fb.eq.1.and.tshock>4e8)then
                   print *, ' TKWARN: told, tnew [1e8]=', 
     &             sngl(tshock/1e8), sngl(tmaxshock/1e8)
                endif
               eint_new = tmaxshock/t2*d(ic,jc,kc)
             endif
             if(consv_e.ne.1) eint_new=eint_old
             te(ic,jc,kc) = eint_new/d(ic,jc,kc)


          else ! PPM ############################## 

             d0   = d(ic,jc,kc)
             uold = u(ic,jc,kc)
             vold = v(ic,jc,kc)
             wold = w(ic,jc,kc)
             u(ic,jc,kc) = (uold*d0+psn_solid_x+dpsn_solid_x)/dsolid
             v(ic,jc,kc) = (vold*d0+psn_solid_y+dpsn_solid_y)/dsolid
             w(ic,jc,kc) = (wold*d0+psn_solid_z+dpsn_solid_z)/dsolid

             ! velocity limit (experimental)
             if(u(ic,jc,kc)*kms>vlimit) u(ic,jc,kc)=vlimit/kms
             if(v(ic,jc,kc)*kms>vlimit) v(ic,jc,kc)=vlimit/kms
             if(w(ic,jc,kc)*kms>vlimit) w(ic,jc,kc)=vlimit/kms
 
             ! deposit mass
             d(ic,jc,kc) = dsolid

             ! for energy conservation
             if(idual.eq.1)then
                eint_old = ge(ic,jc,kc)
             else
                eint_old = te(ic,jc,kc)
             endif
             eint_old = eint_old*d0

             ekin_old = ekin_old                               ! from the entrained gas from the host cell + original motion of a star
     &                + psn_solid*vsn_solid/two                ! from SNe
     &                + (uold**2 + vold**2 + wold**2)/two*d0   ! from the existing cell

             ekin_new =  0.5_RKIND*dsolid*(u(ic,jc,kc)**2
     &                   +v(ic,jc,kc)**2 + w(ic,jc,kc)**2)
             eint_new = max(ekin_old-ekin_new,tiny)

             ! consider shock progation through the central cell
             ! i.e. shock-heated gas is likely to exist in the central cell as well if eint_new is larger
             ebubb = ebubb + eint_new*wei27(indsn)

             eint_new = eint_new + eint_old + eint_host*fsolid 
             eint_new = eint_new/dsolid
             if(consv_e.lt.1) eint_new = eint_old/d0
             !ekin_new = ekin_new/dsolid

c            limit the temperature of a shock if you wish (not physical though)
             tshock   = eint_new*t2
             ekin_new = ekin_new/dsolid

             if(verbose_fb.eq.1)then
                if(tshock>4e8)then
                   print *, ' TKWARN: told, tnew [1e8]=', 
     &             sngl(tshock_old/1e8), sngl(tshock/1e8)
                endif
             endif

             if(tshock>tmaxshock)then
               eint_new = tmaxshock/t2
             endif

             te(ic,jc,kc) = eint_new + ekin_new
             if(idual.eq.1) ge(ic,jc,kc) = eint_new
          endif

          if(verbose_fb.eq.1)then
             uej = psn_solid_x/dsolid*v1/1e5
             vej = psn_solid_y/dsolid*v1/1e5
             wej = psn_solid_z/dsolid*v1/1e5
             if(u(ic,jc,kc)*v1>3e8) 
     &          write(*,'("NOTE: SN Vx=",7(F7.2,1x))') 
     &               uold*v1/1e5,u(ic,jc,kc)*v1/1e5,chi_solid,uej
     &               ,vsn_solid*v1/1e5,d0/dsolid
     &               , psn_solid_x/(chi_solid*mejecta_solid)*v1/1e5
             if(v(ic,jc,kc)*v1>3e8)
     &          write(*,'("NOTE: SN Vy=",7(F7.2,1x))') 
     &               vold*v1/1e5,v(ic,jc,kc)*v1/1e5,chi_solid,vej
     &               ,vsn_solid*v1/1e5,d0/dsolid
     &               , psn_solid_y/(chi_solid*mejecta_solid)*v1/1e5
             if(w(ic,jc,kc)*v1>3e8)
     &          write(*,'("NOTE: SN Vz=",7(F7.2,1x))') 
     &               wold*v1/1e5,w(ic,jc,kc)*v1/1e5,chi_solid,wej
     &               ,vsn_solid*v1/1e5,d0/dsolid
     &               , psn_solid_z/(chi_solid*mejecta_solid)*v1/1e5
          endif


229       continue

        end do  ! indsn=1,nbub 


        ! increase the temperature of the central cell
        Tbubb = ebubb/nsight/d(i,j,k) ! simple average with the geometry
        Tbubb = Tbubb/10
        if(idual.eq.1)then
           te(i,j,k)=min(te(i,j,k)+Tbubb,min(3d8,tmaxshock)/t2) 
           ge(i,j,k)=min(ge(i,j,k)+Tbubb,min(3d8,tmaxshock)/t2)
        else
           te(i,j,k)=min(te(i,j,k)+Tbubb,min(3d8,tmaxshock)/t2)
        endif 


        ! experimental - conserve linear momentum
        if(consv_p.eq.1.and.imethod.ne.2)then 
          do indsn=1,nbub
            ic = i+indsn27(indsn,1)
            jc = j+indsn27(indsn,2)
            kc = k+indsn27(indsn,3)

            dnew1 = d(ic,jc,kc)
            u(ic,jc,kc)=u(ic,jc,kc)-corr_x/dnew1*fsolid
            v(ic,jc,kc)=v(ic,jc,kc)-corr_y/dnew1*fsolid
            w(ic,jc,kc)=w(ic,jc,kc)-corr_z/dnew1*fsolid
          enddo
        endif 

230     continue
      end do
      end do
      end do

      do k=1+ibuff,nz-ibuff
      do j=1+ibuff,ny-ibuff 
      do i=1+ibuff,nx-ibuff
         if(imethod.eq.2)then
         eth1=eth1+te(i,j,k)*d(i,j,k)
         ekk1=ekk1+d(i,j,k)*(u(i,j,k)+u(i-1,j,k))**2/4.
         ekk1=ekk1+d(i,j,k)*(v(i,j,k)+v(i,j-1,k))**2/4.
         ekk1=ekk1+d(i,j,k)*(w(i,j,k)+w(i,j,k-1))**2/4.
         else
         eth1=eth1+ge(i,j,k)*d(i,j,k)
         ekk1=ekk1+d(i,j,k)*u(i,j,k)**2/2.
         ekk1=ekk1+d(i,j,k)*v(i,j,k)**2/2.
         ekk1=ekk1+d(i,j,k)*w(i,j,k)**2/2.
         endif
      end do
      end do
      end do 

300   continue

      !if(z.lt.5.5) print *, 'HJ DEBUG 3', t, xstart, ystart, zstart
      return
      end
c=======================================================================
c////////////////////////  SUBROUTINE  \\\\\\\\\\\\\\\\\\\\\\\
c       Generate the delay times of SN2 for a star particle
c=======================================================================
      subroutine get_number_sn2 (age_yr,idp,zp_star,nsn_tot,
     &                           nsn_sofar,nsn)
      implicit none
#include "fortran_types.def"
      R_PREC age_yr,zp_star
      INTG_PREC idp,nsn,nsn_tot,nsn_sofar,nsn_ok,i
      INTG_PREC localseed
!      SAVE localseed
!      DATA localseed /-340005/
      R_PREC age_star,xdum,ydum,logzpsun
!     These numbers are based on the SN frequency from Starburst99
!     (for the IMF EXPONENT of 1.3 between 8-100 Msun)
      R_PREC co0,co1,co2,coa(1:3),cob(1:3),coc(1:3)
      DATA coa /-2.677292E-01,1.392208E-01,-5.747332E-01/
      DATA cob /4.208666E-02, 2.152643E-02, 7.893866E-02/
      DATA coc /-8.612668E-02,-1.698731E-01,1.867337E-01/
      R_PREC ran1
      EXTERNAL ran1

      nsn=0
      logzpsun=max(min(zp_star,0.05),0.008) ! no extrapolation
      logzpsun=log10(logzpsun/0.02)
      co0 = coa(1)+coa(2)*logzpsun+coa(3)*logzpsun**2
      co1 = cob(1)+cob(2)*logzpsun+cob(3)*logzpsun**2
      co2 = coc(1)+coc(2)*logzpsun+coc(3)*logzpsun**2
     
      if(age_yr/1e6.le.(-co2/co1))return

      nsn_ok    = 0
      localseed = -mod(idp,2147483647)
      do i=1,nsn_tot
!         localseed = -340005 !!!!!!!!!!!!!!!!!!!!!!!!!!!! ONLY FOR DEBUG
         !xdum =  ran0(localseed)
         xdum =  ran1(localseed)
         ! inverse function for y=co0+sqrt(co1*x+co2)
         ydum = ((xdum-co0)**2.-co2)/co1
         if(ydum.le.age_yr/1e6) then
            nsn_ok=nsn_ok+1
         endif
      enddo

      nsn = max(nsn_ok - nsn_sofar,0)
      end

c=======================================================================
c////////////////////////  Function  \\\\\\\\\\\\\\\\\\\\\\\
c       Minimal Standard generatore: Numerical recipe
c=======================================================================
      function ran0(idum)
      implicit none
#include "fortran_types.def"
      INTG_PREC idum, IA, IM, IQ, IR, MASK,k
      R_PREC    ran0, AM
      parameter(IA=16807,IM=2147483647,AM=1./IM)
      parameter(IQ=127773,IR=2836,MASK=123459876)

      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      end
c=======================================================================
      function ran1(idum)
      implicit none
#include "fortran_types.def"
      INTG_PREC idum,IA,IM,IQ,IR,NTAB,NDIV
      R_PREC ran1,AM,EPS,RNMX
      parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
      parameter(NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTG_PREC j,k,iv(NTAB),iy
      save iv,iy
      data iv /NTAB*0/, iy /0/ 
      if (idum.le.0.or.iy.eq.0) then! initialize
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
         end do
         iy=iv(1)
      end if
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      end
c=======================================================================
c////////////////////////  Function  \\\\\\\\\\\\\\\\\\\\\\\
c       Calculate the number of snIa between [t-dt,t]
c=======================================================================
      subroutine get_number_snIa(minitial,age_yr,dt_yr,
     &                           nsnIa,nsnIa_sofar)
      implicit none
#include "fortran_types.def"
      INTG_PREC nsnIa, nsnIa_sofar
      R_PREC coe(1:5),age_yr,dt_yr,logt,minitial,fsn1,fsn2
      DATA   coe /3.586799E-01,-1.551665E-01,2.480943E-02,
     &           -1.741565E-03,4.559893E-05/

c     minitial in Msun
      logt = log10(max(age_yr-dt_yr,1e-6))
      fsn1 = minitial*(coe(1)+coe(2)*logt+coe(3)*logt**2
     &         + coe(4)*logt**3+coe(5)*logt**4)
      if(logt<7.5)fsn1=0.
      logt = log10(age_yr)
      fsn2 = minitial*(coe(1)+coe(2)*logt+coe(3)*logt**2
     &         + coe(4)*logt**3+coe(5)*logt**4)
      if(logt<7.5)fsn2=0.

      nsnIa_sofar = floor(fsn1)
      nsnIa = floor(fsn2) - nsnIa_sofar

      end

