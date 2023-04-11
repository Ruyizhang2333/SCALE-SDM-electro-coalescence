!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coalescence subroutines for the SDM
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-12 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2016-07-20 (S.Shima) [add] Coalescence of the "liqice" attribute added 
!! @li      2019-09-29 (S.Shima) [add] electro coalescence kernels
!! @li      2019-09-30 (S.Shima) [mod] electro_coalescence_efficiency_coulomb
!! @li      2019-09-30 (S.Shima) [mod] to output the collision efficiency
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_coalescence
  use scale_precision

  implicit none
  private
  public :: sdm_coales

  real(RP) :: alpha=0.0D0 ! coefficient to decide the charge amount (Andronache 2004)

contains
  subroutine sdm_coales(sdm_colkrnl,sdm_colbrwn,               &
                        sdm_aslset,sdm_aslrho,                 &
                        sdm_dtcol,            &
                        pres_scale, t_scale,DENS,              &
                        zph_crs,                &
                        ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl, &
                        sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,&
                        sort_id,sort_key,sort_freq,sort_tag,   &
                        sd_rng,sd_rand,                        &
                        sort_tag0,fsort_id,icp,sd_perm,c_rate  )
    use scale_process, only:  &
         PRC_IsMaster,          &
         PRC_MPIstop
    use gadg_algorithm, only: &
         gadg_count_sort
    use rng_uniform_mt, only: &
         c_rng_uniform_mt, gen_rand_array
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_const, only:  &
         t0 => CONST_TEM00, &    ! Gas Constant of vapor [J/K/kg]
         rw => CONST_Rvap,  &    ! Gas Constant of vapor [J/K/kg]
         rd => CONST_Rdry,  &    ! Gas Constant of dry air [J/K/kg]
         cp => CONST_CPdry, &
         p0 => CONST_PRE00       ! Reference Pressure [Pa]
    use m_sdm_common, only: &
         VALID2INVALID,INVALID,knum_sdm, &
         rho_amsul,rho_nacl,ONE_PI,m2micro,r0col,ratcol,ecoll,micro2m,dxiv_sdm,dyiv_sdm,F_THRD,O_THRD,rrst,boltz,mass_air,i2
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort, sdm_getperm
    use m_sdm_motion, only: &
         sdm_getvz_liq
    !  Input variables
    integer, intent(in) :: sdm_colkrnl   ! Kernel type for coalescence process
    integer, intent(in) :: sdm_colbrwn   ! Control flag of Brownian Coagulation and Scavenging process
    integer, intent(in) :: sdm_aslset    ! Control flag to set species and way of chemical material as water-soluble aerosol
    real(RP),intent(in) :: sdm_aslrho(20)! User specified density of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sdm_dtcol   ! tims step of {stochastic coalescence} process
    real(RP), intent(in) :: pres_scale(KA,IA,JA)  ! Pressure
    real(RP), intent(in) :: t_scale(KA,IA,JA)    ! Temperature
    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: zph_crs(KA,IA,JA) ! z physical coordinate
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer, intent(in) :: sd_numasl ! Number of kind of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    ! Input and output variables
    type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
    real(RP),intent(inout) :: sd_rand(1:sd_num) ! random numbers
    integer, intent(inout) :: sort_id(1:sd_num) ! super-droplets sorted by SD-grids
    integer, intent(inout) :: sort_key(1:sd_num) ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
    integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2) ! accumulated number of super-droplets in each SD-grid
    integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Output variables
    integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
    integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    integer, intent(out) :: icp(1:sd_num) ! index of coalescence pair
    integer, intent(out) :: sd_perm(1:sd_num) ! random permutations
    real(RP), intent(out) :: c_rate(1:sd_num) ! coalescence probability
    ! Internal shared variables
    real(RP) :: sd_aslrho(1:22) ! Density of chemical material contained as water-soluble aerosol in super droplets
    integer(RP) :: sd_ncol ! how many times coalescence occurs
    integer :: freq_max ! get the maximum number of super-droplets in each grid
    integer :: hfreq_max ! hfreq_max / 2
    integer :: ipremium  ! premium coef. for coalescence
    ! Work variables
    real(RP) :: dmask(1:22) ! mask for vactorization
    real(RP) :: sd_asl1(1:sd_numasl)
    real(RP) :: sd_asl2(1:sd_numasl) ! aerosol mass of super-droplets with large/small multiplicity
    real(RP) :: sd_cc1  ! slip correction of super-droplets
    real(RP) :: sd_cc2  ! with large/small multiplicity
    real(RP) :: sd_dia1 ! diameter of super-droplets
    real(RP) :: sd_dia2 ! with large/small multiplicity
    real(RP) :: sd_lmd1 ! mean free path of super-droplets
    real(RP) :: sd_lmd2 ! with large/small multiplicity
    real(RP) :: sd_m1   ! mass of super-droplets
    real(RP) :: sd_m2   ! with large/small multiplicity
    real(RP) :: sd_r1   ! radius of super-droplets
    real(RP) :: sd_r2   ! with large/small multiplicity
    real(RP) :: sd_rk1  ! index[k/real] of super-droplets with large multiplicity
    real(RP) :: sd_rw1  ! radius of water parts in super-droplets
    real(RP) :: sd_rw2  ! with large/small multiplicity
    real(RP) :: sd_tmasl1
    real(RP) :: sd_tmasl2 ! total mass of aerosol part in super droplets with large/small multiplicity
    real(RP) :: sd_tvasl1 ! total volume of aerosol part in super
    real(RP) :: sd_tvasl2 ! droplets with large/small multiplicity
    real(RP) :: sd_v1   ! volume of super-droplets
    real(RP) :: sd_v2   ! with large/small multiplicity
    real(RP) :: sd_vz1  ! terminal veolocity of super-droplets
    real(RP) :: sd_vz2  ! 
    real(RP) :: sd_c1   ! temporary
    real(RP) :: sd_c2
    real(RP) :: sd_d1
    real(RP) :: sd_d2
    real(RP) :: sd_g1
    real(RP) :: sd_g2
    real(RP) :: lmd_crs ! air mean free path
    real(RP) :: p_crs   ! pressure
    real(RP) :: pt_crs  ! potential temperarure
    real(RP) :: t_crs   ! temperarure
    real(RP) :: vis_crs ! dynamic viscosity
    real(RP) :: sumdia  ! sum of variables of a pair of droplets
    real(RP) :: sumd
    real(RP) :: sumc
    real(RP) :: sumg
    real(RP) :: sumr
    real(RP) :: k12     ! brownian coagulation coefficient
    real(RP) :: dvz     ! difference in terminal velocity of a pair of super-droplets
    real(RP) :: dtmp    ! temporary
    real(RP) :: frac    ! fraction parts
    real(RP) :: ivvol(KA,IA,JA)   ! inverse of a grid volume
    real(RP) :: tdeg    ! temperature in degree
    real(RP) :: rq      ! radius ratio of a pair of super-droplets
    real(RP) :: ek      ! temporary
    real(RP) :: p       ! temporary
    real(RP) :: q       ! temporary
    real(RP) :: crate_elc   ! electro coalescence kernel
    real(RP) :: eff_elc     ! electro coalesence efficiency [] 

    integer(DP) :: sd_nmax  ! maximum multiplicity
    integer(DP) :: sd_n1    ! multiplicity of super-droplets with large multiplicity
    integer(DP) :: sd_n2    ! multiplicity of super-droplets  with small multiplicity
    integer(i2) :: sd_li1, sd_li2

    integer, allocatable :: fsort_tag(:) ! buffer for sorting
    integer, allocatable :: fsort_freq(:) ! buffer for sorting

    integer :: idx_nasl(1:22)  ! index for vactorization


    integer :: gnum          ! grid number

    integer :: in1, in2, in3, in4 ! index
    integer :: irr, iqq           ! index of coef.
    integer :: i, j, k, m, n, s   ! index
    integer :: t, tc, tp          ! index
    integer :: ix, jy

    integer :: sort_tag0m
    integer :: sort_freqm
    integer :: icptc, icptp

    logical,parameter :: output_kernel = .false. ! Flag to output the kernel value for debug  
    character(len=20) :: fname    ! output filename
    integer           :: stat    ! Runtime status
    integer           :: FILE_SD
    integer,allocatable :: sd_itmp1(:),sd_itmp2(:),sd_itmp3(:)
    
    !--------------------------------------------------------------------
   
    ! open(unit=10,file='kernel.txt')

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_coales",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    ! ni_sdm=IE-IS+1, nj_sdm=JE-JS+1, knum_sdm=floor(rkumax)+1)-KS+1
    gnum = ni_sdm * nj_sdm * knum_sdm

    freq_max = 1

    !### aerosol type ###!
    ! why we need sd_aslrho? Isn't it same to sdm_aslrho?? Check later
    do n=1,22
       sd_aslrho(n) = 1.0_RP
    end do

    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then

          !### numasl=1 @ init : NaCl ###!

          sd_aslrho(1) = real(rho_nacl,kind=RP)

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslrho(1) = real(rho_amsul,kind=RP)
          sd_aslrho(2) = real(rho_nacl,kind=RP)

       end if
    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)
       sd_aslrho(2) = real(rho_nacl,kind=RP)

       !         do n=1,20
       !            call getrname( id_sdm_aslrho + (n-1), sd_aslrho(n+2) )
       !         end do
       do n=1,20
          sd_aslrho(n+2) = sdm_aslrho(n)
       end do

    end if

    ! Sorting super-droplets.
    
    call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
         sort_id,sort_key,sort_freq,sort_tag,'valid')

    ! Initialize
    
    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do

    do n=1,gnum+2
       sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
    end do

    do n=1,sd_num
       c_rate(n) = 0.0_RP
    end do

    ! Get the maximum number of super-droplets in each grid
    
    do m=1,gnum
       freq_max = max( freq_max, sort_freq(m) )
    end do

    hfreq_max = int(freq_max/2)

    ! Sorting the grids by the number of super-droplets in each grid

    allocate( fsort_tag(0:freq_max+1) )
    allocate( fsort_freq(0:freq_max)  )

    call gadg_count_sort( sort_freq(1:gnum), 0, freq_max, &
         fsort_freq, fsort_tag, fsort_id )

    fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

    ! Get random number using random number generator
    
    call gen_rand_array( sd_rng, sd_rand )

    ! Get random permutation layout of super-droplets in each grid

    call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
         &                 sort_tag0,fsort_tag,fsort_id,sd_rand,sd_perm)

    ! Get random number using random number generator

    call gen_rand_array( sd_rng, sd_rand )

    ! Select a pair of super-droples in random permutation layout

!!$      do n=1,hfreq_max
!!$
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle
       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          in1 = 2 * n
          in2 = in1 - 1

          !### select pair of super-droplets in each grid ###!

!!$            in3 = sd_perm( sort_tag0(m) + in1 )
!!$            in4 = sd_perm( sort_tag0(m) + in2 )
          in3 = sd_perm( sort_tag0m + in1 )
          in4 = sd_perm( sort_tag0m + in2 )

          !### set the random index ###!
!!$            tc = sort_tag0(m) + n
          tc = sort_tag0m + n
!!$            tp = tc + int(sort_freq(m)/2)
          tp = tc + int(sort_freqm/2)

!!$            icp(tc) = sort_id( sort_tag0(m) + in3 )
!!$            icp(tp) = sort_id( sort_tag0(m) + in4 )
          icp(tc) = sort_id( sort_tag0m + in3 )
          icp(tp) = sort_id( sort_tag0m + in4 )

       end do

    end do

! To output kenel value for debug 
      if( output_kernel ) then

         ! Get random number using random number generator

         call gen_rand_array( sd_rng, sd_rand )

!OCL NORECURRENCE
         do m=1,gnum
            if( sort_freq(m) <= 1 ) cycle

            sort_tag0m = sort_tag0(m)
            sort_freqm = sort_freq(m)

            do n=1,sort_freqm/2

               tc = sort_tag0m + n
               tp = tc + sort_freqm/2

               sd_r(icp(tc)) = 1.0d-6 ! fix the size as you like [m]
               sd_r(icp(tp)) = exp( log(1.0d-8)+(log(1.0d-3)-log(1.0d-8))*sd_rand(tc) )! uniform sample in log(r) [m]
             
            end do

         end do

         allocate(sd_itmp1(1:sd_num))
         allocate(sd_itmp2(1:sd_num))
         allocate(sd_itmp3(1:sd_num))

         !! evaluate the terminal velocity
         call sdm_getvz_liq(pres_scale,DENS,t_scale,            &
              sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
              sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )

         deallocate(sd_itmp1)
         deallocate(sd_itmp2)
         deallocate(sd_itmp3)
         
         ! Get random number using random number generator (in case it will be used)

         call gen_rand_array( sd_rng, sd_rand )

      end if

    ! Get effective collision probability for "Gravitational Settling" 
    
    if( sdm_colkrnl==0 ) then

       !### Golovin's kernel (m^3/s?) ###!
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = sd_r(icptc)
             sd_r2 = sd_r(icptp)

             c_rate(tc) = 1500.0_RP * 1.333333_RP * ONE_PI              &
                  * ( sd_r1*sd_r1*sd_r1 + sd_r2*sd_r2*sd_r2 )

          end do

       end do

    else if( sdm_colkrnl==1 ) then

       !### Long's kernel (m^2)  ###!

!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = max( sd_r(icptc), sd_r(icptp) )  !! large
             sd_r2 = min( sd_r(icptc), sd_r(icptp) )  !! small

             if( sd_r1 <= 5.E-5_RP ) then
                c_rate(tc) = 4.5E+8_RP * ( sd_r1*sd_r1 )                 &
                     * ( 1.0_RP - 3.E-6_RP/(max(3.01E-6_RP,sd_r1)) )
             else
                c_rate(tc) = 1.d0
             end if

             sumr = sd_r1 + sd_r2

             c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

          end do

       end do

    else if( sdm_colkrnl==2 ) then

       !### Hall's kernel (m^2) ###!

!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = max( sd_r(icptc), sd_r(icptp) )  !! large
             sd_r2 = min( sd_r(icptc), sd_r(icptp) )  !! small

             rq    = sd_r2 / sd_r1
             sd_r1 = sd_r1 * m2micro    !! [m] => [micro-m]
             sd_r2 = sd_r2 * m2micro    !! [m] => [micro-m]

             !! Get index of the array {r0,rat}.

             if( sd_r1 <= r0col(1) ) then
                irr = 1
             else if( sd_r1 <= r0col(2) ) then
                irr = 2
             else if( sd_r1 <= r0col(3) ) then
                irr = 3
             else if( sd_r1 <= r0col(4) ) then
                irr = 4
             else if( sd_r1 <= r0col(5) ) then
                irr = 5
             else if( sd_r1 <= r0col(6) ) then
                irr = 6
             else if( sd_r1 <= r0col(7) ) then
                irr = 7
             else if( sd_r1 <= r0col(8) ) then
                irr = 8
             else if( sd_r1 <= r0col(9) ) then
                irr = 9
             else if( sd_r1 <= r0col(10) ) then
                irr = 10
             else if( sd_r1 <= r0col(11) ) then
                irr = 11
             else if( sd_r1 <= r0col(12) ) then
                irr = 12
             else if( sd_r1 <= r0col(13) ) then
                irr = 13
             else if( sd_r1 <= r0col(14) ) then
                irr = 14
             else if( sd_r1 <= r0col(15) ) then
                irr = 15
             else
                irr = 16
             end if

             if( rq <= ratcol(2) ) then
                iqq = 2
             else if( rq <= ratcol(3) ) then
                iqq = 3
             else if( rq <= ratcol(4) ) then
                iqq = 4
             else if( rq <= ratcol(5) ) then
                iqq = 5
             else if( rq <= ratcol(6) ) then
                iqq = 6
             else if( rq <= ratcol(7) ) then
                iqq = 7
             else if( rq <= ratcol(8) ) then
                iqq = 8
             else if( rq <= ratcol(9) ) then
                iqq = 9
             else if( rq <= ratcol(10) ) then
                iqq = 10
             else if( rq <= ratcol(11) ) then
                iqq = 11
             else if( rq <= ratcol(12) ) then
                iqq = 12
             else if( rq <= ratcol(13) ) then
                iqq = 13
             else if( rq <= ratcol(14) ) then
                iqq = 14
             else if( rq <= ratcol(15) ) then
                iqq = 15
             else if( rq <= ratcol(16) ) then
                iqq = 16
             else if( rq <= ratcol(17) ) then
                iqq = 17
             else if( rq <= ratcol(18) ) then
                iqq = 18
             else if( rq <= ratcol(19) ) then
                iqq = 19
             else if( rq <= ratcol(20) ) then
                iqq = 20
             else
                iqq = 21
             end if
             !! Get c_rate

             if( irr>=16 ) then

                q  = (rq-ratcol(iqq-1)) / (ratcol(iqq)-ratcol(iqq-1))
                ek = (1.0_RP-q)*ecoll(15,iqq-1) + q*ecoll(15,iqq)

                c_rate(tc) = min( ek, 1.0_RP )

             else if( irr>=2 .and. irr<16 ) then

                p = (sd_r1-r0col(irr-1))/(r0col(irr)-r0col(irr-1))
                q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                c_rate(tc) = (1.0_RP-p)*(1.0_RP-q)*ecoll(irr-1,iqq-1)     &
                     + p*(1.0_RP-q)*ecoll(irr,iqq-1)              &
                     + q*(1.0_RP-p)*ecoll(irr-1,iqq)              &
                     + p*q*ecoll(irr,iqq)

             else

                q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                c_rate(tc) = (1.0_RP-q)*ecoll(1,iqq-1) + q*ecoll(1,iqq)

             end if

             sd_r1 = sd_r1 * micro2m    !! [micro-m] => [m]
             sd_r2 = sd_r2 * micro2m    !! [micro-m] => [m]

             !! Get c_rate

             sumr = sd_r1 + sd_r2

             c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

          end do

       end do

    else if( sdm_colkrnl==3 ) then

       !### no coalescence effeciency hydrodynamic kernel (m^2) ###!
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sumr = sd_r(icptc) + sd_r(icptp)

             c_rate(tc) = ONE_PI * (sumr*sumr)

          end do

       end do

    end if

    ! -----

    if( sdm_colkrnl==0 ) then

       !### Golovin's kernel [-] ###!
!!$         do n=1,hfreq_max
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2
             
             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             !! Get location of Super-Droplets
             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do
       end do

    else

       !### Long's kernel, Hall's kernel,       ###!
       !### no col_effi hydrodynamic kernel [-] ###!
       
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2
       
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             dvz = abs( sd_vz(icptc) - sd_vz(icptp) )

             !! Get location of Super-Droplets

             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP)        &
                  * ivvol(k,i,j) * dvz

          end do

       end do

    end if

! -----

! Calculate the electro coalescence probability

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2
             
          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          sd_r1 = sd_r(icp(tc))
          sd_r2 = sd_r(icp(tp))
          sd_vz1= sd_vz(icp(tc))
          sd_vz2= sd_vz(icp(tp))

          !! Calculate the Coulomb force kernel
            
          !### location of super-droplets ###

          ! index in center grid
          i = floor(sd_ri(icptc))+1
          j = floor(sd_rj(icptc))+1
          k = floor(sd_rk(icptc))+1

          p_crs = pres_scale(k,i,j) !! [Pa]
          t_crs = t_scale(k,i,j)    !! [K]

          !### dynamic viscosity [Pa*s]  ###!
          !### (Pruppacher & Klett,1997) ###!

          tdeg = t_crs - t0     !! [K] => [degC]

          if( tdeg>=0 ) then

             vis_crs = ( 1.718d0 + 4.9d-3*tdeg ) * 1.d-5

          else

             vis_crs = ( 1.718d0 + 4.9d-3*tdeg                        &
                  &                                   -1.2d-5*tdeg*tdeg ) * 1.d-5
          end if

          !### air mean free path [m] ###!

          dtmp = dsqrt(8.d0*mass_air*1.d-3/(ONE_PI*rrst*t_crs))

          lmd_crs = (2.d0*vis_crs)/(p_crs*dtmp)

!            call electro_coalescence_efficiency_coulomb(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
!            call electro_coalescence_efficiency_image_charge(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)

          if ( (sd_r1 > (sd_r2*1.0d2)) .or. (sd_r2 > (sd_r1*1.0d2)) ) then
!             call electro_coalescence_efficiency_image_charge(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
          else
!             call electro_coalescence_efficiency_conducting_sphere(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
          end if
          
          if ((sd_r1< 1.0e-7) .or. (sd_r2<1.0e-7)) then 
          !   eff_elc = 0.0
          end if 
          eff_elc = 0.0        

          crate_elc = ONE_PI* (sd_r1+sd_r2)**2 * abs(sd_vz1-sd_vz2) * eff_elc
 
         



          ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
               / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
       !    ivvol = 1.d0 / real(dx_sdm*dy_sdm,kind=r8)               &
       !          &               / real(zph_crs(i,j,k+1)-zph_crs(i,j,k),kind=r8)
  
          
!          crate_elc = 0.0d0

          c_rate(tc) = c_rate(tc) + &
               & crate_elc * real(sdm_dtcol,kind=RP)        &
               &                    * ivvol(k,i,j)

          c_rate(tc) = max(c_rate(tc),0.0d0)
          
         

         end do

      end do

! -----

    ! Get effective collision for "Brownian Coagulation and Scavenging
    ! (Seinfeld & Pandis,2006)" 
    ! This is mechanisim is effective for droplets less than micrometer-size ( below 1um )
    if( sdm_colbrwn>0 ) then
!!$        do n=1,hfreq_max
!!$
!!$          do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)
          
          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             !### information of super-droplets ###

             !! radius of water parts in droplets

             sd_rw1 = sd_r(icptc)    !! [m]
             sd_rw2 = sd_r(icptp)

             !! mass and volume of aerosol parts in droplets

             sd_tmasl1 = 0.0_RP
             sd_tmasl2 = 0.0_RP
             sd_tvasl1 = 0.0_RP
             sd_tvasl2 = 0.0_RP

             do k=1,22

                s = idx_nasl(k)

                sd_tmasl1 = sd_tmasl1 + sd_asl(icptc,s) * dmask(k)
                sd_tmasl2 = sd_tmasl2 + sd_asl(icptp,s) * dmask(k)

                sd_tvasl1 = sd_tvasl1                                    &
                     + sd_asl(icptc,s)/sd_aslrho(s) * dmask(k)
                sd_tvasl2 = sd_tvasl2                                    &
                     + sd_asl(icptp,s)/sd_aslrho(s) * dmask(k)

             end do

             sd_tmasl1 = sd_tmasl1  * 1.E-3_RP    !! [g]=>[kg]
             sd_tmasl2 = sd_tmasl2  * 1.E-3_RP

             !! diameter and mass and volume of droplets

             dtmp = ONE_PI * F_THRD

             sd_v1 = sd_tvasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1)
             sd_v2 = sd_tvasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2)

             sd_m1 = sd_tmasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1) * rw
             sd_m2 = sd_tmasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2) * rw

             sd_dia1 = (6.0_RP*sd_v1/ONE_PI)**O_THRD
             sd_dia2 = (6.0_RP*sd_v2/ONE_PI)**O_THRD

             !### location of super-droplets ###

             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             p_crs = pres_scale(k,i,j) !! [Pa]
             t_crs = t_scale(k,i,j)    !! [K]

             !### dynamic viscosity [Pa*s]  ###!
             !### (Pruppacher & Klett,1997) ###!
             tdeg = t_crs - t0     !! [K] => [degC]
             if( tdeg>=0.0_RP ) then
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg                        &
                     -1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             !### air mean free path [m] ###!
             dtmp = dsqrt(8.0_RP*mass_air*1.E-3_RP/(ONE_PI*rrst*t_crs))
             lmd_crs = (2.0_RP*vis_crs)/(p_crs*dtmp)

             !### slip correction of droplets [-]  ###!
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia1/lmd_crs)
             sd_cc1 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia1)
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia2/lmd_crs)
             sd_cc2 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia2)

             !### diffusion term [m*m/s] ###!
             dtmp = (boltz*t_crs)/(3.0_RP*ONE_PI*vis_crs)
             sd_d1 = dtmp * (sd_cc1/sd_dia1)
             sd_d2 = dtmp * (sd_cc2/sd_dia2)

             !### velocity term [m/s] ###!
             dtmp = (8.0_RP*boltz*t_crs)/ONE_PI
             sd_c1 = dsqrt(dtmp/sd_m1)
             sd_c2 = dsqrt(dtmp/sd_m2)

             !### mean free path of droplets [m] ###!
             dtmp = 8.0_RP/ONE_PI
             sd_lmd1 = dtmp * (sd_d1/sd_c1)
             sd_lmd2 = dtmp * (sd_d2/sd_c2)

             !### length term [m] ###!
             dtmp = (sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)&
                  - (sd_dia1*sd_dia1+sd_lmd1*sd_lmd1)**1.50_RP
             sd_g1 = dtmp/(3.0_RP*sd_dia1*sd_lmd1) - sd_dia1
             dtmp = (sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)&
                  - (sd_dia2*sd_dia2+sd_lmd2*sd_lmd2)**1.50_RP
             sd_g2 = dtmp/(3.0_RP*sd_dia2*sd_lmd2) - sd_dia2

             !### Brownian Coagulation Coefficient K12 [m3/s] ###!
             sumdia = sd_dia1 + sd_dia2
             sumd   = sd_d1   + sd_d2
             sumc   = dsqrt( sd_c1*sd_c1 + sd_c2*sd_c2 )
             sumg   = dsqrt( sd_g1*sd_g1 + sd_g2*sd_g2 )

             dtmp = sumdia/(sumdia+2.0_RP*sumg) + (8.0_RP*sumd)/(sumdia*sumc)
             k12 = 2.0_RP*ONE_PI * sumdia*sumd/dtmp

             !### add effective collision [-] ###!
             c_rate(tc) = c_rate(tc)                                     &
                  + k12 * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do

       end do

    end if

! -----

! Output kenel value for debug then stop
    if( output_kernel  ) then

       if ( PRC_IsMaster ) then
          write(fname,'("coalescence-kernel.txt")')

!          open( FILE_SD, iostat=stat, err=110,                              &
!               &      file=trim(fname), access='sequential', form='formatted',    &
!               &      blank='null', action='write' )
!110       call chkerr(stat)

          open( FILE_SD, iostat=stat,                              &
               &      file=trim(fname), access='sequential', form='formatted',    &
               &      blank='null', action='write' )

!OCL NORECURRENCE
          do m=1,gnum
             if( sort_freq(m) <= 1 ) cycle
             
             sort_tag0m = sort_tag0(m)
             sort_freqm = sort_freq(m)

             do n=1,sort_freqm/2

                tc = sort_tag0m + n
                tp = tc + sort_freqm/2

                ! index in center grid
                i = floor(sd_ri(icptc))+1
                j = floor(sd_rj(icptc))+1
                k = floor(sd_rk(icptc))+1

                ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                     / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

                write(FILE_SD,'(3G15.8)'), &   
                     & sd_r(icp(tc)),sd_r(icp(tp)),c_rate(tc)/(real(sdm_dtcol,kind=RP) * ivvol(k,i,j))

             end do
          end do
         
          close(FILE_SD)      

       end if
       call PRC_MPIstop

    end if

    ! Get total effective collision of droplets
    
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          sd_nmax  = max( sd_n(icptc), sd_n(icptp) )
          ! maximum multiplicity
          ipremium = sort_freqm - 1 + iand(sort_freqm,1)
          ! IAND(sort_freq(i),1) => even:0, odd:1

          c_rate(tc) = c_rate(tc) * real( sd_nmax*ipremium, kind=RP )
       end do
    end do

    ! Stochastic coalescence process.
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          !### set coalescence count ###!
          sd_ncol = int( c_rate(tc), kind=RP )
          frac = c_rate(tc) - real( sd_ncol, kind=RP )

          ! judge coalescence by random number and fractional part

          if( sd_rand(tc) < frac ) then
             sd_ncol =  sd_ncol + 1
          end if

          if( sd_ncol<=0 ) cycle  !! no coalesecense

          !### coalescence procudure ###!

          if( sd_n(icptc) > sd_n(icptp) ) then

             sd_n1  = sd_n( icptc )
             sd_r1  = sd_r( icptc )
             sd_rk1 = sd_rk( icptc )
             sd_m1  = sd_r1 * sd_r1 * sd_r1
             sd_li1 = sd_liqice( icptc )

             sd_n2  = sd_n( icptp )
             sd_r2  = sd_r( icptp )
             sd_m2  = sd_r2 * sd_r2 * sd_r2
             sd_li2 = sd_liqice( icptp )

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptc,s )
                sd_asl2(s) = sd_asl( icptp,s )
             end do

          else

             sd_n1  = sd_n( icptp )
             sd_r1  = sd_r( icptp )
             sd_rk1 = sd_rk( icptp )
             sd_m1  = sd_r1 * sd_r1 * sd_r1
             sd_li1 = sd_liqice( icptp )

             sd_n2  = sd_n( icptc )
             sd_r2  = sd_r( icptc )
             sd_m2  = sd_r2 * sd_r2 * sd_r2
             sd_li2 = sd_liqice( icptc )

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptp,s )
                sd_asl2(s) = sd_asl( icptc,s )
             end do

          end if

          sd_ncol = min( sd_ncol, int(sd_n1/sd_n2,kind=RP) )

          if( sd_n1 > sd_n2*sd_ncol ) then

             sd_n1 = sd_n1 - sd_n2*sd_ncol
             sd_r2 = exp( O_THRD                                      &
                  * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
             !ORG           sd_r2 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD
             sd_li2 = max(sd_li1,sd_li2) ! STAT_LIQ=01, STAT_ICE=10_i2

             do k=1,22
                s = idx_nasl(k)
                dtmp = sd_asl1(s) * real(sd_ncol,kind=RP)
                sd_asl2(s) = sd_asl2(s) + dmask(k) * dtmp
             end do

          else

             !! coalescent SDs with same multiplicity
             !!  - do not change the order of equations for
             !!  - vectorization

             sd_n1 = int( sd_n2/2, kind=RP )
             sd_n2 = sd_n2 - sd_n1

             sd_r1 = exp( O_THRD                                      &
                  * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
             !ORG           sd_r1 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD
             sd_r2 = sd_r1

             sd_li1 = max(sd_li1,sd_li2) ! STAT_LIQ=01, STAT_ICE=10_i2
             sd_li2 = sd_li1

             do k=1,22

                s = idx_nasl(k)
                dtmp = sd_asl1(s)*real(sd_ncol,kind=RP) + sd_asl2(s)
                sd_asl1(s) = sd_asl1(s) + dmask(k)*(dtmp-sd_asl1(s))
                sd_asl2(s) = sd_asl1(s)

             end do

             !! invalid by collisions between SDs with
             !! sd_n1=sd_n2*sd_ncol and sd_n2=1

             if( sd_n1==0 ) then
                sd_rk1 = INVALID
             end if

          end if

          !! This never happens
!!$            !! check muliplicity
!!$
!!$            if( sd_n1>(2.0_RP**63._RP) .or. sd_n2>(2.0_RP**63._RP) ) then
!!$               iexced = -1
!!$               cycle
!!$            end if

          if( sd_n(icptc) > sd_n(icptp) ) then
             sd_n( icptc )  = sd_n1
             sd_r( icptc )  = sd_r1
             sd_rk( icptc ) = sd_rk1
             sd_liqice( icptc ) = sd_li1

             sd_n( icptp )  = sd_n2
             sd_r( icptp )  = sd_r2
             sd_liqice( icptp ) = sd_li2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptc,s ) = sd_asl1(s)
                sd_asl( icptp,s ) = sd_asl2(s)
             end do


          else

             sd_n( icptp )  = sd_n1
             sd_r( icptp )  = sd_r1
             sd_rk( icptp ) = sd_rk1
             sd_liqice( icptp ) = sd_li1

             sd_n( icptc )  = sd_n2
             sd_r( icptc )  = sd_r2
             sd_liqice( icptc ) = sd_li2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptp,s ) = sd_asl1(s)
                sd_asl( icptc,s ) = sd_asl2(s)
             end do

          end if

       end do

    end do

    ! Deallocate
    deallocate( fsort_tag  )
    deallocate( fsort_freq )

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_coales",1,1)
#endif
    return
  end subroutine sdm_coales
!***********************************************************************
      subroutine electro_coalescence_efficiency_coulomb(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
        ! calculate the Coulomb force efficiency
        use m_sdm_common, only: ONE_PI
        real(RP), intent(out) :: eff_elc    ! electro collision efficiency [] 
        real(RP), intent(in)  :: sd_r1      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_r2      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_vz1     ! terminal velocity of particle 1 [m/2] 
        real(RP), intent(in)  :: sd_vz2     ! terminal velocity of particle 2 [m/2] 
        real(RP), intent(in)  :: lmd_crs ! air mean free path [m]
        real(RP), intent(in)  :: vis_crs ! dynamic viscosity [Pa s]
       
        real(RP) :: k_elc ! constant for collection efficiency of electrostatic [Nm^2/C^2]
        real(RP) :: sd_q1 ! electrical charge of particle 1 [C]
        real(RP) :: sd_q2 ! electrical charge of particle 2 [C]
        real(RP) :: sd_rl ! electrical charge of particle 1 [C]
        real(RP) :: sd_rs ! electrical charge of particle 2 [C]
  
        real(RP) :: Cc1  ! Cunningham number of particle 1
        real(RP) :: Cc2  ! Cunningham number of particle 2

           k_elc = 9.0e9
           sd_q1 = alpha*0.83e-6*(2*sd_r1)**2
           sd_q2 = alpha*0.83e-6*(2*sd_r2)**2
          
           if (sd_r1 > sd_r2) then
                sd_rl = sd_r1
                sd_rs = sd_r2
           else
                sd_rl = sd_r2
                sd_rs = sd_r1
           end if

 
           if (sd_q1 < 1.6e-19) then
              sd_q1 = 1.6e-19
           end if

           if (sd_q2 < 1.6e-19) then
              sd_q2 = 1.6e-19
           end if

           Cc1 = 1.0+2.0*lmd_crs/(2.0*sd_r1)*(1.257+0.4*exp(-1.1*(2.0*sd_r1)/2.0/lmd_crs))
           Cc2 = 1.0+2.0*lmd_crs/(2.0*sd_r2)*(1.257+0.4*exp(-1.1*(2.0*sd_r2)/2.0/lmd_crs))

           eff_elc = 1.0*k_elc * sd_q1 * sd_q2/(3.0* ONE_PI * vis_crs *  &
       &             abs(sd_vz1-sd_vz2) * (sd_r1+sd_r2+0.01*sd_rl/2.0)**2.0) *          &
       &             (Cc1/sd_r1+Cc2/sd_r2)
   
        return
      end subroutine electro_coalescence_efficiency_coulomb
!***********************************************************************
      subroutine electro_coalescence_efficiency_conducting_sphere(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
        ! calculate the conducting sphere electro force efficiency
        use m_sdm_common, only: ONE_PI
        real(RP), intent(out) :: eff_elc    ! Conducting sphere electro force coalesence efficiency [] 
        real(RP), intent(in)  :: sd_r1      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_r2      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_vz1     ! terminal velocity of particle 1 [m/2] 
        real(RP), intent(in)  :: sd_vz2     ! terminal velocity of particle 2 [m/2] 
        real(RP), intent(in)  :: lmd_crs    ! air mean free path [m]
        real(RP), intent(in)  :: vis_crs    ! dynamic viscosity [Pa s]
       
        real(RP) :: k_elc=9.0D9     ! constant for collection efficiency of electrostatic [Nm^2/C^2]
        real(RP) :: sd_rl     ! radius of larger particle [m]
        real(RP) :: sd_rs     ! radius of smaller particle [m]
        real(RP) :: sd_ql     ! electrical charge of larger particle [C]
        real(RP) :: sd_qs     ! electrical charge of smaller particle [C]
        real(RP) :: Ccl       ! Cunningham number of larger particle
        real(RP) :: Ccs       ! Cunningham number of smaller particle

        integer ixe,iye,ize,ike,iue,ine,ire   !index
        real(RP) sr_elc,rn_elc,rara,kqda,N_elc,D1R1,D1R2,D2R2,u1_elc,u2_elc,b_elc  !basic
        real(RP) kpj,alf,gam,F5,F6,F7,Fcs    !force calculation
        real(RP) zeta,yae,ybe,yce,y2a,y2b,y2c,tae,tbe,tce,wae,wbe,wce   !iue loop
        real(RP) S0A,S0B,S0C,T1A,T1B,T1C
        real(RP) zeta1,U0A,U1A,U0B,U1B,U0C,U1C    !ire loop
        real(RP) S01,S02,S03,p11,p12,p22    !S and p
        real(RP) K52,K56,K510,K513,K514,K517,K518,K521,K522    !F5
        real(RP) K62,K66,K610,K613,K614,K617,K618,K621,K622    !F6
        real(RP) K72,K76,K710,K713,K714,K717,K718,K721,K722    !F7
        real(RP) T11,T12,T13,U01,U02,U03,U11,U12,U13    !T and U
        real(RP) T11ans,T12ans,T13ans,U01ans,U02ans,U03ans,U11ans,U12ans,U13ans
        real(RP) S01ans,S02ans,S03ans,T11total,T12total,T13total
        real(RP) S01total,S02total,S03total,U01total,U02total,U03total
        real(RP) U11total,U12total,U13total

	T11ans=0.D0
        T12ans=0.D0
        T13ans=0.D0
        U01ans=0.D0
        U02ans=0.D0
        U03ans=0.D0
        U11ans=0.D0
        U12ans=0.D0
        U13ans=0.D0
        S01ans=0.D0
        S02ans=0.D0
        S03ans=0.D0
        T11total=0.D0
        T12total=0.D0
        T13total=0.D0
        S01total=0.D0
        S02total=0.D0
        S03total=0.D0
        U01total=0.D0
        U02total=0.D0
        U03total=0.D0
        U11total=0.D0
        U12total=0.D0
        U13total=0.D0


	if (sd_r1 > sd_r2) then
		sd_rl = sd_r1
             	sd_rs = sd_r2
      	else
             	sd_rl = sd_r2
             	sd_rs = sd_r1
       	end if

	! calculate the charge amount (Andronache 2004)
       	sd_ql = alpha*0.83D-6*(2.D0*sd_rl)**2
       	sd_qs = alpha*0.83D-6*(2.D0*sd_rs)**2

	! limiter to make it larger than 1e
       	if ((sd_rs<50.0D-6) .and. (sd_qs>40.0*1.6e-19)) then
       	   sd_ql = 40.0*1.6D-19
       	end if

       	if ((sd_rs < 50.0D-06) .and. (sd_qs<5.0*1.6e-19)) then
           sd_qs = 5.0*1.6D-19
      	end if

        ! Cunningham correction factor   
       	Ccl = 1.0D0+2.0D0*lmd_crs/(2.0D0*sd_rl)*(1.257D0+0.4D0*exp(-1.1D0*(2.0D0*sd_rl)/2.0D0/lmd_crs))
       	Ccs = 1.0D0+2.0D0*lmd_crs/(2.0D0*sd_rs)*(1.257D0+0.4D0*exp(-1.1D0*(2.0D0*sd_rs)/2.0D0/lmd_crs))

      		sr_elc=0.01D0    !s/sd_rs, parameter to determine the separation distance of two droplets
      !		rn_elc=(sd_rs/sd_rl)*(1.D0+sr_elc*(10**((log10(sd_rl/sd_rs)**0.8)))+sd_rl/sd_rs)

             rn_elc   =(sd_rs/sd_rl)*(1.D0+sr_elc*(sd_rl/sd_rs/2.0)+sd_rl/sd_rs)

             	rara = sd_rl/sd_rs       !droplet/aerosol ratio
             	kqda = sd_ql/sd_qs       !Charge ratio
             	D1R1 = dble((0.5D0)*(rn_elc+(1.D0/rn_elc)-(1.D0/((rara**2)*rn_elc))))
             	D1R2 = dble((0.5D0)*((rn_elc*rara)+(rara/rn_elc)-(1.D0/(rara*rn_elc))))
             	D2R2 = dble((0.5D0)*((rn_elc*rara)+(1.D0/(rara*rn_elc))-(rara/rn_elc)))
             	u1_elc = dble(log(D1R1+(((D1R1**2)-1.D0)**0.5D0)))
             	u2_elc = dble(log(D2R2+(((D2R2**2)-1.D0)**0.5D0)))
             	b_elc  = dble(u1_elc+u2_elc)
             	kpj = dble(1.D0/((D1R2**2)-(rara**2)))
             	alf  = dble(exp(2.D0*u1_elc))
             	gam  = dble((0.5D0)*(exp(2.D0*u2_elc)+1.D0))
             	ike = 24 ! originally 200
             	N_elc = dble(2.D0*ike + 1.D0)        

             	!**ine loop**
             	ine=0
               	do iye=1,ike,1  
               		S02ans=dble(S02ans+((exp(-(2.0D0*ine*u2_elc+u2_elc)))/(1.0D0-exp(-(2.0D0*ine*b_elc+b_elc)))))
               		S03ans=dble(S03ans+((exp(-(2.0D0*ine*u1_elc+u1_elc)))/(1.0D0-exp(-(2.0D0*ine*b_elc+b_elc)))))
               		S01ans=dble(S01ans+(exp(-(2.D0*ine*b_elc+b_elc))/(1.D0-exp(-(2.D0*ine*b_elc+b_elc)))))
               		T11ans=dble(T11ans+(((2.D0*ine+1.D0)*exp((2.D0*ine+1.D0)*(2.D0*u1_elc+u2_elc)-2.D0*(2.D0*ine*b_elc+b_elc)))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))**2)))
               		T12ans=dble(T12ans+(((2.D0*ine+1.D0)*exp((2.D0*ine+1.D0)*(u1_elc+u2_elc)-2.D0*(2.D0*ine*b_elc+b_elc)))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))**2)))
               		T13ans=dble(T13ans+(((2.D0*ine+1.D0)*exp((2.D0*ine*u2_elc+u2_elc)-2.D0*(2.D0*ine*b_elc+b_elc)))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))**2)))
               		U01ans=dble(U01ans+((exp((2.D0*ine+1.D0)*(2.D0*u1_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+(2.D0*ine+3.D0))))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))*(1.D0-exp(-(2.D0*ine*b_elc+3.D0*b_elc))))))
               		U02ans=dble(U02ans+((exp((2.D0*ine+1.D0)*(u1_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+(2.D0*ine+3.D0))))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))*(1.D0-exp(-(2.D0*ine*b_elc+3.D0*b_elc))))))
               		U03ans=dble(U03ans+((exp((2.D0*ine*u2_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+(2.D0*ine+3.D0))))/&
           		&      ((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))*(1.D0-exp(-(2.D0*ine*b_elc+3.D0*b_elc))))))
               		U11ans=dble(U11ans+(((2.D0*ine+1.D0)*exp((2.D0*ine+1.D0)*(2.D0*u1_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+&
           		&      (2.D0*ine+3.D0))))/((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))*(1.D0-exp(-(2.D0*ine*b_elc+3.D0*b_elc))))))
               		U12ans=dble(U12ans+(((2.D0*ine+1.D0)*exp((2.D0*ine+1.D0)*(u1_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+&
           		&      (2.D0*ine+3.D0))))/((1.D0-exp(-(2.D0*ine*b_elc+b_elc)))*(1.D0-exp(-(2.D0*ine*b_elc+3.D0*b_elc))))))
               		U13ans=dble(U13ans+(((2.D0*ine+1.D0)*exp((2.D0*ine*u2_elc+u2_elc)-b_elc*((2.D0*ine+1.D0)+&
           		&      (2.D0*ine+3.D0))))/((1.D0-exp(-(2.D0*ine+1.D0)*b_elc))*(1.D0-exp(-(2.D0*ine+3.D0)*b_elc)))))
               		ine=ine+1
              	enddo

             	!**iue loop**
               	iue=0
               	ire=1
               	do ize=1,2,1  ! 2 was 120 originally
               		zeta=0.
               		yae = dble((2.D0*u1_elc+u2_elc)-(b_elc*iue)-(2.D0*b_elc))
               		ybe = dble((u1_elc+u2_elc)-(b_elc*iue)-(2.D0*b_elc))
               		yce = dble(u2_elc-(b_elc*iue)-(2.D0*b_elc))
               		y2a = dble(0.D0-(b_elc*iue)-b_elc)
               		y2b = dble(u1_elc-(b_elc*iue)-b_elc)
               		y2c = dble(u2_elc-(b_elc*iue)-b_elc)
               		tae = dble(exp(2.D0*yae))
               		tbe = dble(exp(2.D0*ybe))
               		tce = dble(exp(2.D0*yce))
               		wae = dble(1.D0/(1.D0-tae))
               		wbe = dble(1.D0/(1.D0-tbe))
               		wce = dble(1.D0/(1.D0-tce))
               		S0A=dble(((exp((2.D0*ike+1.D0)*(y2a)))/(1.D0-(exp(2.0D0*y2a)))))
               		S0B=dble(((exp((2.D0*ike+1.D0)*(y2b)))/(1.D0-(exp(2.0D0*y2b)))))
               		S0C=dble(((exp((2.D0*ike+1.D0)*(y2c)))/(1.D0-(exp(2.0D0*y2c)))))
               		S01total=dble(S01total+S0A)
               		S02total=dble(S02total+S0B)
               		S03total=dble(S03total+S0C)
               		T1A=dble((iue+1.D0)*((2.D0*tae*wae**2+N_elc*wae)*(exp(N_elc*yae))))
               		T11total=dble(T11total+T1A)
               		T1B=dble((iue+1.D0)*((2.D0*tbe*wbe**2+N_elc*wbe)*(exp(N_elc*ybe))))
               		T12total=dble(T12total+T1B)
               		T1C=dble((iue+1.D0)*((2.D0*tce*wce**2+N_elc*wce)*(exp(N_elc*yce))))
               		T13total=dble(T13total+T1C)

               		!**ire loop**
                 	do ixe=1,iue+1,1
                 		zeta1=dble(exp(-2.0D0*ire*(u1_elc+u2_elc)))
                 		zeta = dble(zeta + zeta1)
                 		U0A=dble((zeta)*(wae*exp(N_elc*yae)))
                 		U01total=dble(U01total+U0A)
                 		U1A=dble((zeta)*(2.D0*tae*wae**2+N_elc*wae)*(exp(N_elc*yae)))
                 		U11total=dble(U11total+U1A)
                 		U0B=dble((zeta)*(wbe*exp(N_elc*ybe)))
                 		U02total=dble(U02total+U0B)
                 		U1B=dble((zeta)*(2.D0*tbe*wbe**2+N_elc*wbe)*(exp(N_elc*ybe)))
                 		U12total=dble(U12total+U1B)
                 		U0C=dble((zeta)*(wce*exp(N_elc*yce)))
                 		U03total=dble(U03total+U0C)
                 		U1C=dble((zeta)*(2.D0*tce*wce**2+N_elc*wce)*(exp(N_elc*yce)))
                 		U13total=dble(U13total+U1C)
                 		ire=ire+1
                 	enddo
               		iue=iue+1
           	enddo

             	!**Calculate S and p values**
             	S01=dble(S01ans+S01total)
             	S02=dble(S02ans+S02total)
             	S03=dble(S03ans+S03total)
             	p11=dble(S02/(2.0D0*((S02*S03)-(S01)**2)))
             	p22=dble(S03/(2.0D0*((S02*S03)-(S01)**2)))
             	p12=dble(S01/(2.0D0*((S02*S03)-(S01)**2)))

             	!**Calculate T,U and K values**
             	U13=dble(U13ans+U13total)
             	K522=dble(-gam*p11*p11)
             	K622=dble(-2.D0*gam*p11*p12)
             	K722=dble(-gam*p12*p12)
             	U03=dble(U03ans+U03total)
             	K521=dble(-gam*p11*p11)
             	K621=dble(-2.D0*gam*p11*p12)
             	K721=dble(-gam*p12*p12)
             	U12=dble(U12ans+U12total)
             	K518=dble(gam*p11*p12*(alf+1.D0))
             	K618=dble(gam*(alf+1.D0)*(p11*p22+p12*p12))
             	K718=dble(gam*p12*p22*(alf+1.D0))
             	U02=dble(U02ans+U02total)
             	K517=dble(gam*p11*p12*(alf+1.D0))
             	K617=dble(gam*(alf+1.D0)*(p11*p22+p12*p12))
             	K717=dble(gam*p12*p22*(alf+1.D0))
             	U11=dble(U11ans+U11total)
             	K514=dble(-gam*alf*p12*p12)
             	K614=dble(-2.D0*gam*alf*p12*p22)
             	K714=dble(-gam*alf*p22*p22)
             	U01=dble(U01ans+U01total)
             	K513=dble(-gam*alf*p12*p12)
             	K613=dble(-2.D0*gam*alf*p12*p22)
             	K713=dble(-gam*alf*p22*p22)
             	T13=dble(T13ans+T13total)
             	K510=dble(p11*p11)
             	K610=dble(2.D0*p11*p12)
             	K710=dble(p12*p12)
             	T12=dble(T12ans+T12total)
             	K56=dble(-2.D0*p11*p12)
             	K66=dble(-2.D0*(p11*p22+p12*p12))
             	K76=dble(-2.D0*p12*p22)
             	T11=dble(T11ans+T11total)
             	K52=dble(p12*p12)
             	K62=dble(2.D0*p12*p22)
             	K72=dble(p22*p22)
 
            	!**Force Coefficient Calculations**

             	F5 = kpj*((K52*T11)+(K56*T12)+(K510*T13)+(K513*U01)+(K514*U11)+(K517*U02)+&
               		(K518*U12)+(K521*U03)+(K522*U13))
             	F6 = kpj*((K62*T11)+(K66*T12)+(K610*T13)+(K613*U01)+(K614*U11)+(K617*U02)+&
               		(K618*U12)+(K621*U03)+(K622*U13))
             	F7 = kpj*((K72*T11)+(K76*T12)+(K710*T13)+(K713*U01)+(K714*U11)+(K717*U02)+&
               		(K718*U12)+(K721*U03)+(K722*U13))

             	Fcs = k_elc * (  sd_ql * sd_qs * F6 + sd_qs**2 * F7 + sd_ql**2 * F5) / sd_rs**2    !Fcs is positive for attractive
    
               if (Fcs < 0) then 
                  Fcs = 0.0
               endif
               
             	eff_elc = 1.0D0 * Fcs /(3.D0*ONE_PI * vis_crs * abs(sd_vz1 - sd_vz2)) * (Ccl/sd_rl + Ccs/sd_rs)
 
        return
      end subroutine electro_coalescence_efficiency_conducting_sphere
!***********************************************************************
      subroutine electro_coalescence_efficiency_image_charge(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
        ! calculate the image charge electro force efficiency
        use m_sdm_common, only: ONE_PI
        real(RP), intent(out) :: eff_elc    ! Conducting sphere electro force coalesence efficiency [] 
        real(RP), intent(in)  :: sd_r1      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_r2      ! radius of particle 1 [m]
        real(RP), intent(in)  :: sd_vz1     ! terminal velocity of particle 1 [m/2] 
        real(RP), intent(in)  :: sd_vz2     ! terminal velocity of particle 2 [m/2] 
        real(RP), intent(in)  :: lmd_crs    ! air mean free path [m]
        real(RP), intent(in)  :: vis_crs    ! dynamic viscosity [Pa s]

        real(RP) :: sd_rl     ! radius of larger particle [m]
        real(RP) :: sd_rs     ! radius of smaller particle [m]
        real(RP) :: sd_ql     ! electrical charge of larger particle [C]
        real(RP) :: sd_qs     ! electrical charge of smaller particle [C]
        real(RP) :: Ccl       ! Cunningham number of larger particle
        real(RP) :: Ccs       ! Cunningham number of smaller particle

        real(RP) sr_elc,rn_elc
        real(RP) Fimg    !force calculation

	if (sd_r1 > sd_r2) then
           sd_rl = sd_r1
           sd_rs = sd_r2
        else
           sd_rl = sd_r2
           sd_rs = sd_r1
        end if

	! calculate the charge amount (Andronache 2004)
       	sd_ql = alpha*0.83D-6*(2.D0*sd_rl)**2
       	sd_qs = alpha*0.83D-6*(2.D0*sd_rs)**2

	! limiter to make it larger than 1e
       	if ((sd_ql < 5.0* 1.6D-19) .and. (sd_rl>1.0e-7)) then
           sd_ql = 1.0*1.6D-19
       	end if

       	if ((sd_qs < 5.0*1.6D-19) .and. (sd_rs>1.0e-7)) then
           sd_qs = 1.0*1.6D-19
      	end if

       ! Cunningham correction factor
       Ccl = 1.0D0+2.0D0*lmd_crs/(2.0D0*sd_rl)*(1.257D0+0.4D0*exp(-1.1D0*(2.0D0*sd_rl)/2.0D0/lmd_crs))
       Ccs = 1.0D0+2.0D0*lmd_crs/(2.0D0*sd_rs)*(1.257D0+0.4D0*exp(-1.1D0*(2.0D0*sd_rs)/2.0D0/lmd_crs))

       sr_elc=0.01D0    !s/sd_rs, parameter to determine the separation distance of two droplets
!       rn_elc=(sd_rs/sd_rl)*(1.D0+sr_elc+sd_rl/sd_rs)
!       rn_elc=(sd_rs/sd_rl)*(1.D0+sr_elc*(10**((log10(sd_rl/sd_rs)**0.8)))+sd_rl/sd_rs)

!       rn_elc   =(sd_rs/sd_rl)*(1.D0+sr_elc*(1e-7/sd_rs)*(sd_rl/sd_rs)+sd_rl/sd_rs)
       rn_elc=(sd_rs/sd_rl)*(1.D0+sr_elc*(sd_rl/sd_rs/2.0)+sd_rl/sd_rs)

       if ((sd_ql> sd_qs*100.) .or. (sd_qs>sd_ql*100.)) then
!          rn_elc=(sd_rs/sd_rl)*(1.D0+sr_elc*(10**((log10(sd_rl/sd_rs)**1.2)))+sd_rl/sd_rs)
       end if

       Fimg = ((4.0d0*9.0d9 * (sd_qs**2)*rn_elc/((rn_elc**2 -1.0d0)**2))&
            &                + 4.0d0*(9.0d9* sd_qs * sd_ql / rn_elc**2) &
            &                - 4.0d0*(9.0d9 * (sd_qs**2)/((rn_elc**2)*rn_elc)))/sd_rl**2
       if (Fimg<0) then
           Fimg = 0
       end if

       eff_elc = 1.0D0 * Fimg /(3.D0*ONE_PI * vis_crs * abs(sd_vz1 - sd_vz2)) * (Ccl/sd_rl + Ccs/sd_rs)
        
       return
     end subroutine electro_coalescence_efficiency_image_charge
!***********************************************************************
end module m_sdm_coalescence
