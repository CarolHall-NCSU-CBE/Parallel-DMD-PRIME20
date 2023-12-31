! =========================================================================
! This module defines global parameters and variables used for simulations
! Last modified on 10/12/2023 by Van Nguyen
! =========================================================================

      module GLOBAL
	use inputreadin

      integer, parameter :: k12=selected_int_kind(12)
      integer, parameter :: k1=selected_int_kind(1)
      real drandm,dtime
      integer iflag,time
      external dtime,drandm,time,srand

!VN: Create allocatable arrays:
	integer(kind=k1) :: ev_code(10000,10000)
	integer, allocatable :: hp1(:),hp2(:)
	integer :: nptnr(10000)
	integer, allocatable :: npt(:),nb(:),na_npt(:),npt_dn(:),dnnab(:),nnabdn(:)
	integer :: tlinks(10000),tlinks2(10000)
	integer :: clinks(10000)
	real, allocatable :: old_rx(:),old_ry(:),old_rz(:)

!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      
      integer ncim1,ncai,ncaj,nnjp1
!LR: Added a third species numbeads variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer res(0:360,0:360), chnnum(10000)
	real bm(10000), bdln(175), tim(10000)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      !integer bptnr((nop1+nop2)),identity((nop1+nop2)),nptnr((nop1+nop2)+3)
	integer identity(10000), bptnr(10000)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
	!integer npt((nop1+nop2)+1),nb(maxnbs*(nop1+nop2)+1),na_npt((nop1+nop2)),npt_dn((nop1+nop2)+1),dnnab(maxnbs*(nop1+nop2)+1),nnabdn((nop1+nop2))
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer bin(0:numbin+1),nbin
      integer, dimension(:), allocatable::cell,wrap_map
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Added firstside3 to hold side-chain data for species 3
      integer num_cell
      integer map(n_nab_cell)
#ifdef equil
      integer(kind=k12) quarter, hb_sum(4),hh_sum(4)
      integer qt_count(4)
#endif
      integer(kind=k12) coll,ncoll
!      integer*8 coll,ncoll     
      logical success
!LR: Added a third species chnln variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real analysis_boxl,boxl,boxl_orig,t,analysis_setemp,setemp,ev_param(3,50)
!LR: Here, I did not add a chnln3 variable, because these are bond-only arrays, and the nanoparticle cannot have bonds.
      real hdelr,rlsq(50)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real sigma(28),welldia(28),epsilon(28)
      real sigma_sq(28,28),sigma_2b(28,28),welldia_sq(28,28),ep_sqrt(28,28),sig_max_all,shlddia_sq(28,28)
      real tfalse,pi,interval,t_fact,interval_max,sortsize,tbin_off,n_forced
      real width, half
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real ep(28,28),bds(28,28),wel(28,28) !! by mookyung
!LR: Here, I did not add a chnln3 variable, because these are bond-only arrays, and the nanoparticle cannot have bonds.
      real del_bdln(175),del_blrn(175),del_blrc(175) !! by mookyung
      real sqz610(5,28) !! by mookyung
      integer noptotal

!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Note - I am not sure who originally coded the -glycine flag, but I found it to be extremely buggy and incomplete.
!LR: Edits I made within -glycine flagged code do not constitute my usage of that feature, I was just trying to be as complete with my edits as possible.
#ifdef glycine 
      real bdln_dummy((nop1+nop2))
#endif

!VN: Arrays that are converted to allocatable arrays
!	integer(kind=k1) ev_code(10000,10000)
!LR: Added hp3 to hold hydrophobic data for species 3
	integer fside1(31),fside2(31), hp(175)
	real bl_rn(175),bl_rc(175)
	real sv(6,10000)
	integer extra_repuls(10000,4), coltype(10000)
      	real shder_dist1,shder_dist2,shder_dist3,shder_dist4		
end module global

