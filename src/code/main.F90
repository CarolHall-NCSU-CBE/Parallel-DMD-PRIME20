!1234567890123456789012345678901234567890123456789012345678901234567890123
!        10        20        30        40        50        60        70
!	discontinuous molecular dynamics main program
 
!     	variables defined in the makefile: 

!     	nop      number of particles
!     	ncoll    number of collisions desired
!     	setemp   simulation temperature 
!     	chnln    # of residues ** nop and chnln depend on each other **
!     	numbeads # of beads per chain 
!     	nphipsi  number of times to take phi-psi data
!     	pwl      well diameter for side chain beads = pwl*(bead diameter)
!     	numsheets  =3 for preset ideal fibril w/ 3 sheets of 4 chains each


!     	flags defined in the makefile:

!     	canon    flag to implement ghost collisions
!     	no_h     quick way to set all side chains to polar regardless
!      	  of hp.inp file
!     	no_hbs   quick way to turn off hydrogen bonds
!     	avs_full to generate avs and phi-psi data over the entire run
!     	runh     to set initial h-bonds for ideal helix or helices
!               ** must also use a starting configuration that's helical
!     	runb     to set initial h-bonds for ideal helix (8 16-mers) or helices
!               ** must also use a starting configuration that's helical
!     	runs     to set initial h-bonds for a specific ideal 12-chain 
!               fibril that is 3 sheets of 4 chains each
!               ** must also use this specific fibril starting config
!     	runsc    to set initial h-bonds for a specific ideal 12-chain
!               fibril that is 3 sheets of 4 chains each and extra loose chains
!     	runr     to start the simulation from a random (previous) config
!     	repuls2  makes chain nbors of hydrogen-bonded pair appear to each
!               other to be larger than they actually are ... this is 
!               based on wolynes method (jcp 110 (23) 1999 pp13096-29)
!               and is used to make h-bonds more stable
!     	write_phipsi  to write out phi and psi bond angle information


!     	parameters set in header:

!     	dnc = c_alpha_i to n_i (covalent bond)
!     	dcc = c_alpha_i to c_i (covalent bond)
!     	dcn = c_i to n_i+1, peptide bond (covalent bond)
!     	dtie = c_alpha_i to n_i+1 (pseudobond)
!     	dtie2 = c_i to c_alpha_i+1 and c_i to n_i (pseudobond)
!     	dcaca = c_alpha_i to c_alpha_i+1 (pseudobond)
!     	del:  bond tolerance; bond collision occurs at (1-del)*bdln,
!            bond extension occurs at (1+del)*bdln
!     	n_b_hydro = required minimum number of residues between 
!                  hydrophobic interactions
!     	n_b_hbond = required minimum number of residues between h bonds
!     	nfsteps: number of false steps to take between actual position 
!               updates and checks for neighbor list updates
!     	also used at very end of code to avoid end configuration overlaps
!     	sqz#     11 different local overlaps defined in make_code to 
!               allow near neighbors to overlap resulting in appropriate
!               allowed phi-psi regions for residues  
!               based on wolynes method (jcp 110 (23) 1999 pp13096-29)
   
#include "def.h"
#include "header.f"
 
      	program dmd
        
	use mpi
      	use global
	use inputreadin

      	implicit none

      	logical over,update,xpulse_del
      	character*64 filename  
      	real*8 tarray_nv(2),delta_nv,extime_nv,colrat_nv
      	real*8 tgho,xr,rating,avegtime,ered,tred,sumvel,e_int
      	real*8 old_tfalse,w,v1,v2,v3,r,fact,ran_non,ran_brk,tstar
      	real*8 rxij,ryij,rzij,vxij,vyij,vzij,bij,vijsq,rijsq,diff,rg_avg,e2e_avg
	integer totalsteps

#ifdef equil
      	real per_hb,per_hh
#endif

#ifndef canon
      	real*8 pre_energy,e_pot,e_kin
#endif

      	integer boundbad, unboundbad, pre_boundbad,pre_unboundbad,hb_partner
      	integer i,j,ii,jj,k,k_j,kk,kk_j,kkk,l,m,mm,num_before_sheet,evcode,writeprop_count
      	integer nbrsum,nl,nm,hb_alpha,hb_ii,hb_ij,hh_ii,hh_ij,n_interval
      	integer*8 nupdates,numghosts,nforcedupdate
      	integer*8 ncoll_max,nevents(30)
      	real*8  ehh_ii,ehh_ij

#ifdef write_phipsi
      	integer nphipsi
#endif

#ifdef __ifc
      integer ifc_time(2),idumb
#endif

! Yiming: mpi variable initialization
		integer cols,row,temp_numsent,n,binsize
		parameter (cols=31)
		parameter (binsize=70)
		real  buffer(cols),ans(cols),x,sort_time_list(binsize),temp_ans(binsize,cols+1)
		integer y,z,kstart,kend, myid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
		integer sort_bead_list(binsize),sort_nptnr_list(binsize)
		integer numsent,numsent_curr, numrcvd, sender, numevents, numevents_old,numevents_new, anstype
		integer  dynamic_numevents,count_round
        	  
		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	
		master=0
		if (myid .eq. master) then !Yiming: the start of master node

!LR: Calculates the total number of particles in the system.  I changed this to be a modular variable,
! as there are a number of places that this streamlines the code, and it will hopefully prevent some of the 
! tediousness that I had in adapting Dave's non-modular 2-species code to a modular 3-species code.

		!VN: get paths for file openning and creating:
	call getcwd(rundir)
	call get_command_argument(0,path)
	realpath = scan (path,'/',back)
	mydir = path(1:realpath)
      	data nevents / 30*0/
      	pi=4.d0*datan(1.d0)
	call readinputs()
	noptotal = nop1+nop2
	call allocatearrays()
      	if (newold == 0) then
		totalsteps = annealingsteps+numsim
	else
		totalsteps = numsim
	endif	
	do stepcount = 1,totalsteps
		if (newold == 0) then
			if (stepcount .le. annealingsteps) then
      				setemp = temps(stepcount)
      				setemp =  setemp*12.d0
      				ncoll = collset(stepcount)
			else
				setemp = simtemp
				setemp =  setemp*12.d0
				ncoll = simcoll
			endif
			if (stepcount == 1) then		
				call genconfig()
			endif
		else
			setemp = simtemp
			setemp =  setemp*12.d0
			ncoll = simcoll
		endif
      	filename = 'results/run'//'1000.lastvel'
      	inquire(file=filename, exist = success)
      	call file_opener()
      	write(fileout,*)'number of collisions requested', ncoll

#ifdef equil
      	ncoll_max=ncoll*4
      	quarter = ncoll/4
      	do k = 1, 4
         	hb_sum(k) = 0
         	hh_sum(k) = 0
         	qt_count(k) = 0
      	enddo
#endif
        numevents = 0
	    t_fact = 0.00005D0
      	n_forced = 150.0d0
      	n_interval = 75000/8
      	interval = t_fact/dsqrt(setemp)
      	interval_max = n_forced*interval
      	sortsize = interval_max/dble(numbin)
      	write(fileout,*) 't_fact=',t_fact
      	write(fileout,*) 'n_forced=',n_forced
      	write(fileout,*) 'n_interval=',n_interval

#ifdef canon
      	write(fileout,*)
      	avegtime=0.00005d0/dsqrt(setemp)
      	write(fileout,*)'avegtime is ',avegtime
      	write(fileout,*)
#endif

!     	to have different random numbers each time, use iflag = time()

#ifdef __ifc
      	call gettimeofday(ifc_time,idumb)
      	iflag =ifc_time(1)
#else
      	iflag =time()
#endif

!     	to have the same random numbers each time, use iflag = constant
     	iflag = 1058472402
      	write(fileout,*)'iflag=',iflag
      	xr=drandm(iflag)
      	call srand(iflag)
!     	each successive drandm call should be of the form drandm(0)
      	!call resultsfile_opener
      	writeprop_count = 0
      	call inputinfo()
      	call make_code() 

#ifdef canon
      	write(fileout,*)'implementing canon code'
#else
      	write(fileout,*)'leave out canon code to check for energy conservation'
#endif

#ifdef write_phipsi
      	nphipsi=0

      	do l=0,360
         	do m=0,360    
            		res(l,m)=0.d0
         	enddo
      	enddo
#endif

      	update=.false.
      	tstar=0.d0
      	t=0.d0
      	tfalse=0.d0
      	nbin = 1
      	tbin_off=0.d0

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
     	 	sv(1,k) = sv(1,k) - dnint(sv(1,k))
         	sv(2,k) = sv(2,k) - dnint(sv(2,k))
         	sv(3,k) = sv(3,k) - dnint(sv(3,k))
         	old_rx(k)=sv(1,k)
         	old_ry(k)=sv(2,k)
         	old_rz(k)=sv(3,k)
         	tim(k)=interval_max+ltstep
         	coltype(k)=-1
         	nptnr(k)=-1
         	npt(k)=(k-1)*maxnbs+1
         	npt_dn(k)=(k-1)*maxnbs+1
         	bptnr(k)=0
	 	do kk = 1, 4
	    		extra_repuls(k,kk)=0 
	    		extra_repuls(k,kk)=0
	    		extra_repuls(k,kk)=0
	 	end do
      	enddo

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	npt(noptotal+1)=(noptotal)*maxnbs+1
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	npt_dn(noptotal+1)=(noptotal)*maxnbs+1

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=(noptotal)+1,(noptotal)+3
         	nptnr(k)=-2
         	coltype(k)=-2
      	enddo 

65   	format(f24.16,i4,i5,i5)
66   	format(f24.16,i6,i5,i5)
      	write(fileout,*)

#ifdef runr
      	!open(unit=7,file='results/run'//fname_digits_pre//'.bptnr',status='old',form='unformatted')

      	do while (.true.)
        	read(prepartner,end=121) coll,bptnr
      	end do   
121   	continue

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k = 1, (noptotal)-1
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
         	do k_j = k+1, (noptotal)
            		if (identity(k)+identity(k_j).eq.5 .and. ev_code(k,k_j).eq.15) then
               		rxij=sv(1,k)-sv(1,k_j)
               		ryij=sv(2,k)-sv(2,k_j)
               		rzij=sv(3,k)-sv(3,k_j)
               		rxij=rxij-dnint(rxij)
               		ryij=ryij-dnint(ryij)
               		rzij=rzij-dnint(rzij)
               		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
               		diff=rijsq-welldia_sq(identity(k),identity(k_j))
               		if (diff .lt. 0.d0) then
               			if (k .le. nop1) then
                  				kk = k-((chnnum(k)-1)*numbeads1)
                  		!LR: Changed an open else statement to a constrained else-if statement, since only species 1 and species 2 can hydrogen bond.
               			elseif(k .le. nop1+nop2) then
                  				kk = k-nop1-((chnnum(k)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
               			endif
               			if (k_j .le. nop1) then
                  				kk_j = k_j-((chnnum(k_j)-1)*numbeads1)
               			!LR: Changed an open else statement to a constrained else-if statement, since only species 1 and species 2 can hydrogen bond.
               			elseif (k_j .le. nop1+nop2) then
                  				kk_j = k_j-nop1-((chnnum(k_j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
               			endif
                  			if ((k .le. nop1) .and. (k_j .le. nop1)) then
                  				if ((kk.ne.chnln1+1).and.(kk.ne.3*chnln1).and.(kk_j.ne.chnln1+1).and.(kk_j.ne.3*chnln1))then
                     				if (identity(k).eq.1) then
                        					call repuls_add(k,k_j)
                     				else
                        					call repuls_add(k_j,k)
                     				end if
                  				endif
                  			elseif ((k .le. nop1) .and. (k_j .gt. nop1)) then
                  				if ((kk.ne.chnln1+1).and.(kk.ne.3*chnln1).and.(kk_j.ne.numbeads1+chnln2+1).and.(kk_j.ne.numbeads1+3*chnln2))then
                     				if (identity(k).eq.1) then
                        					call repuls_add(k,k_j)
                     				else
                        					call repuls_add(k_j,k)
                     				end if
                  				endif
                  			elseif ((k .gt. nop1) .and. (k_j .gt. nop1)) then
                  				if ((kk.ne.numbeads1+chnln2+1).and.(kk.ne.numbeads1+3*chnln2).and.(kk_j.ne.numbeads1+chnln2+1).and.(kk_j.ne.numbeads1+3*chnln2))then
                     				if (identity(k).eq.1) then
                        					call repuls_add(k,k_j)
                     				else
                        					call repuls_add(k_j,k)
                     				end if
                  				endif
                  			endif
               		endif
            		endif
            		if (k_j .eq. bptnr(k)) then
               		if (k .le. nop1) then
                  			kk = k-((chnnum(k)-1)*numbeads1)
               		else
                  			kk = k-nop1-((chnnum(k)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
               		endif
               		if (k_j .le. nop1) then
                  			kk_j = k_j-((chnnum(k_j)-1)*numbeads1)
               		else
                  			kk_j = k_j-nop1-((chnnum(k_j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
               		endif
               		if (identity(k) == 1) then
                  			identity(k)=5
                  			identity(k_j)=8
               		else
                  			identity(k)=8
                  			identity(k_j)=5
               		endif
            		endif
         	enddo
      	enddo

#else
      	write(fileout,*)' '
      	write(fileout,*)'reassigning for helical h-bonds'
      	do m=1,nop/numbeads
         	do k=(m-1)*numbeads+chnln+5,(m-1)*numbeads+2*chnln
            		k_j=k+chnln-4
         	enddo
      	enddo
#endif

!     	initialize event counters
      	numghosts=0
      	nupdates=0
      	nforcedupdate=0

#ifdef debugging
      	boundbad=0
      	unboundbad=0
      	pre_boundbad=0
      	pre_unboundbad=0
#endif
      
!     	check starting configuration for overlaps w/subroutine
!     	and if variable over is true (overlap exists) stop program
      	call checkover(over)

      	if (over) then
         	write(fileout,*)'found overlap/underlap in initial configuration'
!        	call exit(-1)
      	else
         	write(fileout,*)'no overlap/underlap in initial configuration'
      	endif
      
!     	calculate energy with subroutine
      	call energy(ered,tred,sumvel,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij)
      	e_int=ered-0.5d0*sumvel

#ifndef canon
      	pre_energy = ered
      	e_pot = e_int
      	e_kin = 0.5d0*sumvel
#endif
            
      	call radgyr(rg_avg)
      	call end_to_end(e2e_avg)
      	write(fileout,*)'the initial total energy of system is ',ered
      	write(fileout,*)'the initial internal energy of system is ',e_int
      	write(fileout,*)'the initial kinetic energy of system is ',ered-e_int
      	write(fileout,*)'the initial temperature of system is ',tred
      	write(fileout,*)' '
      	write(fileout,*)'the initial number of alpha-helical hb',hb_alpha
      	write(fileout,*)'the initial number of hydrogen bonds',hb_ii+hb_ij
      	write(fileout,*)'the initial number of hydrophobic interactions',ehh_ii+ehh_ij
      	call flush(6)
      	write(rune,22223) 0,(t+tfalse)*dsqrt(setemp)/(sigma(1)*boxl_orig),ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
      	call flush(rune)

#ifdef write_phipsi
!    	call phipsi()
!    	nphipsi = nphipsi + 1
#endif

!     	write initial configuration to .conf file structures taken from full run
!     	write(799) numchains,numbeads,sngl(1.d0),sngl(sigma(2)/boxl_orig/2.d0),0
!     	call config()
!     	set up neighbor list and pointer with subroutine 
      	call nbor_setup()
      	num_cell=int(boxl/(sig_max_all*rl_const)*n_wrap)+2*n_wrap
      	write(fileout,*)
      	write(fileout,*)'at start, number of cells =',num_cell
      	allocate(cell(num_cell**3+1))
      	allocate(wrap_map(num_cell**3))
      	width = boxl/dble(num_cell-2*n_wrap)
      	half = boxl/2.d0
      	call cell_link
      	call nbor()
      	nbrsum=0
      	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)-1
         	nbrsum=nbrsum+na_npt(k)
      	enddo
      	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	write(fileout,*)'at start, avg number of nbors per particle=',dble(nbrsum)/dble(noptotal)
      	write(fileout,*)
      
#ifdef canon
      	tgho=0.d0

      	do while ((tgho .lt. 1.d-18) .or. (tgho .eq. 1.d0))
         	tgho = drandm(0)
      	enddo
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tim((noptotal)+1)=-1.d0*dlog(tgho)*avegtime*.0000001
#else 
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tim((noptotal)+1)=1000000000.d0   
#endif
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tim((noptotal)+2)=interval
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tim((noptotal)+3)=3.3/(dsqrt(setemp))+5
      	call events() 
      	call check_nc_int(boundbad, unboundbad)

#ifdef debugging
      	if ((boundbad .ne. pre_boundbad) .or. (unboundbad .ne. pre_unboundbad)) then
         	print*, 'boundbad and unboundbad',boundbad,unboundbad
         	pre_boundbad = boundbad
         	pre_unboundbad = unboundbad
      	end if
#endif
	call MPI_BCAST(noptotal,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST noptotal',noptotal
	call MPI_BCAST(chnnum,noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST chnnum'
	!call MPI_BCAST(identity,noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)  
	call MPI_BCAST(bm,noptotal,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST bm'
	call MPI_BCAST(ev_param,3*50,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	!call MPI_BCAST(ev_code,noptotal*noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(setemp,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sigma_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
!   add sidechain HP and backbone HB interaction	
	call MPI_BCAST(ep_sqrt,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ep,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(welldia_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fside1,chnln1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fside2,chnln2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST fside'
	call MPI_BCAST(hp,numbeads1+numbeads2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST hp'
	call MPI_BCAST(sqz610,5*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bdln,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bl_rc,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bl_rn,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_bdln,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_blrc,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_blrn,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	write(fileout,*) 'complete sending BCAST bondlength'
	
	call MPI_BCAST(shlddia_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist3,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist4,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(epsilon,28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

	
!**********************************************************************
!                         main loop begins
!**********************************************************************
    
      	delta_nv = dtime(tarray_nv)
      	coll = 0
		n = 0

		do i =1,binsize
        do j =1,cols+1
              temp_ans(i,j) = 0.d0
	    enddo
        enddo 		 
		
		old_tfalse = tfalse

#ifdef equil
6    	do while ((coll .le. ncoll) .and. (ncoll .le. ncoll_max))
#else
777    	do while (coll .le. ncoll)
#endif

        temp_numsent = 0
		count_round = 0
		
	!coll = coll + 1

	!xpulse_del=.false.
!    	finding first collision without tree:
      	old_tfalse = tfalse

      	do while (bin(nbin) == 0)
     		nbin = nbin + 1
     	end do

     	i = bin(nbin)
      	tfalse=tim(i)
		k = tlinks(i)

!Yiming: insertion sorting in the 1st bin
	    !numevents = 1
		numevents_new = 0
		numevents_old = 0
		numsent = 0
		numsent_curr = 0
        do kkk=1,numevents
			sort_bead_list(kkk) = 0
			sort_nptnr_list(kkk)= 0
			sort_time_list(kkk) = 0.d0
		enddo

		    numevents = 1
			sort_bead_list(1) = i 
			sort_time_list(1) = tim(i)	
			sort_nptnr_list(1)= nptnr(i)
		
		do while (k .ne. 0)
			numevents = numevents + 1
			sort_bead_list(numevents) = k 
			sort_nptnr_list(numevents) = nptnr(k)			
			sort_time_list(numevents) = tim(k)
     		k = tlinks(k)
    	end do

			do i = 2, numevents
				x = sort_time_list(i)
				y = sort_bead_list(i)
	            z = sort_nptnr_list(i)
				
				j = i - 1
				do while (j .ge. 1)
					if (sort_time_list(j).le. x) exit
					sort_time_list(j+1) = sort_time_list(j)
					sort_bead_list(j+1) = sort_bead_list(j)
					sort_nptnr_list(j+1) = sort_nptnr_list(j)
					j = j - 1
				end do
					sort_time_list(j+1) = x
					sort_bead_list(j+1) = y
					sort_nptnr_list(j+1) = z
			end do 

#ifdef debugging    
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable     
	if (i .gt. (noptotal)+3) then
	    	write(fileout,*) 'greater than (noptotal)+3', coll, i
	    	call exit(-1)
	end if
#endif

     	if (tfalse .lt. old_tfalse-ltstep) then
75         	format('event ',i10,' btwn ',i4,' and ', i4 ,' (event type ',i2,') is',f19.16, ' less then', f19.16)
            	write(fileout,75) coll,i,nptnr(i),coltype(i), tfalse, old_tfalse
    	endif

678        count_round = count_round + 1
           dynamic_numevents = numevents - (count_round-1)*(numprocs-1)
!   Yiming: send the colliding pair to each of the numprocs-1 of slaves
	do k=1,min(numprocs-1,dynamic_numevents) 
		i = sort_bead_list(k+(count_round-1)*(numprocs-1)) 
	!	write(22,*) 'line 331', numevents,k,sort_bead_list(k),sort_time_list(k),tim(i)!,nptnr(sort_bead_list(k))	
	!	write(fileout,*) 'line 331 2nd', coll,k+1,sort_bead_list(k+1),sort_time_list(k+1),tim(sort_bead_list(k+1)),nptnr(sort_bead_list(k+1))
		if (i .le. (noptotal)) then
			j = nptnr(i) 
			! write(22,*) coll,k,i,j,sv(1,i),sv(1,j),tfalse
			buffer(1) = dble(i)
			buffer(2) = dble(j)
		  !tfalse = sort_time_list(k) !write(fileout,*)'before coll: ',k,i,j,sv(1,i),sv(1,j),tfalse,buffer(2)
			buffer(3) = sort_time_list(k+(count_round-1)*(numprocs-1))	! same as tim(i)
			buffer(4) = sv(1,i)
			buffer(5) = sv(2,i)
			buffer(6) = sv(3,i)
			buffer(7) = sv(4,i)
			buffer(8) = sv(5,i)
			buffer(9) = sv(6,i)
			buffer(10) = sv(1,j)
			buffer(11) = sv(2,j)
			buffer(12) = sv(3,j)
			buffer(13) = sv(4,j)
			buffer(14) = sv(5,j)
			buffer(15) = sv(6,j)
			buffer(16) = dble(coltype(i))
            buffer(17) = dble(ev_code(i,j))
			buffer(18) = dble(identity(i))
			buffer(19) = dble(identity(j))
			buffer(20) = dble(bptnr(i))
			buffer(21) = dble(bptnr(j))
			buffer(22) = dble(extra_repuls(i,1))
			buffer(23) = dble(extra_repuls(i,2))
            buffer(24) = dble(extra_repuls(i,3))
			buffer(25) = dble(extra_repuls(i,4))
			buffer(26) = dble(extra_repuls(j,1))
			buffer(27) = dble(extra_repuls(j,2))
            buffer(28) = dble(extra_repuls(j,3))
			buffer(29) = dble(extra_repuls(j,4))
            !if (bptnr(i).ne.0) buffer(30) = dble(identity(bptnr(i)))
            !if (bptnr(j).ne.0) buffer(31) = dble(identity(bptnr(j)))
			if (extra_repuls(i,4).ne.0)  buffer(30) = dble(identity(extra_repuls(i,4)))
			if (extra_repuls(j,4).ne.0)  buffer(31) = dble(identity(extra_repuls(j,4)))	
			
       else if (i .eq. (noptotal+1)) then
			buffer(1) = dble(i)
			!tfalse = sort_time_list(k) 	!write(fileout,*)'before coll: ',k,i,sv(1,i),tfalse,buffer(2)
            buffer(3) = sort_time_list(k+(count_round-1)*(numprocs-1))
		!	write(fileout,*) 'line 380',i
       else if (i .eq. (noptotal+2)) then
			buffer(1) = dble(i)
			!tfalse = sort_time_list(k)
            buffer(3) = sort_time_list(k+(count_round-1)*(numprocs-1))
       else if (i .eq. (noptotal+3)) then
			buffer(1) = dble(i)
			!tfalse = sort_time_list(k)
            buffer(3) = sort_time_list(k+(count_round-1)*(numprocs-1))
       endif	   

	    call MPI_SEND(buffer,cols,MPI_DOUBLE_PRECISION,k,numsent+1,MPI_COMM_WORLD,ierr)
		!do kk = 1,cols
		!    buffer(kk) = 0.0
		!enddo
   		numsent = numsent + 1
	end do

8888    call MPI_RECV(ans,cols,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        sender = status(MPI_SOURCE)  ! identity of one slave node
		numsent_curr = status(MPI_TAG)    
           do k =1,cols
              temp_ans(numsent_curr,k) = ans(k)
		   enddo
		     temp_ans(numsent_curr,cols+1) = sender
   		!write(fileout,*) 'line403', sender,numsent_curr,temp_ans(numsent_curr,1),ans(1) 
		      temp_numsent = temp_numsent + 1
		   if (temp_numsent .lt. numsent)  goto 8888  ! one round events all received, thus proceed them
		   
      do n = 1+ (count_round-1)*(numprocs-1), temp_numsent
		coll = coll + 1
	   xpulse_del=.false.
		! if (mod(coll,1000000).eq.0) write(fileout,*) 'simulation progress% ', float(coll)/float(ncoll)*100.0
        i = int(temp_ans(n,1))
		tfalse = temp_ans(n,3)
		! write(fileout,*) int(coll),i,n,numsent,numevents 
		
	   if (i .le. (noptotal)) then
			j = int(temp_ans(n,2))				  
	!		evcode = ev_code(i,j)
			sv(1,i) = temp_ans(n,4)
			sv(2,i) = temp_ans(n,5)
			sv(3,i) = temp_ans(n,6)
			sv(4,i) = temp_ans(n,7)
			sv(5,i) = temp_ans(n,8)
			sv(6,i) = temp_ans(n,9)
			sv(1,j) = temp_ans(n,10)
			sv(2,j) = temp_ans(n,11)
			sv(3,j) = temp_ans(n,12)
			sv(4,j) = temp_ans(n,13)
			sv(5,j) = temp_ans(n,14)
		    sv(6,j) = temp_ans(n,15)
			coltype(i) = int(temp_ans(n,16))
			ev_code(i,j) = int(temp_ans(n,28))

	if (coltype(i).eq.14) then 
	   if (temp_ans(n,26).eq.1.0) then
	   if (identity(i).lt.identity(j)) then 
		extra_repuls(i,1) = int(temp_ans(n,17))
		extra_repuls(i,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(i,1)
		nnjp1 = extra_repuls(i,2)

		extra_repuls(ncaj,3)=i
		ev_code(ncaj,i)=xrepuls2
    	ev_code(i,ncaj)=xrepuls1

		extra_repuls(nnjp1,3)=i
		ev_code(i,nnjp1)=xrepuls1
     	ev_code(nnjp1,i)=xrepuls2		

		extra_repuls(j,1) = int(temp_ans(n,19))
		extra_repuls(j,2) = int(temp_ans(n,20))
		ncai = extra_repuls(j,1)
		ncim1 = extra_repuls(j,2)

		extra_repuls(ncai,3)=j
        ev_code(ncai,j)=xrepuls2
    	ev_code(j,ncai)=xrepuls1
		
		extra_repuls(ncim1,3)=j
		ev_code(j,ncim1)=xrepuls1
		ev_code(ncim1,j)=xrepuls2
	   else
		extra_repuls(j,1) = int(temp_ans(n,17))
		extra_repuls(j,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(j,1)
		nnjp1 = extra_repuls(j,2)

		extra_repuls(ncaj,3)=j
		ev_code(ncaj,j)=xrepuls2
    	ev_code(j,ncaj)=xrepuls1

		extra_repuls(nnjp1,3)=j
		ev_code(j,nnjp1)=xrepuls1
     	ev_code(nnjp1,j)=xrepuls2		

		extra_repuls(i,1) = int(temp_ans(n,19))
		extra_repuls(i,2) = int(temp_ans(n,20))
		ncai = extra_repuls(i,1)
		ncim1 = extra_repuls(i,2)

		extra_repuls(ncai,3)=i
        ev_code(ncai,i)=xrepuls2
    	ev_code(i,ncai)=xrepuls1
		
		extra_repuls(ncim1,3)=i
		ev_code(i,ncim1)=xrepuls1
		ev_code(ncim1,i)=xrepuls2	
       endif
	    extra_repuls(i,4)=j
    	extra_repuls(j,4)=i
		endif

	elseif (coltype(i).eq.20) then
	   if (int(temp_ans(n,24)+temp_ans(n,25)).eq.13) then	
			bptnr(i) = int(temp_ans(n,22))
			bptnr(j) = int(temp_ans(n,23))
			identity(i) = int(temp_ans(n,24))
			identity(j) = int(temp_ans(n,25))
	   if (temp_ans(n,26).eq.1.0) then
	   if (identity(i).lt.identity(j)) then 
		extra_repuls(i,1) = int(temp_ans(n,17))
		extra_repuls(i,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(i,1)
		nnjp1 = extra_repuls(i,2)

		extra_repuls(ncaj,3)=i
		ev_code(ncaj,i)=xrepuls2
    	ev_code(i,ncaj)=xrepuls1

		extra_repuls(nnjp1,3)=i
		ev_code(i,nnjp1)=xrepuls1
     	ev_code(nnjp1,i)=xrepuls2		

		extra_repuls(j,1) = int(temp_ans(n,19))
		extra_repuls(j,2) = int(temp_ans(n,20))
		ncai = extra_repuls(j,1)
		ncim1 = extra_repuls(j,2)

		extra_repuls(ncai,3)=j
        ev_code(ncai,j)=xrepuls2
    	ev_code(j,ncai)=xrepuls1
		
		extra_repuls(ncim1,3)=j
		ev_code(j,ncim1)=xrepuls1
		ev_code(ncim1,j)=xrepuls2
	   else
		extra_repuls(j,1) = int(temp_ans(n,17))
		extra_repuls(j,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(j,1)
		nnjp1 = extra_repuls(j,2)

		extra_repuls(ncaj,3)=j
		ev_code(ncaj,j)=xrepuls2
    	ev_code(j,ncaj)=xrepuls1

		extra_repuls(nnjp1,3)=j
		ev_code(j,nnjp1)=xrepuls1
     	ev_code(nnjp1,j)=xrepuls2		

		extra_repuls(i,1) = int(temp_ans(n,19))
		extra_repuls(i,2) = int(temp_ans(n,20))
		ncai = extra_repuls(i,1)
		ncim1 = extra_repuls(i,2)

		extra_repuls(ncai,3)=i
        ev_code(ncai,i)=xrepuls2
    	ev_code(i,ncai)=xrepuls1
		
		extra_repuls(ncim1,3)=i
		ev_code(i,ncim1)=xrepuls1
		ev_code(ncim1,i)=xrepuls2	
       endif
	    extra_repuls(i,4)=j
    	extra_repuls(j,4)=i
		endif	
      endif		
	elseif (coltype(i).eq.16) then 
	   if (temp_ans(n,21) .eq. 1.0)  then
	          xpulse_del = .true.
	    if (identity(i).lt.identity(j)) then 
		extra_repuls(i,1) = int(temp_ans(n,17))
		extra_repuls(i,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(i,1)
		nnjp1 = extra_repuls(i,2)

		ev_code(ncaj,i)=1
    	ev_code(i,ncaj)=1
		ev_code(i,nnjp1)=1
     	ev_code(nnjp1,i)=1

		extra_repuls(j,1) = int(temp_ans(n,19))
		extra_repuls(j,2) = int(temp_ans(n,20))
		ncai = extra_repuls(j,1)
		ncim1 = extra_repuls(j,2)

        ev_code(ncai,j)=1
    	ev_code(j,ncai)=1
		ev_code(j,ncim1)=1
		ev_code(ncim1,j)=1
		else
		extra_repuls(j,1) = int(temp_ans(n,17))
		extra_repuls(j,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(j,1)
		nnjp1 = extra_repuls(j,2)

		ev_code(ncaj,j)=1
    	ev_code(j,ncaj)=1
		ev_code(j,nnjp1)=1
     	ev_code(nnjp1,j)=1

		extra_repuls(i,1) = int(temp_ans(n,19))
		extra_repuls(i,2) = int(temp_ans(n,20))
		ncai = extra_repuls(i,1)
		ncim1 = extra_repuls(i,2)

        ev_code(ncai,i)=1
    	ev_code(i,ncai)=1
		ev_code(i,ncim1)=1
		ev_code(ncim1,i)=1		
		endif
	   endif
	elseif (coltype(i).eq.21) then 

	   if (temp_ans(n,21) .eq. 1.0)  then
	          xpulse_del = .true.
	    if (identity(i).lt.identity(j)) then 
		extra_repuls(i,1) = int(temp_ans(n,17))
		extra_repuls(i,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(i,1)
		nnjp1 = extra_repuls(i,2)

		ev_code(ncaj,i)=1
    	ev_code(i,ncaj)=1
		ev_code(i,nnjp1)=1
     	ev_code(nnjp1,i)=1

		extra_repuls(j,1) = int(temp_ans(n,19))
		extra_repuls(j,2) = int(temp_ans(n,20))
		ncai = extra_repuls(j,1)
		ncim1 = extra_repuls(j,2)

        ev_code(ncai,j)=1
    	ev_code(j,ncai)=1
		ev_code(j,ncim1)=1
		ev_code(ncim1,j)=1
		else
		extra_repuls(j,1) = int(temp_ans(n,17))
		extra_repuls(j,2) = int(temp_ans(n,18))
		ncaj = extra_repuls(j,1)
		nnjp1 = extra_repuls(j,2)

		ev_code(ncaj,j)=1
    	ev_code(j,ncaj)=1
		ev_code(j,nnjp1)=1
     	ev_code(nnjp1,j)=1

		extra_repuls(i,1) = int(temp_ans(n,19))
		extra_repuls(i,2) = int(temp_ans(n,20))
		ncai = extra_repuls(i,1)
		ncim1 = extra_repuls(i,2)

        ev_code(ncai,i)=1
    	ev_code(i,ncai)=1
		ev_code(i,ncim1)=1
		ev_code(ncim1,i)=1		
		endif
	   endif

	    if (temp_ans(n,29).eq.1.0) then	
			bptnr(i) = int(temp_ans(n,22))
			bptnr(j) = int(temp_ans(n,23))
			identity(i) = int(temp_ans(n,24))
			identity(j) = int(temp_ans(n,25))
        endif

	elseif (coltype(i).eq.24) then 	
       		hb_partner = int(temp_ans(n,27))	
        if (ev_code(i,j).lt.45) then
           if (int(temp_ans(n,24)+temp_ans(n,25)) .eq. 13) then
			bptnr(i) = int(temp_ans(n,22))
			bptnr(hb_partner) = int(temp_ans(n,23))
			identity(i) = int(temp_ans(n,24))
			identity(hb_partner) = int(temp_ans(n,25))
			endif
		else
           if (int(temp_ans(n,24)+temp_ans(n,25)).eq. 13) then
			bptnr(j) = int(temp_ans(n,22))
			bptnr(hb_partner) = int(temp_ans(n,23))
			identity(j) = int(temp_ans(n,24))
			identity(hb_partner) = int(temp_ans(n,25))		
		   endif
		endif		
	elseif (coltype(i).eq.25) then 		
       		hb_partner = int(temp_ans(n,27))
        if (ev_code(i,j).lt.45) then
            if (temp_ans(n,29) .eq. 1.0) then
			bptnr(i) = int(temp_ans(n,22))
			bptnr(hb_partner) = int(temp_ans(n,23))
			identity(i) = int(temp_ans(n,24))
			identity(hb_partner) = int(temp_ans(n,25))
			endif
		else
            if (temp_ans(n,29) .eq. 1.0) then
			bptnr(j) = int(temp_ans(n,22))
			bptnr(hb_partner) = int(temp_ans(n,23))
			identity(j) = int(temp_ans(n,24))
			identity(hb_partner) = int(temp_ans(n,25))		
		   endif
		endif

    endif	

!		write(fileout,*) 'line 320',numevents,coll,i,j,nptnr(i),sort_time_list(i),tim(i),sv(1,i),sv(1,j)
!		call eventdyn(i,j,evcode)
   
	    nevents(coltype(i)) = nevents(coltype(i))+1  ! 1=core collision; 2=bond collision; 3=bond stretch 
	
#ifdef debugging                        
#ifndef canon
        call energy(ered,tred,sumvel,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij)
        if (nint(ered) .ne. nint(pre_energy)) then !	
        !if (nint((ered-pre_energy)*100000.d0).gt.1.d0) then
			print*, 'no more conserving energy at coll', coll, i,j,coltype(i)
			print*, pre_energy, ered
			print*, e_pot, ered-0.5d0*sumvel
			print*, e_kin, 0.5d0*sumvel   
			call exit(-1)
        end if
          	e_pot=ered-0.5d0*sumvel
           	e_kin = 0.5d0*sumvel
#endif
#endif            
	    call partial_events(i,j,xpulse_del)  
		      	old_tfalse = tfalse ! Yiming

	!	write(fileout,*)'after coll, I and J:',coll,i,j,sv(1,i),sv(1,j),tim(i),tim(j)!,'bead 13',sv(1,13),sv(4,13),tim(13),nptnr(13),sv(1,nptnr(13)),sv(4,nptnr(13))
	
   ! if (coll.eq.1497) then
     !  do ii=1,noptotal+3	 
	 ! write(21,'(3i4,7f35.30)') coll,ii,nptnr(ii),tim(ii),sv(1,ii),sv(2,ii),sv(3,ii),sv(4,ii),sv(5,ii),sv(6,ii)
	!  write(21,'(3i4,7f35.30)') coll,i,nptnr(i),tim(i),sv(1,i),sv(2,i),sv(3,i),sv(4,i),sv(5,i),sv(6,i)
	 !    enddo
    !  endif
	  
		   numevents_old = numevents_new
		do while (bin(nbin) == 0)
     	   	nbin = nbin + 1
        end do
		ii = bin(nbin)
		k = tlinks(ii)		
		numevents_new = 1
			
     	do while (k .ne. 0)
			numevents_new = numevents_new + 1
     		k = tlinks(k)
    	end do
 		!write(fileout,*)'line 634',coll,numsent,numevents,numevents_new,numevents_old				

!Yiming: if next coll bead is related to this previous one, update their event time   
		   do kkk =n+1,numevents
         ! if ((sort_nptnr_list(kkk).eq.i).or.(sort_bead_list(kkk).eq.j).or.(sort_nptnr_list(kkk).eq.j)) then
		  if ((sort_bead_list(kkk).le.noptotal).and.(tim(sort_bead_list(kkk)).ne.sort_time_list(kkk))) then
		!  write(fileout,*) 'real1',coll,i

			! write(fileout,*)'line 419' !,coll,i,j,kkk,sort_bead_list(kkk),sort_nptnr_list(kkk)
		     go to 777	
		  endif	
		enddo	
		
		if (numevents_old.eq.0) then	
		 if ((numevents_new.ge.numevents).or.(numevents-numevents_new.ge.2)) then
           !    write(fileout,*) 'real2',coll,i

		    go to 777 ! resort 1st event bin 	 
		 endif	
	 		
		else if (numevents_old.gt.1) then
		 if ((numevents_new.ge.numevents_old).or.(numevents_old-numevents_new.ge.2)) then
		!    write(fileout,*) 'real3',coll,i

		    go to 777 ! resort 1st event bin 	 
		 endif
		endif					
		
    !   cycle
#ifdef canon
    elseif (i == (noptotal) + 1) then
3      	i = int(drandm(0)*(noptotal)) + 1
	    	if (i == (noptotal)+1) goto 3
!           t is total simulation time, split is time since last internal energy change
!           tfalse is time since last position update
!           implement ghost collision           	
            numghosts=numghosts+1
           	w=0.d0
!          	advance position of i
           	sv(1,i)=sv(1,i)+sv(4,i)*tfalse
           	sv(2,i)=sv(2,i)+sv(5,i)*tfalse
           	sv(3,i)=sv(3,i)+sv(6,i)*tfalse
4         	v1 = 2.d0*drandm(0)-1.d0
           	v2 = 2.d0*drandm(0)-1.d0
           	r=v1*v1+v2*v2

           	if ((r.eq.0.d0).or.(r.ge.1.d0)) goto 4
		 !if (coll.eq.2) write(22,'(2i4,6f35.30)') coll,i,setemp,bm(i),r,log(r),sqrt(-2.d0*setemp*bm(i)*log(r)/r),fact

           	fact=sqrt(-2.d0*setemp*bm(i)*log(r)/r)
           	sv(4,i) = v1*fact/bm(i)
           	sv(5,i) = v2*fact/bm(i)
       !  write(22,'(2i4,6f35.30)') coll,i,r,log(r),-2.d0*setemp*bm(i)*log(r)/r,fact,sv(4,i),sv(5,i)
			
5        	v1 = 2.d0*drandm(0)-1.d0
           	v2 = 2.d0*drandm(0)-1.d0
           	r=v1*v1+v2*v2

           	if ((r.eq.0.d0).or.(r.ge.1.d0)) goto 5
           	fact=sqrt(-2.d0*setemp*bm(i)*log(r)/r)
			
			! if (coll.eq.15) write(22,'(2i4,5f35.30)') coll,i,v1,r,setemp,fact,bm(i)		
           	sv(6,i) = v1*fact/bm(i)
!          	take care of false positions by rewinding
           	sv(1,i)=sv(1,i)-sv(4,i)*tfalse
           	sv(2,i)=sv(2,i)-sv(5,i)*tfalse
           	sv(3,i)=sv(3,i)-sv(6,i)*tfalse
!    		find new ghost time(s) for the bead(s) involved in the last
!          	event; otherwise, the ghost event time list doesn't change

            if (tlinks2((noptotal)+1) .ne. 0) call del_tbin((noptotal)+1)

	    	tgho=0.d0
	    	do while ((tgho .lt. 1.d-18) .or. (tgho .eq. 1.d0))
               	tgho = drandm(0) 
	    	enddo  
           
            tim((noptotal)+1) = -1.d0*alog(tgho)*avegtime + tfalse
            if (tim((noptotal)+1) .lt. interval_max) call add_tbin((noptotal)+1)
				
            if (tfalse .lt. old_tfalse) tfalse = old_tfalse	

	    	call partial_events(i,0,xpulse_del)

   !  if (coll.le.ncoll) then
      ! do ii=1,noptotal+3	 
	  !write(21,'(4i4,7f35.30)') coll,noptotal+1,ii,nptnr(ii),tim(ii),sv(1,ii),sv(2,ii),sv(3,ii),sv(4,ii),sv(5,ii),sv(6,ii)
  	  !   enddo
	!  write(21,'(4i4,7f35.30)') coll,noptotal+1,i,nptnr(i),tim(i),sv(1,i),sv(2,i),sv(3,i),sv(4,i),sv(5,i),sv(6,i)
    !  endif
	!		write(fileout,*)'after ghost coll:',coll,i,identity(i),sv(1,i),sv(4,i),sv(5,i),sv(6,i),tim(i),nptnr(i),tim(nptnr(i)),tim(noptotal+1) !,tfalse,old_tfalse!,'bead 13',sv(1,13),sv(4,13),tim(13),nptnr(13),sv(1,nptnr(13)),sv(4,nptnr(13))
				
			  old_tfalse = tfalse ! Yiming
					
!Yiming: if next coll bead is related to this previous one, update their event time   		
		   numevents_old = numevents_new
		do while (bin(nbin) == 0)
     	   	nbin = nbin + 1
        end do
		ii = bin(nbin)
		k = tlinks(ii)		
		numevents_new = 1
			
     	do while (k .ne. 0)
			numevents_new = numevents_new + 1
     		k = tlinks(k)
    	end do
 		!write(fileout,*)'line 634',coll,numsent,numevents,numevents_new,numevents_old

		if (n .lt. numevents) then
		   do kkk =n+1,numevents

		  if (tim(i).le.sort_time_list(kkk)) then !.or.(tim(nptnr(i)).le.sort_time_list(kkk))) then  
		!   write(fileout,*) 'ghost1',coll,i

		! write(fileout,*)'line 533' !,coll,i,kkk,sort_bead_list(kkk),sort_nptnr_list(kkk)
		     go to 777			
		 endif

		  if (tim(noptotal+1).le.tim(sort_bead_list(kkk))) then 
		  	!	   write(fileout,*) 'ghost2',coll,i

		! write(fileout,*)'line 522' !, coll,i,kkk,sort_bead_list(kkk),sort_nptnr_list(kkk)
		     go to 777			 
		 endif
	
        !  if ((sort_bead_list(kkk).eq.i).or.(sort_nptnr_list(kkk).eq.i)) then
		!  if ((sort_bead_list(kkk).le.noptotal).and.(tim(sort_bead_list(kkk)).ne.sort_time_list(kkk))) then
		  if (tim(sort_bead_list(kkk)).ne.sort_time_list(kkk)) then
		  	!	   write(fileout,*) 'ghost3',coll,i

		! write(fileout,*)'line 512' !,coll,i,kkk,sort_bead_list(kkk),sort_nptnr_list(kkk)
		     go to 777	
		  endif			 
		 
		enddo
		
      endif		
	  
!Yiming: pale event in 1st bin appears after accepting this collision  
!Yiming: if there is newly generated event into the 1st bin after accepting this collision  
		if (numevents_old.eq.0) then	
		 if ((numevents_new.ge.numevents).or.(numevents-numevents_new.ge.2)) then
		 		!   write(fileout,*) 'ghost4',coll,i

		    go to 777 ! go back to resort 1st event bin 	 
		 endif	
	 		
		else if (numevents_old.gt.1) then
		 if ((numevents_new.ge.numevents_old).or.(numevents_old-numevents_new.ge.2)) then
		  		!   write(fileout,*) 'ghost5',coll,i

		    go to 777 ! go back to resort 1st event bin 	 
		 endif
		endif
	  
       !     cycle
#endif
    
     else if (i == (noptotal) + 2) then
	  !      tfalse = ans(2)
!       	advance all particles to real positions and check if nbor list should be updated
	    	t = t + tfalse            
		!	write(fileout,*)'line 523', coll, t, tfalse
!           reduce collision time list by time of current collision

            do k=1,(noptotal)+3
              	tim(k)=tim(k)-tfalse
            enddo

	    	interval_max = interval_max - tfalse
	    	tbin_off = tbin_off + tfalse
!           	advance false positions to real positions
            	do k=1,(noptotal)
               	sv(1,k)=sv(1,k)+sv(4,k)*tfalse
               	sv(2,k)=sv(2,k)+sv(5,k)*tfalse
               	sv(3,k)=sv(3,k)+sv(6,k)*tfalse
            	enddo
	    	tfalse = 0.d0
!           	calculate displacement vector for each mc w/subroutine and
!           	set logical update to know whether to update neighbor lists
            	call displ(update)
!           	if variable update is true, update neighbor lists
            	if ((update) .or. (interval .gt. interval_max)) then
	       	if (update == .false.) then
                  		nforcedupdate=nforcedupdate+1
	          		n_forced = n_forced * 1.01d0
                  		print*, 'forced updated', t, 'n_forced', n_forced
             endif
               	interval_max = interval * n_forced
				sortsize = interval_max/dble(numbin)
               	tbin_off = 0.d0
               	update=.false.
               	nupdates=nupdates+1
					do k=1,(noptotal)
                  		sv(1,k) = sv(1,k) - anint(sv(1,k))
                  		sv(2,k) = sv(2,k) - anint(sv(2,k))
                  		sv(3,k) = sv(3,k) - anint(sv(3,k))
						old_rx(k)=sv(1,k)
                  		old_ry(k)=sv(2,k)
                  		old_rz(k)=sv(3,k)	  
					enddo	
		
               	call nbor()
!              	reset all times to interval_max+ltstep
					do k=1,(noptotal)
                  		tim(k)=interval_max+ltstep
                  		coltype(k)=-1
                  		nptnr(k)=-1
					enddo
               	call events()
				nbin = 1      !reset to check from the 1st bin
           endif
				if (tlinks2(i) .ne. 0) call del_tbin(i)            
					tim(i) = interval*0.999d0 
				call add_tbin(i)
				
      	old_tfalse = tfalse ! Yiming
		!write(12,*) coll, old_tfalse

		    go to 777 ! go back to resort 1st event bin 	
	
      !      	cycle
	
	else if (i == (noptotal) + 3) then
		!        tfalse = ans(2)
#ifdef debugging
		call checkover(over)
            	if (over) then
               	write(fileout,*)'found overlap/underlap at coll', coll
               	call exit(-1)
            	endif
            	call check_nc_int(boundbad, unboundbad)
            	if ((boundbad .ne. pre_boundbad) .or. (unboundbad .ne. pre_unboundbad)) then
               	print*, 'boundbad and unboundbad',boundbad,unboundbad
               	pre_boundbad = boundbad
               	pre_unboundbad = unboundbad
            	end if
#endif
	    	call flush(6)
           	call energy(ered,tred,sumvel,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij)
	    	call radgyr(rg_avg)
	    	call end_to_end(e2e_avg)
	    	write(rune,22223) coll,(t+tfalse)*sqrt(setemp)/(sigma(1)*boxl_orig),ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
            call flush(rune)

#ifndef canon
            if (nint(ered) .ne. nint(pre_energy)) then
            	print*, 'no more conserving energy', coll
               	print*, pre_energy, ered
               	print*, e_pot, ered-0.5d0*sumvel
               	print*, e_kin, 0.5d0*sumvel
               	call exit(-1)
            end if
				e_pot=ered-0.5d0*sumvel
            	e_kin = 0.5d0*sumvel
#endif
!           write each config. during the entire simul
            call config()
			
#ifdef write_phipsi	
	    	if (coll .gt. ncoll/2) then
               	call phipsi()
	       	nphipsi = nphipsi + 1
	    	endif
#endif

#ifdef equil
!	    	summation of hb and hb to check for equilibrium
	    	if (coll .ge. quarter) call sum(hb_ii+hb_ij,hh_ii+hh_ij)
#endif

			if (tlinks2(i) .ne. 0) call del_tbin(i)
			tim(i)=3.3d0/(sqrt(setemp)) + 5 + tfalse
			if (tim(i) .lt. interval_max) call add_tbin(i)
		!	write(fileout,*)'after noptotal+3 coll:',coll,1,sv(1,1),sv(4,1)	
			old_tfalse = tfalse
 		!write(12,*) coll, old_tfalse
	
	 end if	   


	if (coll.ge.ncoll) then 
	        call  MPI_SEND(1.0,0,MPI_DOUBLE_PRECISION,sender,0, MPI_COMM_WORLD, ierr)	
	 go to 999	
     endif
    enddo

     if (numsent .lt. numevents)  goto 678 ! events left in the bin, continue sendind to slaves
    !if (temp_numsent .lt. numsent) goto 8888        ! more events in the bin sent to slave nodes
	  
    Enddo
!  Yiming: simulation output and analysis

#ifdef equil
!     	check for equilibrium
      	if (coll-1 .eq. ncoll) then
	 	write(fileout,*) ' '
	 	write(fileout,*) 'at coll', ncoll
	 	call average(per_hb,per_hh)
	 	call flush(6)
	 	if ((per_hb .gt. 5.0) .or. (per_hh .gt. 5.0)) then
            		ncoll = ncoll*2
	    		quarter = ncoll/4
	    		call doubling
	    		goto 6
	 	end if
      	end if	    
#endif

!**********************************************************************
!                     main loop ends
!**********************************************************************

999    	delta_nv = dtime(tarray_nv)
!     	delta is in seconds - change to hours (or minutes)
      	extime_nv = delta_nv/3600.
      	colrat_nv = float(coll-1)/extime_nv/1.0e6
!     	advance false positions to real positions

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sv(1,k)=sv(1,k)+sv(4,k)*tfalse
         	sv(2,k)=sv(2,k)+sv(5,k)*tfalse
         	sv(3,k)=sv(3,k)+sv(6,k)*tfalse
         	sv(1,k) = sv(1,k) - dnint(sv(1,k))
         	sv(2,k) = sv(2,k) - dnint(sv(2,k))
         	sv(3,k) = sv(3,k) - dnint(sv(3,k))
      	enddo

      	tfalse=0.0     
!     	calculate and print out properties of interest
      	call energy(ered,tred,sumvel,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij)
      	e_int=ered-0.5d0*sumvel
      	call radgyr(rg_avg)
      	call end_to_end(e2e_avg)
      	write(rune,22223) coll,t*dsqrt(setemp)/(sigma(1)*boxl_orig),ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
      	write(fileout,*)' '
      	write(fileout,*)'the final total energy of system (e) is ',ered
     	write(fileout,*)'the final potential energy of system is ',e_int
      	write(fileout,*)'the final kinetic energy of system is ',ered-e_int
      	write(fileout,*)'the final temperature of system (kt) is ',tred
      	write(fileout,*)
      	write(fileout,*)'the final number of alpha-helical hb',hb_alpha
      	write(fileout,*)'the final number of hydrogen bonds',hb_ii+hb_ij
      	write(fileout,*)'the final number of hydrophobic interactions',ehh_ij+ehh_ij
      	write(fileout,*)' '

#ifdef write_phipsi
      	call phipsi()
      	nphipsi = nphipsi + 1
#endif

      	call config()
!     	check end configuration for overlaps w/subroutine
      	call checkover(over)
      	call scale_up()
!LR: I separated the analysis code from the main code, but I wanted to leave these here for reference.
	!analysis_setemp=setemp
	!analysis_boxl=boxl
	!call analysis_main(analysis_setemp,analysis_boxl)

#ifdef write_phipsi
!     	print out phi-psi data as fort.residue number plus 1000
      	nphipsi = nphipsi*((chnln1-2)*(nop1/numbeads1)+(chnln2-2)*(nop2/numbeads2))
      	do l=0,360
         	do m=0,360
	    		if (res(l,m) .ne. 0) then
               		nl=l-180
               		nm=m-180
               		write(runphipsi,*) nl,nm,float(res(l,m))/nphipsi*100.d0
	    		endif 
        	enddo
	enddo
#endif

      	filename = 'results/run'//fname_digits//'.pdb'
      	open(unit=11113, file=filename)
      	call write_rasmol
      	close(unit=11113)      
!     	********** time ********** 
      	write(fileout,*)'t=',t + tfalse
      	tstar=(t+tfalse)*dsqrt(setemp)/sigma(1)
!     	this is gulati's calculation:  tdl=t*dsqrt(tred)/sigma(2)
      	write(fileout,*)'total simulation time is ',t+tfalse
      	write(fileout,*)'dimensionless time is ',tstar
      	write(fileout,*)'collision rate is ',ncoll/tstar
      	write(fileout,*)
!     	********** event tally ********** 
      	do k = 1, 30
	 	if (nevents(k) .ne. 0) write(fileout,*)'event', k, nevents(k)/dble(ncoll)*100.0,'%'
      	enddo
      	write(fileout,*)  
     	write(fileout,*)'number of ghost events = ',numghosts
      	write(fileout,*)'    = ',dble(numghosts)/dble(ncoll)*100.,'%'
      	write(fileout,*)'number of  neighbor list updates= ',nupdates - nforcedupdate
      	write(fileout,*)'number of  forced neighbor list updates= ',nforcedupdate
      	write(fileout,*)
     
#ifdef debugging
      	write(*,*) 'number of times hbing particles caught w/ bad angle', boundbad
      	write(*,*) 'number of times non-hbing particles should have been hbing', unboundbad
#endif

      	write(fileout,*)
     	write(fileout,*)'details on the time in the main loop ...'
      	write(fileout,'('' execution time (cpu hours)  '',f15.4)') extime_nv
      	write(fileout,'('' millions of coll/cpu hour   '',f15.3)') colrat_nv

      	deallocate(cell)
      	deallocate(wrap_map)

22222 	format(i15,3f12.4,5i8,6f12.4)
22223 	format(i15,3f12.4,3i8,4f12.4)

      	go to 888
	enddo
	Else

!   Yiming: slave receives broadcasted constants from master
	call MPI_BCAST(noptotal,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(chnnum,noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	!call MPI_BCAST(identity,noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)  
	call MPI_BCAST(bm,noptotal,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ev_param,3*50,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	!call MPI_BCAST(ev_code,noptotal*noptotal,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(setemp,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sigma_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
!   add sidechain HP and backbone HB interaction	
	call MPI_BCAST(ep_sqrt,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ep,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(welldia_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fside1,chnln1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fside2,chnln2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(hp,numbeads1+numbeads2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sqz610,5*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bdln,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bl_rc,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bl_rn,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_bdln,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_blrc,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(del_blrn,chnln1+chnln2,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	
	call MPI_BCAST(shlddia_sq,28*28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist3,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(shder_dist4,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(epsilon,28,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
	
	!call MPI_BCAST(xrepuls1,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
	!call MPI_BCAST(xrepuls2,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

!   slave receive one event in the bin from master
90   call MPI_RECV(buffer,cols,MPI_DOUBLE_PRECISION, master,MPI_ANY_TAG, MPI_COMM_WORLD, status,ierr)

	if (status(MPI_TAG) .ne. 0) then
	     xpulse_del = .false.
		 hb_partner = 0 
		row=status(MPI_TAG) ! row is numsent 
	   do ii = 1,cols	
		ans(ii) = 0.d0
       enddo
  		i = int(buffer(1))
	!write(fileout,*) 'slave, line 1196', i,row 
!   Yiming: judge collision type & calculate  Δt, V, x
		if  (i .le. (noptotal)) then
		    j = int(buffer(2))
			tfalse = buffer(3)	
			sv(1,i) = buffer(4)
			sv(2,i) = buffer(5)
			sv(3,i) = buffer(6)
			sv(4,i) = buffer(7)
			sv(5,i) = buffer(8)
			sv(6,i) = buffer(9)	
			sv(1,j) = buffer(10)
			sv(2,j) = buffer(11)
			sv(3,j) = buffer(12)
			sv(4,j) = buffer(13)
			sv(5,j) = buffer(14)
			sv(6,j) = buffer(15)
			coltype(i)=int(buffer(16))
			evcode = int(buffer(17)) !ev_code(i,j)
			identity(i) = int(buffer(18))
			identity(j) = int(buffer(19))
			bptnr(i) = int(buffer(20))
			bptnr(j) = int(buffer(21))
			extra_repuls(i,1) = int(buffer(22))
			extra_repuls(i,2) = int(buffer(23))
            extra_repuls(i,3) = int(buffer(24))
			extra_repuls(i,4) = int(buffer(25))
			extra_repuls(j,1) = int(buffer(26))
			extra_repuls(j,2) = int(buffer(27))
            extra_repuls(j,3) = int(buffer(28))
			extra_repuls(j,4) = int(buffer(29))
        !    if (bptnr(i).ne.0) identity(bptnr(i)) = int(buffer(30))
        !    if (bptnr(j).ne.0) identity(bptnr(j)) = int(buffer(31))
			if (extra_repuls(i,4).ne.0) identity(extra_repuls(i,4)) = int(buffer(30))
			if (extra_repuls(j,4).ne.0) identity(extra_repuls(j,4)) = int(buffer(31))		
		!	write(fileout,*)'slave receive normal coll:',i,j,tfalse,evcode,coltype(i),sv(1,i),sv(4,i),bm(i),bm(j)
		
!!!!!HP+HB !!!!!
          	if (i .le. nop1) then
                  	ii = i-((chnnum(i)-1)*numbeads1)
                	 !LR: Changed an open else statement to a constrained else-if statement, since only species 1 
                	 ! and species 2 need the math for later "end-bead" work
            	elseif (i .le. nop1+nop2) then
                  	ii = i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
           	endif
           	if (j .le. nop1) then
                  	jj = j-((chnnum(j)-1)*numbeads1)
            	 !LR: Changed an open else statement to a constrained else-if statement, since only species 1 
            	 ! and species 2 need the math for later "end-bead" work
          	elseif (j .le. nop1+nop2) then
                  	jj = j-nop1-((chnnum(j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
         	endif
      		!evcode = ev_code(i,j)
			
!           	1=core collision; 2=bond collision; 3=bond stretch 
!           	if next event is square-well and involves n and c, must check 
!           	if the square-shoulder interaction should be turned on
            	if (coltype(i).eq.7) then
!             	event is either a capture or no event, tij has already been calc'd
	         	if ((extra_repuls(i,4) .eq. 0) .and. (extra_repuls(j,4) .eq. 0)) then
                    if ((i .le. nop1) .and. (j .le. nop1)) then
                        		if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.chnln1+1).and.(jj.ne.3*chnln1)) then
!                           		capture didn't involve an end bead
                            		if (identity(i) < identity(j)) then
!                                			i is n, j is c
		                   			call repuls_check(i,j,rating)
                            		else
!                               			i is c, j is n
	                           			call repuls_check(j,i,rating)
                            		endif
		              		       if (rating .le. 10.d0) then
!                                   		capture event occurs if rating < 10 (rating of 1 is the best)
                                        		coltype(i)=4
                            		else
!                                    		no event, only calcs will be to push i,j off discontinuity
                                       		coltype(i)=14
                            		endif
                        		else
!                 				n or c or both are terminal beads randomly decide if it should be a 
!		  				non-event based on non-event rate with non-end beads
                                  		ran_non= drandm(0)
                             		if (ran_non .le. 0.2d0) then
		        				coltype(i)=4
                               		else
                        				coltype(i)=9
                             		endif
                        		endif
              		elseif ((i .le. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    			capture didn't involve an end bead
                     			if (identity(i) < identity(j)) then
!                       				i is n, j is c
		        				call repuls_check(i,j,rating)
                     			else
!                       				i is c, j is n
	                				call repuls_check(j,i,rating)
                     			endif
		     				if (rating .le. 10.d0) then
!                    				capture event occurs if rating < 10 (rating of 1 is the best)
                        				coltype(i)=4
                     			else
!                       				no event, only calcs will be to push i,j off discontinuity
                        				coltype(i)=14
                     			endif
                  			else
!                 				n or c or both are terminal beads randomly decide if it should be a 
!		  				non-event based on non-event rate with non-end beads
                     			ran_non= drandm(0)
                     			if (ran_non .le. 0.2d0) then
		        				coltype(i)=4
                     			else
                        				coltype(i)=9
                     			endif
                  			endif
              		elseif ((i .gt. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    			capture didn't involve an end bead
                     			if (identity(i) < identity(j)) then
!                       				i is n, j is c
		        				call repuls_check(i,j,rating)
                     			else
!                       				i is c, j is n
	                				call repuls_check(j,i,rating)
                     			endif
		     				if (rating .le. 10.d0) then
!                    				capture event occurs if rating < 10 (rating of 1 is the best)
                        				coltype(i)=4
                     			else
!                       				no event, only calcs will be to push i,j off discontinuity
                        				coltype(i)=14
                     			endif
                  			else
!                 				n or c or both are terminal beads randomly decide if it should be a 
!		  				non-event based on non-event rate with non-end beads
                     			ran_non= drandm(0)
                     			if (ran_non .le. 0.2d0) then
		        				coltype(i)=4
                     			else
                        				coltype(i)=9
                     			endif
              			   endif
              		endif
	       	else
		  		coltype(i)=9
	       	endif
            	elseif (coltype(i).eq.10) then
!              	must decide if atoms are moving away by dissociation via square square shoulder (=5) or non-event (=15)
               	if (evcode < 45) then
		  		hb_partner = extra_repuls(i,4)
                  		if (identity(i) < identity(hb_partner)) then
!                    		i is n, hb_partner is c
                     		call repuls_check_3(i,hb_partner,j,rating)
                  		else
!                    		i is c, hb_partner is n
                     		call repuls_check_3(hb_partner,i,j,rating)
                  		endif
                  		if (rating .le. 10.d0) then
		     			coltype(i) = 5
 	          		else
 	             			coltype(i) = 15
                  		endif
               	else
                  		hb_partner = extra_repuls(j,4)
						!write(fileout,*) i,j,hb_partner
                  		if (identity(j) < identity(hb_partner)) then
!                    		j is n, hb_partner is c
                     		call repuls_check_3(j,hb_partner,i,rating)
                  		else
!                    		j is c, hb_partner is n
                     		call repuls_check_3(hb_partner,j,i,rating)
                  		endif
		  		       if (rating .le. 10.d0) then
                     		coltype(i) = 5
                  		else
                     		coltype(i) = 15
                   		endif
	         	endif
		elseif (coltype(i).eq.12) then
!              	must decide if atoms are moving toward each other via square square shoulder (=6) or non-event (=15) 
               	if (evcode < 45) then
                  		hb_partner = extra_repuls(i,4)
		  		    if (bptnr(i) == hb_partner) then
                     		call check_sigma(i,hb_partner,rating)
		     			if (rating .gt. 10.d0) then
                        			coltype(i) = 6
	             			else
	                			coltype(i) = 13
		     			endif
	          		else
                     		coltype(i) = 15
		  		    endif
                 else
                  		hb_partner = extra_repuls(j,4)
		  		        if (bptnr(j) == hb_partner) then
	             			call check_sigma(j,hb_partner,rating)  
                     		if (rating .gt. 10.d0) then
                        			coltype(i) = 6
	             			else
		        			coltype(i) = 13
		     			    endif
		   		       else
		     			    coltype(i) = 15
                       endif
	          	endif
          endif
!!!!!HP+HB !!!!!
		if (coltype(i) .lt. 14)  call eventdyn(i,j,evcode)
!!!!!HP+HB !!!!!
            	if (coltype(i).eq.20) then
               	   if (identity(i)+identity(j) .eq. 5) then
!                 		it's a hydrogen-bond capture
                  		bptnr(i)=j
                  		bptnr(j)=i
                  		identity(i)=identity(i)+4
                  		identity(j)=identity(j)+4
			ans(22) = dble(bptnr(i))
			ans(23) = dble(bptnr(j))
			ans(24) = dble(identity(i))
			ans(25) = dble(identity(j))
                  		if ((i .le. nop1) .and. (j .le. nop1)) then
                  			if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.chnln1+1).and.(jj.ne.3*chnln1)) then
!                 				capture didn't involve an end bead
                     			if (identity(i).lt.identity(j)) then
!                       				i is n, j is c
  		        				call repuls_add(i,j)
                     			else
!                       				i is c, j is n
		        				call repuls_add(j,i)
                     			endif
		  			       endif
                  		elseif ((i .le. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                 				capture didn't involve an end bead
                     			if (identity(i).lt.identity(j)) then
!                       				i is n, j is c
  		        				call repuls_add(i,j)
                     			else
!                       				i is c, j is n
		        				call repuls_add(j,i)
                     			endif
		  			        endif
                  		elseif ((i .gt. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                 				capture didn't involve an end bead
                     			if (identity(i).lt.identity(j)) then
!                       				i is n, j is c
  		        				call repuls_add(i,j)
                     			else
!                       				i is c, j is n
		        				call repuls_add(j,i)
                     			endif
		  			       endif
                	   endif
	            	end if
             
			 if (identity(i).lt.identity(j)) then
	        ans(17) = dble(extra_repuls(i,1))
	        ans(18) = dble(extra_repuls(i,2))
	        ans(19) = dble(extra_repuls(j,1))
	        ans(20) = dble(extra_repuls(j,2))
             else
	        ans(17) = dble(extra_repuls(j,1))
	        ans(18) = dble(extra_repuls(j,2))
	        ans(19) = dble(extra_repuls(i,1))
	        ans(20) = dble(extra_repuls(i,2))
             endif			 
			if (ans(17).ne.0.0) ans(26)=1.0
			
            	elseif (coltype(i).eq.21) then
!           		dissociation
	            	if (identity(i).le.8) then
                  		if ((i .le. nop1) .and. (j .le. nop1)) then
                  			if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.chnln1+1).and.(jj.ne.3*chnln1)) then
!                    			didn't involve an end bead
                     			if (identity(i).lt.(identity(j))) then
!                       				i is n, j is c
                        				call repuls_del_a(i,j)
                     			else
!                       				i is c, j is n
                        				call repuls_del_a(j,i)
                     			endif
                     			xpulse_del = .true.
                  			endif
                  		elseif ((i .le. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    			didn't involve an end bead
                     			if (identity(i).lt.(identity(j))) then
!                       				i is n, j is c
                        				call repuls_del_a(i,j)
                    				else
!                       				i is c, j is n
                        			       call repuls_del_a(j,i)
                     			endif
                     			xpulse_del = .true.
                  			endif
                  		elseif ((i .gt. nop1) .and. (j .gt. nop1)) then
                  			if ((ii.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    			didn't involve an end bead
                     			if (identity(i).lt.(identity(j))) then
!                       				i is n, j is c
                        				call repuls_del_a(i,j)
                     			else
!                       				i is c, j is n
                        				call repuls_del_a(j,i)
                     			endif
                     			xpulse_del = .true.
                  			endif
                  		endif
				      if (bptnr(i) == j) then
!                    		it's a hydrogen bond dissociation or break
                     		bptnr(i)=0
                     		bptnr(j)=0
                     		identity(i)=identity(i)-4
                     		identity(j)=identity(j)-4
			ans(22) = dble(bptnr(i))
			ans(23) = dble(bptnr(j))
			ans(24) = dble(identity(i))
			ans(25) = dble(identity(j))
			ans(29) = 1.0
	          		   end if
	         	endif

       			if (identity(i).lt.identity(j)) then
	        ans(17) = dble(extra_repuls(i,1))
	        ans(18) = dble(extra_repuls(i,2))
	        ans(19) = dble(extra_repuls(j,1))
	        ans(20) = dble(extra_repuls(j,2))
                else
	        ans(17) = dble(extra_repuls(j,1))
	        ans(18) = dble(extra_repuls(j,2))
	        ans(19) = dble(extra_repuls(i,1))
	        ans(20) = dble(extra_repuls(i,2))
                endif					
			if (xpulse_del) ans(21) = 1.0

			
		    elseif (coltype(i).eq.24) then
               	if (evcode < 45) then
                  		hb_partner = extra_repuls(i,4)
                  		if (identity(i) + identity(hb_partner) .eq. 5) then
                     		bptnr(i)=hb_partner
                     		bptnr(hb_partner)=i
                     		identity(i)=identity(i)+4
                     		identity(hb_partner)=identity(hb_partner)+4

			ans(22) = dble(bptnr(i))
			ans(23) = dble(bptnr(hb_partner))
			ans(24) = dble(identity(i))			
			ans(25) = dble(identity(hb_partner))						
                  		endif
               	else
                  		hb_partner = extra_repuls(j,4)
                  		if (identity(j) + identity(hb_partner) .eq. 5) then
                     		bptnr(j)=hb_partner
                     		bptnr(hb_partner)=j
                     		identity(j)=identity(j)+4
                     		identity(hb_partner)=identity(hb_partner)+4

			ans(22) = dble(bptnr(j))
			ans(23) = dble(bptnr(hb_partner))
			ans(24) = dble(identity(j))
			ans(25) = dble(identity(hb_partner))						
                  		endif

               	endif
			ans(27) = dble(hb_partner)
			
           	elseif (coltype(i).eq.25) then
               	if (evcode < 45) then   
                  		hb_partner = extra_repuls(i,4)
                  		if (identity(i) .ge. 5) then
                     		bptnr(i)=0
                     		bptnr(hb_partner)=0
                     		identity(i)=identity(i)-4
                     		identity(hb_partner)=identity(hb_partner)-4
			ans(22) = dble(bptnr(i))
			ans(23) = dble(bptnr(hb_partner))
			ans(24) = dble(identity(i))
			ans(25) = dble(identity(hb_partner))
			ans(29) = 1.0
                  		endif
               	else
                  		hb_partner = extra_repuls(j,4)
                  		if (identity(j) .ge. 5) then
                     		bptnr(j)=0
                     		bptnr(hb_partner)=0
                     		identity(j)=identity(j)-4
                     		identity(hb_partner)=identity(hb_partner)-4
			ans(22) = dble(bptnr(j))
			ans(23) = dble(bptnr(hb_partner))
			ans(24) = dble(identity(j))
			ans(25) = dble(identity(hb_partner))
            ans(29) = 1.0			
                  		endif
               	end if  	 
			ans(27) = dble(hb_partner)
			
            elseif (coltype(i).eq.14) then
!             	advance to get off square-well this is the only result of this type of collision
	       	     call bumpoff(i,j,evcode)
!	       	turning on the square-shoulder potential for the hb psuedobonds
                  	if ((i .le. nop1) .and. (j .le. nop1)) then
                  		if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.chnln1+1).and.(jj.ne.3*chnln1)) then
!                 			capture didn't involve an end bead
                     		if (identity(i).lt.identity(j)) then
!                       			i is n, j is c
  		        			call repuls_add(i,j)
                     		else
!                       			i is c, j is n
		        			call repuls_add(j,i)
                     		endif
		  		        endif
                  	elseif ((i .le. nop1) .and. (j .gt. nop1)) then
                  		if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                 			capture didn't involve an end bead
                     		if (identity(i).lt.identity(j)) then
!                       			i is n, j is c
  		        			call repuls_add(i,j)
                     		else
!                       			i is c, j is n
		        			call repuls_add(j,i)
                     		endif
		  		        endif
                  	elseif ((i .gt. nop1) .and. (j .gt. nop1)) then
                  		if ((ii.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                 			capture didn't involve an end bead
                     		if (identity(i).lt.identity(j)) then
!                       			i is n, j is c
  		        			call repuls_add(i,j)
                     		else
!                       			i is c, j is n
		        			call repuls_add(j,i)
                     		endif
		  		        endif
                	endif
            
			if (identity(i).lt.identity(j)) then			
			ans(17) = dble(extra_repuls(i,1))
	        ans(18) = dble(extra_repuls(i,2))
	        ans(19) = dble(extra_repuls(j,1))
	        ans(20) = dble(extra_repuls(j,2))
			else
			ans(17) = dble(extra_repuls(j,1))
	        ans(18) = dble(extra_repuls(j,2))
	        ans(19) = dble(extra_repuls(i,1))
	        ans(20) = dble(extra_repuls(i,2))			
			endif
			if (ans(17).ne.0.0) ans(26)=1.0
			
            	elseif (coltype(i).eq.15) then
!              	advance to get off square-shoulder this is the only result of this type of collision
	       	call bumpoff(i,j,evcode)
            	elseif (coltype(i).eq.16) then
!             	advance to get off square-well this is the only result of this type of collision
	       	call bumpoff(i,j,evcode)
!              	turning off the square-shoulder potential for the hb psuedobonds
                  	if ((i .le. nop1) .and. (j .le. nop1)) then
                  		if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.chnln1+1).and.(jj.ne.3*chnln1)) then
!                    		didn't involve an end bead
                     		if (identity(i).lt.(identity(j))) then
!                       			i is n, j is c
                        			call repuls_del_a(i,j)
                     		else
!                       			i is c, j is n
                        			call repuls_del_a(j,i)
                     		endif
                     		xpulse_del = .true.
                  		endif
                  	elseif ((i .le. nop1) .and. (j .gt. nop1)) then
                  		if ((ii.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    		didn't involve an end bead
                     		if (identity(i).lt.(identity(j))) then
!                       			i is n, j is c
                        			call repuls_del_a(i,j)
                     		else
!                       			i is c, j is n
                        			call repuls_del_a(j,i)
                     		endif
                     		xpulse_del = .true.
                  		endif
                  	elseif ((i .gt. nop1) .and. (j .gt. nop1)) then
                  		if ((ii.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+3*chnln2)) then
!                    		didn't involve an end bead
                     		if (identity(i).lt.(identity(j))) then
!                       			i is n, j is c
                        			call repuls_del_a(i,j)
                     		else
!                       			i is c, j is n
                        			call repuls_del_a(j,i)
                     		endif
                     		xpulse_del = .true.
                  		endif
                  	endif
					
            if (identity(i).lt.identity(j)) then
			ans(17) = dble(extra_repuls(i,1))
	        ans(18) = dble(extra_repuls(i,2))
	        ans(19) = dble(extra_repuls(j,1))
	        ans(20) = dble(extra_repuls(j,2))
			else
			ans(17) = dble(extra_repuls(j,1))
	        ans(18) = dble(extra_repuls(j,2))
	        ans(19) = dble(extra_repuls(i,1))
	        ans(20) = dble(extra_repuls(i,2))
			endif
     	    if (xpulse_del) ans(21) = 1.0 

            	endif
	
!!!!!HP+HB !!!!!
            ans(1) = dble(i)
			ans(2) = dble(j)
			ans(3) = tfalse 
			ans(4) = sv(1,i)
			ans(5) = sv(2,i)
			ans(6) = sv(3,i)
			ans(7) = sv(4,i)
			ans(8) = sv(5,i)
			ans(9) = sv(6,i)
			ans(10) = sv(1,j)
			ans(11) = sv(2,j)
			ans(12) = sv(3,j)
			ans(13) = sv(4,j)
			ans(14) = sv(5,j)
			ans(15) = sv(6,j)
			ans(16) = dble(coltype(i))	
            ans(28)	= dble(evcode)	

		!	write(fileout,*)'slave calculate normal coll:',i,j,tfalse,evcode,coltype(i),sv(1,i),sv(4,i),bm(i),bm(j)
         call MPI_SEND(ans,cols,MPI_DOUBLE_PRECISION,master,row,MPI_COMM_WORLD,ierr)
		 	    goto 90	
#ifdef canon
    	elseif (i == (noptotal) + 1) then
            ans(1) = dble(i)
			!tfalse = buffer(3)	
			ans(3) = buffer(3)
         call MPI_SEND(ans,cols,MPI_DOUBLE_PRECISION,master,row,MPI_COMM_WORLD,ierr)
		 	goto 90	
#endif  
      	elseif (i == (noptotal) + 2) then
            ans(1) = dble(i)
			!tfalse = buffer(3)	
			ans(3) = buffer(3)
			call MPI_SEND(ans,cols,MPI_DOUBLE_PRECISION,master,row,MPI_COMM_WORLD,ierr)
         	goto 90	
	!	endif
		elseif (i == (noptotal) + 3) then
	        ans(1) = dble(i)
			!tfalse = buffer(3)	
			ans(3) = buffer(3)
			call MPI_SEND(ans,cols,MPI_DOUBLE_PRECISION,master,row,MPI_COMM_WORLD,ierr)
	        goto 90	
		end if
		
	 end if
	     	write(fileout,*) 'terminate slave',status(MPI_TAG),MPI_TAG
	End if
!	Yiming: enddo for the major loop in serial code   
888		write(fileout,*) 'mpi finalized' 
    call MPI_FINALIZE(ierr)
    End	
#include "readinputs.f"
#include "allocatearrays.f"
#include "genconfig.f"
#include "bumped.f"      
#include "sqshlder.f"      
#include "check_nc_int.f"
#include "add_tbin.f"
#include "bond.f"
#include "cell_add.f"
#include "cell_link.f"
#include "check_sigma.f" 
#include "checkover.f"
#include "config.f"
#include "core.f"
#include "displ.f"
#include "energy.f"
#include "eventdyn.f"
#include "eventredo_up.f"
#include "eventredo_down.f"
#include "events.f"
#include "del_tbin.f"
#include "inputinfo.f"
#include "make_code.f"
#include "nbor.f"
#include "nbor_setup.f"
#include "nc_sqwel.f"
#include "partial_events.f"
#include "files_opn.f"
#include "write_rasmol-YM.f"
#include "repuls_add.f"
#include "repuls_check.f"
#include "repuls_check_3.f"
#include "repuls_del_a.f"
#include "repuls_del_b.f"
#include "scale_down.f"
#include "scale_up.f"
#include "sqwel.f"
#include "files_close.f"
!LR: I separated the analysis code out of the main code.
#include "radgyr.f"  
#include "end2end.f"

#ifdef equil
#include "equil.f"
#endif

#ifdef write_phipsi
#include "phipsi.f"
#endif
