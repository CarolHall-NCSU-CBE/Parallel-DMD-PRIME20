! ============================================================================
! Subroutine allocatearrays defines allocatable arrays corresponding to each system 	  
! Last modified on 10/12/2023 by Van Nguyen
! ============================================================================

	subroutine allocatearrays()
   
#include "def.h"

	use global
	use inputreadin
 
	implicit none
	
	allocate (old_rx(noptotal),old_ry(noptotal),old_rz(noptotal))
	allocate (hp1(numbeads1),hp2(numbeads2))
	allocate (npt(noptotal+1),nb(maxnbs*noptotal+1),na_npt(noptotal))
	allocate (npt_dn(noptotal+1),dnnab(maxnbs*noptotal+1),nnabdn(noptotal))
	
      	return

      	end




