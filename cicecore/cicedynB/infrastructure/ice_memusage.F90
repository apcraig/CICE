! Provides methods for querying memory use

MODULE ice_memory

!-------------------------------------------------------------------------------
! PURPOSE: memory use query methods
!    Should call ice_memory_init once before calling other interfaces
!-------------------------------------------------------------------------------

   use ice_kinds_mod, only : dbl_kind

   implicit none
   private
    
! PUBLIC: Public interfaces

   public ::  ice_memory_getusage, &
	      ice_memory_init, &
              ice_memory_print

   logical(log_kind), public :: memory_stats

! PRIVATE DATA:

   real(dbl_kind) :: mb_blk = 1.0_dbl_kind
   logical        :: initset = .false.

!===============================================================================

contains

!===============================================================================
! Initialize memory conversion to MB

subroutine ice_memory_init(iunit)

   implicit none

   !----- arguments -----

   integer, optional :: iunit   !< output unit number for optional writes
     
   !----- local -----

   ! --- Memory stats --- 
   integer :: msize                   ! memory size (high water)
   integer :: mrss                    ! resident size (current memory use)
   integer :: msize0,msize1           ! temporary size
   integer :: mrss0,mrss1,mrss2       ! temporary rss
   integer :: mshare,mtext,mdatastack
   integer :: ierr
 
   integer :: GPTLget_memusage

   real(dbl_kind),allocatable :: mem_tmp(:)
   character(*),parameter  :: subname = '(ice_memory_init)'
    
   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   ierr = GPTLget_memusage (msize, mrss0, mshare, mtext, mdatastack)
   allocate(mem_tmp(1024*1024))    ! 1 MWord, 8 MB
   mem_tmp = -1.0
   ierr = GPTLget_memusage (msize, mrss1, mshare, mtext, mdatastack)
   deallocate(mem_tmp)
   ierr = GPTLget_memusage (msize, mrss2, mshare, mtext, mdatastack)
   mb_blk = 1.0_dbl_kind
   if (mrss1 - mrss0 > 0) then
      mb_blk = (8.0_dbl_kind)/((mrss1-mrss0)*1.0_dbl_kind)
      initset = .true.
   endif

   if (present(iunit)) then
      write(iunit,'(A,l4)')    subname//' Initset conversion flag is ',initset
      write(iunit,'(A,f16.2)') subname//' 8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
      write(iunit,'(A,f16.2)') subname//' 8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
      write(iunit,'(A,f16.2)') subname//' Memory block size conversion in bytes is ',mb_blk*1024_dbl_kind*1024.0_dbl_kind
   endif

end subroutine ice_memory_init

!===============================================================================
! Determine memory use

subroutine ice_memory_getusage(r_msize,r_mrss)

   implicit none

   !----- arguments ---
   real(dbl_kind),intent(out) :: r_msize  !< memory usage value
   real(dbl_kind),intent(out) :: r_mrss   !< memory usage value

   !----- local ---
   integer :: msize,mrss
   integer :: mshare,mtext,mdatastack
   integer :: ierr
   integer :: GPTLget_memusage
   character(*),parameter  :: subname = '(ice_memory_getusage)'

   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   ierr = GPTLget_memusage (msize, mrss, mshare, mtext, mdatastack)
   r_msize = msize*mb_blk
   r_mrss  = mrss*mb_blk

end subroutine ice_memory_getusage

!===============================================================================
! Print memory use

subroutine ice_memory_print(iunit,string)

   implicit none

   !----- arguments ---
   integer, intent(in) :: iunit    !< unit number to write to
   character(len=*),optional, intent(in) :: string  !< optional string

   !----- local ---   
   real(dbl_kind)     :: msize,mrss
   character(len=128) :: lstring
   character(*),parameter  :: subname = '(ice_memory_print)'

   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   lstring = ' '
   if (present(string)) then
      lstring = string
   endif

   call ice_memory_getusage(msize,mrss)

   if (initset) then
      write(iunit,'(2a,2f14.4,1x,a)') subname,' memory use (MB) = ',msize,mrss,trim(lstring)
   else
      write(iunit,'(2a,2f14.4,1x,a)') subname,' memory use (??) = ',msize,mrss,trim(lstring)
   endif

end subroutine ice_memory_print

!===============================================================================

END MODULE ice_memory
