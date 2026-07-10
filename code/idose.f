!=======================================================================
!  Program: IDOSE version 0.1-alpha
!  Purpose: Estimate inhalation dose coefficients for selected
!           nuclide/age/intake-type with aerosol size adjustment
!=======================================================================
      program idose
      implicit none

!-------------------------------
! Declarations and parameters
!-------------------------------
      integer, parameter :: NREGION=9, NORG=35, NAGE=6, NTARGET=10, NDIAM=50, NNUCS=1252
      real(8), parameter :: CONVFACTOR = 3.1536D7
      
      integer, dimension(25) :: indices = [1,2,3,4,5,6,7,8,9,10,11,20,21,22,23,24,25,26,27,28,29,30,31,32,33] !CAP88 organ list

      ! Selection inputs
      character(len=10) :: XNUC
      character(len=10) :: CHEMFORM
      character(len=1)  :: DINTAKE, XINTAKE

      ! Names and labels
      character(len=9)  :: ORGANNAME(NORG)
      character(len=9)  :: TARGET(NTARGET)
      character(len=32) :: AGEGROUPNAME(NAGE)
      character(len=9)  :: NUCNAME(NNUCS)

      ! Weights, masses, factors
      real(8) :: REGIONCOEFFS(NREGION,NORG)
      real(8) :: SIZECOEFFS(NREGION)
      real(8) :: REGWEIGHTS(NORG)
      real(8) :: REMMASS(NORG)
      real(8) :: ORGCOEFFS(NORG)
      real(8) :: ADIAMFACTOR(NAGE,NTARGET,NDIAM)
      real(8) :: ARRAYWEIGHT(NORG,7) !7 age groups (2 adult groups)
      real(8) :: TEDE(3) !3 lung absoption types
      character(len=10) :: CTEMP(NORG)

      ! Other ancillary information
      real(8) :: QAERO(NNUCS)
      character(len=2) :: QINTAKE(NNUCS)
      character(len=3) :: REMQA
      integer :: QAGE(NNUCS)

      ! Dose calc variables
      real(8) :: HREM, HEFF
      real(8) :: SUMREG, MAXREGC
      real(8) :: SUMREM, SUMREMW, MAXREMC, SUMREM2, SUMREMW2
      real(8) :: WGONOV, WGONTE
      real(8) :: best, diff
      integer :: i, j, k, m, n, xref, idxov, idxte, idxmaxrem, number

      ! I/O helpers
      character(len=500) :: inline
      character(len=10)  :: nuc
      integer            :: xage, age, ios

!-------------------------------
! User selections and defaults
!-------------------------------

      CHEMFORM = 'INORGANIC'

!-------------------------------
! Read in index of nuclides
!-------------------------------

!          OPEN(9,FILE='INPUT.NDX',status='old',form='formatted')
          OPEN(9,FILE='allnucslist.txt',status='old',form='formatted')

          read(9, '(i)') NUMBER
           
          do i = 1, NUMBER 
          READ(9, '(A7,f9.2,A2,I5)') NUCNAME(i),qaero(i),qintake(i),qage(i)
          enddo !end loop over nuclides on this file

          close(9)

! Open output files

          open(8,file='output.complyinhalationdc.csv',status='unknown',form='formatted')
  412 FORMAT('Nuclide,Default_Type,Age,AerosolDiam,F,M,S,QAFLAG')
          write(8,412)

          open(11,file='output.cap88v4inhalationdc.csv',status='unknown',form='formatted')
  413 FORMAT('IDOSE v0.1-alpha : dose coefficients provided as Sv/Bq')
          write(11,413)

! Start loop over all nuclides from index list

      do n = 1, NUMBER

! Default intake type by nuclide (simple heuristic)

      XNUC = NUCNAME(n)

      if (XNUC(1:2) == 'Cs' .or. XNUC(1:2) == 'H-') then
       DINTAKE = 'F'
      elseif (XNUC(1:2) == 'Se' .or. XNUC(1:2) == 'I-') then
       DINTAKE = 'F'
      elseif (XNUC(1:2) == 'Th') then
       DINTAKE = 'S'
      else
       DINTAKE = 'M'
      end if

      ! Default intake type more complex situations
      if (XNUC(1:2) == 'Po' .and. CHEMFORM == 'ORGANIC') then
       DINTAKE = 'F'
      elseif (XNUC(1:2) == 'Po' .and. CHEMFORM /= 'INORGANIC') then
       DINTAKE = 'M'
      endif 

      if (XNUC(1:2) == 'Hg' .and. CHEMFORM == 'METHYL') then
       DINTAKE = 'M'
      endif


      REMQA = ""

      do m = 1, 3

       if(m == 1) XINTAKE = 'F'
       if(m == 2) XINTAKE = 'M'
       if(m == 3) XINTAKE = 'S'

!-------------------------------
! Initialize arrays to zero
!-------------------------------
      REGIONCOEFFS = 0.0D0
      SIZECOEFFS   = 0.0D0
      REGWEIGHTS   = 0.0D0
      REMMASS      = 0.0D0
      ORGCOEFFS    = 0.0D0
      ADIAMFACTOR  = 0.0D0
      ARRAYWEIGHT  = 0.0D0
      ORGANNAME    = ' '

!-------------------------------
! Read tissue weights (regular)
!-------------------------------
      open(unit=9, file='namelist_regular.txt', status='old', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then 
       write(*,*) 'Cannot open namelist_regular.txt'
       stop 
      endif

      read(9,'(A)', iostat=ios) inline
      if (ios == 0) read(inline,'(A9,7I10)', iostat=ios) nuc, (age, j=1,7)
      if (ios /= 0) then
       write(*,*) 'Header read error in namelist_regular.txt'
       stop
      endif

      do k = 1, 33
       read(9,'(A9,7F10.2)', iostat=ios) ORGANNAME(k), (ARRAYWEIGHT(k,j), j=1,7)
       if (ios /= 0) then
       write(*,*) 'Read error in namelist_regular.txt body'
       stop
       else 
       REGWEIGHTS(k) = ARRAYWEIGHT(k,7)
       endif
      end do
      close(9)

      ORGANNAME(34) = 'REM'
      ORGANNAME(35) = 'E50'
      REGWEIGHTS(34:35) = 0.0D0

      if(m==1 .AND. n==1) then
       write(11,313) (ORGANNAME(indices(i)),i=1,size(indices))
      endif

  313 FORMAT(25(A9),'E_50')

!-------------------------------
! Read remainder masses
!-------------------------------
      open(unit=9, file='namelist_remainder.txt', status='old', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then
      write(*,*) 'Cannot open namelist_remainder.txt'
      stop
      endif

      read(9,'(A)', iostat=ios) inline
      if (ios == 0) read(inline,'(A9,7I10)', iostat=ios) nuc, (age, j=1,7)
      if (ios /= 0) then
      write(*,*) 'Header read error in namelist_remainder.txt'
      stop 
      endif

      do k = 1, 33
       read(9,'(A9,7F10.2)', iostat=ios) inline, (ARRAYWEIGHT(k,j), j=1,7)
       if (ios /= 0) then 
        write(*,*) 'Read error in namelist_remainder.txt body'
        stop
       else
        REMMASS(k) = ARRAYWEIGHT(k,7)
       endif
      end do
      close(9)
      REMMASS(34:35) = 0.0D0

!-------------------------------
! Read aerosol size factors
!-------------------------------
      open(unit=9, file='DC_PAK3.DEP', status='old', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then
       write(*,*) 'Cannot open DC_PAK3.DEP'
       stop 
      endif

      do i = 1, NAGE
       read(9,'(A)', iostat=ios) AGEGROUPNAME(i)
       if (ios /= 0) then
        write(*,*) 'Read error: AGEGROUPNAME in DC_PAK3.DEP'
        stop
       endif

       do j = 1, NTARGET
        read(9,'(A9,50E10.2)', iostat=ios) TARGET(j), (ADIAMFACTOR(i,j,k), k=1, NDIAM)
        if (ios /= 0) then
         write(*,*) 'Read error: target lines in DC_PAK3.DEP'
         stop
        endif
       end do
      end do
      close(9)

!-------------------------------
! Select nearest aerosol bin
!-------------------------------
      xref = -1
      best = 1.0D30
      do k = 1, NDIAM
       diff = dabs(ADIAMFACTOR(6,1,k) - QAERO(n)) !QAERO is the user specified aerosol diameter
       if (diff < best) then
         best = diff
         xref = k
       end if
      end do
      if (xref < 1) stop 'No matching aerosol diameter bin found'
!      write(*,'(A,I4,2A,1PE12.4)') 'Aerosol bin: ', xref, '  AD= ', '', ADIAMFACTOR(6,1,xref)

!-------------------------------
! Build size factors for 9 regions
!-------------------------------
      do k = 1, NREGION
       SIZECOEFFS(k) = ADIAMFACTOR(6, k+1, xref)  ! adult age index=6, targets 2..10
      end do

!-------------------------------
! Read 9 region HDB files
!-------------------------------
       call read_region_hdb('AI.HDB',      1, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBE-GEL.HDB', 2, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBE-SOL.HDB', 3, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBE-SEQ.HDB', 4, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBI-GEL.HDB', 5, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBI-SOL.HDB', 6, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('BBI-SEQ.HDB', 7, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('ET1.HDB',     8, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)
       call read_region_hdb('ET2.HDB',     9, XNUC, QAGE(n), XINTAKE, REGIONCOEFFS)

!-------------------------------
! Aggregate per-organ coefficients
!-------------------------------
      do i = 1, NORG
       ORGCOEFFS(i) = 0.0D0
        do k = 1, NREGION
         ORGCOEFFS(i) = ORGCOEFFS(i) + REGIONCOEFFS(k,i) * SIZECOEFFS(k)
        end do
      end do

!-------------------------------
! Gonad weighting selection
!-------------------------------
      idxov = 21
      idxte = 26
      WGONOV = 0.0D0
      WGONTE = 0.0D0
      if (ORGCOEFFS(idxov) >= ORGCOEFFS(idxte)) then
       WGONOV = 0.2D0
      else
       WGONTE = 0.2D0
      end if

!-------------------------------
! Sum regular component
!-------------------------------
      SUMREG = 0.0D0
      do i = 1, NORG
       if (i == idxov) then
        SUMREG = SUMREG + ORGCOEFFS(i) * WGONOV
       elseif (i == idxte) then
        SUMREG = SUMREG + ORGCOEFFS(i) * WGONTE
       else
        SUMREG = SUMREG + ORGCOEFFS(i) * REGWEIGHTS(i)
       end if
      end do

! Max regular coefficient (for exceptional remainder case)
      MAXREGC = 0.0D0
      do i = 1, NORG
        if (REGWEIGHTS(i) > 0.0D0) then
        if (ORGCOEFFS(i) > MAXREGC) MAXREGC = ORGCOEFFS(i)
        end if
      end do
      if (WGONOV > 0.0D0) then
       if (ORGCOEFFS(idxov) > MAXREGC) MAXREGC = ORGCOEFFS(idxov)
      end if
       if (WGONTE > 0.0D0) then
       if (ORGCOEFFS(idxte) > MAXREGC) MAXREGC = ORGCOEFFS(idxte)
      end if

!-------------------------------
! Remainder dose (HREM)
!-------------------------------
      SUMREM   = 0.0D0
      SUMREMW  = 0.0D0
      MAXREMC  = 0.0D0
      idxmaxrem = 0

      do i = 1, NORG
       if (REMMASS(i) > 0.0D0) then
        SUMREM  = SUMREM  + ORGCOEFFS(i) * REMMASS(i)
        SUMREMW = SUMREMW + REMMASS(i)
        if (ORGCOEFFS(i) > MAXREMC) then
        MAXREMC   = ORGCOEFFS(i)
        idxmaxrem = i
        end if
       end if
      end do

      if (SUMREMW <= 0.0D0) then !this never happens with the current list of tissue weights
        HREM = 0.0D0
        REMQA = trim(REMQA) // "Z"
      else
       if (MAXREGC > MAXREMC) then !max regular tissue coef < max remainder tissue coef
        HREM = SUMREM / SUMREMW
        REMQA = trim(REMQA) // "A" 
       else !max remainder tissue coef > max regular tissue coef
        SUMREM2  = 0.0D0
        SUMREMW2 = 0.0D0
        do i = 1, NORG
         if (REMMASS(i) > 0.0D0 .and. i /= idxmaxrem) then
           SUMREM2  = SUMREM2  + ORGCOEFFS(i) * REMMASS(i)
           SUMREMW2 = SUMREMW2 + REMMASS(i)
         end if
        end do
         if (SUMREMW2 > 0.0D0) then !remainder tissue weights excluding the max value > 0
          HREM = 0.5D0 * ((SUMREM2 / SUMREMW2) + MAXREMC)
          REMQA = trim(REMQA) // "B"
         else !remainder tissue weights excluding max value <= 0.0 which would never happen with the current list
          HREM = MAXREMC
          REMQA = trim(REMQA) // "C"
         end if
       end if
      end if

!-------------------------------
! Effective dose
!-------------------------------

      HEFF = SUMREG + 0.05D0 * HREM
      ORGCOEFFS(34) = HREM
      ORGCOEFFS(35) = HEFF
      TEDE(m) = HEFF

!-------------------------------
! Output results
!-------------------------------
!      write(*,*) 
!      write(*,*) 'Per-organ inhalation dose coefficients'
!      write(*,*) 'Index,Organ,Dose (Sv/Bq), Dose (mrem/pCi)'
!      do i = 1, NORG
!       write(*,'(I2,2X,A9,2X,1PE12.4,2X,1PE12.4)') i, ORGANNAME(i), ORGCOEFFS(i), ORGCOEFFS(i)*3700.0D0

       if(XINTAKE == QINTAKE(n)) then !output CAP88 format 
       write(11,314) NUCNAME(n),QINTAKE(n),QAGE(n),QAERO(n)
       do i = 1, NORG
       write(CTEMP(i), '(1PE10.3)') ORGCOEFFS(i)
       write(CTEMP(i), '(A)') ADJUSTL(CTEMP(i))
       enddo
       write(11,315) (CTEMP(indices(i)),i=1,size(indices)),HEFF
       endif

  314 FORMAT(A9,',',A2,',',I4,',',F10.2)
  315 FORMAT(25(A10),1PE10.3)

       enddo !lung absorption types

       !output that COMPLY will need (adult age, 1 um aerosol)
       write(*,'(A9,2X,A2,2X,I4,2X,F10.2,2X,3(2X,1PE12.4),1X,A3)') NUCNAME(n),DINTAKE,QAGE(n),QAERO(n),(TEDE(i),i=1,3),REMQA
       write(8,414) NUCNAME(n),DINTAKE,QAGE(n),QAERO(n),(TEDE(i),i=1,3),REMQA
       enddo !NNUCS

  414 FORMAT(A9,',',A2,',',I4,',',F10.2,3(',',1PE12.4),',',A3)

      contains

!-------------------------------
! Convert a string to uppercase
!-------------------------------
      subroutine upcase(str)
       character(len=*), intent(inout) :: str
       integer :: i, ia
      do i = 1, len_trim(str)
       ia = iachar(str(i:i))
       if (ia >= iachar('a') .and. ia <= iachar('z')) str(i:i) = achar(ia-32)
      end do
      end subroutine upcase

!---------------------------------------------------------
! Read one region HDB file and fill REGIONCOEFFS(kreg,:)
! Match on (NUC, AGE, TYPE). Stops after first match.
!---------------------------------------------------------
      subroutine read_region_hdb(fname, kreg, XNUC, XAGE, XINTAKE, regcoef)
      character(len=*), intent(in)  :: fname
      integer,          intent(in)  :: kreg
      character(len=*), intent(in)  :: XNUC, XINTAKE
      integer,          intent(in)  :: XAGE
      real(8),          intent(inout) :: regcoef(NREGION,NORG)

      integer :: lu, ios, j, age
      character(len=10) :: nuc
      character(len=1)  :: itype
      real(8) :: aerodiam, tmp(NORG)
      character(len=256) :: line

      lu = 20 + kreg
      open(unit=lu, file=fname, status='old', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then
       write(*,*) 'Warning: cannot open ', trim(fname)
       return
      end if

      ! skip two header lines
      read(lu,'(A)', iostat=ios) line
      read(lu,'(A)', iostat=ios) line

      do
       read(lu,'(A7,I5,A2,E8.1,35E10.2)', iostat=ios) nuc, age, itype, aerodiam, (tmp(j), j=1, NORG)
       if (ios /= 0) exit
       if (age > 7300) age = 7300
       if (nuc == XNUC .and. age == XAGE .and. itype == XINTAKE) then
         regcoef(kreg,1:NORG) = tmp(1:NORG)
         exit
       end if
       end do

       close(lu)
       end subroutine read_region_hdb

       end program idose

