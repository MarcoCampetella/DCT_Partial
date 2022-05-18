MODULE globale
! atomlist = list of selected atoms ; natom = number of atoms of the list 
INTEGER, ALLOCATABLE :: atomlist(:)
INTEGER              :: natom
! atom_num_vec() = vector of atomic number ; shell_type() = type of each shell
! npXshell() = number of primitive for each shell ; shell_atom_map() = owing of each shell
! coordinates(:,:) = n*3 coordinates matrix ; prim_exp() = exp of each primitive
! prim_coeffic() = coeff of each primitive ; prim_coeffic_sp = coeff of each p primitive of type SP
! shell_center() = center of each shell ; nshellXatom = number of shell for each atom
! nexpXatom() = number of exponent for each atom ; nbasisXatom = number of atomic basis function x atom
! selectedAO() = vector of the selected atomic orbital ; selectedEXP() = vector of the selected EXP
! selectedEXP() = selectedCOEFF() ; selectedSHELL() = vector of the selected SHELL
! mapAO(,) = map containing the AO label, first and last exp index, the type and the related atom
! for each AO
! densGS(:,:) = density matrix coefficients of GS , densEX(:,:) = density matrix coefficients of EX
INTEGER, ALLOCATABLE            :: atom_num_vec(:),shell_type(:),npXshell(:)
INTEGER, ALLOCATABLE            :: shell_atom_map(:),nshellXatom(:),nexpXatom(:)
INTEGER, ALLOCATABLE            :: nbasisXatom(:),selectedAO(:),selectedEXP(:)
INTEGER, ALLOCATABLE            :: selectedSHELL(:),mapAO(:,:)
REAL*8, ALLOCATABLE             :: coordinates(:,:),prim_exp(:),mo_coeff(:)
REAL*8, ALLOCATABLE             :: prim_coeffic(:),prim_coeffic_sp(:),shell_center(:)
REAL*8, ALLOCATABLE             :: gs_dens(:),ci_dens(:),ci_udens(:),densGS(:,:),densEX(:,:)
REAL*8, ALLOCATABLE             :: densUEX(:,:) 
CHARACTER(LEN=10), ALLOCATABLE  :: AOtype(:)
INTEGER                         :: natomtot,nshell,nexp,nbasis,nbasisi,nbasisSelected,nEXPSelected 
INTEGER                         :: nSHELLselected,nbasis2,nbasisD
! grid(:,:,:) = it's the grid point where the density is evaluated
INTEGER                         :: nptot,npgrid,npextx,npexty,npextz
REAL*8, ALLOCATABLE             :: grid(:,:,:,:),x(:),y(:),z(:)
REAL*8, DIMENSION(3)            :: origin
! densSELGS(:) = GS selected density ; densSELEX(:) = EX selected density 
REAL*8, ALLOCATABLE             :: densSELGS(:),densSELEX(:),densUSELEX(:),matrixAO(:,:)
END MODULE globale


PROGRAM SelectedDct
USE globale
IMPLICIT NONE
! filein = .fchk file to parse; selection = list of atoms
INTEGER             :: IOut,istatus,i,j,k,tempcount,gridpoint
INTEGER             :: AllocateStatus,STAT,t1,t2,clock_rate,clock_max
CHARACTER(LEN=100)  :: filein,selection,cubext
LOGICAL             :: verbosity,parallel,debug,gridext,selectionV
REAL*8              :: DX,DY,DZ,gridfactor 

1001 FORMAT(/,' The input file is                    : ',(A),/)
1011 FORMAT(' The Grid is extracted from           : ',(A),/)
1002 FORMAT(' The number of Grid Points are        : ',I10,/)
1003 FORMAT(' The Total Number of Atoms is         : ',I10,/)
1004 FORMAT(' The Total Number of A.O. is          : ',I10,/)
1005 FORMAT(' The Total Number of Exp is           : ',I10,/)
1006 FORMAT(' The Total Number of Shell is         : ',I10,/)
1007 FORMAT(' The Number of Selected Atoms is      : ',I10,/)
1008 FORMAT(' The Number of Selected A.O. is       : ',I10,/)
1009 FORMAT(' The Number of Selected Exp is        : ',I10,/)
1010 FORMAT(' The Number of Selected Shell is      : ',I10,/)

INTERFACE
   SUBROUTINE Parselist(selection)
      CHARACTER(LEN=50), INTENT(IN) :: selection
   END SUBROUTINE Parselist
END INTERFACE

INTERFACE
   SUBROUTINE readfchk(filein,IOut,verbosity,debug,selectionV)
      CHARACTER(LEN=50), INTENT(IN) :: filein
      INTEGER                       :: IOut
      LOGICAL                       :: verbosity,debug,selectionV
   END SUBROUTINE readfchk
END INTERFACE

INTERFACE
   SUBROUTINE MakeGrid(gridfactor,verbosity,DX,DY,DZ)
      LOGICAL                       :: verbosity
      REAL*8                        :: DX,DY,DZ,gridfactor
   END SUBROUTINE MakeGrid
END INTERFACE

INTERFACE
   SUBROUTINE GridFromExt(verbosity,DX,DY,DZ,cubext)
      LOGICAL                       :: verbosity
      REAL*8                        :: DX,DY,DZ
      CHARACTER(LEN=100)            :: cubext
   END SUBROUTINE GridFromExt
END INTERFACE

INTERFACE
   SUBROUTINE calcDensPar(IOut,verbosity,debug,DX,DY,DZ)
      INTEGER                       :: IOut
      LOGICAL                       :: verbosity,debug
      REAL*8                        :: DX,DY,DZ 
   END SUBROUTINE calcDensPar
END INTERFACE

INTERFACE
   SUBROUTINE calcDensParNew(IOut,verbosity,debug,DX,DY,DZ)
      INTEGER                       :: IOut
      LOGICAL                       :: verbosity,debug
      REAL*8                        :: DX,DY,DZ 
   END SUBROUTINE calcDensParNew
END INTERFACE

INTERFACE
   SUBROUTINE calcDens(IOut,verbosity,debug,DX,DY,DZ)
      INTEGER                       :: IOut
      LOGICAL                       :: verbosity,debug
      REAL*8                        :: DX,DY,DZ 
   END SUBROUTINE calcDens
END INTERFACE

INTERFACE
   SUBROUTINE calcDensNew(IOut,verbosity,debug,DX,DY,DZ)
      INTEGER                       :: IOut
      LOGICAL                       :: verbosity,debug
      REAL*8                        :: DX,DY,DZ 
   END SUBROUTINE calcDensNew
END INTERFACE

CALL SYSTEM_CLOCK ( t1, clock_rate, clock_max )

CALL RdOpts(filein,selection,gridfactor,verbosity,debug,parallel,cubext,& 
            gridext,selectionV)  ! RdOpts reads the options from keyboard
IOut = 750
OPEN (UNIT=750, FILE='resume.dat')
CALL PrtHdr(IOut)  ! write the banner of the program

!IF (selectionL .EQV. .TRUE.) THEN
!END IF
IF (selectionV) THEN
   CALL Parselist(selection) ! count the selectioned atoms and expands the selection in a vector
END IF

!WRITE(*,*) atomlist
!WRITE(*,*) natomtot
!WRITE(*,*) TRIM(filein)

CALL readfchk(filein,IOut,verbosity,debug,selectionV) ! read the data from the .fchk file and 

DO i = 1,natom
   IF (atomlist(i) > natomtot) THEN
      WRITE(IOut,*) "SELECTION LIST WRONG : SELECTION OUT OF RANGE"
      STOP
   END IF
END DO

! create the matrix used for the subsequently calculations 
IF (gridext) THEN
   CALL GridFromExt(verbosity,DX,DY,DZ,cubext)
ELSE
   CALL MakeGrid(gridfactor,verbosity,DX,DY,DZ)
END IF

WRITE(IOut,1001) TRIM(filein)
IF (gridext) THEN
   WRITE(IOut,1011) TRIM(cubext)
END IF
WRITE(IOut,1002) nptot
WRITE(IOut,1003) natomtot 
WRITE(IOut,1004) nbasis 
WRITE(IOut,1005) nexp 
WRITE(IOut,1006) nshell 
IF (verbosity) THEN
   WRITE(IOut,1007) natom    
   WRITE(IOut,1008) nbasisSelected
   WRITE(IOut,1009) nEXPSelected
   WRITE(IOut,1010) nSHELLselected
   WRITE(Iout,*) "The Selected Atoms are               : "
   WRITE(Iout,*) 
   WRITE(Iout,'(33X,8I5)') (atomlist(i), i=1,natom)
   WRITE(Iout,*)
   WRITE(Iout,*) "The Selected AO are                  : "
   WRITE(Iout,*) 
   WRITE(Iout,'(33X,8I5)') (selectedAO(i), i=1,nbasisSelected)
   WRITE(Iout,*)
   WRITE(Iout,*) "The Selected EXP are                 : "
   WRITE(Iout,*) 
   WRITE(Iout,'(33X,8I5)') (selectedEXP(i), i=1,nEXPSelected)
   WRITE(Iout,*)
   WRITE(Iout,*) "The Selected SHELL are               : "
   WRITE(Iout,*) 
   WRITE(Iout,'(33X,8I5)') (selectedSHELL(i), i=1,nSHELLselected)
   WRITE(Iout,*) 
END IF
!WRITE(Iout,*)
!WRITE(*,*) npgrid


IF (parallel) THEN
   CALL calcDensParNew(IOut,verbosity,debug,DX,DY,DZ)
ELSE
   CALL calcDensNew(IOut,verbosity,debug,DX,DY,DZ)
END IF

CALL calcDCT(DX,DY,DZ,IOut)

IF (verbosity) THEN
   IF (gridext) THEN
      CALL writeDensEXT(DX,DY,DZ)
   ELSE
      CALL writeDens(DX,DY,DZ)
   END IF
END IF

CALL SYSTEM_CLOCK ( t2, clock_rate, clock_max )

write (IOut,'(A,F12.4,A)') ' Elapsed time                           : ', REAL(t2-t1) / (REAL(clock_rate)*60),' min'

CLOSE(750)




END PROGRAM SelectedDct
