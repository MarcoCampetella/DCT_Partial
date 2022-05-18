!--------------------------------------------------------------------
SUBROUTINE PrtHdr(IOut)
INTEGER IOut
!
! Print header
!
1000 format(/,1x,45('-'),/,21x,'DctSelected',20x,/,1x,45('-'),/ &
            '                                     CNRS tool',/,&
            ' M.Campetella                          ver 1.0',/,&
            1x,45('-'),/)
  Write(IOut,1000)
  Return
END SUBROUTINE
!--------------------------------------------------------------------

SUBROUTINE RdOpts(filein,selection,gridfactor,verbosity,debug,parallel,cubext,gridext,selectionV)
USE globale 
!
!
INTEGER                                      :: IGet,IArg,num_args,istatus
REAL*8                                       :: gridfactor  
CHARACTER(LEN=100)                           :: filein,selection,cubext
CHARACTER(LEN=50), DIMENSION(:), ALLOCATABLE :: args
CHARACTER(LEN=6)                             :: npgridCHAR,gridfactorCHAR
LOGICAL                                      :: verbosity,parallel,debug,gridext,selectionV


1010 FORMAT(' Error file missing : ','"'(A)'"')
1020 FORMAT(' The input file is  : ','"'(A)'"')
1000 FORMAT(' Error in input stream: nothing followed option ','"'(A)'"')
1001 FORMAT(' Error in input stream: option ',(A),' unknown')
1009 FORMAT(' The total number of arguments are ',I5,/)
1200 FORMAT('                                                                                     ',/, &
            ' Use the following options to run the program:                                       ',/, &
            '   -h      (optional) ... Get this help message                                      ',/, &
            '   -grid   (optional) ... Number of Grid Points ; DEF. = 80*80*80                    ',/, &
            '   -rescal (optional) ... Rescal factor for the Box Edge ; DEF. = 2                  ',/, &
            '   -verb   (optional) ... Increases the verbosity ; DEF. = .FALSE.                   ',/, &
            '   -db     (optional) ... Debug mode ; DEF. = .FALSE.                                ',/, &
            '   -ext    (optional) ... External Cube File ; DEF. = .FALSE.                        ',/, &
            '   -par    (optional) ... Active the OMP protocol ; DEF. = .FALSE.                   ',/, &
            '   -sel    (optional) ... atom selection, syntax : ...1-12,25..., DEF.: all          ',/, & 
            '   -f      (required) ... .fchk file name                                            ',/ )

  num_args = command_argument_count()
  ALLOCATE(args(num_args))
  IGet         = 0
  filein       = ''
  npgrid       = 80
  gridfactor   = 2
  verbosity    = .FALSE.
  debug        = .FALSE.
  parallel     = .FALSE.
  gridext      = .FALSE.
  selectionV   = .FALSE.
!  NoSelection  = .TRUE.
!  selectionL   = .FALSE.
!  fileout      = 'last'

  IF (num_args.lt.2) THEN
    WRITE(*,1200)
    STOP
  END IF

  DO IArg = 1, num_args
    CALL get_command_argument(IArg,args(IArg))
  END DO
!  write(Iout,1009) num_args

  DO IArg = 1, num_args
     IF (IGet.NE.0) THEN
        IF (IGet.EQ.1) THEN
           filein = args(IArg)
        ELSEIF (IGet.EQ.2) THEN
           selection = args(IArg)
        ELSEIF (IGet.EQ.3) THEN
           npgridCHAR = args(IArg)
           READ(npgridCHAR,*) npgrid
        ELSEIF (IGet.EQ.4) THEN
           gridfactorCHAR = args(IArg)
           READ(gridfactorCHAR,*) gridfactor
        ELSEIF (IGet.EQ.5) THEN
           cubext  = args(IArg)
           gridext = .TRUE.
        ENDIF
        IGet = 0
     ELSEIF (args(IArg).EQ.'-h') THEN
        WRITE(IOut,1200)
        STOP
     ELSEIF (args(IArg).EQ.'-f') THEN
        IGet = 1
        IF (IArg.EQ.num_args) GOTO 800
     ELSEIF (args(IArg).EQ.'-sel') THEN
        IGet = 2
        selectionV = .TRUE.
        IF (IArg.EQ.num_args) GOTO 800
     ELSEIF (args(IArg).EQ.'-grid') THEN
        IGet = 3
        IF (IArg.EQ.num_args) GOTO 800
     ELSEIF (args(IArg).EQ.'-rescal') THEN
        IGet = 4
        IF (IArg.EQ.num_args) GOTO 800
     ELSEIF (args(IArg).EQ.'-ext') THEN
        IGet = 5
        IF (IArg.EQ.num_args) GOTO 800
     ELSEIF (args(IArg).EQ.'-verb') THEN
        verbosity = .TRUE.
     ELSEIF (args(IArg).EQ.'-db') THEN
        debug = .TRUE.
     ELSEIF (args(IArg).EQ.'-par') THEN
        parallel = .TRUE.
     ELSE
        GOTO 801
     ENDIF
  ENDDO

  GOTO 900

800 WRITE(*,1000) TRIM(args(IArg))
  STOP

801 WRITE(*,1001) args(IArg)
  WRITE(*,*) IArg
  STOP

805 WRITE(*,*) 'wrong integer data'
  STOP

900 CONTINUE

!!!!!!!!!!!!! CHECK THE INPUT FILE !!!!!!!!!!!!!!!!
OPEN (UNIT=9, FILE=filein, ACTION='READ',STATUS='OLD', iostat=istatus)
IF (istatus.NE.0) THEN
   WRITE(*,1010) TRIM(filein)
   STOP
END IF
CLOSE(9)
IF (gridext) THEN
   OPEN (UNIT=9, FILE=cubext, ACTION='READ',STATUS='OLD', iostat=istatus)
   IF (istatus.NE.0) THEN
      WRITE(*,1010) TRIM(cubext)
      STOP
   END IF
   CLOSE(9)
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RETURN

End Subroutine RdOpts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Parselist(list)
USE globale
IMPLICIT NONE
INTEGER :: i,j,k,countpartial,indice,atom1int,atom2int,AllocateStatus,STAT
CHARACTER(LEN=50) :: list
CHARACTER(LEN=6) :: partial
! partialMatrix CONTAINS EVERY FRAGMENT OF THE LIST
CHARACTER(LEN=20), ALLOCATABLE :: partialMatrix(:)
CHARACTER(LEN=4) :: atom1string,atom2string

countpartial=0

DO i=1,LEN_TRIM(list) ! COUNT THE NUMBER OF COMMA IN LIST
   IF (list(i:i).EQ.',') THEN
      countpartial=countpartial+1
   END IF
END DO
! EVERY ELEMENT OF partialMatrix ARE THE ELEMENTS OF list SEPARATED BY ,
ALLOCATE(partialMatrix(countpartial+1))

IF (SIZE(partialMatrix).EQ.1) THEN
   partialMatrix = list
ELSE
   indice=1 ! TAKE INTO ACCOUNT THE RUNNING OVER THE LIST
   j=1 ! COUNT THE NUMBER OF partialMatrix ELEMENT
   DO i=1,LEN_TRIM(list)
      IF (list(i:i).EQ.',') THEN
         partialMatrix(j) = TRIM(list(indice:i-1))
         indice=i+1
         j=j+1
         IF (j.EQ.SIZE(partialMatrix)) THEN
            partialMatrix(j) = TRIM(list(indice:LEN_TRIM(list)))
         END IF
      END IF
   END DO
END IF

!WRITE(*,*) partialMatrix
! COUNT THE TOTAL NUMBER OF ATOMS
j=0
DO i=1,SIZE(partialMatrix)
   IF (SCAN(partialMatrix(i),"-").EQ.0) THEN
      j=j+1
   ELSE
      READ(partialMatrix(i)(1:(SCAN(partialMatrix(i),"-")-1)),*) atom1int 
      READ(partialMatrix(i)((SCAN(partialMatrix(i),"-")+1):LEN(partialMatrix(i))),*) atom2int
      j=j+ABS((atom2int-atom1int)+1)
   END IF
END DO
natom = j

ALLOCATE(atomlist(natom),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

! DEFINES THE ATOMLIST VECTOR
j=0
DO i=1,SIZE(partialMatrix)
   IF (SCAN(partialMatrix(i),"-").EQ.0) THEN
      j=j+1
      READ(partialMatrix(i),*) atomlist(j)
   ELSE
      READ(partialMatrix(i)(1:(SCAN(partialMatrix(i),"-")-1)),*) atom1int
      READ(partialMatrix(i)((SCAN(partialMatrix(i),"-")+1):LEN(partialMatrix(i))),*) atom2int
      DO k = atom1int,atom2int
         j=j+1
         atomlist(j) = k
      END DO
   END IF
END DO

END SUBROUTINE Parselist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MakeGrid(gridfactor,verbosity,deltax,deltay,deltaz)
USE globale

INTEGER                                 :: i,j,k,l,tempcount,AllocateStatus,tempCub
REAL*8, DIMENSION(3)                    :: com,coord_max,r12
REAL*8, DIMENSION(npgrid,npgrid,npgrid) :: cube
REAL*8                                  :: r,temp_r,delta,xmin,xmax,ymin,ymax,zmin,zmax
REAL*8                                  :: deltax,deltay,deltaz,rx,ry,rz,tempDelta,gridfactor
REAL*8, DIMENSION(natomtot)             :: xg,yg,zg,charge
LOGICAL                                 :: verbosity
REAL*8, PARAMETER                       :: au2ang = 0.52917

nptot   = npgrid**3

!WRITE(*,*) npgrid
!WRITE(*,*) nptot

ALLOCATE(grid(npgrid,npgrid,npgrid,3),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(x(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(y(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(z(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

! calculation of the center of mass, it will be the new origin
com = 0.
DO i = 1,natomtot
   com(:) = com(:) + coordinates(i,:)*atom_num_vec(i)
END DO
com = com/SUM(atom_num_vec)

xg(:) = coordinates(:,1)
yg(:) = coordinates(:,2)
zg(:) = coordinates(:,3)
xmin = MINVAL(xg)
ymin = MINVAL(yg)
zmin = MINVAL(zg)
xmax = MAXVAL(xg)
ymax = MAXVAL(yg)
zmax = MAXVAL(zg)
rx = ( xmax - xmin ) 
IF ( rx < 4. ) rx = 6.*(gridfactor/2)
ry = ( ymax - ymin ) 
IF ( ry < 4. ) ry = 6.*(gridfactor/2)
rz = ( zmax - zmin ) 
IF ( rz < 4. ) rz = 6.*(gridfactor/2)
deltax = gridfactor*2.0 * rx / npgrid
!WRITE(*,*) gridfactor , rx , npgrid
!WRITE(*,*) gridfactor*2.0 * rx / npgrid
deltay = gridfactor*2.0 * ry / npgrid
deltaz = gridfactor*2.0 * rz / npgrid
!WRITE(*,'(9F8.4)') xmin , xmax , deltax, ymin , ymax , deltay, zmin , zmax , deltaz
!WRITE(*,*)
!WRITE(*,'(9F8.4)') xmin- (rx/2) , xmin- (rx/2) + npgrid*deltax , rx, ymin- (rz/2) , ymin- (rz/2)+ npgrid*deltay , &
!                   ry, zmin - (rz/2), zmin - (rz/2) + npgrid*deltaz , rz
!WRITE(*,*) coordinates(1,:)
!DO i = 1,natomtot
   
!   r12 = coordinates(i,:) - com(:)
!   temp_r = dsqrt(r12(1)**2.0+r12(2)**2.0+r12(3)**2.0)
!   IF (temp_r.GT.r) THEN
!      r = temp_r
!      coord_max(:) = coordinates(i,:)
!   END IF 
!END DO

! DEFINE THE GRID
DO i = 1,npgrid
   DO j = 1,npgrid
      DO k = 1,npgrid
         grid(k,j,i,1) = com(1) - (gridfactor*rx) + (i-1)*deltax 
         grid(k,j,i,2) = com(2) - (gridfactor*ry) + (j-1)*deltay 
         grid(k,j,i,3) = com(3) - (gridfactor*rz) + (k-1)*deltaz 
      END DO
   END DO
END DO

!DO i = 1,npgrid
!   DO j = 1,npgrid
!         WRITE(*,'(10F12.4)') (grid(k,j,i,1), k=1,npgrid) 
!   END DO
!END DO
!WRITE(*,*) npgrid,nptot
!x = 1.
!WRITE(*,'(6F8.4)') (x(i), i = 1,nptot)

!WRITE(*,'(6F10.4)') (x(i), i=1,ntot)
tempcount = 1
DO i = 1,npgrid
   DO j = 1,npgrid
      DO k = 1,npgrid
         x(tempcount) = grid(k,j,i,1)
         y(tempcount) = grid(k,j,i,2)
         z(tempcount) = grid(k,j,i,3)
         tempcount    = tempcount + 1
      END DO
   END DO
END DO
!WRITE(*,'(6F8.4)') (x(i), i = 1,nptot)
! IF VERBOSITY IS ON, WRITE THE GRID
IF (verbosity) THEN
   OPEN(UNIT=20,FILE='grid.xyz',ACTION='WRITE')
   WRITE(20,*) (npgrid)**3
   WRITE(20,*)
   DO i = 1,npgrid
      DO j = 1,npgrid
         DO k = 1,npgrid
            WRITE(20,'(A6,3F18.10)') "H",(grid(k,j,i,l)*au2ang,l=1,3) 
         END DO
      END DO
   END DO
   CLOSE(20)
END IF

END SUBROUTINE MakeGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GridFromExt(verbosity,deltax,deltay,deltaz,filein)
USE globale

INTEGER                                 :: i,j,k,l,tempcount,AllocateStatus,tempCub,reason,tempint
INTEGER                                 :: istatus
REAL*8                                  :: deltax,deltay,deltaz,rx,ry,rz,tempDelta,tempreal
REAL*8, DIMENSION(natomtot)             :: xg,yg,zg,charge
LOGICAL                                 :: verbosity
CHARACTER(LEN=100)                      :: filein,linea
REAL*8, PARAMETER                       :: au2ang = 0.52917

100 format (A)

OPEN(UNIT=10,FILE=trim(filein),STATUS='OLD',ACTION='READ',IOSTAT= istatus)
   READ(10,*)
   READ(10,*)
   READ(10,*) tempreal,origin(1),origin(2),origin(3),tempreal
   READ(10,*) npextx,deltax,tempreal,tempreal
   READ(10,*) npexty,tempreal,deltay,tempreal
   READ(10,*) npextz,tempreal,tempreal,deltaz
CLOSE(10)

!WRITE(*,*) npextx,npexty,npextz,npextx*npexty*npextz
nptot   = npextx*npexty*npextz

!WRITE(*,*) npgrid
!WRITE(*,*) nptot

ALLOCATE(grid(npextz,npexty,npextx,3),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(x(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(y(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
ALLOCATE(z(nptot),STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"


! DEFINE THE GRID
DO i = 1,npextx
   DO j = 1,npexty
      DO k = 1,npextz
         grid(k,j,i,1) = origin(1) + (i-1)*deltax 
         grid(k,j,i,2) = origin(2) + (j-1)*deltay 
         grid(k,j,i,3) = origin(3) + (k-1)*deltaz 
      END DO
   END DO
END DO

tempcount = 1
DO i = 1,npextx
   DO j = 1,npexty
      DO k = 1,npextz
         x(tempcount) = grid(k,j,i,1)
         y(tempcount) = grid(k,j,i,2)
         z(tempcount) = grid(k,j,i,3)
         tempcount    = tempcount + 1
      END DO
   END DO
END DO
!WRITE(*,'(6F8.4)') (x(i), i = 1,nptot)
! IF VERBOSITY IS ON, WRITE THE GRID
IF (verbosity) THEN
   OPEN(UNIT=20,FILE='grid.xyz',ACTION='WRITE')
   WRITE(20,*) nptot
   WRITE(20,*)
   DO i = 1,npgrid
      DO j = 1,npgrid
         DO k = 1,npgrid
            WRITE(20,'(A6,3F18.10)') "H",(grid(k,j,i,l)*au2ang,l=1,3) 
         END DO
      END DO
   END DO
   CLOSE(20)
END IF

END SUBROUTINE GridFromExt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE readfchk(filein,IOut,verbosity,debug,selectionV)
USE globale

IMPLICIT NONE
INTEGER              :: istatus,reason,i,j,k,IOut,tempexp,tempbasis
INTEGER              :: tempatom,tempcount,tempshell,factor,tempindex
INTEGER              :: tempcountEXP,tempindexEXP,tempindexSHELL,tempcountSHELL
INTEGER              :: tempexpF,tempexpL
REAL*8, ALLOCATABLE  :: temp_coord(:) 
CHARACTER(LEN=100)   :: filein
CHARACTER(LEN=200)   :: linea
CHARACTER(LEN=100)   :: temp_char
REAL*8               :: TraceGS,TraceEX,TraceUEX
LOGICAL              :: verbosity,debug,selectionV

100 format (A)
1001 FORMAT(' GS Density Matrix Trace              : ',F12.6,/)
1002 FORMAT(' EX Rel Density Matrix Trace          : ',F12.6,/)
1003 FORMAT(' EX Unrel Density Matrix Trace        : ',F12.6,/)

OPEN(UNIT=10,FILE=trim(filein),STATUS='OLD',ACTION='READ',IOSTAT= istatus)
DO
   READ(10,100,IOSTAT=reason) linea
   IF (reason < 0 ) THEN
      EXIT
   ELSE IF (reason == 0) THEN
      IF (INDEX(linea,"Number of atoms").NE.0) THEN
         READ(linea,'(A55,I10)') temp_char , natomtot
         !WRITE(*,*) 'natom :',natomtot
         IF (selectionV.EQV..FALSE.) THEN
            natom = natomtot
            ALLOCATE(atomlist(natom))
            DO i=1,natom
               atomlist(i) = i
            END DO
            !WRITE(*,*) (atomlist(i), i =1,natom)
         END IF
         DO i = 1,natom
            IF (atomlist(i) > natomtot) THEN
               WRITE(IOut,*) "SELECTION LIST WRONG : SELECTION OUT OF RANGE"
               STOP
            END IF
         END DO
         ALLOCATE(nbasisXatom(natomtot),nshellXatom(natomtot),nexpXatom(natomtot),&
                  temp_coord(natomtot*3),coordinates(natomtot,3),atom_num_vec(natomtot))
         !WRITE(*,*) INDEX(linea,"Number of atoms")
         !WRITE(*,*) temp_char
      END IF
      IF (INDEX(linea,"Number of basis functions").NE.0) THEN
         READ(linea,'(A55,I10)') temp_char , nbasis
         nbasisD = nbasis*(nbasis+1)/2
         !WRITE(*,*) 'nbasis,nbasis2,nbasisD : ',nbasis,nbasis2,nbasisD
         ALLOCATE(gs_dens(nbasisD),ci_dens(nbasisD),ci_udens(nbasisD),&
                  mapAO(nbasis,4),AOtype(nbasis),densGS(nbasis,nbasis),densEX(nbasis,nbasis),&
                  densUEX(nbasis,nbasis))
      END IF
      IF (INDEX(linea,"Number of independent functions").NE.0) THEN
         READ(linea,'(A55,I10)') temp_char , nbasisi
         nbasis2 = nbasis*nbasisi
         !WRITE(*,*) 'nbasis,nbasis2,nbasisD : ',nbasis,nbasis2,nbasisD
         ALLOCATE(mo_coeff(nbasis2))
      END IF
      IF (INDEX(linea,"Number of contracted shells").NE.0) THEN
         READ(linea,'(A55,I10)') temp_char , nshell
         !WRITE(*,*) 'nshell : ',nshell
!         ALLOCATE(shell_type(nshell),npXshell(nshell),shell_atom_map(nshell),shell_center(nshell*3))
         ALLOCATE(shell_type(nshell),npXshell(nshell),shell_atom_map(nshell))
      END IF
      IF (INDEX(linea,"Number of primitives per shell").NE.0) THEN
         READ(10,*) (npXshell(i), i=1,nshell)
         nexp = SUM(npXshell)
         !WRITE(*,*) 'next tot : ', nexp
         ALLOCATE(prim_exp(nexp),prim_coeffic(nexp),prim_coeffic_sp(nexp))
      END IF
      IF (INDEX(linea,"Atomic numbers").NE.0) THEN
         READ(10,'(8X,6I12)') (atom_num_vec(i), i=1,natomtot)
      END IF
      IF (INDEX(linea,"Current cartesian coordinates").NE.0) THEN
         tempcount = natomtot*3
         READ(10,*) (temp_coord(i), i=1,natomtot*3)
         j=1
         DO i = 1,natomtot
            coordinates(i,1) = temp_coord(j)
            coordinates(i,2) = temp_coord(j+1)
            coordinates(i,3) = temp_coord(j+2)
            j=j+3
            !WRITE(*,'(3F15.8)') (coordinates(i,k),k=1,3)
         END DO
      END IF
      IF (INDEX(linea,"Alpha MO coefficients").NE.0) THEN
         !READ(10,*) (mo_coeff(i), i=1,nbasis2)
         READ(10,'(5E16.8)') (mo_coeff(i), i=1,nbasis2)
         !WRITE(*,'(5E16.8)') (mo_coeff(i), i=1,nbasis*nbasis)
      END IF
      IF (INDEX(linea,"Total SCF Density").NE.0) THEN
         READ(10,'(5E16.8)') (gs_dens(i), i=1,nbasisD)
         !WRITE(*,'(5E16.8)') (mo_coeff(i), i=1,nbasis*nbasis)
      END IF
      IF (INDEX(linea,"Total CI Rho").NE.0) THEN
         READ(10,'(5E16.8)') (ci_udens(i), i=1,nbasisD)
         !WRITE(*,'(5E16.8)') (ci_dens(i), i=1,nbasis*(nbasis+1)/2)
      END IF
      IF (INDEX(linea,"Total CI Density").NE.0) THEN
         READ(10,'(5E16.8)') (ci_dens(i), i=1,nbasisD)
         !WRITE(*,'(5E16.8)') (ci_dens(i), i=1,nbasis*(nbasis+1)/2)
      END IF
      IF (INDEX(linea,"Shell types").NE.0) THEN
         READ(10,*) (shell_type(i), i=1,nshell)
      END IF
      IF (INDEX(linea,"Shell to atom map").NE.0) THEN
         READ(10,*) (shell_atom_map(i), i=1,nshell)
      END IF
      IF (INDEX(linea,"Primitive exponents").NE.0) THEN
         READ(10,'(5E16.8)') (prim_exp(i), i=1,nexp)
      END IF 
      IF (linea(1:24).EQ."Contraction coefficients") THEN
         READ(10,'(5E16.8)') (prim_coeffic(i), i=1,nexp)
         !READ(10,*) (prim_coeff(i), i=1,nexp)
         !WRITE(11,'(5E16.8)') (prim_coeff(i), i=1,nexp)
         !DO i=1,nexp
         !    WRITE(11,'(E16.8)') prim_coeffic(i)
         !END DO
      END IF 
      IF (INDEX(linea,"P(S=P) Contraction coefficients").NE.0) THEN
         READ(10,'(5E16.8)') (prim_coeffic_sp(i), i=1,nexp)
         !READ(10,*) (prim_coeff(i), i=1,nexp)
         !WRITE(11,'(5E16.8)') (prim_coeff(i), i=1,nexp)
         !DO i=1,nexp
         !    WRITE(11,'(E16.8)') prim_coeffic(i)
         !END DO
      END IF 
!      IF (INDEX(linea,"Coordinates of each shell").NE.0) THEN
!         READ(10,'(5E16.8)') (shell_center(i), i=1,nshell*3)
!         !WRITE(*,'(5E16.8)') (shell_center(i), i=1,nshell*3)
!      END IF 
   END IF
END DO
CLOSE(10)
!WRITE(*,*) "temp_coord"
!WRITE(*,'(6F8.4)') temp_coord
!WRITE(*,*) "mo_ceff" 
!WRITE(*,'(6F8.4)') mo_coeff
!WRITE(*,*) "gs_dens"
!WRITE(*,'(6F8.4)') gs_dens
!WRITE(*,*) "ci_dens"
!WRITE(*,'(6F8.4)') ci_dens
!WRITE(*,*) "shell_type"
!WRITE(*,'(6I8)') shell_type
!WRITE(*,*) "shell_atom_map"
!WRITE(*,'(6I8)') shell_atom_map
!WRITE(*,*) "prim_exp"
!WRITE(*,'(6F12.4)') prim_exp
!WRITE(*,*) "prim_coeffic"
!WRITE(*,'(6F8.4)') prim_coeffic
! TRANSFORM THE DENSITY VECTOR IN MATRIX
tempcount = 1
!DO i = 1,nbasis
!   DO j = i,nbasis
!      densGS(j,i) = gs_dens(tempcount)
!      densEX(j,i) = ci_dens(tempcount)
!      tempcount = tempcount + 1 
!   END DO
!END DO
!WRITE(*,*) "densGS"
DO i = 1,nbasis
   DO j = 1,i
      densGS(i,j)  = gs_dens(tempcount)
      densEX(i,j)  = ci_dens(tempcount)
      densUEX(i,j) = ci_udens(tempcount)
      !densGS(j,i) = densGS(i,j)
      !densEX(j,i) = densEX(i,j)
      tempcount = tempcount + 1 
   END DO
END DO

IF (verbosity) THEN
   TraceGS  = 0.
   TraceEX  = 0.
   TraceUEX = 0.
   DO i = 1,nbasis
      TraceGS  = TraceGS  + densGS(i,i)
      TraceEX  = TraceEX  + densEX(i,i)
      TraceUEX = TraceUEX + densUEX(i,i)
   !   WRITE(*,'(6F8.4)') (densGS(i,j), j=1,i)
   END DO
   
   WRITE(IOut,1001) traceGS 
   WRITE(IOut,1002) traceEX 
   WRITE(IOut,1003) traceUEX 
END IF

!CREATE THE nshellXatom VECTOR
tempatom = 1
tempcount = 1
tempshell = shell_atom_map(1)
DO i = 2,nshell
   IF (tempshell.NE.shell_atom_map(i)) THEN
      nshellXatom(tempatom) = tempcount
      tempcount = 1
      tempatom = tempatom + 1
      tempshell = shell_atom_map(i)
   ELSE
      tempcount = tempcount + 1
   END IF
   IF (i.EQ.nshell) THEN
      nshellXatom(tempatom) = tempcount 
   END IF
END DO
!CREATE THE nbasisXatom and nexpXatom VECTORS
tempbasis = 0
tempcount = 1
tempexp   = 0
DO i = 1,natomtot
   DO j = 1,nshellXatom(i)
      IF (shell_type(tempcount).EQ.0) THEN
         factor = 1
      ELSE IF (shell_type(tempcount).EQ.1) THEN
         factor = 3
      ELSE IF (shell_type(tempcount).EQ.-1) THEN
         factor = 4
      ELSE IF (shell_type(tempcount).EQ.2) THEN
         factor = 6
      ELSE IF (shell_type(tempcount).EQ.-2) THEN
         factor = 5
      ELSE IF (shell_type(tempcount).EQ.3) THEN
         factor = 10
      ELSE IF (shell_type(tempcount).EQ.-3) THEN
         factor = 7
      END IF
      tempbasis = tempbasis + factor
      tempexp   = tempexp + npXshell(tempcount)
      tempcount = tempcount + 1
   END DO
   nexpXatom(i)    = tempexp
   nbasisXatom(i)  = tempbasis
   tempexp   = 0
   tempbasis = 0
END DO
!WRITE(*,*) nbasisXatom
!WRITE(*,*)
!WRITE(*,*) nexpXatom  
!WRITE(*,*)
!WRITE(*,*) nshellXatom
!!WRITE(*,*) SUM(nbasisXatom)
!WRITE(*,*)
!WRITE(*,*) atomlist
!WRITE(*,*)

! DEFINE THE NUMBER OF SELECTED EXP, BASIS AND SHELL
tempbasis = 0
tempEXP   = 0
tempshell = 0
DO i = 1,natom
   tempbasis = tempbasis + nbasisXatom(atomlist(i))
   tempEXP   = tempEXP   + nexpXatom(atomlist(i))
   tempshell = tempshell + nshellXatom(atomlist(i))
END DO
nbasisSelected = tempbasis
nEXPSelected   = tempEXP
nSHELLselected = tempshell
!WRITE(*,*) 'ATOMLIST : '
!!DO i = 1,natom
!WRITE(*,'(6I6)') atomlist
!!END DO
!WRITE(*,*) "nEXPSelected"
!WRITE(*,'(6I6)') nEXPSelected
!WRITE(*,*) "nexpXatom"
!WRITE(*,'(6I6)') nexpXatom
ALLOCATE(selectedAO(nbasisSelected),selectedEXP(nEXPSelected),selectedSHELL(nSHELLselected))
!WRITE(*,*) nbasisSelected
!WRITE(*,*) nSHELLselected
!WRITE(*,*) nbasisXatom

!CREATE THE selectedAO, selectedSHELL and selectedEXP VECTORS
tempcount      = 1
tempcountEXP   = 1
tempcountSHELL = 1
tempindex      = 0
tempindexEXP   = 0
tempindexSHELL = 0
DO i = 1,natom
!   IF (atomlist(i).GT.1) THEN
   DO j = 1,atomlist(i)-1
         tempindex      = tempindex      + nbasisXatom(j)
         tempindexEXP   = tempindexEXP   + nexpXatom(j)
         tempindexSHELL = tempindexSHELL + nshellXatom(j)
   END DO
!   END IF
   !WRITE(*,*) tempindex
   DO j = 1,nbasisXatom(atomlist(i))
      selectedAO(tempcount) = tempindex + j
      tempcount = tempcount + 1
   END DO
   DO j = 1,nexpXatom(atomlist(i))
      selectedEXP(tempcountEXP) = tempindexEXP + j
      tempcountEXP = tempcountEXP + 1
   END DO
   DO j = 1,nshellXatom(atomlist(i))
      selectedSHELL(tempcountSHELL) = tempindexSHELL + j
      tempcountSHELL = tempcountSHELL + 1
   END DO
   tempindex      = 0
   tempindexEXP   = 0
   tempindexSHELL = 0
END DO
!WRITE(*,*) 'selectedAO'
!!DO i = 1,natom
!WRITE(*,'(6I6)') selectedAO
!!END DO
!WRITE(*,*) "selectedEXP"
!WRITE(*,'(6I6)') selectedEXP
!WRITE(*,*) "selectedSHELL"
!WRITE(*,'(6I6)') selectedSHELL

!CREATE THE mapAO matrix AND AOType vector
tempcount  = 0
tempexp    = 0
tempexpF   = 0
tempexpL   = 0
DO i = 1,nshell
   IF (shell_type(i).EQ.0) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      mapAO(tempcount,2) = tempexp
      tempexp = tempexp + npXshell(i) - 1
      mapAO(tempcount,3) = tempexp
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "S"
   ELSE IF (shell_type(i).EQ.1) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "PX"
      DO j = 1,2
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "PY"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "PZ"
         END IF
      END DO
   ELSE IF (shell_type(i).EQ.-1) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "S"
      DO j = 1,3
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "PXS"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "PYS"
         ELSE IF (j .EQ. 3) THEN
            AOtype(tempcount) = "PZS"
         END IF
      END DO
   ELSE IF (shell_type(i).EQ.2) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "DX2"
      DO j = 1,5
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "DY2"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "DZ2"
         ELSE IF (j .EQ. 3) THEN
            AOtype(tempcount) = "DXY"
         ELSE IF (j .EQ. 4) THEN
            AOtype(tempcount) = "DXZ"
         ELSE IF (j .EQ. 5) THEN
            AOtype(tempcount) = "DYZ"
         END IF
      END DO
   ELSE IF (shell_type(i).EQ.-2) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "D3Z2R2"
      DO j = 1,4
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "DXZ"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "DYZ"
         ELSE IF (j .EQ. 3) THEN
            AOtype(tempcount) = "DX2Y2"
         ELSE IF (j .EQ. 4) THEN
            AOtype(tempcount) = "DXY"
         END IF
      END DO
   ELSE IF (shell_type(tempcount).EQ.3) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "FX3"
      DO j = 1,9
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "FY3"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "FZ3"
         ELSE IF (j .EQ. 3) THEN
            AOtype(tempcount) = "FXY2"
         ELSE IF (j .EQ. 4) THEN
            AOtype(tempcount) = "FX2Y"
         ELSE IF (j .EQ. 5) THEN
            AOtype(tempcount) = "FX2Z"
         ELSE IF (j .EQ. 6) THEN
            AOtype(tempcount) = "FXZ2"
         ELSE IF (j .EQ. 7) THEN
            AOtype(tempcount) = "FYZ2"
         ELSE IF (j .EQ. 8) THEN
            AOtype(tempcount) = "FY2Z"
         ELSE IF (j .EQ. 9) THEN
            AOtype(tempcount) = "FXYZ"
         END IF
      END DO
   ELSE IF (shell_type(tempcount).EQ.-3) THEN
      tempcount = tempcount + 1 
      tempexp = tempexp + 1
      mapAO(tempcount,1) = tempcount
      tempexpF = tempexp
      mapAO(tempcount,2) = tempexpF
      tempexp = tempexp + npXshell(i) - 1
      tempexpL = tempexp
      mapAO(tempcount,3) = tempexpL
      mapAO(tempcount,4) = shell_atom_map(i)
      AOtype(tempcount) = "FZ3ZR2"
      DO j = 1,6
         tempcount = tempcount + 1
         mapAO(tempcount,1) = tempcount
         mapAO(tempcount,2) = tempexpF
         mapAO(tempcount,3) = tempexpL
         mapAO(tempcount,4) = shell_atom_map(i)
         IF (j .EQ. 1) THEN 
            AOtype(tempcount) = "FXZ2XR2"
         ELSE IF (j .EQ. 2) THEN
            AOtype(tempcount) = "FYZ2YR2"
         ELSE IF (j .EQ. 3) THEN
            AOtype(tempcount) = "FX2ZY2Z"
         ELSE IF (j .EQ. 4) THEN
            AOtype(tempcount) = "FXYZ"
         ELSE IF (j .EQ. 5) THEN
            AOtype(tempcount) = "FX3XY2"
         ELSE IF (j .EQ. 6) THEN
            AOtype(tempcount) = "FX2YY3"
         END IF
      END DO
   END IF
END DO

IF (verbosity) THEN
   WRITE(IOut,*) "AO Index | First Exp. | Last Exp. | Atom Index | AO Type"
   DO i = 1,nbasis
      WRITE(IOut,'(I5,I13,I12,I13,A20)') mapAO(i,1),mapAO(i,2),mapAO(i,3),mapAO(i,4),AOtype(i)
   END DO
      !WRITE(IOut,*)
END IF
!WRITE(*,*) tempcount
!DO I = 1,nbasis
!   WRITE(*,*) mapAO(i,1),mapAO(i,2),mapAO(i,3),mapAO(i,4),AOtype(i)
!END DO
!WRITE(*,*) nbasisXatom 
!WRITE(*,*) selectedAO 

RETURN

END SUBROUTINE readfchk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcAOPAR(centro,coeff,tipoAO,expo,n,funz)
USE globale

INTEGER                      :: tipo,n,i,j
REAL*8, DIMENSION(3)         :: centro
REAL*8, DIMENSION(20)        :: coeff,expo
REAL*8, DIMENSION(nptot)     :: funz
REAL*8                       :: pref,pref2,temppref,pref3
REAL*8, PARAMETER            :: pi=3.1415927
CHARACTER(LEN=10)            :: tipoAO

funz = 0.
! N = (2*alpha/pi)**(3/4) * {[8*alpha]**(l+m+n)*l!*m!*n!/(2l!*2m!*2n!)}**1/2



!$omp parallel private(j,pref,pref2,pref1,pref3,i)
!$omp do
DO j = 1,nptot
   DO i = 1,n
      IF ( (TRIM(tipoAO).EQ."S") ) THEN
            pref = (2*expo(i)/pi)**(0.75)
            funz(j) = funz(j) + coeff(i)*pref*EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+&
                      (z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PX") .OR. (TRIM(tipoAO).EQ."PXS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PY") .OR. (TRIM(tipoAO).EQ."PYS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PZ") .OR. (TRIM(tipoAO).EQ."PZS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DX2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(x(j)-centro(1))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DY2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DZ2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(z(j)-centro(3))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DXY") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DXZ") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DYZ") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DX2Y2") ) THEN
            pref1 = ((2048*(expo(i)**7))/(9*(pi**3)))**(-0.5)
            pref2 = ((2048*(expo(i)**7))/(pi**3))**(-0.5)
            pref3 = (1./(2*(pref1-pref2)))**(0.5)
            !pref = (1./2.)*((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref3*((x(j)-centro(1))*(x(j)-centro(1))- &
                      (y(j)-centro(2))*(y(j)-centro(2)))*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."D3Z2R2") ) THEN
            pref1 = ((2048*(expo(i)**7))/(9*(pi**3)))**(-0.5)
            pref2 = ((2048*(expo(i)**7))/(pi**3))**(-0.5)
            pref3 = (1./(6*(pref1-pref2)))**(0.5)
            funz(j) = funz(j) + coeff(i)*pref3*((2*(z(j)-centro(3))*(z(j)-centro(3)))-&
                      (((x(j)-centro(1))**2)+(y(j)-centro(2))**2))*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FY3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((y(j)-centro(2))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FZ3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((z(j)-centro(3))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXY2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*((y(j)-centro(2))**2)*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2Y") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**2)*(y(j)-centro(2))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2Z") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**2)*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXZ2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*((z(j)-centro(3))**2)*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FYZ2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*((z(j)-centro(3))**2)*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FY2Z") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((y(j)-centro(2))**2)*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXYZ") ) THEN
            pref = ((524288*(expo(i)**9))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(y(j)-centro(2))*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FZ3ZR2") ) THEN
            pref1 = ((131072*(expo(i)**9))/((pi**3)))**(-0.5)
            pref2 = ((524288*(expo(i)**9))/(9*(pi**3)))**(-0.5)
            pref3 = (1./(2*(pref1+pref2)))**(0.5)
            funz(j) = funz(j) + coeff(i)*pref3*(-((x(j)-centro(1))**2)*(z(j)-centro(3))-& 
                      ((y(j)-centro(2))**2)*(z(j)-centro(3)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXZ2XR2") ) THEN
            pref1  = ((32768*(expo(i)**9))/(225*(pi**3)))**(-0.5)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(-0.5)
            pref3 = (1./(pref1+2.*pref2))**(0.5)
            funz(j) = funz(j) + coeff(i)*(-1)*pref3*((((x(j)-centro(1))**3))+&
                      ((x(j)-centro(1))*((y(j)-centro(2))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FYZ2YR2") ) THEN
            pref1  = ((32768*(expo(i)**9))/(225*(pi**3)))**(-0.5)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(-0.5)
            pref3 = (1./(pref1+2.*pref2))**(0.5)
            funz(j) = funz(j) + coeff(i)*(-1)*pref3*((((y(j)-centro(2))**3))+&
                      ((y(j)-centro(2))*((x(j)-centro(1))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX3XY2") ) THEN
            pref1  = ((32768*(expo(i)**9))/(225*(pi**3)))**(-0.5)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(-0.5)
            pref3 = (1./(pref1-pref2))**(0.5)
            funz(j) = funz(j) + coeff(i)*pref3*((((x(j)-centro(1))**3))-&
                      ((x(j)-centro(1))*((y(j)-centro(2))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2YY3") ) THEN
            pref1  = ((32768*(expo(i)**9))/(225*(pi**3)))**(-0.5)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(-0.5)
            pref3 = (1./(pref1-pref2))**(0.5)
            funz(j) = funz(j) + coeff(i)*pref3*((-((y(j)-centro(2))**3))+&
                      (((x(j)-centro(1))**2)*(y(j)-centro(2))))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      END IF  
   END DO
END DO
!$omp end do


!$omp end parallel

END SUBROUTINE calcAOPAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcAO(centro,coeff,tipoAO,expo,n,funz)
USE globale

INTEGER                        :: tipo,n,i,j
REAL*8, DIMENSION(3)           :: centro
REAL*8, DIMENSION(20)          :: coeff,expo
REAL*8, DIMENSION(nptot)       :: funz
REAL*8                         :: pref,pref2,temppref
REAL*8, PARAMETER              :: pi=3.1415927
CHARACTER(LEN=10)              :: tipoAO

funz = 0.

DO j = 1,nptot
   DO i = 1,n
      IF ( (TRIM(tipoAO).EQ."S") ) THEN
            pref = (2*expo(i)/pi)**(0.75)
            funz(j) = funz(j) + coeff(i)*pref*EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+&
                      (z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PX") .OR. (TRIM(tipoAO).EQ."PXS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PY") .OR. (TRIM(tipoAO).EQ."PYS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ((TRIM(tipoAO).EQ."PZ") .OR. (TRIM(tipoAO).EQ."PZS")) THEN
            pref = ((128*(expo(i)**5))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DX2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(x(j)-centro(1))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DY2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DZ2") ) THEN
            pref = ((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(z(j)-centro(3))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DXY") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(y(j)-centro(2))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DXZ") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DYZ") ) THEN
            pref = ((2048*(expo(i)**7))/(pi**3))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*(z(j)-centro(3))*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."DX2Y2") ) THEN
            pref = (1/SQRT(2.))*((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))*(x(j)-centro(1))- &
                      (y(j)-centro(2))*(y(j)-centro(2)))*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."D3Z2R2") ) THEN
            pref = (1/SQRT(6.))*((2048*(expo(i)**7))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((2*(z(j)-centro(3))*(z(j)-centro(3)))-&
                      (((x(j)-centro(1))**2)+(y(j)-centro(2))**2))*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+(y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FY3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((y(j)-centro(2))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FZ3") ) THEN
            pref = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((z(j)-centro(3))**3)*EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXY2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*((y(j)-centro(2))**2)*&
                      EXP(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2Y") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**2)*(y(j)-centro(2))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2Z") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((x(j)-centro(1))**2)*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXZ2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*((z(j)-centro(3))**2)*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FYZ2") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(y(j)-centro(2))*((z(j)-centro(3))**2)*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FY2Z") ) THEN
            pref = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*((y(j)-centro(2))**2)*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXYZ") ) THEN
            pref = ((524288*(expo(i)**9))/(9*(pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(x(j)-centro(1))*(y(j)-centro(2))*(z(j)-centro(3))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FZ3ZR2") ) THEN
            pref = (1/SQRT(2.))*((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*pref*(-((x(j)-centro(1))**2)*(z(j)-centro(3))-& 
                      ((y(j)-centro(2))**2)*(z(j)-centro(3)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FXZ2XR2") ) THEN
            pref  = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*(-1)*(1/SQRT(2.))*((pref*((x(j)-centro(1))**3))+&
                      (pref2*(x(j)-centro(1))*((y(j)-centro(2))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FYZ2YR2") ) THEN
            pref  = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*(-1)*(1/SQRT(2.))*((pref*((y(j)-centro(2))**3))+&
                      (pref2*(y(j)-centro(2))*((x(j)-centro(1))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX3XY2") ) THEN
            pref  = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*(1/SQRT(2.))*((pref*((x(j)-centro(1))**3))-&
                      (pref2*(x(j)-centro(1))*((y(j)-centro(2))**2)))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      ELSE IF ( (TRIM(tipoAO).EQ."FX2YY3") ) THEN
            pref  = ((32768*(expo(i)**9))/(225*(pi**3)))**(0.25)
            pref2 = ((131072*(expo(i)**9))/((pi**3)))**(0.25)
            funz(j) = funz(j) + coeff(i)*(1/SQRT(2.))*((-pref*((y(j)-centro(2))**3))+&
                      (pref2*((x(j)-centro(1))**2)*(y(j)-centro(2))))*&
                      exp(-expo(i)*((x(j)-centro(1))**2+&
                      (y(j)-centro(2))**2+(z(j)-centro(3))**2))
      END IF  
   END DO
END DO

END SUBROUTINE calcAO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcDensParNew(IOut,verbosity,debug,deltax,deltay,deltaz)
USE globale

INTEGER                                      :: i,j,k,m
REAL*8, DIMENSION(nptot)                     :: funcI
REAL*8, DIMENSION(20,nbasisSelected)         :: expoI,coeffI
REAL*8, DIMENSION(3,nbasisSelected)          :: centroI
INTEGER, DIMENSION(nbasisSelected)           :: nexpI
REAL*8                                       :: deltax,deltay,deltaz,area,DV
LOGICAL                                      :: verbosity,debug

ALLOCATE(densSELGS(nptot),densSELEX(nptot),densUSELEX(nptot),matrixAO(nptot,nbasisSelected))

countI      = 0
countY      = 0
countAOext  = 0
countAOint  = 0
countEXP    = 0

expoI  = 0.
coeffI = 0.
nexpI  = 0.

densSELGS  = 0.
densSELEX  = 0.
densUSELEX = 0.
matrixAO   = 0.

!DO i = 1,nbasisSelected 
!! EXTRACT THE COEFF AND EXP FOR ORBITAL I
!   nexpI = (mapAO(selectedAO(i),3)-mapAO(selectedAO(i),2)) + 1
!   DO k = 1,nexpI
!      expoI(k)  = prim_exp(mapAO(selectedAO(i),2)+k-1)
!   END DO
!   IF ( (AOtype(selectedAO(i))(1:3).EQ.'PXS') .OR. ((AOtype(selectedAO(i))(1:3).EQ.'PYS')).OR.&
!      ((AOtype(selectedAO(i))(1:3).EQ.'PZS')) ) THEN
!      DO k = 1,nexpI
!         coeffI(k) = prim_coeffic_sp(mapAO(selectedAO(i),2)+k-1)
!      END DO
!   ELSE
!      DO k = 1,nexpI
!         coeffI(k) = prim_coeffic(mapAO(selectedAO(i),2)+k-1)
!      END DO
!   END IF
!! SET THE CENTER OF ORBITAL I
!   DO k = 1,3
!      centroI(k) = coordinates(mapAO(selectedAO(i),4),k)
!   END DO
!! CALCULATE THE AO I
!   CALL calcAOPAR(centroI,coeffI,AOtype(selectedAO(i)),expoI,nexpI,matrixAO(:,i))
!!   IF (debug) THEN
!!      WRITE(IOut,'(A,3I6)') " ITERATION, AO AND NEXPI : ",I,selectedAO(i),nexpI 
!!      WRITE(IOut,*) "EXPONENTS"
!!      WRITE(IOut,'(5F12.4)') expoI
!!      WRITE(IOut,*) "COEFFICIENTS"
!!      WRITE(IOut,'(5F12.4)') coeffI
!!      WRITE(IOut,*) "The RELATIVE CENTER IS : " 
!!      WRITE(IOut,'(3F12.4)') centroI
!!      area = 0.
!!      DV = deltax*deltay*deltaz
!!      DO m = 1,nptot
!!         area = area + funcI(m)*DV
!!      END DO
!!      WRITE(IOut,'(A,F12.4)') " THE AREA OF THE ORBITAL I IS : ",area
!!      WRITE(IOut,*) 
!!   END IF
!   expoI  = 0.
!   coeffI = 0.
!END DO


DO i = 1,nbasisSelected 
! EXTRACT THE COEFF AND EXP FOR ORBITAL I
   nexpI(i) = (mapAO(selectedAO(i),3)-mapAO(selectedAO(i),2)) + 1
   DO k = 1,nexpI(i)
      expoI(k,i)  = prim_exp(mapAO(selectedAO(i),2)+k-1)
   END DO
   IF ( (AOtype(selectedAO(i))(1:3).EQ.'PXS') .OR. ((AOtype(selectedAO(i))(1:3).EQ.'PYS')).OR.&
      ((AOtype(selectedAO(i))(1:3).EQ.'PZS')) ) THEN
      DO k = 1,nexpI(i)
         coeffI(k,i) = prim_coeffic_sp(mapAO(selectedAO(i),2)+k-1)
      END DO
   ELSE
      DO k = 1,nexpI(i)
         coeffI(k,i) = prim_coeffic(mapAO(selectedAO(i),2)+k-1)
      END DO
   END IF
! SET THE CENTER OF ORBITAL I
   DO k = 1,3
      centroI(k,i) = coordinates(mapAO(selectedAO(i),4),k)
   END DO
   IF (debug) THEN
      WRITE(IOut,'(A,3I6)') " ITERATION, AO AND NEXPI : ",I,selectedAO(i),nexpI 
      WRITE(IOut,*) "EXPONENTS"
      WRITE(IOut,'(5F12.4)') (expoI(k,i),k=1,nexpI(i))
      WRITE(IOut,*) "COEFFICIENTS"
      WRITE(IOut,'(5F12.4)') (coeffI(k,i),k=1,nexpI(i))
      WRITE(IOut,*) "The RELATIVE CENTER IS : " 
      WRITE(IOut,'(3F12.4)') (centroI(k,i),k=1,3)
      area = 0.
      DV = deltax*deltay*deltaz
      DO m = 1,nptot
         area = area + funcI(m)*DV
      END DO
      WRITE(IOut,'(A,F12.4)') " THE AREA OF THE ORBITAL I IS : ",area
      WRITE(IOut,*) 
   END IF
END DO

!!$omp parallel private(i)
!!$omp do
DO i = 1,nbasisSelected
! CALCULATE THE AO I
   CALL calcAOPAR(centroI(:,i),coeffI(:,i),AOtype(selectedAO(i)),expoI(:,i),nexpI(i),matrixAO(:,i))
END DO
!!$omp end do
!!$omp end parallel


   !WRITE(IOut,*) "These are the EXP for each I"
DO i = 1,nbasisSelected 
   DO j = i,nbasisSelected
      IF (selectedAO(i) .EQ. selectedAO(j)) THEN
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + densGS(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
            densSELEX(m)  = densSELEX(m)  + densEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j) 
            densUSELEX(m) = densUSELEX(m) + densUEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
         END DO
         !$omp end do
         !$omp end parallel
      ELSE
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + 2.*densGS(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j) 
            densSELEX(m)  = densSELEX(m)  + 2.*densEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
            densUSELEX(m) = densUSELEX(m) + 2.*densUEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
         END DO
         !$omp end do
         !$omp end parallel
      END IF
   END DO
END DO

END SUBROUTINE calcDensParNew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcDensPar(IOut,verbosity,debug,deltax,deltay,deltaz)
USE globale

INTEGER                       :: i,j,k,m,nexpI,nexpJ
REAL*8, DIMENSION(nptot)      :: funcI,funcJ
REAL*8, DIMENSION(20)         :: expoI,coeffI,expoJ,coeffJ
REAL*8, DIMENSION(3)          :: centroI,centroJ
REAL*8                        :: deltax,deltay,deltaz,area,DV
LOGICAL                       :: verbosity,debug

ALLOCATE(densSELGS(nptot),densSELEX(nptot),densUSELEX(nptot))

countI      = 0
countY      = 0
countAOext  = 0
countAOint  = 0
countEXP    = 0

expoI  = 0.
coeffI = 0.
expoJ  = 0.
coeffJ = 0.

densSELGS  = 0.
densSELEX  = 0.
densUSELEX = 0.

   !WRITE(IOut,*) "These are the EXP for each I"
DO i = 1,nbasisSelected 
! EXTRACT THE COEFF AND EXP FOR ORBITAL I
   nexpI = (mapAO(selectedAO(i),3)-mapAO(selectedAO(i),2)) + 1
   DO k = 1,nexpI
      expoI(k)  = prim_exp(mapAO(selectedAO(i),2)+k-1)
   END DO
   IF ( (AOtype(selectedAO(i))(1:3).EQ.'PXS') .OR. ((AOtype(selectedAO(i))(1:3).EQ.'PYS')).OR.&
      ((AOtype(selectedAO(i))(1:3).EQ.'PZS')) ) THEN
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic_sp(mapAO(selectedAO(i),2)+k-1)
      END DO
   ELSE
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic(mapAO(selectedAO(i),2)+k-1)
      END DO
   END IF
! SET THE CENTER OF ORBITAL I
   DO k = 1,3
      centroI(k) = coordinates(mapAO(selectedAO(i),4),k)
   END DO
! CALCULATE THE AO I
   CALL calcAOPAR(centroI,coeffI,AOtype(selectedAO(i)),expoI,nexpI,funcI)
   IF (debug) THEN
      WRITE(IOut,'(A,3I6)') " ITERATION, AO AND NEXPI : ",I,selectedAO(i),nexpI 
      WRITE(IOut,*) "EXPONENTS"
      WRITE(IOut,'(5F12.4)') expoI
      WRITE(IOut,*) "COEFFICIENTS"
      WRITE(IOut,'(5F12.4)') coeffI
      WRITE(IOut,*) "The RELATIVE CENTER IS : " 
      WRITE(IOut,'(3F12.4)') centroI
      area = 0.
      DV = deltax*deltay*deltaz
      DO m = 1,nptot
         area = area + funcI(m)*DV
      END DO
      WRITE(IOut,'(A,F12.4)') " THE AREA OF THE ORBITAL I IS : ",area
      WRITE(IOut,*) 
   END IF
! EXTRACT THE COEFF AND EXP FOR ORBITAL J
   DO j = i,nbasisSelected
!   DO j = 1,i
      nexpJ = (mapAO(selectedAO(j),3)-mapAO(selectedAO(j),2)) + 1
      DO k = 1,nexpJ
         expoJ(k)  = prim_exp(mapAO(selectedAO(j),2)+k-1)
      END DO
      IF ( (AOtype(selectedAO(j))(1:3).EQ.'PXS').OR.(AOtype(selectedAO(j))(1:3).EQ.'PYS').OR.&
         (AOtype(selectedAO(j))(1:3).EQ.'PZS') ) THEN
         DO k = 1,nexpJ
            coeffJ(k) = prim_coeffic_sp(mapAO(selectedAO(j),2)+k-1)
         END DO
      ELSE
         DO k = 1,nexpJ
            coeffJ(k) = prim_coeffic(mapAO(selectedAO(j),2)+k-1)
         END DO
      END IF
! SET THE CENTER OF ORBITAL J
      DO k = 1,3
         centroJ(k) = coordinates(mapAO(selectedAO(j),4),k)
      END DO
! CALCULATE THE AO J
      CALL calcAOPAR(centroJ,coeffJ,AOtype(selectedAO(j)),expoJ,nexpJ,funcJ)
      IF (selectedAO(i) .EQ. selectedAO(j)) THEN
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + densGS(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densSELEX(m)  = densSELEX(m)  + densEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densUSELEX(m) = densUSELEX(m) + densUEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
         END DO
         !$omp end do
         !$omp end parallel
      ELSE
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + 2.*densGS(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densSELEX(m)  = densSELEX(m)  + 2.*densEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densUSELEX(m) = densUSELEX(m) + 2.*densUEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
         END DO
         !$omp end do
         !$omp end parallel
      END IF
      expoJ  = 0.
      coeffJ = 0.
   END DO
   expoI  = 0.
   coeffI = 0.
END DO

END SUBROUTINE calcDensPar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcDensNew(IOut,verbosity,debug,deltax,deltay,deltaz)
USE globale

INTEGER                       :: i,j,k,m,nexpI,nexpJ
REAL*8, DIMENSION(nptot)      :: funcI,funcJ
REAL*8, DIMENSION(20)         :: expoI,coeffI,expoJ,coeffJ
REAL*8, DIMENSION(3)          :: centroI,centroJ
REAL*8                        :: deltax,deltay,deltaz,area,DV
LOGICAL                       :: verbosity,debug

ALLOCATE(densSELGS(nptot),densSELEX(nptot),densUSELEX(nptot),matrixAO(nptot,nbasisSelected))

countI      = 0
countY      = 0
countAOext  = 0
countAOint  = 0
countEXP    = 0

expoI  = 0.
coeffI = 0.
expoJ  = 0.
coeffJ = 0.

densSELGS  = 0.
densSELEX  = 0.
densUSELEX = 0.
matrixAO   = 0.

DO i = 1,nbasisSelected 
! EXTRACT THE COEFF AND EXP FOR ORBITAL I
   nexpI = (mapAO(selectedAO(i),3)-mapAO(selectedAO(i),2)) + 1
   DO k = 1,nexpI
      expoI(k)  = prim_exp(mapAO(selectedAO(i),2)+k-1)
   END DO
   IF ( (AOtype(selectedAO(i))(1:3).EQ.'PXS') .OR. ((AOtype(selectedAO(i))(1:3).EQ.'PYS')).OR.&
      ((AOtype(selectedAO(i))(1:3).EQ.'PZS')) ) THEN
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic_sp(mapAO(selectedAO(i),2)+k-1)
      END DO
   ELSE
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic(mapAO(selectedAO(i),2)+k-1)
      END DO
   END IF
! SET THE CENTER OF ORBITAL I
   DO k = 1,3
      centroI(k) = coordinates(mapAO(selectedAO(i),4),k)
   END DO
! CALCULATE THE AO I
   CALL calcAO(centroI,coeffI,AOtype(selectedAO(i)),expoI,nexpI,matrixAO(:,i))
   IF (verbosity) THEN
      WRITE(IOut,'(A,3I6)') " ITERATION, AO AND NEXPI : ",I,selectedAO(i),nexpI 
      WRITE(IOut,*) "EXPONENTS"
      WRITE(IOut,'(5F12.4)') expoI
      WRITE(IOut,*) "COEFFICIENTS"
      WRITE(IOut,'(5F12.4)') coeffI
      WRITE(IOut,*) "The RELATIVE CENTER IS : " 
      WRITE(IOut,'(3F12.4)') centroI
      area = 0.
      DV = deltax*deltay*deltaz
      DO m = 1,nptot
         area = area + funcI(m)*DV
      END DO
      WRITE(IOut,'(A,F12.4)') " THE AREA OF THE ORBITAL I IS : ",area
      WRITE(IOut,*) 
   END IF
   expoI  = 0.
   coeffI = 0.
END DO

   !WRITE(IOut,*) "These are the EXP for each I"
DO i = 1,nbasisSelected 
   DO j = i,nbasisSelected
      IF (selectedAO(i) .EQ. selectedAO(j)) THEN
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + densGS(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
            densSELEX(m)  = densSELEX(m)  + densEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j) 
            densUSELEX(m) = densUSELEX(m) + densUEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
         END DO
         !$omp end do
         !$omp end parallel
      ELSE
         !$omp parallel private(m)
         !$omp do
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + 2.*densGS(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j) 
            densSELEX(m)  = densSELEX(m)  + 2.*densEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
            densUSELEX(m) = densUSELEX(m) + 2.*densUEX(selectedAO(j),selectedAO(i))*matrixAO(m,i)*matrixAO(m,j)
         END DO
         !$omp end do
         !$omp end parallel
      END IF
   END DO
END DO

END SUBROUTINE calcDensNew


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcDens(IOut,verbosity,debug,deltax,deltay,deltaz)
USE globale

INTEGER                        :: i,j,k,m,nexpI,nexpJ
REAL*8, DIMENSION(nptot)       :: funcI,funcJ
REAL*8, DIMENSION(20)          :: expoI,coeffI,expoJ,coeffJ
REAL*8, DIMENSION(3)           :: centroI,centroJ
REAL*8                         :: deltax,deltay,deltaz,area,DV
LOGICAL                        :: verbosity,debug

ALLOCATE(densSELGS(nptot),densSELEX(nptot),densUSELEX(nptot))

countI      = 0
countY      = 0
countAOext  = 0
countAOint  = 0
countEXP    = 0

expoI  = 0.
coeffI = 0.
expoJ  = 0.
coeffJ = 0.

densSELGS = 0.
densSELEX = 0.
densUSELEX = 0.

   !WRITE(IOut,*) "These are the EXP for each I"
DO i = 1,nbasisSelected 
! EXTRACT THE COEFF AND EXP FOR ORBITAL I
   nexpI = (mapAO(selectedAO(i),3)-mapAO(selectedAO(i),2)) + 1
   DO k = 1,nexpI
      expoI(k)  = prim_exp(mapAO(selectedAO(i),2)+k-1)
   END DO
   IF ( (AOtype(selectedAO(i))(1:3).EQ.'PXS') .OR. ((AOtype(selectedAO(i))(1:3).EQ.'PYS')).OR.&
      ((AOtype(selectedAO(i))(1:3).EQ.'PZS')) ) THEN
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic_sp(mapAO(selectedAO(i),2)+k-1)
      END DO
   ELSE
      DO k = 1,nexpI
         coeffI(k) = prim_coeffic(mapAO(selectedAO(i),2)+k-1)
      END DO
   END IF
! SET THE CENTER OF ORBITAL I
   DO k = 1,3
      centroI(k) = coordinates(mapAO(selectedAO(i),4),k)
   END DO
! CALCULATE THE AO I
   CALL calcAO(centroI,coeffI,AOtype(selectedAO(i)),expoI,nexpI,funcI)
   IF (verbosity) THEN
      WRITE(IOut,'(A,3I6)') " ITERATION, AO AND NEXPI : ",I,selectedAO(i),nexpI 
      WRITE(IOut,*) "EXPONENTS"
      WRITE(IOut,'(5F12.4)') expoI
      WRITE(IOut,*) "COEFFICIENTS"
      WRITE(IOut,'(5F12.4)') coeffI
      WRITE(IOut,*) "The RELATIVE CENTER IS : " 
      WRITE(IOut,'(3F12.4)') centroI
      area = 0.
      DV = deltax*deltay*deltaz
      DO m = 1,nptot
         area = area + funcI(m)*DV
      END DO
      WRITE(IOut,'(A,F12.4)') " THE AREA OF THE ORBITAL I IS : ",area
      WRITE(IOut,*) 
   END IF
! EXTRACT THE COEFF AND EXP FOR ORBITAL J
   DO j = i,nbasisSelected
      nexpJ = (mapAO(selectedAO(j),3)-mapAO(selectedAO(j),2)) + 1
      DO k = 1,nexpJ
         expoJ(k)  = prim_exp(mapAO(selectedAO(j),2)+k-1)
      END DO
      IF ( (AOtype(selectedAO(j))(1:3).EQ.'PXS').OR.(AOtype(selectedAO(j))(1:3).EQ.'PYS').OR.&
         (AOtype(selectedAO(j))(1:3).EQ.'PZS') ) THEN
         DO k = 1,nexpJ
            coeffJ(k) = prim_coeffic_sp(mapAO(selectedAO(j),2)+k-1)
         END DO
      ELSE
         DO k = 1,nexpJ
            coeffJ(k) = prim_coeffic(mapAO(selectedAO(j),2)+k-1)
         END DO
      END IF
! SET THE CENTER OF ORBITAL J
      DO k = 1,3
         centroJ(k) = coordinates(mapAO(selectedAO(j),4),k)
      END DO
! CALCULATE THE AO J
      CALL calcAO(centroJ,coeffJ,AOtype(selectedAO(j)),expoJ,nexpJ,funcJ)
      IF (selectedAO(i) .EQ. selectedAO(j)) THEN
!         DO k = 1,npgrid*npgrid*npgrid
!            densSELGS(k) = densSELGS(k) + densGS(selectedAO(j),selectedAO(i))*funcI(k)*funcJ(k) 
!            densSELEX(k) = densSELEX(k) + densEX(selectedAO(j),selectedAO(i))*funcI(k)*funcJ(k) 
!         END DO
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + densGS(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densSELEX(m)  = densSELEX(m)  + densEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densUSELEX(m) = densUSELEX(m) + densUEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
         END DO
      ELSE
!         DO k = 1,npgrid*npgrid*npgrid
!            densSELGS(k) = densSELGS(k) + 2.*densGS(selectedAO(j),selectedAO(i))*funcI(k)*funcJ(k) 
!            densSELEX(k) = densSELEX(k) + 2.*densEX(selectedAO(j),selectedAO(i))*funcI(k)*funcJ(k) 
!         END DO
         DO m=1,nptot
            densSELGS(m)  = densSELGS(m)  + 2.*densGS(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densSELEX(m)  = densSELEX(m)  + 2.*densEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
            densUSELEX(m) = densUSELEX(m) + 2.*densUEX(selectedAO(j),selectedAO(i))*funcI(m)*funcJ(m) 
         END DO
      END IF
      expoJ  = 0.
      coeffJ = 0.
   END DO
   expoI  = 0.
   coeffI = 0.
END DO

END SUBROUTINE calcDens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calcDCT(DX,DY,DZ,IOut)
USE globale

INTEGER                  :: i,IOut
REAL*8, DIMENSION(nptot) :: DeltaDens,DeltaUDens
REAL*8, DIMENSION(3)     :: Bpos,Bneg,BposU,BnegU,DipoleVersor,DipoleVersorU
REAL*8, DIMENSION(3)     :: EndCilinder,EndCilinderU
REAL*8                   :: d1,d2,dist,charge,DX,DY,DZ,d1U,d2U,distU,chargeU
REAL*8, PARAMETER        :: au2ang = 0.52917 

1000 FORMAT(' The Dct for the Selected Fragment is   : ',F12.8,' ang.',/)
1001 FORMAT(' The Transferred Charge is              : ',F12.8,' au. ',/)
1009 FORMAT(' The UDct for the Selected Fragment is  : ',F12.8,' ang.',/)
1010 FORMAT(' The Transferred UCharge is             : ',F12.8,' au.',/)
1002 FORMAT("draw cylinder {"3F8.4"} {"3F8.4"} radius 0.1")
1003 FORMAT("draw cone {"3F8.4"} {"3F8.4"} radius 0.2")
1004 FORMAT(' The Bar. of Positive Rel. Density is   : ',3F12.8,/)
1005 FORMAT(' The Bar. of Negative Rel. Density is   : ',3F12.8,/)
1020 FORMAT(' The Dct vector of Rel. Density is      : ',3F12.8,/)
1006 FORMAT(' The Bar. of Positive Urel. Density is  : ',3F12.8,/)
1007 FORMAT(' The Bar. of Negative Urel. Density is  : ',3F12.8,/)
1021 FORMAT(' The Dct vector of URel. Density is     : ',3F12.8,/)

!write(*,*) ' Calculation of the variation of the electronic...............'
DO i=1,nptot
   DeltaDens(i)  = densSELEX(i)  - densSELGS(i)    !!!! Definition of the Delta Matrix!!!!!
   DeltaUDens(i) = densUSELEX(i) - densSELGS(i)    !!!! Definition of the Delta Matrix!!!!!
END DO

!OPEN (UNIT=900, FILE='buttami.dat')
!WRITE(900,'(6E13.5)') (densSELGS(i),i=1,nptot)
!CLOSE(900)

Bpos = 0.0d0   ! Baricenter for the positive values
Bneg = 0.0d0   ! Baricenter for the negative values
d1   = 0.0d0   ! Total positive density
d2   = 0.0d0   ! Total negative density

BposU = 0.0d0   ! Baricenter for the positive values
BnegU = 0.0d0   ! Baricenter for the negative values
d1U   = 0.0d0   ! Total positive density
d2U   = 0.0d0   ! Total negative density

!write(*,*) ' Calculation of the two centroids positions...................'
DO i=1,nptot
        IF(DeltaDens(i) > 0) THEN
                Bpos(1)=Bpos(1)+DeltaDens(i)*x(i)
                Bpos(2)=Bpos(2)+DeltaDens(i)*y(i)
                Bpos(3)=Bpos(3)+DeltaDens(i)*z(i)
                d1 = d1+DeltaDens(i)
                BposU(1)=BposU(1)+DeltaUDens(i)*x(i)
                BposU(2)=BposU(2)+DeltaUDens(i)*y(i)
                BposU(3)=BposU(3)+DeltaUDens(i)*z(i)
                d1U = d1U+DeltaUDens(i)
                !charge=charge+Cube3(i) 
        ELSE IF(DeltaDens(i) < 0) THEN
                Bneg(1)=Bneg(1)+DeltaDens(i)*x(i)
                Bneg(2)=Bneg(2)+DeltaDens(i)*y(i)
                Bneg(3)=Bneg(3)+DeltaDens(i)*z(i)
                d2 = d2+DeltaDens(i)
                BnegU(1)=BnegU(1)+DeltaUDens(i)*x(i)
                BnegU(2)=BnegU(2)+DeltaUDens(i)*y(i)
                BnegU(3)=BnegU(3)+DeltaUDens(i)*z(i)
                d2U = d2U+DeltaUDens(i)
        END IF
END DO

Bpos(1)  =  Bpos(1)/d1
Bpos(2)  =  Bpos(2)/d1
Bpos(3)  =  Bpos(3)/d1
BposU(1) =  BposU(1)/d1
BposU(2) =  BposU(2)/d1
BposU(3) =  BposU(3)/d1

Bneg(1)  =  Bneg(1)/d2
Bneg(2)  =  Bneg(2)/d2
Bneg(3)  =  Bneg(3)/d2
BnegU(1) =  BnegU(1)/d2
BnegU(2) =  BnegU(2)/d2
BnegU(3) =  BnegU(3)/d2

dist    = (SQRT((Bpos(1)-Bneg(1))**2+(Bpos(2)-Bneg(2))**2+(Bpos(3)-Bneg(3))**2))*au2ang
distU   = (SQRT((BposU(1)-BnegU(1))**2+(BposU(2)-BnegU(2))**2+(BposU(3)-BnegU(3))**2))*au2ang
charge  = ((d1-d2)*DX*DY*DZ)/2
chargeU = ((d1U-d2U)*DX*DY*DZ)/2
!write(*,*)
Bpos  = Bpos*au2ang
Bneg  = Bneg*au2ang
BposU = BposU*au2ang
BnegU = BnegU*au2ang

WRITE(IOut,1000)  dist
WRITE(IOut,1001)  charge
WRITE(IOut,1009)  distU
WRITE(IOut,1010)  chargeU

DO i = 1,3
   DipoleVersor(i) = (Bpos(i)-Bneg(i))/dist
   EndCilinder(i)  = Bneg(i) + DipoleVersor(i)*dist*0.8
   DipoleVersorU(i) = (BposU(i)-BnegU(i))/distU
   EndCilinderU(i)  = BnegU(i) + DipoleVersorU(i)*distU*0.8
END DO

OPEN (UNIT=850, FILE='arrow_rel.tcl')
WRITE(850,'(A)') "draw color red"
WRITE(850,1002) Bneg(1),Bneg(2),Bneg(3),EndCilinder(1),EndCilinder(2),EndCilinder(3)
WRITE(850,1003) EndCilinder(1),EndCilinder(2),EndCilinder(3),Bpos(1),Bpos(2),Bpos(3)
CLOSE(850)

OPEN (UNIT=850, FILE='arrow_urel.tcl')
WRITE(850,'(A)') "draw color green"
WRITE(850,1002) BnegU(1),BnegU(2),BnegU(3),EndCilinderU(1),EndCilinderU(2),EndCilinderU(3)
WRITE(850,1003) EndCilinderU(1),EndCilinderU(2),EndCilinderU(3),BposU(1),BposU(2),BposU(3)
CLOSE(850)

OPEN (UNIT=850, FILE='coord.xyz')
WRITE(850,'(I6)') natomtot
WRITE(850,*)
DO i = 1,natomtot
   WRITE(850,'(I6,3F18.10)') atom_num_vec(i),(coordinates(i,j)*au2ang, j=1,3)
END DO
CLOSE(850)


WRITE(IOut,1004) Bpos
WRITE(IOut,1005) Bneg
WRITE(IOut,1020) Bpos(1)-Bneg(1),Bpos(2)-Bneg(2),Bpos(3)-Bneg(3)
WRITE(IOut,1006) BposU
WRITE(IOut,1007) BnegU
WRITE(IOut,1021) BposU(1)-BnegU(1),BposU(2)-BnegU(2),BposU(3)-BnegU(3)

END SUBROUTINE calcDCT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE writeDens(deltax,deltay,deltaz)
USE globale

REAL*8                                  :: deltax,deltay,deltaz,tempDelta
REAL*8, DIMENSION(natomtot)             :: charge
REAL*8, DIMENSION(npgrid,npgrid,npgrid) :: tempGSdensity,tempEXdensity,tempUEXdensity
INTEGER                                 :: tempCub,tempcount,i,j,k

tempcount = 1
DO i = 1,npgrid
   DO j = 1,npgrid
      DO k = 1,npgrid
         tempGSdensity(k,j,i)  = densSELGS(tempcount)
         tempEXdensity(k,j,i)  = densSELEX(tempcount)
         tempUEXdensity(k,j,i) = densUSELEX(tempcount)
         tempcount             = tempcount + 1
      END DO
   END DO
END DO
tempDelta = 0.
tempCub   = 1
OPEN(UNIT=10,FILE='gs_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npgrid,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F8.4,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npgrid
   DO j = 1 , npgrid
      WRITE(10,'(6E13.5)') (tempGSdensity(k,j,i),k=1,npgrid)
   END DO
END DO 
CLOSE(10)

OPEN(UNIT=10,FILE='rel_ex_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Relaxed Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npgrid,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F8.4,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npgrid
   DO j = 1 , npgrid
      WRITE(10,'(6E13.5)') (tempEXdensity(k,j,i),k=1,npgrid)
   END DO
END DO 
CLOSE(10)

OPEN(UNIT=10,FILE='urel_ex_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Unrelaxed Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npgrid,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npgrid,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F8.4,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npgrid
   DO j = 1 , npgrid
      WRITE(10,'(6E13.5)') (tempUEXdensity(k,j,i),k=1,npgrid)
   END DO
END DO 
CLOSE(10)

END SUBROUTINE writeDens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE writeDensEXT(deltax,deltay,deltaz)
USE globale

REAL*8                                  :: deltax,deltay,deltaz,tempDelta
REAL*8, DIMENSION(natomtot)             :: charge
REAL*8, DIMENSION(npextz,npexty,npextx) :: tempGSdensity,tempEXdensity,tempUEXdensity
INTEGER                                 :: tempCub,tempcount,i,j,k

tempcount = 1
DO i = 1,npextx
   DO j = 1,npexty
      DO k = 1,npextz
         tempGSdensity(k,j,i)  = densSELGS(tempcount)
         tempEXdensity(k,j,i)  = densSELEX(tempcount)
         tempUEXdensity(k,j,i) = densUSELEX(tempcount)
         tempcount             = tempcount + 1
      END DO
   END DO
END DO
tempDelta = 0.
tempCub   = 1
OPEN(UNIT=10,FILE='gs_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npextx,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npexty,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npextz,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F18.10,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npextx
   DO j = 1 , npexty
      WRITE(10,'(6E13.5)') (tempGSdensity(k,j,i),k=1,npextz)
   END DO
END DO 
CLOSE(10)

OPEN(UNIT=10,FILE='rel_ex_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Relaxed Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npextx,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npexty,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npextz,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F18.10,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npextx
   DO j = 1 , npexty
      WRITE(10,'(6E13.5)') (tempEXdensity(k,j,i),k=1,npextz)
   END DO
END DO 
CLOSE(10)

OPEN(UNIT=10,FILE='urel_ex_selected.cube',ACTION='WRITE')
!WRITE(10,'(A)') TRIM(filein)  
WRITE(10,'(A)') "Test"
WRITE(10,'(A)') "This is the Grid for the Unrelaxed Density Evaluation" 
WRITE(10,'(I6,3F18.10,I6)') natomtot,grid(1,1,1,1),grid(1,1,1,2),grid(1,1,1,3),tempCub
WRITE(10,'(I6,3F18.10)') npextx,deltax,tempDelta,tempDelta
WRITE(10,'(I6,3F18.10)') npexty,tempDelta,deltay,tempDelta
WRITE(10,'(I6,3F18.10)') npextz,tempDelta,tempDelta,deltaz
DO i = 1,natomtot
   charge(i) = atom_num_vec(i)
   WRITE(10,'(I6,F18.10,3F18.10)') atom_num_vec(i),charge(i),(coordinates(i,j),j=1,3)
END DO
DO i = 1 , npextx
   DO j = 1 , npexty
      WRITE(10,'(6E13.5)') (tempUEXdensity(k,j,i),k=1,npextz)
   END DO
END DO 
CLOSE(10)

END SUBROUTINE writeDensEXT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
