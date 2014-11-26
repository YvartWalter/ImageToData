!------------------------------------------------------------------!
! PGM2SPECTRUM                                                     !
! Convertie une image PGM en un spectre (x,y)                      !
! version 1.0.2                                                    !
! version 1.0.5 : bouclage                                         !
! version 1.1.0 post traitement pour l'extraction des nuages de    !
!               points                                             !
!------------------------------------------------------------------!
MODULE POST
 IMPLICIT NONE

 CONTAINS

 SUBROUTINE LOADXY(ofile,ncol,Ymin)
  IMPLICIT NONE
  CHARACTER(LEN=255), INTENT(In) :: ofile ! nom.out
  INTEGER, INTENT(In) :: ncol
  REAL(4), INTENT(In) :: Ymin
  REAL(4), DIMENSION(:), ALLOCATABLE :: AY, AX, PX, SX, SY
  INTEGER :: ufile=300,outfile=400,i,nncol,j,w,cn
  REAL(4) :: x,y,mx,my,xo,maxi,loc,k
  OPEN(unit=ufile,file=trim(ofile)//".out", status="OLD",action="READ")
  OPEN(unit=outfile,file="SP"//trim(adjustl(trim(ofile)))//".out", status="UNKNOWN",action="WRITE")
  write(*,*) "1"
  ALLOCATE(AY(1:ncol))
  ALLOCATE(AX(1:ncol))
  ALLOCATE(PX(1:ncol))
  write(*,*) "2"
  PX(:) = 0.0
  AX(:) = 0.0
  AY(:) = 0.0
  maxi=0.0
  mx=0.0
  my=0.0
  xo=0.0
  nncol=1
  write(*,*) "3"
  DO i=1,ncol
   READ(ufile,*) x,y
   AX(i) = x
   AY(i) = y
  END DO
  write(*,*) "4"
  w=0
  j=1
  nncol=0
  write(*,*) "5"
  DO i=1,ncol !recherche des pics, ils seront encadré par la valeur PX- et PX+
   if(AY(i).gt.Ymin) then
    if(w.eq.0) then
      PX(j) = AX(i)
      j=j+1
    end if
    w=1
   else
    if(w.eq.1) then
      PX(j) = AX(i-1)
      j=j+1
      nncol = nncol + 1
    end if
    w = 0
   end if
  END DO
  nncol=nncol
  write(*,*) nncol
  write(*,*) "6"
  ALLOCATE(SY(1:nncol))
  ALLOCATE(SX(1:nncol))
  SY(:) = 0.0
  SX(:) = 0.0
  write(*,*) "7"
  i=1
  k=1
  loc=PX(i) !est le point PX-
  DO WHILE(loc.gt.Ymin)!crible, recherche le point moyen entre PX- et PX+
    y=0.0
    x=0.0
    maxi=0.0
    DO j=0,ncol
      IF((AX(j).ge.loc).AND.(AX(j).le.PX(i+1))) THEN
        x=x+AX(j)
        y=y+AY(j)
        maxi=maxi+1.0
      END IF
    END DO
    SX(k) = x/maxi
    SY(k) = y/maxi
    i=i+2 !aller au point PX- suivant
    loc=PX(i)
    k=k+1
  END DO
  write(*,*) "Scatter plot : ",nncol,"points"
  CLOSE(unit=ufile)
  DO i=1,nncol
    WRITE(outfile,*) SX(i),SY(i)
  END DO
  CLOSE(unit=outfile)
  DEALLOCATE(SX)
  DEALLOCATE(SY)
  DEALLOCATE(PX)
  DEALLOCATE(AX)
  DEALLOCATE(AY)
 END SUBROUTINE LOADXY

END MODULE POST

PROGRAM PGM2SPECTRUM
 USE POST
 Implicit None
 INTEGER, PARAMETER :: ufile = 100, ofile = 200
 CHARACTER(LEN=255) :: cfile,oofile
 CHARACTER(len=2) :: cmt
 CHARACTER(len=1) :: ISnp
 INTEGER :: ncol = -1, nline = -1, i, nmcl,cc, cl, &
            ios, ii, ib
 INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIC
 REAL(4), DIMENSION(:), ALLOCATABLE :: SPECTRE
 REAL(4) :: pixel, mcl, &
            deltaX = 1.0,Xmin = 1.0, deltaY = 1.0, Ymin = 1.0, &
            Multf = 1.0, &
            Minspec = 0.0
 LOGICAL :: isfile = .TRUE.

 ISnp="F"
 ib = 1
 READ(5,*,iostat=ios) ib
 if(ios.lt.0) ib = 1
 WRITE(*,*) "START LOOP UP TO ",ib
 !----------------Debut de la boucle-------------------!
 DO ii = 1,ib
 WRITE(*,*) "CONVERT PGM",ii," on ",ib
 !-----------------------------------------------------!
 !Ouverture du fichier
 READ(5,*,iostat=ios) cfile
 IF(ios.lt.0) THEN
   WRITE(*,*) "I/Error : file name"
   STOP
 END IF
!  READ(5,*,iostat=ios) oofile
!  IF(ios.lt.0) WRITE(*,*) "I/warning : OUT = PGM2SPEC.out"
 WRITE(oofile,*) cfile !la sortie a le meme nom que la source
 !Lecture de la largeur totale en unité en x et de l'unité de départ
 READ(5,*,iostat=ios) deltaX,Xmin
 IF(ios.lt.0) WRITE(*,*) "I/warning : X = 1"
 !de même pour les unités en ordonnées
 READ(5,*,iostat=ios) deltaY, Ymin
 IF(ios.lt.0) WRITE(*,*) "I/warning : Y = 1"
 READ(5,'(A1)',ADVANCE='NO',iostat=ios) ISnp
 IF(ios.lt.0) WRITE(*,*) "I/warning : NP = False"
 READ(5,*,iostat=ios) Multf !facteur multiplicatif apparaissant sur les courbes
 IF(ios.lt.0) WRITE(*,*) "I/warning : Mfactor = 1"
 WRITE(*,*) "X",deltaX,Xmin
 WRITE(*,*) "Y",deltaY,Ymin
 WRITE(*,*) "Reading file ["//trim(cfile)//".pgm"//"]"
 !Lecture et initialisation
 INQUIRE(file=trim(cfile)//".pgm",exist = isfile)
 IF(isfile) THEN
    OPEN(unit=ufile, file=trim(cfile)//".pgm", action="READ", &
         access="SEQUENTIAL", status="OLD", form="FORMATTED")
 ELSE
   WRITE(*,*) "I/Error : file : "//trim(cfile)//".pgm does not exist !"
   STOP
 END IF
 ! Ecrite du spectre obtenu dans un fichier de sortie
 OPEN(unit=ofile, file=trim(oofile)//".out", status="UNKNOWN", &
      action="WRITE", access="SEQUENTIAL", form="FORMATTED")
 READ(ufile,*) cmt
  WRITE(*,*) "pgm magic number = "//cmt
 READ(ufile,*) cmt
 READ(ufile,*) ncol, nline
 Write(*,*) "Width = ",ncol, "Height = ",nline
 IF((nline*ncol).gt.0) THEN
    ALLOCATE(PIC(1:ncol,1:nline)) !largeur hauteur
   ALLOCATE(SPECTRE(1:ncol)) !spectre final, 1 valeur de y=moy(des nline (== y) noir) pour chaque valeur de ncol = x
 ELSE
   STOP
 END IF
  PIC(:,:) = 0.0
 SPECTRE(:) = 0.0
 !On commence à lire l'image par le haut
 ! # L'image est codée ligne par ligne en partant du haut
 ! # Chaque ligne est codée de gauche à droite
 ! Mais ici on fait partir notre référence du point en BAS à GAUCHE, ce qui détermine le zéro du repère
 ! L'objectif est de convertir le tracé noir en courbe (x,y) traçable par gnuplot
  !  cl = ligne courante
  !  cc = colonne courante
 cl = 0
 DO i = 1, nline
  cl = cl + 1
  DO cc = 1,ncol
   READ(ufile,*) pixel
   !Cette méthode ne selectionne que les pixels noirs
   IF(pixel.gt.0.0) THEN
    PIC(cc,cl) = 0
   ELSE
    PIC(cc,cl) = 1
   END IF
  END DO
 END DO
 !Maintenant latable PIC contient des 1 pour un pixel noir
 !et des 0 pour les pixels blancs

!Verification du code
!  DO cl=1, nline
! !   cl = cl + 1
!   DO cc = 1,ncol
!     WRITE(ofile,'(I1)',advance='NO') PIC(cc,cl)
!   END DO
!   WRITE(ofile,'(A)',advance='YES')
!  END DO

 DO cc=1,ncol
   !Parcour des colonnes : 1 valeur de x
   ! Il faut s'attendre à avoir plusieurs valeurs de y pour cette valeur de x, dans ce cas
   ! on peut en faire la moyenne pour n'avoir qu'un seul point
  mcl = 0.0
  nmcl = 0
   cl = 0
   DO i = 1,nline
     ! On part de la ligne en bas
     cl = cl + 1
     pixel = PIC(cc,cl)
     mcl = mcl + real(cl)*real(pixel)
     nmcl = nmcl + pixel
   END DO
   IF(nmcl.eq.0) THEN
    mcl = 0.0
   ELSE
    mcl = mcl / real(nmcl)
   END IF
   SPECTRE(cc) = mcl
 END DO
 !On place la plus faible valeur de spectre au 0
 Minspec = MINVAL(SPECTRE)
 DO i=1,ncol
  IF(SPECTRE(i).gt.0.0) SPECTRE(i) = real(nline) - SPECTRE(i)
  SPECTRE(i) = SPECTRE(i) - Minspec
  WRITE(ofile,*) (real(i) / real(ncol))*deltaX +Xmin, &
                 (SPECTRE(i)/real(nline))*deltaY/Multf + Ymin/Multf
 END DO
 !Fin
 DEALLOCATE(PIC)
 DEALLOCATE(SPECTRE)
 CLOSE(ufile)
 CLOSE(ofile)
 !----------------Fin de la boucle-------------------!
 READ(5,*) cmt
 WRITE(*,*) "Output send to "//trim(oofile)//".out"
 IF(ISnp.eq."T") THEN
   CALL LOADXY(oofile,ncol,Ymin)
   WRITE(*,*) "Scatter plot post treated and send to SP"//trim(adjustl(trim(oofile)))//".out"
 END IF
 END DO
 !-----------------------------------------------------!
 WRITE(*,*) "--------END--------"
END PROGRAM PGM2SPECTRUM
