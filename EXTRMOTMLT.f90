MODULE CONTOUR_DETECTION
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  CONTAINS

  SUBROUTINE CONTOUR(sx,sy,img,listecontour,nbpoints)
   IMPLICIT NONE
   INTEGER,INTENT(In) :: sx,sy
   REAL(KIND=DP),INTENT(In) :: img(sx,sy)
   REAL(KIND=DP),INTENT(Out) :: listecontour
   INTEGER, INTENT(Out) :: nbpoints

   REAL(KIND=DP), DIMENSION(1:sx,1:sy) :: imgcontour
   INTEGER,DIMENSION(1:2) :: point_depart,p,q
   INTEGER :: i,j,direction,nbpoints,px,py, &
              direction_courante,boucle,voisin
  ! 1. Initialisation
   point_depart(:) = 0
   p(:) = 0
   q(:) = 0
   direction = 1; ! vers la droite
   nbpoints = 0;  ! Nombre de points sur le contour
   DO i=1,sy
     DO j=1,sx
       imgcontour(j,i) = img(j,i)
     END DO
   END DO
  !2. Recherche du premier pixel noir - on suppose qu'il y en a au moins 1
   i=1
   j=1
   DO WHILE(img(i,j).gt.0.0_DP)
     i=i+1
     IF(i.gt.sx) THEN
      i=1
      j=j+1
     END IF
   END DO
  ! 3. Reperage du point_depart
   p(1) = i
   p(2) = j
   point_depart(1) = i
   point_depart(2) = j
  ! 4. direction_courante = direction
   direction_courante = direction
  ! 5. Sélectionner le voisin q de p dans la direction_courante indiquée
  !    par le code de freeman - on suppose qu'il existe
   boucle=1
   DO WHILE (boucle.eq.1)
     voisin=1
     DO WHILE (voisin.eq.1)
        voisin=1
        q(1) = p(1);
        q(2) = p(2);
        SELECT CASE(direction_courante)
          CASE(0)
           q(2) = q(2)+1
          CASE(1)
           q(1) = q(1)-1
           q(2) = q(2)+1
          CASE(2)
           q(1) = q(1)-1
          CASE(3)
           q(1) = q(1)-1
           q(2) = q(2)-1
          CASE(4)
           q(2) = q(2)-1
          CASE(5)
           q(1) = q(1)+1
           q(2) = q(2)-1
          CASE(6)
           q(1) = q(1)+1
          CASE(7)
           q(1) = q(1)+1
           q(2) = q(2)+1
          CASE DEFAULT
        END SELECT
       ! 6. si q n'est pas un point de la région ou s'il n'est pas voisin
       !    d'un pixel blanc (connexité 4 V) ou s'il a deja ete identifie
       !    comme un point du contour alors
        IF( (dble(img(q(1)-1,q(2)))+&
             dble(img(q(1)+1,q(2)))+&
             dble(img(q(1),q(2)-1))+&
              dble(img(q(1),q(2)+1)).eq.0).OR.&
            (.NOT.(img(q(1),q(2)).eq.0)).OR.&
            (imgcontour(q(1),q(2)).eq.128) ) THEN
           ! Si ce n'est pas le point de depart alors changement de direction et retour en 5
           IF ( .NOT.((q(1).eq.point_depart(1)).AND.(q(2)==point_depart(2))) ) THEN
              ! On tourne dans le sens des aiguilles d'une montre
              direction_courante = mod(direction_courante + 7 , 8)
              voisin=0
           END IF
        END IF
     END DO
     ! 7. Marquer p comme etant un point du contour
     nbpoints=nbpoints+1
     listecontour(nbpoints,1) = p(1)
     listecontour(nbpoints,2) = p(2)
     imgcontour(listecontour(nbpoints,1),listecontour(nbpoints,2)) = 128
     ! 8. direction = direction_courante
     direction = direction_courante
      ! On cherchera le voisin suivant le plus "a gauche"
     direction_courante=mod(direction_courante+3,8)
     ! 9. p = q
     p(1) = q(1)
     p(2) = q(2)
     ! 10. tant que p est different de point_depart alors retourner en 5
     IF((p(1).eq.point_depart(1)).AND.(p(2).eq.point_depart(2))) THEN
       boucle = 1
     ELSE
       boucle = 0
     END IF
   END DO
END SUBROUTINE CONTOUR

END MODULE CONTOUR_DETECTION
