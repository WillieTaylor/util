C     catomain adds peptides to CA file by (modified) method of Levitt (1977)
C     Original by Willie Taylor
C     Modified 14.11.96 By Andrew C.R. Martin, UCL to remove main program
C     and use a C wrapper so that it can be used as a filter
      
      SUBROUTINE ATOMG(FNM)
C     block for CA coordinates and data
      COMMON  /CAS/   ARES(1000),ASEQ(1000),AINS(1000),NAT,
     -     CA(3,1000),FRAME(3,3,1000)
      COMMON /ATOMS/  CB(3,1000),N(3,1000),C(3,1000),O(3,1000),H(3,1000)
      character*160 fnm
      REAL N
      INTEGER ASEQ,NAT
      CHARACTER*1 AINS
      CHARACTER*3 ARES
      CHARACTER*132 TEXT
      CHARACTER*1 CHAIN
      LOGICAL FIRST

      if(fnm(1:1).eq.' ') then
         ibrk = 5
      else
         open(ibrk,file=fnm,status='old')
      endif

      SCALE = 1.0
      NAT=0
      FIRST=.TRUE.
C     BEGIN loop over file
 300  CONTINUE
      READ(IBRK,'(A)',END=301) TEXT
      IF(TEXT(1:4).NE.'ATOM') GOTO 300
      IF (TEXT(14:16).NE.'CA') GOTO 300
      IF(FIRST) THEN
C     set the chain id
         FIRST=.FALSE.
         CHAIN=TEXT(22:22)
      ELSE
C     check for the same chain
         IF(TEXT(22:22).NE.CHAIN) THEN
C     different chain so exit
            PRINT*,'*NB* next chain'
            RETURN
         END IF
      END IF
      NAT = NAT + 1
C      PRINT*,NAT
      FRAME(1,1,NAT) = 1.0
      FRAME(2,2,NAT) = 1.0
      FRAME(3,3,NAT) = 1.0
      ARES(NAT) = TEXT(18:20)
      READ(TEXT(23:26),9000) ASEQ(NAT)
      AINS(NAT) = TEXT(27:27)
      READ(TEXT(31:38),9010) CA(1,NAT)
      READ(TEXT(39:46),9010) CA(2,NAT)
      READ(TEXT(47:54),9010) CA(3,NAT)
      IF(NAT.EQ.1) THEN
         D = DIST3(CA(1,NAT))
         DO I=1,3
            N(I,NAT) = CA(I,NAT)*(D-1.5)/D
            H(I,NAT) = CA(I,NAT)*(D-2.5)/D
         END DO
      END IF
      IF(NAT.LT.3) GOTO 300
      CALL ADDMAIN(CA(1,NAT-2),CA(1,NAT-1),CA(1,NAT),
     -     C(1,NAT-2), N(1,NAT-1),
     -     O(1,NAT-2), H(1,NAT-1),
     -     CB(1,NAT-1), 
     -     SCALE, ARES(NAT-1))
      CALL ADDCB(N(1,NAT-2),CA(1,NAT-2),C(1,NAT-2),CB(1,NAT-2),1.538)
      GOTO 300
C     FINISH of loop over file
 301  CONTINUE
      CALL ADDMAIC(CA(1,NAT-2),CA(1,NAT-1),CA(1,NAT),
     -     C(1,NAT-1), N(1,NAT),
     -                            O(1,NAT-1), H(1,NAT),
     -     CB(1,NAT-1), SCALE)
      CALL ADDCB(CA(1,NAT-1),N(1,NAT-1),C(1,NAT-1),
     -     CB(1,NAT-1), 1.538)
      D = DIST3(CA(1,NAT))
      DO I=1,3
         C(I,NAT) = CA(I,NAT)*(D-1.5)/D
         O(I,NAT) = CA(I,NAT)*(D-2.0)/D
      END DO
      CALL ADDCB(CA(1,NAT),N(1,NAT),C(1,NAT),
     -     CB(1,NAT), 1.538)
 9000 FORMAT(I4)
 9010 FORMAT(F8.3)
 9020 FORMAT(30X,3F8.3)
      RETURN
      END
      
      subroutine addmain(nca,ca,cca,c,n,o,h,cb,scale,res)
C     adds a dummy mainchain atoms BEFORE Ca
      real ca(3), c(3), n(3), o(3), h(3), cb(3)
      real nca(3), cca(3), mid(3)
      real a(3), b(3), aXb(3)
      character*3 res
      obond = 1.8
      cbond = 1.5
      hbond = 1.0
      do i=1,3
         mid(i) = (ca(i)+nca(i))/2.0
         a(i) = nca(i)-ca(i)
         b(i) = cca(i)-ca(i)
      end do
      call addcb(nca,ca,cca,cb,cbond)
      if (res.eq.'PRO') then
         do i=1,3
            aXb(i) = cb(i) - ca(i)
         end do
         dab = scale/cbond
         hbond = obond
      else
         call cross (a,b,aXb)
         dab = scale/dist3(aXb)
      end if
      do i=1,3
         n(i) = mid(i) - a(i)*0.125 + aXb(i)*0.25*dab
         c(i) = mid(i) + a(i)*0.125 - aXb(i)*0.50*dab
         o(i) = mid(i) - aXb(i)*obond*dab
         h(i) = mid(i) + aXb(i)*hbond*dab
      end do
      return
      end                     
      
      subroutine addmaiC(nca,ca,cca,c,n,o,h,cb,bond)
C     adds a dummy mainchain atoms AFTER Ca
      real ca(3), c(3), n(3), o(3), h(3), cb(3)
      real nca(3), cca(3), mid(3)
      real a(3), b(3), aXb(3)
      
      do i=1,3
         mid(i) = (ca(i)+cca(i))/2.0
         a(i) = nca(i)-ca(i)
         b(i) = cca(i)-ca(i)
         cb(i) = a(i) + b(i)
      end do
      call cross (a,b,aXb)
      dab = bond/dist3(aXb)
      dcb = bond*3.0/(2.0*dist3(cb)) 
      do i=1,3
         n(i) = mid(i) + b(i)*0.125
         c(i) = mid(i) - b(i)*0.125 + aXb(i)*0.5*dab
         o(i) = mid(i) + aXb(i)*2.0*dab
         h(i) = mid(i) - aXb(i)*dab
         cb(i) = ca(i) - cb(i)*dcb
      end do
      return
      end                     
      
      SUBROUTINE CAORD(fnm)
C     writes the .BRK file 
C     14.11.96 Corrected resnum printing to use ASEQ(I) not I By: ACRM
      
C     block for CA coordinates and data
      COMMON  /CAS/   ARES(1000),ASEQ(1000),AINS(1000),NAT,
     -     CA(3,1000),FRAME(3,3,1000)
      COMMON /ATOMS/  CB(3,1000),N(3,1000),C(3,1000),O(3,1000),H(3,1000)
      REAL N, CG(3)
      INTEGER ASEQ,NAT
      CHARACTER*1 AINS
      CHARACTER*3 ARES
      CHARACTER*160 fnm

      if(fnm(1:1).eq.' ') then
         ibrk = 6
      else
         open(ibrk,file=fnm,status='unknown')
      endif
      
      NA=0
      DO 100 I=1,NAT
         WRITE(IBRK,2000) NA+1,'  N   ',ARES(I),ASEQ(I),( N(J,I),J=1,3)
         WRITE(IBRK,2000) NA+2,'  CA  ',ARES(I),ASEQ(I),(CA(J,I),J=1,3)
         WRITE(IBRK,2000) NA+3,'  O   ',ARES(I),ASEQ(I),( O(J,I),J=1,3)
         WRITE(IBRK,2000) NA+4,'  C   ',ARES(I),ASEQ(I),( C(J,I),J=1,3)
         NA=NA+4
         IF(ARES(I).EQ.'GLY') GOTO 100
         NA=NA+1
         WRITE(IBRK,2000) NA,'  CB  ',ARES(I),ASEQ(I),(CB(J,I),J=1,3)
         IF(ARES(I).NE.'PRO') GOTO 100
         CALL ADDCG(CA(1,I),CB(1,I),N(1,I),H(1,I),CG)
         WRITE(IBRK,2000) NA+1,'  CG  ',ARES(I),ASEQ(I),(CG(J),J=1,3)
         WRITE(IBRK,2000) NA+2,'  CD  ',ARES(I),ASEQ(I),(H(J,I),J=1,3)
         NA=NA+2
 100  CONTINUE
      RETURN
 2000 FORMAT('ATOM',I7,2A,'  ',I4,4X,3F8.3,'  1.00  9.99')
      END
      
      subroutine addcg(ca,cb,n,cd,cg)
      real ca(3),cb(3),n(3),cd(3),cg(3)
      do i=1,3
         cg(i) = (cb(i)+cd(i))*0.5 + (cb(i)-ca(i)+cd(i)-n(i))*0.2
      end do
      return
      end
      
      subroutine addcb(n,ca,c,cb,cabond)
C     adds a dummy C_ to a GLY
      real ca(3), n(3), c(3), cb(3)
      real nca(3), cca(3)
      real x(3), y(3)
      
      do i=1,3
         nca(i) = ca(i)-n(i)
         cca(i) = ca(i)-c(i)
         x(i) = nca(i)+cca(i)
      end do 
      call cross (nca,cca,y)
C      ang = acos(-1.0)/2.0 - asin(1.0/sqrt(3.0))
      ang = 0.9128
      sx = cabond*cos(ang)/dist3(x)
      sy = cabond*sin(ang)/dist3(y)
      do i=1,3
         cb(i) = ca(i) + x(i)*sx + y(i)*sy
      end do
      return
      end                     
      SUBROUTINE CROSS(A,B,CR)
C     
C     routine to find cross product between two
C     vectors A and B
C     
      REAL A(3), B(3), CR(3)
C     
      CR(1) = A(2)*B(3) - B(2)*A(3)
      CR(2) = B(1)*A(3) - A(1)*B(3)
      CR(3) = A(1)*B(2) - B(1)*A(2)
      RETURN
      END
      function dist3 (a)
C     returns the length of vector a
      real    a(3)
      dist3 = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return 
      end
      
