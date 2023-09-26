      PROGRAM SUMMARY
C   PREPARES SUMMARY OF WAVEFUNCTIONS AND ENERGIES
C   MODIFIED VERSION TO CREATE INPUT FILES FOR SPINV ROUTINES.
C   MODIFIED TO ACCEPT COMMAND LINE ARGUMENT FOR NUMBER OF DIGITS.
C      IF NO ARGS ARE PROVIDED, DEFAULT TO 18 DIGITS AFTER DECIMAL POINT
C      IF VALUE OF ARG > 99, CUTS OFF LAST DIGIT (I.E. 275 BECOMES 27)
c      LAST UPDATED 05/02/2023  ERIC ENE, EVAN PETRIMOULX
      IMPLICIT REAL*16 (A-H,O-Z)
      DIMENSION YA(2,36000),WA(2),LA(2,36000),PA(2,36000),SA(36000),
     1   NBLKA(20)
      DIMENSION ISIZE(50),DIFF(50),RAT(50),NN(50),E(50),FYO(50),FXO(50)
     1   ,FYE(50),FXE(50),NSIZE(50),EE(50),N0(2),XN0(2),RN0(2),POW(2)
     2   ,NNOUT(50)
      DIMENSION XNUM(12),WCHK(25)
      INTEGER PA,PB,SA,SB, NUM_FIGS
      CHARACTER FNWVA*18,DATE*14,FN(25)*1,MATMAT*12,MATDAT*12,NAME*15,
     1   Q*1,CWA1*26,line*80,ANUM(12)*8,FORMAT*8,FMT115*32,
     1   FMT112*32,FMT23*32,CWA2(25)*26,MATPOW*12,NSH*13,CN*2,
     1   STATUS(50)*18,STATS*12,DATEW*8,FMT3(3)*12,
     1   FIGS_ARG*2, ARG_FMT1*64, ARG_FMT2*64, ARG_FMT3*64,
     1   RM_COMMAND*64
      CHARACTER*65 TITLE
      LOGICAL LEX,LFMT,LRYD
      DATA RAT(1),RAT(2),DIFF(1)/3*0.0/
      EQUIVALENCE (FN(1),FNWVA)
      OPEN(5,FILE='summary.dat',STATUS='OLD')
  333 FORMAT(1X,A26,3D26.19)
    3 FORMAT(1X,3D26.19)
    4 FORMAT(1X,3D30.23)
      KX3 = 0
C   REMOVE EXISTING FILES
      RM_COMMAND = 'rm 111SPOW.?AT'
      call system(RM_COMMAND)
   
C   RIGHT STATE READING ROUTINE.
C
      OPEN(2,FILE='DATE.DAT',STATUS='OLD')
      READ(2,'(A14)') DATE
C      WRITE(*,*)DATE
      CLOSE(2)
      READ(5,74) TITLE
C      WRITE(*,74) TITLE
      READ(5,*) LRYD
   74 FORMAT(A60,A14)
  161 ITER = 0
   31 FNWVA = '????'
      READ(5,'(A18)',END=160) FNWVA
      IF(FNWVA.EQ.'EXIT'.OR.FNWVA.EQ.'NEXT') GO TO 160
C      WRITE(*,'(A18,$)') FNWVA
      OPEN(1,FILE=FNWVA,STATUS='OLD',ERR=80)
      GO TO 81
   80 WRITE(*,'(/2A/)') FNWVA,' not found by summarv.'
      STOP
   81 READ(1,6) IZA,LRGLA,NSPNA,NW,NEIGA,NAME,TITLE
    6 FORMAT(3I3,I5,I3,1X,A12,A52)
      DATEW = TITLE(37:44)
      Q = 'Q'
C      IF(IZA.GE.10) TITLE = TITLE(2:65)
      IF(TITLE(1:1).NE.'Q'.AND.TITLE(2:2).NE.'Q') Q = ' '
      Z = IZA
      IH = 2
      IF(IZA.EQ.1) IH = 1
      NEIG0 = 0
      ITER = ITER + 1
      NWA = NW + NEIG0
      AMASS = 0.
C   FIND THE RIGHT FORMAT
      FMT3(1) = '(1X,3D26.19)'
      FMT3(2) = '(1X,3D30.23)'
      FMT3(3) = '(1X,3D38.31)'
      DO 40 IFF=1,3
      IF = IFF
      IF(IF.GT.1) BACKSPACE(1)
      READ(1,FMT3(IF),ERR=40) WA(1),WA(2),AMASS
      GO TO 41
   40 CONTINUE
C
   41 SCALE = 1
      IF(LRYD) SCALE = 1.Q0 - AMASS
    8 READ(1,1) NBLA,(NBLKA(I+1),I=1,NBLA),(NN(I+1),I=1,NBLA-IH+1),MAR12
     1   ,MAR1,KONO
    1 FORMAT(16I5)
      do 402 i=1,nbla-ih+1
  402 nnout(i+1) = nn(i+1)
      NBLKA(1) = 0
      DO 12 I=1,NBLA
      NBLKA(I+1) = NBLKA(I+1) + NEIG0
      READ(1,2) NB1,NB2,LA(1,NB1+NEIG0),LA(2,NB1+NEIG0),YA(1,NB1+NEIG0)
     1   ,YA(2,NB1+NEIG0)
C   SKIP A LINE IF THE WAVE FUNDTION FILE CONTAINS FURTHER DATA AFTER KONO
      IF(I.EQ.1.AND.YA(1,NB1+NEIG0).EQ.0.0) THEN
         READ(1,2) NB1,NB2,LA(1,NB1+NEIG0),LA(2,NB1+NEIG0),
     1   YA(1,NB1+NEIG0),YA(2,NB1+NEIG0)
      ENDIF
C   CHECK FOR NO SCREENED HYDROGENIC TERM.
      IF(I.EQ.1.AND.IH.EQ.2.AND.NB2.GT.NB1) THEN
         IH = 1
         BACKSPACE(1)
         BACKSPACE(1)
         GO TO 8
      ENDIF
C      WRITE(*,2) NB1,NB2,LA(1,NB1+NEIG0),LA(2,NB1+NEIG0),YA(1,NB1+NEIG0)
C     1   ,YA(2,NB1+NEIG0)
    2 FORMAT(3I5,I2,2F20.12)
      NB1 = NB1 + NEIG0
      NB2 = NB2 + NEIG0
      DO 13 K=1,2
      DO 13 J=NB1,NB2
      LA(K,J) = LA(K,NB1)
      YA(K,J) = YA(K,NB1)
   13 CONTINUE
      IF(LRGLA.GT.0.AND.I.GT.1.AND.LA(2,NB1).EQ.LA(1,NB1-1).AND.
     1   LA(1,NB1).EQ.LA(2,NB1-1)) KX3 = 1
   12 CONTINUE
C      read(*,*)
C      READ(1,113)(PA(1,K),PA(2,K),SA(K),K=1,NWA)
C  113 FORMAT(10(I3,2I2))
      CLOSE(1,STATUS='KEEP')
c      WRITE(*,996) IZA,LRGLA,NSPNA,NW,NEIGA,TITLE
  996 FORMAT(1H ,3HZ =,I3,6H   L =,I3,6H   S =,I3,6H   N =,I4,
     1   '   NEIGA =',I3/1X,A60)
C   CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
      IF(NEIGA.EQ.0) GO TO 32
      N = LRGLA + NEIGA
      C = 2.Q0*(IZA-1)/(N*IZA)
      YA(1,1) = 1.Q0
      YA(2,1) = C/2.
    7 FORMAT(6I3,A60)
   32 DO 102 K=1,NBLA
      K1 = NBLKA(K+1)
      K2 = K1 + 1
      IF(K2.GT.NWA) K2 = NWA
      YA(1,K) = YA(1,K1)
      YA(2,K) = YA(2,K1)
      LA(1,K) = LA(1,K1)
      LA(2,K) = LA(2,K1)
  102 CONTINUE
      NBLKA(IH) = 0
      E(ITER) = Z*Z*0.5*WA(1)
   30 IF(ITER.GT.1) GO TO 21
C
C   OPEN OUTPUT FILES WITH NAMES ZnSLPOW.MAT FOR SUMMARY OF RESULTS AND
C                                ZnSLPOW.DAT FOR REGENERATED INPUT FILE.
C   'POW' IS THE EXTENSION FOUND ON THE WAVE FUNCTION FILE, AND IS REPLACED
C   BY 'POL' FOR THE FINITE NUCLEAR MASS CASE.
C
      MATMAT =FNWVA(1:4)//FNWVA(9:11)//'.MAT'
      MATDAT =FNWVA(1:4)//FNWVA(9:11)//'.DAT'
      IF(IZA.GE.10) MATMAT=FNWVA(1:5)//FNWVA(10:12)//'.MAT'
      IF(IZA.GE.10) MATDAT=FNWVA(1:5)//FNWVA(10:12)//'.DAT'
  203 OPEN (2,FILE=MATMAT,STATUS='NEW')
      OPEN (3,FILE=MATDAT,STATUS='NEW')
C
      KX1 = 0
      KX2 = 0
      kxx = 0
      IF(YA(1,2).EQ.1.) KX1 = 1
      IF(YA(1,3).EQ.YA(1,2).AND.YA(2,3).EQ.YA(2,2)) KX2 = 1
      NBET = NBLA - KX2 + 1 - IH
C    WRITE REGENERATED INPUT FILE.
      WRITE(3,20) NBET,KX1,KX3,AMASS,0.01,DATE
   20 FORMAT(3(I1,','),F20.18,',',F6.4,',,,,,,',4X,A14)
C
      MATPOW = 'dpold.pow'
      IF(AMASS.NE.0) MATPOW = 'dpold.pol'
      OPEN(8,FILE=MATPOW,STATUS='UNKNOWN')
      CLOSE(8,STATUS='DELETE')
      OPEN(8,FILE=MATPOW,STATUS='NEW')
      WRITE(8,200) NBET,KX1,KX2,AMASS,0.010,
     1',.FALSE.,.false.,0,.FALSE.,0.,0                ||'
  200 FORMAT(3(I1,','),F20.18,',',F5.3,A)
C
   21 KX1 = 0
      KX2 = 0
      kxx = 0
      IF(YA(1,2).EQ.1.) KX1 = 1
      IF(YA(1,3).EQ.YA(1,2).AND.YA(2,3).EQ.YA(2,2)) KX2 = 1
      NBET = NBLA - KX2 + 1 - IH
      KX3 = KX2
      IF(LRGLA.GE.1) KX3 = 1
      WRITE(3,22) WA(1),15,NEIGA
   22 FORMAT(F15.13,',',I2,',',I2)
      IF(LRGLA.GE.1.AND.KX2.EQ.0.and.kx2.eq.1) THEN
        DO 401 II=2,NBET
        I = NBET - II + 2
        NNOUT(I+2) = NN(I+1)
  401   CONTINUE
        kxx = 1
        NNOUT(3) = 0
      ENDIF
      DO 35 IT=35,48
      ITTS = IT
   35 IF(TITLE(IT:IT).EQ.',') GO TO 36
   36 STATS = TITLE(ITTS:ITTS+11)
      DATEW = TITLE(ITTS-8:ITTS-1)
C      write(*,*) 'stats =',stats
      write(*,*) title(35:)
C      write(*,*) itts,datew
      TITLE(ITTS+1:) = ' '
      IF(ITTS.EQ.48) ITTS = 45
      DO 26 IT=1,10
      ITT = IT
   26 IF(TITLE(IT:IT).NE.' ') GO TO 28
   28 IF(TITLE(ITT:ITT).EQ.'Q') ITT = ITT + 1
      DO 33 IT=ITT,ITT+5 
      ITT = IT
   33 IF(TITLE(IT:IT).NE.' ') GO TO 34
   34 IF(TITLE(ITT:ITT+1).NE.'Z=') THEN
        N = LRGLA + NEIGA
        WRITE(CN,'(I2)') N
        LINE = TITLE(ITT:)
        TITLE(ITT:) = 'Z='//FNWVA(1:1)//' '//FNWVA(2:2)//' '//
     1     FNWVA(3:4)//'   '//LINE
        IF(IZA.GT.9) TITLE(ITT:) = 'Z='//FNWVA(1:2)//' '//FNWVA(3:3)
     1     //' '//FNWVA(4:5)//'  '//LINE
        IF(N.GT.9) TITLE(ITT:) = 'Z='//FNWVA(1:1)//' '//FNWVA(2:2)
     1     //' '//CN//FNWVA(4:4)//'  '//LINE
        IF(IZA.GT.9.AND.N.GT.9) TITLE(ITT:) = 'Z='//FNWVA(1:2)//' '//
     1     FNWVA(3:3)//' '//CN//FNWVA(5:5)//' '//LINE
      ENDIF
      WRITE(3,24) IZA,LRGLA,NSPNA,LA(1,NBLA),FNWVA,
     1   ' '//TITLE(ITT:ITT+34)//DATEW//STATS
   24 FORMAT(I2,',',3(I1,','),'''',A12,''',''',A52,'''')
C      STATUS(ITER) = TITLE(ITTS-12:ITTS-1)//STATS
      STATUS(ITER) = ' '//DATEW//STATS
C   ADJUST OUTPUT FORMAT
      LFMT = .FALSE.
      DO 29 I=KX2+IH+1,NBLA
   29 IF(YA(1,I).GE.10.OR.YA(2,I).GE.10) LFMT = .TRUE.
      FMT23 = '(F7.5,'','',F7.5,8(('','',F7.5)))'
      IF(LFMT) FMT23 = '(F7.5,'','',F7.5,8(('','',F8.5)))'
      WRITE(3,FMT23) YA(1,IH),YA(2,IH),(YA(1,I),YA(2,I),I=KX2+IH+1,NBLA)
   23 FORMAT(F6.4,',',F7.5,(8(',',F7.5)))
      WRITE(3,25) (NNOUT(I),I=2,NBLA+2-IH+kxx),0,MAR12,MAR1,KONO
   25 FORMAT(10(I2,','))
C
      WRITE(8,22) WA(1), 0,NEIGA
c      write(*,*) itt,title(itt:)
      WRITE(8,24) IZA,LRGLA,NSPNA,LA(1,NBLA),FNWVA
     1   ,' '//TITLE(ITT:ITT+34)//DATEW//STATS
      WRITE(8,FMT23) YA(1,IH),YA(2,IH),(YA(1,I),YA(2,I),I=KX2+IH+1,NBLA)
      WRITE(8,25) (NNOUT(I),I=2,NBLA+2-IH+kxx),0,MAR12,MAR1,KONO
      IF(ITER.EQ.1) THEN
        WRITE(8,200) NBET,KX1,KX2,AMASS,0.010,
     1  ',.false.,.false., 3,.TRUE. ,0. ,0 ,1 3P states ||'
        WRITE(8,'(A)') '=============================LGENB===LPOW==KDIV=
     1LSHORT===LINC=S============   '
        WRITE(8,22) WA(1), 0,NEIGA
        WRITE(8,24) IZA,LRGLA,NSPNA,LA(1,NBLA),FNWVA,TITLE(ITT:)
        WRITE(8,FMT23)
     1    YA(1,IH),YA(2,IH),(YA(1,I),YA(2,I),I=KX2+IH+1,NBLA)
        WRITE(8,25) (NNOUT(I),I=2,NBLA+2-IH+kxx),0,MAR12,MAR1,KONO
      ENDIF
C
C
C   WRITE SUMMARY FILE.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NFIG = 18
      FORMAT = '(F21.19)'
      WRITE(FORMAT(3:4),'(I2)') NFIG+2
      WRITE(FORMAT(6:7),'(I2)') NFIG
      FMT115 = '(6X,I4,1X,A24,'//FORMAT(2:7)//',F8.2)'
      WRITE(FMT115(12:13),'(I2)') NFIG+6
      FMT112 = '(6X,F22.15,6F9.5)'
      WRITE(FMT112(9:10),'(I2)') NFIG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(ITER.EQ.1) WRITE(2,114) TITLE,DATE,AMASS,(LA(1,I),
     1   LA(2,I),I=IH,NBLA)
      IF(ITER.GT.1) DIFF(ITER) = -(E(ITER) - E(ITER-1))
      RAT(ITER) = 0.
      IF(ITER.GT.2.AND.DIFF(ITER).NE.0.) RAT(ITER) = 
     1   DIFF(ITER-1)/DIFF(ITER)
      IF(RAT(ITER).LT.1..AND.RAT(ITER).NE.0.) THEN
         WRITE(*,'(I4,''  RAT ='',F10.3)') ITER,RAT(ITER)
         RAT(ITER) = 1.1
C         READ(*,*)
      ENDIF
      ISIZE(ITER) = NWA
      NSIZE(ITER) = NN(2) - NEIGA
  114 FORMAT(//6X,A47,3X,A14/
     1   6X,' mu/M =',F15.13,6(:4X,'(',I1,',',I1,')'))
      NSH = ' '
      IF(IH.EQ.1) NSH = 'No scr. hyd. '
      IF(ITER.EQ.1) WRITE(2,118) NSH
  118 FORMAT(6X,A13,'R12 R1 KONO')
      IF(ITER.GT.1) WRITE(2,111)
      WRITE(2,111) FNWVA,Q,MAR12,MAR1,KONO,(NBLKA(I+1)-NBLKA(I),
     1   NN(2+I-IH),I=IH,NBLA)
  111 FORMAT(6X,A12,A1,3I3,6(I5,'(',I2,')'))
      WRITE(2,FMT112) WA(1),(YA(1,K),K=IH,NBLA)
      WRITE(2,117) (YA(2,K),K=IH,NBLA)
  112 FORMAT(6X,F22.15,6F9.5)
  117 FORMAT(6X,22X,6F9.5)
  101 FORMAT(I2,4X,5I4,2F9.5)
      GO TO 31
  160 IF(FNWVA.EQ.'????') STOP 160
      WRITE(*,116)
      WRITE(2,116)
  116 FORMAT(/14X,'ENERGIES',15X,'DIFFERENCES',7X,'RATIOS',5X,'STATUS')
      ICO = 0
      ICE = 0
      DO 61 I=1,ITER
      ET = E(I) - Z**2/2 - (Z-1)**2/(2*N**2)
      OUT1 = ET*SCALE
      OUT2 = DIFF(I)*SCALE

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     get command line argument for number of digits on right side of decimal point

      IF (COMMAND_ARGUMENT_COUNT() .ne. 1) THEN
              NUM_FIGS = 18
      ELSE 
              CALL GET_COMMAND_ARGUMENT(1, FIGS_ARG)
              READ(FIGS_ARG, *) NUM_FIGS
      END IF  

      ARG_FMT1 = '(I4, 1X, F26.21, 33X, A18)'
      WRITE(ARG_FMT1(11:15), '(I2, A1, I2)') NUM_FIGS + 7, '.', NUM_FIGS

      ARG_FMT2 = '(I4, 1X, F26.21, F24.21, 8X, A18)'
      WRITE(ARG_FMT2(11:15), '(I2, A1, I2)') NUM_FIGS + 7, '.', NUM_FIGS
      WRITE(ARG_FMT2(19:23), '(I2, A1, I2)') NUM_FIGS + 7, '.', NUM_FIGS

      ARG_FMT3 = '(I4, 1X, F26.21, F24.21, F7.2, 1X, A18)'
      WRITE(ARG_FMT3(11:15), '(I2, A1, I2)') NUM_FIGS + 7, '.', NUM_FIGS
      WRITE(ARG_FMT3(19:23), '(I2, A1, I2)') NUM_FIGS + 7, '.', NUM_FIGS

  119 IF(I.EQ.1) WRITE(2,ARG_FMT1) ISIZE(I),OUT1,STATUS(I)
      IF(I.EQ.2) WRITE(2,ARG_FMT2) ISIZE(I),OUT1,OUT2,STATUS(I)
      IF(I.GE.3) WRITE(2,ARG_FMT3) ISIZE(I),OUT1,OUT2,RAT(I),STATUS(I)
      IF(I.EQ.1) WRITE(*,ARG_FMT1) ISIZE(I),OUT1,STATUS(I)
      IF(I.EQ.2) WRITE(*,ARG_FMT2) ISIZE(I),OUT1,OUT2,STATUS(I)
      IF(I.GE.3) WRITE(*,ARG_FMT3) ISIZE(I),OUT1,OUT2,RAT(I),STATUS(I)
  115 FORMAT(I4,1X,F26.21,F24.21,F7.2)

c      PRINT *, ARG_FMT1

!   135 FORMAT(I4,1X,F26.21,A18)
!   136 FORMAT(I4,1X,F26.21,F24.21,8X,A18)
!   137 FORMAT(I4,1X,F26.21,F24.21,F7.2,1X,A18)
      IF(I.LT.3) GO TO 61
      IF(IAND(I,1).EQ.0) GO TO 62
      ICO = ICO + 1
      FYO(ICO) = LOG(RAT(I)-1.)
C   CORRECT FOR INCREASING RATIOS
      IF(I.GE.5.AND.RAT(I).GT.RAT(I-2)) THEN
        AVRAT = (RAT(I) + RAT(I-2))/2
        FYO(ICO) = LOG(AVRAT-1.)
        FYO(ICO-1) = LOG(AVRAT-1.)
      ENDIF
      FXO(ICO) = LOG(1.Q0/NSIZE(I))
      NSZO = NSIZE(I)
CC      WRITE(6,FMT115) NSIZE(I),FXO(ICO),FYO(ICO)
      GO TO 61
   62 ICE = ICE + 1
      FYE(ICE) = LOG(RAT(I)-1.)
C   CORRECT FOR INCREASING RATIOS
      IF(I.GE.6.AND.RAT(I).GT.RAT(I-2)) THEN
        AVRAT = (RAT(I) + RAT(I-2))/2
        FYE(ICO) = LOG(AVRAT-1.)
        FYE(ICO-1) = LOG(AVRAT-1.)
        ENDIF
      FXE(ICE) = LOG(1.Q0/NSIZE(I))
      NSZE = NSIZE(I)
CC      WRITE(6,FMT115) NSIZE(I),FXE(ICE),FYE(ICE)
   61 CONTINUE
C
      IF(ICO.GT.1)CALL LSQF(ICO,FYO,FXO,AO,BO,SIGAO,SIGBO,XBARO,YBARO)
      IF(ICE.GT.1)CALL LSQF(ICE,FYE,FXE,AE,BE,SIGAE,SIGBE,XBARE,YBARE)
      IF(ICO.GT.1) GO TO 120
      AO = AE
      BO = 0.5
      SIGAO = SIGAE
      SIGBO = 0
      XBARO = FXO(1)
      YBARO = FYO(1)
  120 IF(ICE.GT.1) GO TO 121
      AE = AO
      BE = 0.5
      SIGAE = SIGAO
      SIGBE = 0
      XBARE = FXE(1)
      YBARE = FYE(1)
  121 I1 = 1
      IF(NSZE.GT.NSZO) I1 = 2
      I2 = 3 - I1
      N0(I1) = NSZO
      XN0(I1) = FXO(ICO)
      RN0(I1) = EXP(YBARO + BO*(XN0(I1)-XBARO)) + 1.
      POW(I1) = BO
      N0(I2) = NSZE
      XN0(I2) = FXE(ICE)
      RN0(I2) = EXP(YBARE + BE*(XN0(I2)-XBARE)) + 1.
      POW(I2) = BE
      CALL EXTRAP(N0,RN0,POW,FAC)
      DEX1 = -DIFF(ITER)*FAC
      EX1 = E(ITER) + DEX1
C
      RN0(I1) = EXP(YBARO + (BO+SIGBO)*(XN0(I1)-XBARO)) + 1.
      POW(I1) = BO + SIGBO
      RN0(I2) = EXP(YBARE + (BE+SIGBE)*(XN0(I2)-XBARE)) + 1.
      POW(I2) = BE + SIGBE
      CALL EXTRAP(N0,RN0,POW,FAC)
      DEX2 = -DIFF(ITER)*FAC
      EX2 = E(ITER) + DEX2
C
      RN0(I1) = EXP(YBARO + (BO-SIGBO)*(XN0(I1)-XBARO)) + 1.
      POW(I1) = BO - SIGBO
      RN0(I2) = EXP(YBARE + (BE-SIGBE)*(XN0(I2)-XBARE)) + 1.
      POW(I2) = BE - SIGBE
      CALL EXTRAP(N0,RN0,POW,FAC)
      DEX3 = -DIFF(ITER)*FAC
      EX3 = E(ITER) + DEX3
      ERROR = ABS((DEX2 - DEX3)/2.)
      IF(SIGBE.EQ.0.OR.SIGBO.EQ.0) ERROR = ABS(DEX1)
      ET = EX1 - Z**2/2 - (Z-1)**2/(2*N**2)
      OUT1 = ET*SCALE
      WRITE(*,122) 'Extp ',OUT1,ERROR*SCALE,-DIFF(ITER)/DEX1
      WRITE(2,122) 'Extp ',OUT1,ERROR*SCALE,-DIFF(ITER)/DEX1
  122 FORMAT(A5,F26.21,'+-',F22.21,F7.2)
      IF(LRYD) THEN
        WRITE(*,'(A,F24.20)') 'RYDBERG MASS SCALING FACTOR =',SCALE
        WRITE(2,'(A,F24.20)') 'RYDBERG MASS SCALING FACTOR =',SCALE
      ENDIF
      CLOSE(2,STATUS='KEEP')
      WRITE(3,'(''0/'')')
      WRITE(3,'(''0/'')')
      CLOSE(3,STATUS='KEEP')
      WRITE(8,201) NBET,KX1,KX2,AMASS,0.010,
     1  ',.false.,.true. , 1,.TRUE. ,0. ,0 ,0 1P states ||'
CCC      IF(AMASS.NE.0) GO TO 202
      WRITE(8,201) NBET,KX1,KX2,AMASS,0.010,
     1  ',.false.,.true. , 3,.false.,0. ,1 ,1 3D states ||'
      WRITE(8,201) NBET,KX1,KX2,AMASS,0.010,
     1  ',.false.,.true. , 3,.false.,0. ,1 ,0 1D states ||'
      WRITE(8,201) NBET+1,KX1,KX2,AMASS,0.010,
     1  ',.false.,.true. , 2,.false.,0. ,0 ,1 3F states ||'
  202 WRITE(8,'(A78)') '0/0/    THE END                                 
     1                            ||'
  201 FORMAT('0/',3(I1,','),F15.13,',',F5.3,A)
      CLOSE(8,STATUS='KEEP')
      WRITE(*,27) MATMAT,MATDAT,MATPOW
   27 FORMAT('OUTPUT FILES ARE',3(2X,A12))
      IF(FNWVA.EQ.'NEXT') READ(*,*)
      IF(FNWVA.EQ.'NEXT') GO TO 161
  162 STOP
      END
      SUBROUTINE LSQF(N,Y,X,A,B,SIGA,SIGB,XBAR,YBAR)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION Y(1),X(1)
      S = N
      SX = 0.
      SY = 0.
      SXX = 0.
      SXY = 0.
      DO 10 I=1,N
      SX = SX + X(I)
      SY = SY + Y(I)
      SXX = SXX + X(I)**2
      SXY = SXY + X(I)*Y(I)
   11 FORMAT(I4,4D15.7)
   10 CONTINUE
      DELTA = S*SXX - SX**2
      A = (SXX*SY - SX*SXY)/DELTA
      B = (S*SXY - SX*SY)/DELTA
      CHI2 = 0.
      DO 12 I=1,N
   12 CHI2 = CHI2 + (Y(I) - A - B*X(I))**2
      NN2 = 1
      IF(N.GT.3) NN2 = N - 2
      SIGA = SQRT(SXX/DELTA*CHI2/(NN2))
      SIGB = SQRT(S/DELTA*CHI2/(NN2))
C      XBAR = SIGN(SQRT(SXX/N),SX)
C      XBAR = SX/N
      XBAR = X(1)
      YBAR = A + B*XBAR
CC      WRITE(*,'(''A,B,SIGA,SIGB='',6F10.3)')A,B,SIGA,SIGB,XBAR,YBAR
CC      WRITE(2,'(''A,B,SIGA,SIGB='',6F10.3)')A,B,SIGA,SIGB,XBAR,YBAR
      RETURN
      END
      SUBROUTINE EXTRAP(N0,RN0,POW,FAC)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION N0(2),RN0(2),POW(2),C(2)
      C(1) = N0(1)**POW(1)*(RN0(1)-1)
      C(2) = N0(2)**POW(2)*(RN0(2)-1)
      SUM = 0.
      PROD = 1.
      DO 10 N=1,100
      RN = 1. + C(2)/(N0(2)+2*N)**POW(2)
      PROD = PROD/RN
      SUM = SUM + PROD
      RN = 1. + C(1)/(N0(1)+2*N)**POW(1)
      PROD = PROD/RN
      SUM = SUM + PROD
C      IF(10*(N/10).EQ.N) WRITE(*,*) N,RN,SUM
   10 CONTINUE
      FAC = SUM
c      WRITE(*,'(''N0,RN0,POW='',I4,2F15.7)') N0(2),RN0(2),POW(2)
C      WRITE(*,'(''N0,RN0,POW='',I4,2F15.7)') N0(1),RN0(1),POW(1)
C      WRITE(*,'(''FAC ='',F15.7)') FAC
C      WRITE(6,'(''FAC ='',F15.7)') FAC
      RETURN
      END
C      FUNCTION QEXP(X)
C      IMPLICIT REAL*16(A-H,O-Z)
C      QEXP = EXP(X)
C      RETURN
C      END

