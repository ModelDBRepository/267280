c 17 March 2021: adapt bask.f to simulate piriform multipolar
c cell.  Increase dendritic areas to allow for spines.
c Data from Tseng and Haberly, e.g. on periodic intrinsic bursts

c 20 March 2001, model of cortex basket cell, taken from nrt.f.
c Dendrites need to be made shorter and surface area reduced.
c May speed up c-current kinetics, reduce d(i), etc.

c 28 Dec. 2000, begin converting interneuron program to nRT cell.
c Soma will be comp. 1.  4 equivalent dendrites, each with 13 comps.
c (so 53 SD compartments).  Branching axon with 6 compartments - 59
c compartments in all.  Try one integration program for whole structure.
c Currents: leak, fast Na (naf), persistent Na (nap), fast DR (kdr),
c A-current (ka), K2 current, M-current (km), C current (kc), AHP
c (kahp), T-current (cat), high-thresh. Ca (CAL), h-current = anomalous
c rectifier (ar).
        PROGRAM multipolar 

       INTEGER J1, I, J, K, L, O, ISEED, K1
       REAL*8  TIMTOT, Z, Z1, Z2, curr(59), c(59), DT, time

c CINV is 1/C, i.e. inverse capacitance
       real*8 v(59), chi(59), cinv(59), mnaf(59), hnaf(59), mkdr(59),
     x mka(59),hka(59),mk2(59),hk2(59),mkm(59),mkc(59),mkahp(59),
     x mcat(59),hcat(59),mcal(59),mar(59),jacob(59,59),betchi(59),
     x gam(0:59,0:59),gL(59),gnaf(59),gnap(59),gkdr(59),gka(59),
     x gk2(59),gkm(59),gkc(59),gkahp(59),gcat(59),gcaL(59),gar(59),
     x gampa(59),gnmda(59),ggaba_a(59),ggaba_b(59),cafor(59),
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL,vk,vna,var,vca,vgaba_a
       real*8 tg1, tg2, tg3, tg4

        INTEGER NEIGH(59,5), NNUM(59)

c In initial version, setup from bask program, but should
c use one from, say, suppyrRS?
c      CALL   multipolar_SETUP
       CALL scort_setup_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL multipolarMAJ (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)

          do i = 1, 59
             cinv(i) = 1.d0 / c(i)
          end do

        MG = 1.0d0
C  IN MILLIMOLAR

        VL = -65.d0
c       VK =  -100.d0
        VK =  -85.d0
        VNA = 50.d0
        VCA = 125.d0
        VAR = -40.d0
        VGABA_A = -75.d0

c       DT = .00125d0
        DT = .00200d0

        TIMTOT = 2500.d0

c ? initialize membrane state variables?
        v = VL

        k1 = idnint (4.d0 * (v(1) + 120.d0))

      hnaf = alphah_naf(k1)/(alphah_naf(k1)+betah_naf(k1))
      hka = alphah_ka(k1)/(alphah_ka(k1)+betah_ka(k1))
      hk2 = alphah_k2(k1)/(alphah_k2(k1)+betah_k2(k1))
      hcat=alphah_cat(k1)/(alphah_cat(k1)+betah_cat(k1))

       O = 0
       TIME = 0.d0

2000      TIME = TIME + DT
          O = O + 1
          IF (TIME .GT. TIMTOT) GO TO 2001

          IF  (MOD(O,200).EQ. 0 )                      THEN
c         IF  ( i.eq.i          )                      THEN
            WRITE (6,904) TIME, v(1), v(11), v(27), v(40),
     X V(53), v(54), v(57), v(59)
          ENDIF
904       FORMAT(2X,F7.2,2X,8f7.2)

c         curr(1) =  0.05d0

c Current pulse
c         if ((time.le.150.d0).or.(time.gt.425.d0)) then
          if  (time.le.150.d0) then
           curr(1) = -0.000d0
          else
c          curr(1) = 0.35d0
           curr(1) = 0.50d0
          endif

c         if (time.le.100.d0) then
c            curr(1) = -0.4d0
c         else
c            curr(1) = 0.40d0
c         endif

c Brief periodic current pulses
c         if (mod(O,5000).le.500) then
c               curr(1) = 1.5d0
c         else
c               curr(1) = 0.d0
c         endif

c Triangular current
c         if (time.le.1000.d0) then
c      curr(1) = 0.65d0 * time / 1000.d0
c         else
c      curr(1) = 0.65d0 - (time-1000.d0)*0.65d0/1000.d0
c         endif

c Antidromic stimulation
c         if (mod(O,10000).le. 200) then
c               curr(57) = 0.4d0
c         else
c               curr(57) = 0.d0
c         endif

c Dendritic stimulation
c         if (mod(O,10000).le. 500) then
c               curr( 8) = 0.5d0
c         else
c               curr( 8) = 0.d0
c         endif

c Dendritic (longer) current pulse
c         if ((time.le.50.d0).or.(time.ge.100.d0)) then
c               curr(11) = 0.0d0
c         else
c               curr(11) = 0.02d0
c         endif

c GABA-A input
             goto 909
              tg1 =  80.d0
              tg2 =  83.d0
              tg3 =  86.d0
              tg4 =  89.d0
                if (time.le.tg1) then
                  ggaba_a(8) = 0.d0
                else
                  z = time - tg1
        ggaba_a(8) = 0.85d-3 * (0.56d0 * dexp( - z / 18.d0)
     X                        +0.44d0 * dexp (-z / 89.d0))
                endif

               if (time.gt.tg2) then
                 z = time - tg2
      ggaba_a(8)=ggaba_a(8)+0.85d-3 * (0.56d0 * dexp( - z / 18.d0)
     X                        +0.44d0 * dexp (-z / 89.d0))
               endif

               if (time.gt.tg3) then
                 z = time - tg3
      ggaba_a(8)=ggaba_a(8)+0.85d-3 * (0.56d0 * dexp( - z / 18.d0)
     X                        +0.44d0 * dexp (-z / 89.d0))
               endif

               if (time.gt.tg4) then
                 z = time - tg4
      ggaba_a(8)=ggaba_a(8)+0.85d-3 * (0.56d0 * dexp( - z / 18.d0)
     X                        +0.44d0 * dexp (-z / 89.d0))
               endif

               ggaba_a(21) = ggaba_a(8)
               ggaba_a(34) = ggaba_a(8)
               ggaba_a(47) = ggaba_a(8)

               ggaba_a( 5) = ggaba_a(8)
               ggaba_a(18) = ggaba_a(8)
               ggaba_a(31) = ggaba_a(8)
               ggaba_a(44) = ggaba_a(8)
909                  CONTINUE

c                 gnaf = 0.d0
c                 gnap = 0.d0
c                 gkdr = 0.d0
c                 gka = 0.d0
                  gk2 = 0.d0
c                 gkm = 0.d0
c                 gkc = 0.d0
c                 gkahp = 0.d0
c                 gcat = 0.d0
c                 gcaL = 0.d0
c                 gar = 0.d0

c below 4 statements disconnect axon from soma.
c Note that axon is rather "leaky"
c            gam(1,54) = 0.d0
c            gam(54,1) = 0.d0
c            jacob(1,54) = 0.d0
c            jacob(54,1) = 0.d0

        CALL   multipolarINT (V,CHI,CINV,mnaf,hnaf,mkdr,mka,hka,mk2,hk2,
     x mkm,mkc,mkahp,mcat,hcat,mcal,mar,dt,neigh,nnum,jacob,mg,
     x vL,vk,vna,var,vca,vgaba_a,betchi,gam,gL,gnaf,gnap,gkdr,gka,
     x gk2,gkm,gkc,gkahp,gcat,gcaL,gar,gampa,gnmda,ggaba_a,
     x ggaba_b,O,time,
     X    alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,
     X    cafor,curr)


              GO TO 2000
2001          CONTINUE


1000    CONTINUE
        END


C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END

C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE multipolar_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = I
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .0125d0 + .1525d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from interneuron data, Martina & Jonas 1997, tau x 0.5
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 58.3d0)/6.7d0))
        tauh = 0.225d0 + 1.125d0/(1.d0+dexp((v+37.d0)/15.d0))
c from interneuron data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-27.d0)/11.5d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif
c Speed-up of C kinetics here.
          alpham_kc(i) = 2.d0 * alpham_kc(i)
           betam_kc(i) = 2.d0 *  betam_kc(i)

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe et al., 1996, pg. 170
         minf = 1.d0/(1.d0 + dexp((-v-52.d0)/7.4d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/5.d0))
         taum = 1.d0 + .33d0/(dexp((v+27.d0)/10.d0) +
     x                  dexp((-v-102.d0)/15.d0))
         tauh = 28.3d0 +.33d0/(dexp((v+48.d0)/4.d0) +
     x                     dexp((-v-407.d0)/50.d0))
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do 2, i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
2      CONTINUE
       END

C       DEBUG SUBCHK
C                   END DEBUG
        SUBROUTINE multipolarMAJ
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 4 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

        REAL*8 C(59),GL(59),GAM(0:59,0:59),GNAF(59),GCAT(59)
        REAL*8 GKDR(59),GKA(59),GKC(59),GKAHP(59),GCAL(59)
        REAL*8 GK2(59),GKM(59),GNAP(59),GAR(59)
        REAL*8 JACOB(59,59),RI_SD,RI_AXON,RM_SD,RM_AXON,CDENS
        INTEGER LEVEL(59)
        REAL*8 GNAF_DENS(0:9), GCAT_DENS(0:9), GKDR_DENS(0:9)
        REAL*8 GKA_DENS(0:9), GKC_DENS(0:9), GKAHP_DENS(0:9)
        REAL*8 GCAL_DENS(0:9), GK2_DENS(0:9), GKM_DENS(0:9)
        REAL*8 GNAP_DENS(0:9), GAR_DENS(0:9)
        REAL*8 RES, RINPUT
        REAL*8 RSOMA, PI, BETCHI(59), CAFOR(59)
        REAL*8 RAD(59), LEN(59), GAM1, GAM2, ELEN(59)
        REAL*8 RIN, D(59), AREA(59), RI
        INTEGER NEIGH(59,5), NNUM(59)
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS

c       RI_SD = 200.d0
        RI_SD = 250.d0
        RM_SD = 50000.d0
c       RM_SD = 25000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

        gnaf_dens(0) = 400.d0
        gnaf_dens(1) =  60.d0
        gnaf_dens(2) =  60.d0
        gnaf_dens(3) =  60.d0
        do i = 4, 9
          gnaf_dens(i) = 30.d0
c         gnaf_dens(i) = 10.d0
        end do

        gkdr_dens(0) = 400.d0
        gkdr_dens(1) = 100.d0
        gkdr_dens(2) = 100.d0
        gkdr_dens(3) = 100.d0
        do i = 4, 9
         gkdr_dens(i) = 20.d0
c        gkdr_dens(i) = 60.d0
        end do

        do i = 1, 9
c         gnap_dens(i) = 0.005d0 * gnaf_dens(i)
          gnap_dens(i) = 0.500d0 * gnaf_dens(i) ! as piriformENDO100
        end do

        do i = 1, 3
          gcat_dens(i) = 0.05d0
        end do
        do i = 4, 9
          gcat_dens(i) = 0.5d0
        end do

        do i = 1, 3
          gcal_dens(i) = 0.5d0
c         gcal_dens(i) = 0.1d0
        end do
        do i = 4, 9
          gcal_dens(i) = 2.5d0
c         gcal_dens(i) = 0.2d0
        end do

        gka_dens(0) = 1.d0
        gka_dens(1) =  2.d0
        gka_dens(2) =  1.d0
        gka_dens(3) =  1.d0
        do i = 4, 9
         gka_dens(i) = 1.0d0
        end do

        do i = 1, 9
c        gkc_dens(i) = 10.00d0
c        gkc_dens(i) = 25.00d0
!        gkc_dens(i) =  5.00d0
         gkc_dens(i) =  0.00d0
        end do
        gkc_dens(1) = 10.d0

        gkm_dens(0) = 8.00d0
        do i = 1, 9
c        gkm_dens(i) = 0.50d0
         gkm_dens(i) = 6.00d0
        end do

        gk2_dens(0) = .0d0
        do i = 1, 9
         gk2_dens(i) = 0.00d0
        end do

        do i = 1, 9
c        gkahp_dens(i) = 0.10d0
c        gkahp_dens(i) = 0.15d0
         gkahp_dens(i) = 0.12d0
        end do

        do i = 1, 9
c        gar_dens(i) = 0.025d0
         gar_dens(i) = 0.02d0
        end do

        WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(L)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 9
          WRITE (6,9990) I, gnaf_dens(i), gcaL_dens(i), gkdr_dens(i),
     X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 41, 13
         level(i) = 2
        end do
        do i = 3, 42, 13
           level(i) = 3
           level(i+1) = 3
        end do
        do i = 5, 44, 13
           level(i) = 4
           level(i+1) = 4
           level(i+2) = 4
        end do
        do i = 8, 47, 13
           level(i) = 5
           level(i+1) = 5
           level(i+2) = 5
        end do
        do i = 11, 50, 13
           level(i) = 6
           level(i+1) = 7
           level(i+2) = 8
           level(i+3) = 9
        end do

        do i = 54, 59
         level(i) = 0
        end do

c connectivity of axon
        nnum(54) = 2
        nnum(55) = 3
        nnum(56) = 3
        nnum(58) = 3
        nnum(57) = 1
        nnum(59) = 1
         neigh(54,1) =  1
         neigh(54,2) = 55
         neigh(55,1) = 54
         neigh(55,2) = 56
         neigh(55,3) = 58
         neigh(56,1) = 55
         neigh(56,2) = 57
         neigh(56,3) = 58
         neigh(58,1) = 55
         neigh(58,2) = 56
         neigh(58,3) = 59
         neigh(57,1) = 56
         neigh(59,1) = 58

c connectivity of SD part
          nnum(1) = 5
          neigh(1,1) = 54
          neigh(1,2) =  2
          neigh(1,3) = 15
          neigh(1,4) = 28
          neigh(1,5) = 41

          do i = 2, 41, 13
           nnum(i) = 3
           neigh(i,1) = 1
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
          end do

          do i = 3, 42, 13
           nnum(i) = 4
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
           neigh(i,4) = i + 3
          end do

          do i = 4, 43, 13
           nnum(i) = 3
           neigh(i,1) = i - 2
           neigh(i,2) = i - 1
           neigh(i,3) = i + 3
          end do

          do i = 5, 44, 13
           nnum(i) = 3
           neigh(i,1) = i - 2
           neigh(i,2) = i + 1
           neigh(i,3) = i + 3
          end do

          do i = 6, 45, 13
           nnum(i) = 3
            neigh(i,1) = i - 3
            neigh(i,2) = i - 1
            neigh(i,3) = i + 3
          end do

          do i = 7, 46, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 3
          end do

          do i = 8, 47, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 3
          end do

          do i = 9, 48, 13
           nnum(i) = 1
           neigh(i,1) = i - 3
          end do

          do i = 10, 49, 13
           nnum(i) = 1
           neigh(i,1) = i - 3
          end do

          do i = 11, 50, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 1
          end do

          do i = 12, 51, 13
           nnum(i) = 2
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
          end do

          do i = 13, 52, 13
           nnum(i) = 2
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
          end do

          do i = 14, 53, 13
           nnum(i) = 1
           neigh(i,1) = i - 1
          end do

         DO 332, I = 1, 59
           WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
     X NEIGH(I,5)
3330     FORMAT(2X,I5,I5,I5,I5,I5,I5)
332      CONTINUE
          DO 858, I = 1, 59
           DO 858, J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
            DO 859, L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
859         CONTINUE
             IF (IT.EQ.0) THEN
              WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
             ENDIF
858       CONTINUE

c length and radius of axonal compartments
          do i = 54, 59
            len(i) = 50.d0
          end do
c         rad(54) = 0.80d0
c         rad(55) = 0.7d0
          rad(54) = 0.70d0
          rad(55) = 0.6d0
          do i = 56, 59
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 20.d0
          rad(1) = 7.5d0

          do i = 2, 53
c          len(i) = 40.d0
           len(i) = 80.d0
          end do

          rad(2) =   1.06d0
          rad(3) =   rad(2) / 1.59d0
          rad(4) =   rad(2) / 1.59d0
          rad(5) =   rad(2) / 2.53d0
          rad(6) =   rad(2) / 2.53d0
          rad(7) =   rad(2) / 1.59d0
          rad(8) =   rad(2) / 2.53d0
          rad(9) =   rad(2) / 2.53d0
          rad(10) =  rad(2) / 1.59d0
          rad(11) =  rad(2) / 2.53d0
          rad(12) =  rad(2) / 2.53d0
          rad(13) =  rad(2) / 2.53d0
          rad(14) =  rad(2) / 2.53d0

          do i = 15, 53
           rad(i) = rad(i-13)
          end do

        WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
        DO 920, I = 1, 59
920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 59
            if (level(i).le.1) then
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
            else
          AREA(I) = 4.d0 * PI * RAD(I) * LEN(I)
            endif
C NO CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
C IN ORIGINAL bask.f, but change that here
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
         DO 1019, I = 2, 53
           Z = Z + AREA(I)
1019     CONTINUE
         WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

        DO 140, I = 1, 59
        DO 140, K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )

         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
140     CONTINUE
c gam computed in microS

        DO 299, I = 1, 59
299       BETCHI(I) = .05d0
        BETCHI( 1) =  .02d0

        DO 300, I = 1, 59
c300     D(I) = 2.D-4
300     D(I) = 1.D-4
        DO 301, I = 1, 59
c        IF (LEVEL(I).EQ.1) D(I) = 5.D-3
         IF (LEVEL(I).EQ.1) D(I) = 2.D-4
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


       DO 160, I = 1, 59
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
C     NOTE CORRECTION

        do 200, i = 1, 59
200     C(I) = 1000.d0 * C(I)
C     TO GO FROM MICROF TO NF.

      DO 909, I = 1, 59
       JACOB(I,I) = - GL(I)
      DO 909, J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
             WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, 59
          do j = 1, 59
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

       DO 500, I = 1, 59
        WRITE (6,501) I,C(I)
501     FORMAT(1X,I2,' C(I) = ',F7.4)
500     CONTINUE
        END


      SUBROUTINE multipolarINT(V,CHI,CINV,mnaf,hnaf,mkdr,mka,hka,mk2,hk2
     x ,mkm,mkc,mkahp,mcat,hcat,mcal,mar,dt,neigh,nnum,jacob,mg,
     x vL,vk,vna,var,vca,vgaba_a,betchi,gam,gL,gnaf,gnap,gkdr,gka,
     x gk2,gkm,gkc,gkahp,gcat,gcaL,gar,gampa,gnmda,ggaba_a,
     x ggaba_b,O,time,
     X    alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,
     X    cafor,curr)

c CINV is 1/C, i.e. inverse capacitance
       real*8 v(59), chi(59), cinv(59), mnaf(59), hnaf(59), mkdr(59),
     x mka(59),hka(59),mk2(59),hk2(59),mkm(59),mkc(59),mkahp(59),
     x mcat(59),hcat(59),mcal(59),mar(59),jacob(59,59),betchi(59),
     x gam(0:59,0:59),gL(59),gnaf(59),gnap(59),gkdr(59),gka(59),
     x gk2(59),gkm(59),gkc(59),gkahp(59),gcat(59),gcaL(59),gar(59),
     x gampa(59),gnmda(59),ggaba_a(59),ggaba_b(59),cafor(59),
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)

        real*8 fastna_shift

c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(59), fchi(59),fmnaf(59),fhnaf(59),fmkdr(59),
     x fmka(59),fhka(59),fmk2(59),fhk2(59),
     x fmkm(59),fmkc(59),fmkahp(59),
     x fmcat(59),fhcat(59),fmcal(59),fmar(59)

c below are for calculating the partial derivatives
       real*8 dfv_dv(59,59), dfv_dchi(59), dfv_dmnaf(59),
     x  dfv_dhnaf(59),dfv_dmkdr(59),dfv_dmka(59),dfv_dhka(59),
     x  dfv_dmk2(59),dfv_dhk2(59),dfv_dmkm(59),dfv_dmkc(59),
     xdfv_dmkahp(59),dfv_dmcat(59),dfv_dhcat(59),dfv_dmcal(59),
     x  dfv_dmar(59)

        real*8 dfchi_dv(59), dfchi_dchi(59),
     x dfmnaf_dmnaf(59), dfmnaf_dv(59),dfhnaf_dhnaf(59),
     x dfhnaf_dv(59),dfmkdr_dmkdr(59),dfmkdr_dv(59),
     x dfmka_dmka(59),dfmka_dv(59),dfhka_dhka(59),dfhka_dv(59),
     x dfmk2_dmk2(59),dfmk2_dv(59),dfhk2_dhk2(59),dfhk2_dv(59),
     x dfmkm_dmkm(59),dfmkm_dv(59),dfmkc_dmkc(59),dfmkc_dv(59),
     x dfmcat_dmcat(59),dfmcat_dv(59),dfhcat_dhcat(59),
     x dfhcat_dv(59),dfmcal_dmcal(59),dfmcal_dv(59),
     x dfmar_dmar(59),dfmar_dv(59),dfmkahp_dchi(59),
     x dfmkahp_dmkahp(59), dt2, outrcd(20), time

         REAL*8 dt,mg,vL,vk,vna,var,vca,vgaba_a,curr(59),Z,Z1,Z2
         INTEGER O, K0, K1, NEIGH(59,5), NNUM(59)
       REAL*8 OPEN(59),gamma(59),gamma_prime(59)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(59), alpham_ahp_prime(59)
       REAL*8 gna_tot(59),gk_tot(59),gca_tot(59),gar_tot(59)
       REAL*8 gca_high(59)
c this will be gCa conductance corresponding to high-thresh channels

       DO 301, I = 1, 59
          FV(I) = -GL(I) * (V(I) - VL) * cinv(i)
          DO 302, J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K) - V(I)) * cinv(i)
301    CONTINUE


c       CALL FNMDA (V, OPEN, MG)

      DO 421, I = 1, 59
421    FV(I) = FV(I) + ( CURR(I)
     X   - ggaba_b(I)*(V(I)-VK)
     X   - (gampa(I) + open(i) * gnmda(I))*V(I)
     X   - ggaba_a(I)*(V(I)-Vgaba_a) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV

       do i = 1, 59
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i))
        if (chi(i).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
       end do

      DO 88, I = 1, 59
       gna_tot(i) = gnaf(i) * (mnaf(i)**3) * hnaf(i) +
     x     gnap(i) * (mnaf(i)**3)
       gk_tot(i) = gkdr(i) * (mkdr(i)**4) +
     x             gka(i)  * (mka(i)**4) * hka(i) +
     x             gk2(i)  * mk2(i) * hk2(i) +
     x             gkm(i)  * mkm(i) +
     x             gkc(i)  * mkc(i) * gamma(i) +
     x             gkahp(i)* mkahp(i)
       gca_tot(i) = gcat(i) * (mcat(i)**2) * hcat(i) +
     x              gcaL(i) * (mcaL(i)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i)**2)
       gar_tot(i) = gar(i) * mar(i)


88     FV(I) = FV(I) - ( gna_tot(i) * (v(i) - vna)
     X  + gk_tot(i) * (v(i) - vK)
     X  + gca_tot(i) * (v(i) - vCa)
     X  + gar_tot(i) * (v(i) - var) ) * cinv(i)

         do i = 1, 59
         do j = 1, 59
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_b(I) + ggaba_a(i) + gampa(i)
     X   + open(i) * gnmda(I) )
          endif
         end do
         end do

          do i = 1, 59
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i) *
     x                     gamma_prime(i) * (v(i)-vK)
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i)**2) *
     X    (gnaf(i) * hnaf(i) + gnap(i)) * (v(i) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i)**3) *
     X                    (v(i) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i)**3)
     X                   * (v(i) - vK)
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i)**3) *
     X                   hka(i) * (v(i) - vK)
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i)**4) *
     X                    (v(i) - vK)
        dfv_dmk2(i)  = - cinv(i) * gk2(i) * hk2(i) * (v(i)-vK)
        dfv_dhk2(i)  = - cinv(i) * gk2(i) * mk2(i) * (v(i)-vK)
        dfv_dmkm(i)  = - cinv(i) * gkm(i) * (v(i) - vK)
        dfv_dmkc(i)  = - cinv(i) * gkc(i) * gamma(i) * (v(i)-vK)
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i) - vK)
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i) *
     X                    hcat(i) * (v(i) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i)**2) *
     X                  (v(i) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i) *
     X                      (v(i) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i) - var)
          end do

         do i = 1, 59
c         fchi(i) = - cafor(i) * gca_tot(i) * (v(i) - vca)
c    x       - betchi(i) * chi(i)
c         dfchi_dv(i) = - cafor(i) * gca_tot(i)
c         dfchi_dchi(i) = - betchi(i)
c Note 2 possibilities: chi can depend on total ICa, or only on
c  ICa through high-thresh. channels
          fchi(i) = - cafor(i) * gca_high(i) * (v(i) - vca)
     x       - betchi(i) * chi(i)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, 59
        alpham_ahp(i) = dmin1(0.2d-4 * chi(i),0.01d0)
        if (chi(i).le.500.d0) then
          alpham_ahp_prime(i) = 0.2d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, 59
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i))
     x                  -.001d0 * mkahp(i)
c    x                  -.0025d0 * mkahp(i) ! NOTE
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i))
       end do

          do i = 1, 59

       K1 = IDNINT ( 4.d0 * (V(I) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

             fastNa_shift = -2.5d0
       K0 = IDNINT ( 4.d0 * (V(I)+      fastNa_shift+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0


        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i)) -
     X              betam_naf(k0) * mnaf(i)
        fhnaf(i) = alphah_naf(k1) * (1.d0 - hnaf(i)) -
     X              betah_naf(k1) * hnaf(i)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i)) -
     X              betam_kdr(k1) * mkdr(i)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i)) -
     X              betam_ka (k1) * mka(i)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i)) -
     X              betah_ka (k1) * hka(i)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i)) -
     X              betam_k2 (k1) * mk2(i)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i)) -
     X              betah_k2 (k1) * hk2(i)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i)) -
     X              betam_km (k1) * mkm(i)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i)) -
     X              betam_kc (k1) * mkc(i)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i)) -
     X              betam_cat(k1) * mcat(i)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i)) -
     X              betah_cat(k1) * hcat(i)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i)) -
     X              betam_caL(k1) * mcaL(i)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i)) -
     X              betam_ar (k1) * mar(i)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i)) -
     X                  dbetam_naf(k0) * mnaf(i)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i)) -
     X                  dbetah_naf(k1) * hnaf(i)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i)) -
     X                  dbetam_kdr(k1) * mkdr(i)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i)) -
     X                  dbetam_ka(k1) * mka(i)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i)) -
     X                  dbetah_ka(k1) * hka(i)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i)) -
     X                  dbetam_k2(k1) * mk2(i)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i)) -
     X                  dbetah_k2(k1) * hk2(i)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i)) -
     X                  dbetam_km(k1) * mkm(i)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i)) -
     X                  dbetam_kc(k1) * mkc(i)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i)) -
     X                  dbetam_cat(k1) * mcat(i)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i)) -
     X                  dbetah_cat(k1) * hcat(i)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i)) -
     X                  dbetam_caL(k1) * mcaL(i)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i)) -
     X                  dbetam_ar(k1) * mar(i)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, 59
          v(i) = v(i) + dt * fv(i)
           do j = 1, 59
        v(i) = v(i) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i) = v(i) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i) = chi(i) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i) = mnaf(i) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        hnaf(i) = hnaf(i) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i) = mkdr(i) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i) =  mka(i) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i) =  hka(i) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i) =  mk2(i) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i) =  hk2(i) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i) =  mkm(i) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i) =  mkc(i) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i) = mkahp(i) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i) =  mcat(i) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i) =  hcat(i) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i) =  mcaL(i) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i) =   mar(i) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
         end do


c      IF (MOD(O,75).EQ.0) THEN
       IF (MOD(O,150).EQ.0) THEN
           OUTRCD(1) = TIME
           OUTRCD(2) = v(1)
           outrcd(3) = v(8)
           outrcd(4) = v(54)
           outrcd(5) = v(57)
           outrcd(6) = chi(1)
           outrcd(7) = chi(8)
           outrcd(8) = gk_tot(1) * (v(1) - vK)
       outrcd(9) = gkdr(1) * (mkdr(1)**4) * (v(1)-vK)
       outrcd(10)= gka(1)  * (mka(1)**4) * hka(1)*(v(1)-vK)
       outrcd(11)= gk2(1)  * mk2(1) * hk2(1) * (V(1)-vK)
       outrcd(12)= gkm(1)  * mkm(1) * (v(1)-vK)
       outrcd(13)= gkc(1)  * mkc(1) * gamma(1) * (v(1)-vK)
       outrcd(14)= gkahp(1)* mkahp(1) * (v(1) - vK)
c          outrcd(9) = gca_tot(1) * (v(1) - vca)
c          outrcd(10) = gna_tot(1) * (v(1) - vna)
c          outrcd(11) = gna_tot(4) * (v(4) - vna)
c          outrcd(12) = v(3)
c          outrcd(13) = v(4)
c          outrcd(14) = v(55)
           outrcd(15) = v(56)
           outrcd(16) = gca_tot(8) * (v(8) - vca)
           outrcd(17) = curr(1)
           outrcd(18) = 1.d3 * ggaba_a(8)
c        OPEN(11,FILE='BASK1.BU')
         OPEN(11,FILE='multipolar25.test')
         WRITE (11,FMT='(20F10.3)') (OUTRCD(I),I=1,20)
C        WRITE (11,900) OUTRCD
C900      FORMAT (100A4)
       ENDIF


              END

               SUBROUTINE FNMDA (VSTOR,OPEN,MG)
               REAL*8 VSTOR( 59), OPEN(59)
               REAL*8 A, BB1, BB2, VM, A1, A2, B1, B2, MG
c modify so that potential is absolute and not relative to
c  "rest"
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)
C  TO DETERMINE VOLTAGE-DEPENDENCE OF NMDA CHANNELS
           DO 1, I = 1, 59
c          VM = VSTOR(I) - 60.
           VM = VSTOR(I)
           A1 = DEXP(-.016d0*VM - 2.91d0)
           A2 = 1000.d0 * MG * DEXP (-.045d0 * VM - 6.97d0)
           B1 = DEXP(.009d0*VM + 1.22d0)
           B2 = DEXP(.017d0*VM + 0.96d0)
        OPEN(I)     = 1.d0/(1.d0 + (A1+A2)*(A1*BB1 + A2*BB2) /
     X   (A*A1*(B1+BB1) + A*A2*(B2+BB2))  )
C  FROM JAHR & STEVENS, EQ. 4A
C               DO 124, J = 1, 19
C          OPEN(J) = 1./(1.+.667* EXP(-0.07*(VSTOR(J)-60.)) )
C  FROM CHUCK STEVENS
1               CONTINUE
                        END

