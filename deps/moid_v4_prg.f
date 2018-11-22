c ==================================================================================================================================
c current version: 4.0
c version 1.0(July 31th, 2013)
c version 2.0(August 23th,2013)   works faster, some unnecessary caclulations removed
c version 3.0(September 5th,2014) handles extremely eccentric orbits, corrections
c                                 done thanks to comments submitted by Zong-Fu,Sie NCU Taiwan
c version 4.0(July 30th,2017)     removes important bug injected in version 3.0,
c                                 corrections done thanks to comments submitted by
c                                 Robert Jedicke,University of Hawaii
c--------------------------------
c This program calculates the MOID between two asteroids - Lutetia and Magdalena.
c It uses all ideas and solutions described with details in the paper by
c T.Wisniowski and H.Rickman "A Fast, Geometric Method for Calculating Accurate
c Minimum Orbit Intersection Distances (MOIDs)" published in 2013 in Acta Astronomica.
c The program is free and may be used without limits as the core of any other program.
c The authors will appreciate for mentioning them, when the program will appear to be useful.
c ==================================================================================================================================
      implicit real*8 (a-h,o-z), integer*4 (i-n)
      real*8 longit,longit_m,longit_o,moid,incliA,incliB
      logical aleft,aright,bleft,bright,calc1,calc2,calc3,calc4
c tables
      real*8, dimension(3)::rAt,rBt,Axt,Ayt,Bxt,Byt,Bzt
      real*8, dimension(10)::tmpmoid,tmptrueB,tmplongit

c....parameters of the program
ccc  (these steps are optimized with respect to speed and reliability)
      cstep=0.12d0    ! scanning step of true anomaly/longitude in radians
      stepini=0.07d0  !initial step of first tuning in radians
      steptresh=1E-5  !final step of first tunig (to choose the MOID) in radians

ccc   (this step depends on expected final accuracy)
      stepmin=1E-14   !threshold step of second tuning in radians

c....constants and given values
      pi=3.141592653589793d0
      twopi=2d0*pi
      degrad=pi/180d0
c....orbital parameters of body A - asteroid Magdalena, angles in [deg]
      saxisA=3.1924186d0  !semiaxis [AU]
      eccenA=0.0843230d0  !eccenticity [deg]
      argpeA=298.98851d0  !argument of perihelion [deg]
      omegaA=161.64831d0  !longitude of ascending node [deg]
      incliA=10.64148d0   !inclination [deg]
c....orbital parameters of body B - asteroid Lutetia, angles in [deg]
      saxisB=2.4354656d0  !semiaxis [AU]
      eccenB=0.1629385d0  !eccenticity [deg]
      argpeB=250.36793d0  !argument of perihelion [deg]
      omegaB=80.89426d0  !longitude of ascending node [deg]
      incliB=3.06405d0   !inclination [deg]

ccccccccc START OF PREPARING THE ORBITS ccccccccccccccc
c....converting angles to [rad]
      argpeA=argpeA*degrad
      omegaA=omegaA*degrad
      incliA=incliA*degrad
      argpeB=argpeB*degrad
      omegaB=omegaB*degrad
      incliB=incliB*degrad
c....computing parameters of transition matrix c11...c33
      c11=dcos(omegaA)*dcos(argpeA)-dsin(omegaA)*dcos(incliA)
     a *dsin(argpeA)
      c12=dsin(omegaA)*dcos(argpeA)+dcos(omegaA)*dcos(incliA)
     a *dsin(argpeA)
      c13=dsin(incliA)*dsin(argpeA)
      c21=-dcos(omegaA)*dsin(argpeA)-dsin(omegaA)*dcos(incliA)
     a *dcos(argpeA)
      c22=-dsin(omegaA)*dsin(argpeA)+dcos(omegaA)*dcos(incliA)
     a *dcos(argpeA)
      c23=dsin(incliA)*dcos(argpeA)
      c31=dsin(incliA)*dsin(omegaA)
      c32=-dsin(incliA)*dcos(omegaA)
      c33=dcos(incliA)
c....calculating new values of Euler angles using transition matrix
        sintmpi=dsin(incliB)
        costmpi=dcos(incliB)
        costmpo=dcos(omegaB)
        sintmpo=dsin(omegaB)
        costmpa=dcos(argpeB)
        sintmpa=dsin(argpeB)
        x1=costmpo*costmpa-sintmpo*costmpi*sintmpa
        x2=sintmpo*costmpa+costmpo*costmpi*sintmpa
        x3=sintmpi*sintmpa
        y1=-costmpo*sintmpa-sintmpo*costmpi*costmpa
        y2=-sintmpo*sintmpa+costmpo*costmpi*costmpa
        y3=sintmpi*costmpa
        z1=sintmpi*sintmpo
        z2=-sintmpi*costmpo
        z3=costmpi
        z1n=c11*z1+c12*z2+c13*z3
        z2n=c21*z1+c22*z2+c23*z3
        z3n=c31*z1+c32*z2+c33*z3
        y3n=c31*y1+c32*y2+c33*y3
        x3n=c31*x1+c32*x2+c33*x3
        incliB=datan2(dsqrt(z1n*z1n+z2n*z2n),z3n)
        omegaB=-datan2(z1n,-z2n)
        argpeB=-datan2(x3n,y3n)
c....helpful precalculated values
        costmpo=dcos(omegaB)
        sintmpo=dsin(omegaB)
        sintmpi=dsin(incliB) !=z1n/sintmpo
        costmpi=z3n          !=dcos(incliB)
        sint=sintmpo*costmpi
        cost=costmpo*costmpi
        radA=saxisA*(1d0-eccenA*eccenA)
        radB=saxisB*(1d0-eccenB*eccenB)
ccccccccc END OF PREPARING THE ORBITS ccccccccccccccc

ccccccccc START OF SCANNING ccccccccccccccccccccccccc
ccc This tool yields a preliminary approach to the minima of the distance function.
ccc By scanning one full revolution of meridional plane we look for the local minima.
c......initial parameters
        trueB=-2d0*cstep
        moid=1E6;dist_o=1E6  !something big
        tmpmoid(1)=1E6;tmpmoid(2)=1E6
        tmpmoid(3)=1E6;tmpmoid(4)=1E6
        iii1=0;jjj1=0
c.....Looking for the minima with rotating meridional plane
c.......a)at first we calculate the coordinates of two additional positions of the plane to create first triplet
          do 307 iii=1,2
            rB=radB/(1d0+eccenB*dcos(trueB))!compute the radius for B
            sintmp=dsin(trueB+argpeB)
            costmp=dcos(trueB+argpeB)
            Bz_sq=sintmpi*sintmp
            Bz_sq=Bz_sq*Bz_sq !square of Z-coordinate for B
            longit=datan2(sintmpo*costmp+sintmp*cost,
     a       costmpo*costmp-sintmp*sint)  !compute the longitude for A
            tmp2=eccenA*dcos(longit)      !temporary value
            rA=radA/(1d0+tmp2)            !compute the radius for A (two possibilities)
            rA2=radA/(1d0-tmp2)
            tmp1=rB*dsqrt(1d0-Bz_sq)      !temporary value
            if (dabs(tmp1-rA)>dabs(tmp1+rA2)) then
              rA=rA2
              longit=longit-pi            !the second possibility gives smaller distance
              tmp1=tmp1+rA2
            else
              tmp1=tmp1-rA
            endif
            dist=rB*rB*Bz_sq+tmp1*tmp1 !square of the distance A-B
            if (iii.eq.1) then
              dist_oo=dist
            else
              dist_o=dist
              trueB_o=trueB
              longit_o=longit
            endif
            trueB=trueB+cstep
  307     continue
c.......b)now we scan doing one full revolution of the meridional plane
          nmax=0 !counts the minima
          dist_min=dist
          do while (trueB<(twopi+cstep))  !loop for true anomaly of B
            rB=radB/(1d0+eccenB*dcos(trueB))!compute the radius for B
            sintmp=dsin(trueB+argpeB)
            costmp=dcos(trueB+argpeB)
            Bz_sq=sintmpi*sintmp
            Bz_sq=Bz_sq*Bz_sq !square of Z-coordinate for B
            longit=datan2(sintmpo*costmp+sintmp*cost,
     a       costmpo*costmp-sintmp*sint) !compute the longitude for A
            tmp2=eccenA*dcos(longit)      !temporary value
            rA=radA/(1d0+tmp2)            !compute the radius for A (two possibilities)
            rA2=radA/(1d0-tmp2)
            tmp1=rB*dsqrt(1d0-Bz_sq)      !temporary value
            if (dabs(tmp1-rA)>dabs(tmp1+rA2)) then
              rA=rA2
              longit=longit-pi            !the second possibility gives smaller distance
              tmp1=tmp1+rA2
            else
              tmp1=tmp1-rA
            endif
            dist=rB*rB*Bz_sq+tmp1*tmp1 !square of the distance A-B
           if ((dist_o<=dist).and.(dist_o<=dist_oo)) then  !the minimum was found
              nmax=nmax+1
              tmptrueB(nmax)=trueB_o
              tmplongit(nmax)=longit_o
              tmpmoid(nmax)=dist_o
            endif
            if (dist_min>dist) dist_min=dist
            dist_oo=dist_o
            trueB_o=trueB
            longit_o=longit
            dist_o=dist
            trueB=trueB+cstep
         end do
ccccccccc END OF SCANNING ccccccccccccccccccccccccc

ccc "WATER" PROCEDURE cccccccccccccccccccccccccccccccccc
ccc In case only one minimum is detected we take a special care
ccc to avoid the risk of missing the minima.
ccc Instead of starting the tuning with one detected minimum,
ccc we start it with four positions of the meridional plane, evenly
ccc distributed along the inclined orbit. The later tuning procedure
ccc moves the points along the orbits similarly as water droplets
ccc are leading by the gravity to the points of minimum height � so
ccc we called this phase �water procedure�. It slightly slows down
ccc the calculations, but minimizes the risk of missing the MOID.
ccc With "water procedure" the speed is 9-12 seconds per 100,000 MOIDs, risk of missing <1E-6
ccc Without "water procedure" the speed is <9 seconds per 100,000 MOIDs, risk of missing about 3E-5
ccc (speed measured on fast single CPU core)

       !goto 510   !if this jump is active - the water procedure is switched off

  405   if (nmax<2) then  !only one minimum was detected
         nmax=4
         do 407 iii=1,4
            tmptrueB(iii)=(.25+.5*iii)*pi  !evenly distributed points
            sintmp=dsin(tmptrueB(iii)+argpeB)
            costmp=dcos(tmptrueB(iii)+argpeB)
            tmplongit(iii)=datan2(sintmpo*costmp+sintmp*cost,
     a       costmpo*costmp-sintmp*sint) !compute the longitude for A
            tmpmoid(iii)=1E6 !something big
  407     continue
        endif
cccccc END OF "WATER" PROCEDURE ccccccccccccccccccccccccc

cccccc START OF PARALLEL TUNING ccccccccccccccccccccccc
ccc After the scanning phase, we typically have a few minima on a meridional plane.
ccc The goal of the tuning procedure is to move objects separately along their orbits
ccc in order to find the smallest possible distance between them, which is no longer a meridional distance.
  510    do 615 jjj=1,nmax+1
         if (jjj.le.nmax) then
            moid=tmpmoid(jjj)
			trueB_m=tmptrueB(jjj)
		    longit_m=tmplongit(jjj)
            step=stepini
            threshold=steptresh
         else
           if (nmax.eq.2) then
             if (dabs(tmpmoid(1)-tmpmoid(2)).lt.1E-4) then
ccc in case of two minima are very close to each other(<1E-4 a.u.) go to "water procedure"
               nmax=1
               goto 405
             else
               if (tmpmoid(1)<moid) then
                  moid=tmpmoid(1)
                  trueB_m=tmptrueB(1)
                  longit_m=tmplongit(1)
               endif
             endif
           else
             do 634 iii=1,nmax-1  !the choice of moids for final tuning
              if (tmpmoid(iii)<moid) then
                moid=tmpmoid(iii)
                trueB_m=tmptrueB(iii)
                longit_m=tmplongit(iii)
               endif
  634        continue
           endif
            step=2d0*stepini  !initial step for final tuning
            threshold=stepmin !terminal step for final tuning
         endif
        rBt(2)=radB/(1d0+eccenB*dcos(trueB_m))
        sintmp=dsin(trueB_m+argpeB)
        costmp=dcos(trueB_m+argpeB)
        Bxt(2)=costmpo*costmp-sintmp*sint
        Byt(2)=sintmpo*costmp+sintmp*cost
        Bzt(2)=sintmpi*sintmp
        rAt(2)=radA/(1d0+eccenA*dcos(longit_m))
        Axt(2)=dcos(longit_m)
        Ayt(2)=dsin(longit_m)
        aleft=.true.;aright=.true.;bleft=.true.;bright=.true.
        do while (step>=threshold)
          lpoints=0
          j1min=1;j1max=3
          i1min=1;i1max=3
          calc1=.false.;calc2=.false.
          calc3=.false.;calc4=.false.
          if (bleft) then
            rBt(1)=radB/(1d0+eccenB*dcos(trueB_m-step))
            sintmp=dsin(trueB_m-step+argpeB)
            costmp=dcos(trueB_m-step+argpeB)
            Bxt(1)=costmpo*costmp-sintmp*sint
            Byt(1)=sintmpo*costmp+sintmp*cost
            Bzt(1)=sintmpi*sintmp
            lpoints=lpoints+1
          endif
          if (bright) then
            rBt(3)=radB/(1d0+eccenB*dcos(trueB_m+step))
            sintmp=dsin(trueB_m+step+argpeB)
            costmp=dcos(trueB_m+step+argpeB)
            Bxt(3)=costmpo*costmp-sintmp*sint
            Byt(3)=sintmpo*costmp+sintmp*cost
            Bzt(3)=sintmpi*sintmp
            lpoints=lpoints+1
          endif
          if (aleft) then
            rAt(1)=radA/(1d0+eccenA*dcos(longit_m-step))
            Axt(1)=dcos(longit_m-step)
            Ayt(1)=dsin(longit_m-step)
            lpoints=lpoints+1
          endif
          if (aright) then
            rAt(3)=radA/(1d0+eccenA*dcos(longit_m+step))
            Axt(3)=dcos(longit_m+step)
            Ayt(3)=dsin(longit_m+step)
            lpoints=lpoints+1
          endif
          j1_t=2;i1_t=2
          if (lpoints.eq.1) then
            if (aleft) i1max=1
            if (aright) i1min=3
            if (bleft) j1max=1
            if (bright) j1min=3
          endif
          if (lpoints.eq.2) then
            if (aleft.and.bright) calc1=.true.
            if (aleft.and.bleft) calc2=.true.
            if (aright.and.bright) calc3=.true.
            if (aright.and.bleft) calc4=.true.
          endif
          do 557 j1=j1min,j1max
            do 555 i1=i1min,i1max
              if (lpoints.eq.2) then
                if (i1.ne.1) then
                  if (((j1.ne.3).and.calc1)
     a             .or.((j1.ne.1).and.calc2)) goto 555
                endif
                if (i1.ne.3) then
                  if (((j1.ne.3).and.calc3)
     a             .or.((j1.ne.1).and.calc4)) goto 555
                endif
              endif
              if ((i1.eq.2).and.(j1.eq.2)) goto 555
              Dx=rBt(j1)*Bxt(j1)-rAt(i1)*Axt(i1)
              Dy=rBt(j1)*Byt(j1)-rAt(i1)*Ayt(i1)
              Dz=rBt(j1)*Bzt(j1)
              dist=(Dx*Dx+Dy*Dy+Dz*Dz)
              if (dist<moid) then
                 moid=dist
                 j1_t=j1;i1_t=i1
              endif
  555       continue
  557     continue
         if ((j1_t.ne.2).or.(i1_t.ne.2)) then
          aleft=.false.;aright=.false.;bleft=.false.;bright=.false.
          if (i1_t.ne.2) then
            if (i1_t.eq.1) then
              aleft=.true.
              longit_m=longit_m-step
              rAt(3)=rAt(2);Axt(3)=Axt(2);Ayt(3)=Ayt(2)
              rAt(2)=rAt(1);Axt(2)=Axt(1);Ayt(2)=Ayt(1)
            else
              aright=.true.
              longit_m=longit_m+step
              rAt(1)=rAt(2);Axt(1)=Axt(2);Ayt(1)=Ayt(2)
              rAt(2)=rAt(3);Axt(2)=Axt(3);Ayt(2)=Ayt(3)
            endif
          endif
          if (j1_t.ne.2) then
            if (j1_t.eq.1) then
              bleft=.true.
              trueB_m=trueB_m-step
              rBt(3)=rBt(2);Bxt(3)=Bxt(2);Byt(3)=Byt(2);Bzt(3)=Bzt(2)
              rBt(2)=rBt(1);Bxt(2)=Bxt(1);Byt(2)=Byt(1);Bzt(2)=Bzt(1)
            else
              bright=.true.
              trueB_m=trueB_m+step
              rBt(1)=rBt(2);Bxt(1)=Bxt(2);Byt(1)=Byt(2);Bzt(1)=Bzt(2)
              rBt(2)=rBt(3);Bxt(2)=Bxt(3);Byt(2)=Byt(3);Bzt(2)=Bzt(3)
            endif
          endif
         else
            aleft=.true.;aright=.true.;bleft=.true.;bright=.true.
            step=step*.15 !.15 is optimal value
         endif
        end do
         if (jjj.le.nmax) then
            tmpmoid(jjj)=moid
            tmptrueB(jjj)=trueB_m
            tmplongit(jjj)=longit_m
         endif
  615   continue
cccccc END OF PARALLEL TUNING ccccccccccccccccccccccc
        moid=dsqrt(moid) !we dealed with squares
        write (*,*) 'MOID [AU] =', moid
      END
