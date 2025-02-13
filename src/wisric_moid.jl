"""
    wisric_moid(saxisA, eccenA, argpeA, omegaA, incliA, saxisB, eccenB, argpeB, omegaB, incliB)

Compute the Minimum Orbit Intersection Distance (MOID) between two orbits A and B for which the
Keplerian elements (a,e,ω,Ω,i) are given. Angular arguments are in [rad].

## Return
   MOID in same units like semi-major axis a (usually au)

## Reference
Based on the paper by T.Wisniowski and H.Rickman "A Fast, Geometric Method for Calculating Accurate Minimum Orbit
Intersection Distances (MOIDs)", Acta Astronomica, Vol. 63 (2013) pp. 293–307 and the Fortran program version
v4.0 (July 30th,2017) written by the authors of the paper cited above.

See also: http://moid.cbk.waw.pl/orbity/default/index
"""
function wisric_moid(saxisA, eccenA, argpeA, omegaA, incliA, saxisB, eccenB, argpeB, omegaB, incliB)

# Tables
    rAt = zeros(3)
    rBt = zeros(3)
    Axt = zeros(3)
    Ayt = zeros(3)
    Bxt = zeros(3)
    Byt = zeros(3)
    Bzt = zeros(3)
    tmpmoid = zeros(10)
    tmptrueB = zeros(10)
    tmplongit = zeros(10)

# Parameters of the program. These steps are optimized with respect to speed and reliability
    cstep = 0.12        # scanning step of true anomaly/longitude in radians
    stepini = 0.07      # initial step of first tuning in radians
    steptresh = 1E-5    # final step of first tunig (to choose the MOID) in radians

# This step depends on the expected final accuracy
    stepmin = 1E-14     # threshold step of second tuning in radians

# Constants and given values
    twopi = 2*pi

# START OF PREPARING THE ORBITS

# Computing parameters of transition matrix C11...C33
    c11 = cos(omegaA)*cos(argpeA) - sin(omegaA)*cos(incliA)*sin(argpeA)
    c12 = sin(omegaA)*cos(argpeA) + cos(omegaA)*cos(incliA)*sin(argpeA)
    c13 = sin(incliA)*sin(argpeA)
    c21 =-cos(omegaA)*sin(argpeA) - sin(omegaA)*cos(incliA)*cos(argpeA)
    c22 =-sin(omegaA)*sin(argpeA) + cos(omegaA)*cos(incliA)*cos(argpeA)
    c23 = sin(incliA)*cos(argpeA)
    c31 = sin(incliA)*sin(omegaA)
    c32 =-sin(incliA)*cos(omegaA)
    c33 = cos(incliA)

# Calculating new values of Euler angles using transition matrix
    sintmpi = sin(incliB)
    costmpi = cos(incliB)
    costmpo = cos(omegaB)
    sintmpo = sin(omegaB)
    costmpa = cos(argpeB)
    sintmpa = sin(argpeB)
    x1 = costmpo*costmpa - sintmpo*costmpi*sintmpa
    x2 = sintmpo*costmpa + costmpo*costmpi*sintmpa
    x3 = sintmpi*sintmpa
    y1 =-costmpo*sintmpa - sintmpo*costmpi*costmpa
    y2 =-sintmpo*sintmpa + costmpo*costmpi*costmpa
    y3 = sintmpi*costmpa
    z1 = sintmpi*sintmpo
    z2 = -sintmpi*costmpo
    z3 = costmpi
    z1n = c11*z1 + c12*z2 + c13*z3
    z2n = c21*z1 + c22*z2 + c23*z3
    z3n = c31*z1 + c32*z2 + c33*z3
    y3n = c31*y1 + c32*y2 + c33*y3
    x3n = c31*x1 + c32*x2 + c33*x3
    incliB = atan(sqrt(z1n*z1n + z2n*z2n),z3n)
    omegaB =-atan(z1n,-z2n)
    argpeB =-atan(x3n,y3n)

# Helpful precalculated values
    costmpo = cos(omegaB)
    sintmpo = sin(omegaB)
    sintmpi = sin(incliB)   # = z1n/sintmpo
    costmpi = z3n           # = cos(incliB)
    sint = sintmpo*costmpi
    cost = costmpo*costmpi
    radA = saxisA*(1.0 - eccenA*eccenA)
    radB = saxisB*(1.0 - eccenB*eccenB)

# END OF PREPARING THE ORBITS

# START OF SCANNING
# This tool yields a preliminary approach to the minima of the distance function.
# By scanning one full revolution of meridional plane we look for the local minima.

# Initial parameters
        iii1 = 0
        jjj1 = 0
        longit_m = 0
        longit_o = 0
        trueB = -2.0*cstep
        trueB_m = 0
        trueB_o = 0
        dist = 0
        dist_oo = 0
        dist_o = 1E6                # something big
        moid = 1E6                  # dito
        tmpmoid = [1E6 for n=1:4]   # dito

# Looking for the minima with rotating meridional plane
# a) At first we calculate the coordinates of two additional positions of the plane to create first triplet

    for iii in 1:2
        rB = radB/(1.0 + eccenB*cos(trueB))     # compute the radius for B
        sintmp = sin(trueB + argpeB)
        costmp = cos(trueB + argpeB)
        Bz_sq = sintmpi*sintmp
        Bz_sq = Bz_sq*Bz_sq                     # square of Z-coordinate for B
        longit = atan(sintmpo*costmp+sintmp*cost,costmpo*costmp-sintmp*sint) # compute the longitude for A
        tmp2 = eccenA*cos(longit)               # temporary value
        rA = radA/(1.0 + tmp2)                  # compute the radius for A (two possibilities)
        rA2 = radA/(1.0 - tmp2)
        tmp1 = rB*sqrt(1.0 - Bz_sq)             # temporary value
        if (abs(tmp1 - rA) > abs(tmp1 + rA2))
            rA = rA2
            longit = longit - pi                 # the second possibility gives smaller distance
            tmp1 = tmp1 + rA2
        else
            tmp1 = tmp1 - rA
        end
        dist = rB*rB*Bz_sq + tmp1*tmp1          # square of the distance A-B
        if (iii == 1)
            dist_oo = dist
        else
            dist_o = dist
            trueB_o = trueB
            longit_o = longit
        end
        trueB = trueB + cstep
    end

# b) Now we scan doing one full revolution of the meridional plane

    nmax = 0                                # counts the minima
    dist_min = dist
    while (trueB < (twopi + cstep))         # loop for true anomaly of B
        rB = radB/(1.0 + eccenB*cos(trueB)) # compute the radius for B
        sintmp = sin(trueB + argpeB)
        costmp = cos(trueB + argpeB)
        Bz_sq = sintmpi*sintmp
        Bz_sq = Bz_sq*Bz_sq                 # square of Z-coordinate for B
        longit = atan(sintmpo*costmp+sintmp*cost,costmpo*costmp-sintmp*sint) # compute the longitude for A
        tmp2 = eccenA*cos(longit)           # temporary value
        rA = radA/(1.0 + tmp2)              # compute the radius for A (two possibilities)
        rA2 = radA/(1.0 - tmp2)
        tmp1 = rB*sqrt(1.0 - Bz_sq)         # temporary value

        if (abs(tmp1-rA) > abs(tmp1+rA2))
            rA = rA2
            longit = longit - pi              # the second possibility gives smaller distance
            tmp1 = tmp1 + rA2
        else
            tmp1 = tmp1 - rA
        end

        dist = rB*rB*Bz_sq + tmp1*tmp1      # square of the distance A-B

        if ((dist_o <= dist) && (dist_o <= dist_oo))  # the minimum was found
            nmax = nmax + 1
            tmptrueB[nmax] = trueB_o
            tmplongit[nmax] = longit_o
            tmpmoid[nmax] = dist_o
        end

        if (dist_min > dist) dist_min = dist end

        dist_oo = dist_o
        trueB_o = trueB
        longit_o = longit
        dist_o = dist
        trueB = trueB + cstep
    end

    # END OF SCANNING

    # "WATER" PROCEDURE (can be skipped to increase computing speed)
    water_procedure(nmax, argpeB, sintmpo, costmpo, sint, cost, tmptrueB, tmplongit, tmpmoid)

    # START OF PARALLEL TUNING
    # After the scanning phase, we typically have a few minima on a meridional plane.
    # The goal of the tuning procedure is to move objects separately along their orbits
    # in order to find the smallest possible distance between them, which is no longer a meridional distance.

    for jjj in 1:nmax+1

        if (jjj <= nmax)
            moid = tmpmoid[jjj]
            trueB_m = tmptrueB[jjj]
            longit_m = tmplongit[jjj]
            step = stepini
            threshold = steptresh
        else
            if (nmax == 2)
                # in case of two minima are very close to each other(<1E-4 a.u.) go to "water procedure"
                if (abs(tmpmoid[1] - tmpmoid[2]) < 1E-4)
                    nmax = 1
                    water_procedure(nmax, argpeB, sintmpo, costmpo, sint, cost, tmptrueB, tmplongit, tmpmoid)
                else
                    if (tmpmoid[1] < moid)
                        moid = tmpmoid[1]
                        trueB_m = tmptrueB[1]
                        longit_m = tmplongit[1]
                    end
                end
            else
                for iii in 1:nmax-1  # the choice of moids for final tuning
                    if (tmpmoid[iii] < moid)
                        moid = tmpmoid[iii]
                        trueB_m = tmptrueB[iii]
                        longit_m = tmplongit[iii]
                    end
                end
        end
        step = 2.0*stepini  # initial step for final tuning
        threshold = stepmin # terminal step for final tuning
    end

    rBt[2] = radB/(1.0 + eccenB*cos(trueB_m))
    sintmp = sin(trueB_m + argpeB)
    costmp = cos(trueB_m + argpeB)
    Bxt[2] = costmpo*costmp - sintmp*sint
    Byt[2] = sintmpo*costmp + sintmp*cost
    Bzt[2] = sintmpi*sintmp
    rAt[2] = radA/(1.0 + eccenA*cos(longit_m))
    Axt[2] = cos(longit_m)
    Ayt[2] = sin(longit_m)
    aleft = true ; aright = true
    bleft = true ; bright = true

    while (step >= threshold)
        lpoints = 0
        j1min = 1 ; j1max = 3
        i1min = 1 ; i1max = 3
        calc1 = false ; calc2 = false
        calc3 = false ; calc4 = false
        if (bleft)
            rBt[1] = radB/(1.0 + eccenB*cos(trueB_m - step))
            sintmp = sin(trueB_m-step + argpeB)
            costmp = cos(trueB_m-step + argpeB)
            Bxt[1] = costmpo*costmp - sintmp*sint
            Byt[1] = sintmpo*costmp + sintmp*cost
            Bzt[1] = sintmpi*sintmp
            lpoints = lpoints + 1
        end
        if (bright)
            rBt[3] = radB/(1.0 + eccenB*cos(trueB_m + step))
            sintmp = sin(trueB_m + step+argpeB)
            costmp = cos(trueB_m + step+argpeB)
            Bxt[3] = costmpo*costmp - sintmp*sint
            Byt[3] = sintmpo*costmp + sintmp*cost
            Bzt[3] = sintmpi*sintmp
            lpoints = lpoints + 1
        end
        if (aleft)
            rAt[1] = radA/(1.0 + eccenA*cos(longit_m - step))
            Axt[1] = cos(longit_m - step)
            Ayt[1] = sin(longit_m - step)
            lpoints = lpoints + 1
        end
        if (aright)
            rAt[3] = radA/(1.0 + eccenA*cos(longit_m + step))
            Axt[3] = cos(longit_m + step)
            Ayt[3] = sin(longit_m + step)
            lpoints = lpoints + 1
        end
        j1_t = 2 ; i1_t = 2
        if (lpoints == 1)
            if (aleft) i1max = 1 end
            if (aright) i1min = 3 end
            if (bleft) j1max = 1 end
            if (bright) j1min = 3 end
        end
        if (lpoints == 2)
            if (aleft && bright) calc1 = true end
            if (aleft && bleft) calc2 = true end
            if (aright && bright) calc3 = true end
            if (aright && bleft) calc4 = true end
        end

        for j1 in j1min:j1max
            for i1 in i1min:i1max
                if (lpoints == 2)
                if (i1 != 1)
                    if (((j1 != 3) && calc1) || ((j1 != 1) && calc2)) continue end
                end
                if (i1 != 3)
                    if (((j1 != 3) && calc3) || ((j1 != 1) && calc4)) continue end
                end
                end
                if ((i1 == 2) && (j1 == 2)) continue end
                Dx = rBt[j1]*Bxt[j1] - rAt[i1]*Axt[i1]
                Dy = rBt[j1]*Byt[j1] - rAt[i1]*Ayt[i1]
                Dz = rBt[j1]*Bzt[j1]
                dist = (Dx*Dx+Dy*Dy+Dz*Dz)
                if (dist < moid)
                    moid = dist
                    j1_t = j1
                    i1_t = i1
                end
            end
        end

        if ((j1_t != 2) || (i1_t != 2))
            aleft = false ; aright = false
            bleft = false ; bright = false
            if (i1_t != 2)
                if (i1_t == 1)
                    aleft = true
                    longit_m = longit_m - step
                    rAt[3] = rAt[2] ; Axt[3] = Axt[2] ; Ayt[3] = Ayt[2]
                    rAt[2] = rAt[1] ; Axt[2] = Axt[1] ; Ayt[2] = Ayt[1]
                else
                    aright = true
                    longit_m=longit_m + step
                    rAt[1] = rAt[2] ; Axt[1] = Axt[2] ; Ayt[1] = Ayt[2]
                    rAt[2] = rAt[3] ; Axt[2] = Axt[3] ; Ayt[2] = Ayt[3]
                end
            end
            if (j1_t != 2)
                if (j1_t == 1)
                    bleft = true
                    trueB_m = trueB_m - step
                    rBt[3] = rBt[2] ; Bxt[3] = Bxt[2] ; Byt[3] = Byt[2] ; Bzt[3] = Bzt[2]
                    rBt[2] = rBt[1] ; Bxt[2] = Bxt[1] ; Byt[2] = Byt[1] ; Bzt[2] = Bzt[1]
                else
                    bright = true
                    trueB_m = trueB_m + step
                    rBt[1] = rBt[2] ; Bxt[1] = Bxt[2] ; Byt[1] = Byt[2] ; Bzt[1] = Bzt[2]
                    rBt[2] = rBt[3] ; Bxt[2] = Bxt[3] ; Byt[2] = Byt[3] ; Bzt[2] = Bzt[3]
                end
            end
            else
                aleft = true ; aright = true ; bleft = true ; bright = true
                step=step*0.15  # 0.15 is optimal value
            end
        end
        if (jjj <= nmax)
            tmpmoid[jjj] = moid
            tmptrueB[jjj] = trueB_m
            tmplongit[jjj] = longit_m
        end
    end

    # END OF PARALLEL TUNING
    return sqrt(moid) # we dealed with squares
end


"""
In case only one minimum is detected we take a special care to avoid the risk of missing the minima.
Instead of starting the tuning with one detected minimum, we start it with four positions of the meridional plane, evenly
distributed along the inclined orbit. The later tuning procedure moves the points along the orbits similarly as water droplets
are leading by the gravity to the points of minimum height, so we called this phase "water procedure". It slightly slows down
the calculations, but minimizes the risk of missing the MOID. With "water procedure" the risk of missing is < 1E-6,
without "water procedure" the risk of missing is about 3E-5.
"""
function water_procedure(nmax, argpeB, sintmpo, costmpo, sint, cost, tmptrueB, tmplongit, tmpmoid)
    if (nmax < 2)    # only one minimum was detected
        nmax = 4
        for iii in 1:4
            tmptrueB[iii] = (0.25 + 0.5*iii)*pi  # evenly distributed points
            sintmp = sin(tmptrueB[iii] + argpeB)
            costmp = cos(tmptrueB[iii] + argpeB)
            tmplongit[iii] = atan(sintmpo*costmp + sintmp*cost,costmpo*costmp - sintmp*sint) # compute the longitude for A
            tmpmoid[iii] = 1E6 # something big
        end
    end
    return nothing
end
