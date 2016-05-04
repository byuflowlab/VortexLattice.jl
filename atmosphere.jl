function atmosphere(altitude::Float64)

    # assumes english units
    # [rho, mu, a, T, P] = atmosphere(altitude)

    # ----------- constants ---------------
    aT = [-6.5 0 1 2.8 0 -2.8 -2]*0.00054864 # temperature gradient (R/ft)
    h = [0 11 20 32 47 51 71 84.852]*3280.8399 # altitude (ft)
    g = 32.174 # gravitational acceleration
    R = 1716.5 # specific gas constant
    Tsl = 518.67 # sea level temperature
    Psl = 2116.21662 # sea level pressure
    musl = 3.73719712e-7 # sea level viscosity
    S = 1.8*110.4 # constant in Sutherlands formula
    gamma = 1.4
    # ---------------------------------------

    if altitude > h[end]
      println("Altitude exceeds standard atmosphere data")
    end

    # ---- find temperature and pressure at defined points ----
    Tpts = zeros(9)
    Ppts = zeros(9)
    Tpts[1] = Tsl
    Ppts[1] = Psl

    for i = 2:8
        Tpts[i] = Tpts[i-1] + aT[i-1]*(h[i]-h[i-1])

        if aT[i-1] == 0
            Ppts[i] = Ppts[i-1]*exp(-g*(h[i]-h[i-1])/R/Tpts[i-1])
        else
            Ppts[i] = Ppts[i-1]*(Tpts[i-1]/Tpts[i])^(g/R/aT[i-1])
        end
    end
    # ---------------------------------------------

    # ------ find values at altitude ---------
    hidx = find(altitude .>= h)[end]

    T = Tpts[hidx] + aT[hidx]*(altitude-h[hidx])

    if aT[hidx] == 0
        P = Ppts[hidx]*exp(-g*(altitude-h[hidx])/R/Tpts[hidx])
    else
        P = Ppts[hidx]*(Tpts[hidx]/T)^(g/R/aT[hidx])
    end

    rho = P/R/T
    # ----------------------------------------------

    # ------- Sutherlands Law ------------
    mu = musl*(T/Tsl)^(3.0/2)*(Tsl+S)/(T+S)
    # -----------------------------------

    # --------------- speed of sound -------------------
    a = sqrt(gamma*R*T)

    return rho, mu, a, T, P

end
