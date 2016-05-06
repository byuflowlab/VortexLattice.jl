function Pdrag(alt, mach, xt, mac, sweep, tc)

    # ------------- constants ---------------
    rho, mu, a, T = atmosphere(alt)
    Re = rho*mach*a*mac/mu
    # -----------------------------------------

    # -------------- Cf ---------------------
    # compute Cf based on transition location
    Rex = Re*xt
    if Rex <= 0
      Rex = 0.0001
    end
    xeff = 38.7*xt*Rex^(-3.0/8) # effective boundary layer length
    Rext = Re*(1-xt+xeff)

    Cfturb = 0.455/(log10(Rext))^2.58
    Cflam = 1.328/sqrt(Rex)
    Cfstart = 0.455/(log10(Re*xeff))^2.58
    Cf_inc = Cflam*xt + Cfturb*(1-xt+xeff) - Cfstart*xeff

    # roughness increment
    Cf_inc = 1.07*Cf_inc

    # effect of mach number
    Tw = 1 + 0.178*mach^2
    Tp = 1 + 0.035*mach^2 + 0.45*(Tw-1)
    mup = Tp^1.5*(T+216)/(Tp*T+216)
    Rp = 1/mup/Tp
    Cf = Cf_inc/Tp/Rp^0.2
    # ---------------------------------------

    # ------------ form factor ----------------------------
    cossw = cos(sweep)
    m0 = 0.5
    z = (2-m0^2)*cossw/sqrt(1-(m0*cossw)^2)
    k = 1 + tc*z + tc^4*100
    # -----------------------------------------------------

    # ---------- wetted area / S ----------------------
    SwetS = 2*(1+0.2*tc)
    # --------------------------------------------

    # parasite drag
    CDp = Cf*k*SwetS

    return CDp
end



function Cdrag(CL, Lambda, tc, mach, supercrit)

    cosL = cos(Lambda)
    clp = CL/cosL^2
    tcp = tc/cosL

    # compute Mcc
    Mcc = 0.954-0.235*clp+0.0259*clp^2
    Mcc = Mcc - (1.963-1.078*clp+0.350*clp^2)*tcp
    Mcc = Mcc + (2.969-2.738*clp+1.469*clp^2)*tcp.^2
    Mcc = Mcc + supercrit*.06
    Mcc = Mcc/cosL

    # compute Cdc
    rm = mach / Mcc
    dm = rm-1

    if rm < .5
      cdc = 0.0
    elseif (rm >= .5 && rm < .8)
      cdc = 1.3889e-4+5.5556e-4*dm+5.5556e-4*dm*dm
    elseif (rm >= .8 && rm < .95)
      cdc = 7.093e-4+6.733e-3*dm+.01956*dm*dm+.01185*dm*dm*dm
    elseif (rm >= .95 && rm < 1.0)
      cdc = .001000+.02727*dm+.4920*dm*dm+3.573*dm*dm*dm
    elseif (rm >= 1.0 && rm < 1.075)
      cdc = .001000+.02727*dm-.1952*dm*dm+19.09*dm*dm*dm
    else
      cdc = 0.01 + 0.33477*(dm-0.075) # linear extension
    end

    cdc = cdc*cosL^3
    
    return cdc
end
