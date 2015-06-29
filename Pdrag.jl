function Pdrag(alt,mach,xt,mac,sweep,tc)

  # ------------- constants ---------------
  (rho,mu,a,T) = atmosphere(alt)
  Re = rho*mach*a*mac/mu
  # -----------------------------------------

  # -------------- Cf ---------------------
  # compute Cf based on transition location
  Rex = Re*xt
  if Rex <=0
      Rex = 0.0001
  end
  xeff = 38.7*xt*Rex^(-3/8) # effective boundary layer length
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

