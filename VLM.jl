function VLM(wing,fs,ref,pdrag,mvr,plots, QC, TE, LE, CP)
  # freestream properties
  if fs.method == "CL"
      CL_mthd = true
  else
    CL_mthd = false
  end

  if CL_mthd # specify CL
      CLref = fs.CL
      alpha = 0
  else
      alpha = fs.alpha
      CLref = 0
  end

  mach = fs.mach
  rho = 7.382e-4
  U = 700
  q = 1/2*rho*U^2

  # reference quantities
  Sref = ref.S
  cref = ref.c

  # viscous drag (2 possible methods)
  if pdrag.method == "pass"
    pass_mthd = true
  else
    pass_mthd = false
  end

  if pass_mthd # Reynolds number dependent method
      alt = pdrag.alt
      xt = pdrag.xt

      # not used
      cd0 = 0
      cd1 = 0
      cd2 = 0
  else # quadratic variation with section lift coefficient
      cd0 = pdrag.polar(1)
      cd1 = pdrag.polar(2)
      cd2 = pdrag.polar(3)

      # not used
      alt = 0
      xt = 0
  end

  # structures
  qmvrN = mvr.qN # ratio of maneuver dynamic pressure to cruise dynamic pressure
  n = mvr.n # load factor
  kbar = mvr.kbar # coefficient used in computing area dependent weight (3.57e4 or 0)

  # stall
  CLmax = ref.CLmax
  # -------------------------------------

  # ------------- geometry ------------------
  (QC,TE,CP,LE) = gen_geometry(wing, QC, TE, CP, LE)

  # compute wing area
  S = 2*sum(CP.chord.*CP.ds)
  # ------------------------------------------

  # ---------- force influence coefficients --------------------
  # aerodynamic influence coefficients
  AIC = getAIC(QC,TE,CP)

  # lift
  LIC = getLIC(CP,rho,U)

  # induced drag
  DIC = getDIC(TE.y,TE.z,rho, TE.y, TE.z)

  # viscous drag
  (D1, D2) = getViscous(cd1,cd2,CP,rho,U)

  # weight
  (WIC, BMM) = getWIC(CP,rho,U)
  # -----------------------------------------------

  # ---------- compute angle of attack necessary to match CL -----
  if CL_mthd
      alpha = getAlpha(CLref*q*Sref,CP,LIC,AIC,U)
  end
  # ------------------------------------------------------------------

  # -------- compute circulation -----------------
  Vn = -U*(cos(alpha)*sin(CP.twist) + sin(alpha)*cos(CP.twist).*cos(CP.dihedral)')
  gamma = (AIC\Vn)
  # ----------------------------------------------

  # --------- aerodynamic forces ------------------
  L = LIC'*gamma
  Di = gamma'*DIC*gamma
  Dp = cd0*q*S + D1*gamma + D2*gamma.^2

  CL = L/q/Sref
  CDi = Di/q/Sref
  CDp = Dp/q/Sref
  # ------------------------------------------------

  # --------- weight (integrated bending moment over thickness --------
  # circulation at maneuver load
  bbb = AIC\cos([CP.dihedral])'
  LL = LIC'*bbb
  gamma_mvr = gamma + ((n/qmvrN-1)*(LIC'*gamma)/LL*bbb')'

  # compute weight
  W = qmvrN*WIC*gamma_mvr

  # add area dependent weight
  W = W + kbar*S

  # weight coefficient
  CW = W/q/Sref/cref
  # -----------------------------------------------------------------

  # ----- estimate compressibility and parasite drag (strip theory + PASS method) -----
  supercrit = 1
  P = round(wing.span/sum(wing.span)*wing.N) # number of panels in each section
  cdc = zeros(1,length(wing.span))
  cdp = zeros(1,length(wing.span))
  area = zeros(1,length(wing.span))
  start = 1

  for i = 1:length(wing.span)
      finish = start + P[i] - 1

      # rename for convenience
      cr = wing.chord[i]
      ct = wing.chord[i+1]
      cbar = 1/2*(cr + ct)
      tcbar = (wing.tc[i]*wing.chord[i] + wing.tc[i+1]*wing.chord[i+1])/(wing.chord[i]+wing.chord[i+1])
      area[i] = cbar*wing.span[i]
      mac = 2/3*(cr + ct - cr*ct/(cr+ct))
      CL_local = sum(gamma[start:finish]'.*CP.ds[start:finish])*2/U/area[i]

      # compressibility drag
      cdc[i] = Cdrag(CL_local,wing.sweep[i],tcbar,mach,supercrit)

      # parasite drag
      cdp[i] = Pdrag(alt,mach,xt,mac,wing.sweep[i],tcbar)

      start = finish + 1
  end

  # compressibility drag - area weighted average
  CDc = 2*sum(cdc.*area)/Sref[1]

  if pass_mthd
      # parasite drag - area weighted average
      CDp = 2*sum(cdp.*area)/Sref[1]

      # add viscous dependent induced drag
      Lambda_bar = area*wing.sweep'/sum(area)
      CDi = CDi + 0.38*CDp*CL^2/cos(Lambda_bar)^2
  end
  # ----------------------------------------------------------------------

  # ----------- cl distribution at CLmax ----------------
  cl = 2/U*gamma./CP.chord
  clmax_dist = cl + (rho*U*(CLmax*Sref-L/q)/LL)[1]*(bbb)./CP.chord #.'

  # clmax as a function of thickness - polynomial fit
  p_cl = Poly([-2.7634e-5, 0.0022307, -0.06567, 0.8013, -1.748])
  clmax = polyval(p_cl,CP.tc*100)

  cl_margin = clmax - clmax_dist
  # ----------------------------------------

  # --------- pitching moment about a.c. --------------------
  # pitching moment about quarter chord
  MIC = getMIC(CP,rho,U,QC.x[1])

  # find aerodynamic center
  dRHS = -sin(alpha)*sin(CP.twist) + cos(alpha)*cos(CP.twist).*cos(CP.dihedral')
  dbc = AIC\dRHS
  dMda = MIC'*dbc
  dLda = LIC'*dbc
  xac = -dMda/dLda + QC.x[1]

  # find Moment about a.c.
  MICac = getMIC(CP,rho,U,xac[1])
  Mac = MICac'*gamma

  Cmac = Mac/q/Sref/cref
  # ----------------------------------------------------------

  # --------------- structures --------------------------
  # bending moment distribution
  Mb = qmvrN*BMM*gamma

  # distance along structural span
  ds_str = CP.ds'./cos(CP.sweep)
  eta_str = [0 cumsum(ds_str)]
  eta_str = 0.5*(eta_str[1:end-1] + eta_str[2:end])
  # --------------------------------------------------------------

  if (plots)
      # ------------- plots --------------------
      # plot wings
#       plot_wing(LE,QC,TE,CP)

      N = length(QC.x)
#       axis equal
      for i = 1:N-1
        PyPlot.plot(-[LE.x[i] LE.x[i+1]], [LE.y[i] LE.y[i+1]])
        PyPlot.plot([TE.y[i] TE.y[i+1]],-[TE.x[i] TE.x[i+1]])
      end
      PyPlot.plot([LE.y[1] TE.y[1]],-[LE.x[1] TE.x[1]])
      PyPlot.plot([LE.y[end] TE.y[end]],-[LE.x[end] TE.x[end]])

      # plot lift
#       figure(50) hold on
      PyPlot.figure()
      l = 2*gamma/U/(ref.c)
      l_mvr = 2*gamma_mvr/U/(ref.c)
      eta = linspace(0,0.5,length(l))
      PyPlot.plot(eta,l,"b")
      PyPlot.plot(eta,l_mvr,"r")
#       xlabel('\xi / b')
#       ylabel('$$\frac{c_l c}{c_{ref}}$$','interpreter','latex','rotation',0)

      # plot cl
      PyPlot.figure()
      PyPlot.plot(eta,cl,"b")
      PyPlot.plot(eta,clmax_dist,"r")
      PyPlot.plot(eta,clmax,"r--")
#       xlabel('\xi / b')
#       ylabel('c_l')

      # plot bending over thickness
      PyPlot.figure()
      PyPlot.plot(eta_str/eta_str[end]*0.5,Mb./(CP.tc.*CP.chord)')
#       xlabel('\xi_{str}/b_{str}')
#       ylabel('M_b/t')
      PyPlot.show()
      # -------------------------------------------------
  end

  return [CL, CDi, CDp, CDc, CW, Cmac, cl_margin, gamma, CP]
end

