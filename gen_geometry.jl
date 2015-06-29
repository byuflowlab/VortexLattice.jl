function gen_geometry(wing, QC, TE, CP, LE)
  # ----------------- rename for convenience --------------------
  b = wing.span
  Lambda = wing.sweep
  chord = wing.chord
  twist = wing.twist
  tc = wing.tc
  phi = wing.dihedral
  N = wing.N
  # --------------------------------------------------------------

  P = int(round(b/sum(b)*N)) # divide up panels
  # -------------  Quarter Chord Locations --------------------------
  QC.x[1] = 0
  QC.y[1] = 0
  QC.z[1] = 0
  last = 1
  c = 0
  t = 0
  thickness = 0
  for i = 1:length(b)
      first = last
      last = first + P[i]
      eta = linspace(0, b[i], P[i]+1)

      QC.x = QC.x[first] + eta*tan(Lambda[i])
      QC.y = QC.y[first] + eta*cos(phi[i])
      QC.z = QC.z[first] + eta*sin(phi[i])
      c = chord[i] + eta*(chord[i+1]-chord[i])/b[i] # chord
      t = twist[i] + eta*(twist[i+1]-twist[i])/b[i] # twist
      thickness = tc[i]*chord[i] + eta*(tc[i+1]*chord[i+1]-tc[i]*chord[i])/b[i] # thickness
  end
  # --------------------------------------------------------------

  # ------------ Trailing Edge Locations --------------------------
  TE.x = QC.x + 3/4*c
  TE.y = QC.y
  TE.z = QC.z
  # ----------------------------------------------------------------

  # -------------- Leading Edge locations ------------------------
  LE.x = QC.x - c/4
  LE.y = QC.y
  LE.z = QC.z
  # -----------------------------------------------------------------

  # ------------- Control Point Locations --------------------------
  CP.chord = 1/2*(c[1:N] + c[2:N+1])
  CP.twist = 1/2*(t[1:N] + t[2:N+1])
  CP.tc = 1/2*(thickness[1:N] + thickness[2:N+1])./CP.chord
  CP.x = 1/2*(QC.x[1:N] + QC.x[2:N+1]) + 1/2*CP.chord
  CP.y = 1/2*(QC.y[1:N] + QC.y[2:N+1])
  CP.z = 1/2*(QC.z[1:N] + QC.z[2:N+1])

  last = 0
  for i = 1:length(b)
      first = last + 1
      last = first + P[i] - 1
      CP.dihedral = phi[i]*ones(1,P[i])
      CP.sweep = Lambda[i]*ones(1,P[i])
  end
  CP.ds = sqrt((QC.y[2:N+1]-QC.y[1:N]).^2 + (QC.z[2:N+1]-QC.z[1:N]).^2)
  return QC, TE, CP, LE
end

