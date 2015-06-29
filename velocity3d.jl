function velocity3d(QC, TE, CP)
  # first calculate induced velocities from vorticies on right side of wing.
  m = length(CP.y)
  n = length(QC.y)-1

  # u = zeros(m,n)
  v = zeros(m,n)
  w = zeros(m,n)

  # Induced velocities due to bound vortex (BC)
  dv, dw = vortex(QC.x[1:n],QC.y[1:n],QC.z[1:n],QC.x[2:n+1],QC.y[2:n+1],QC.z[2:n+1],CP.x,CP.y,CP.z)
  # u = u + du
  v += dv
  w += dw

  # Induced velocities due to AB - left vortex bound to wing
  dv, dw = vortex(TE.x[1:n],TE.y[1:n],TE.z[1:n],QC.x[1:n],QC.y[1:n],QC.z[1:n],CP.x,CP.y,CP.z)
  # u = u + du
  v += dv
  w += dw

  # Induced velocities due to CD - right vortex bound to wing
  dv, dw = vortex(QC.x[2:n+1],QC.y[2:n+1],QC.z[2:n+1],TE.x[2:n+1],TE.y[2:n+1],TE.z[2:n+1],CP.x,CP.y,CP.z)
  # u = u + du
  v += dv
  w += dw

  # Induced velocities due to trailing vortex from A
  l = 1e9  # some large number to represent infinity
  # xINF = TE.x + l*cos(alpha)*ones(size(TE.x))
  # yINF = TE.y
  # zINF = TE.z + l*sin(alpha)*ones(size(TE.z))
  xINF = TE.x + l
  yINF = TE.y
  zINF = TE.z

  dv, dw = vortex(xINF[1:n],yINF[1:n],zINF[1:n],TE.x[1:n],TE.y[1:n],TE.z[1:n],CP.x,CP.y,CP.z)
  # u = u + du
  v += dv
  w += dw

  # Induced velocities due to trailing vortex from D
  dv, dw = vortex(TE.x[2:n+1],TE.y[2:n+1],TE.z[2:n+1],xINF[2:n+1],yINF[2:n+1],zINF[2:n+1],CP.x,CP.y,CP.z)
  # u = u + du
  v += dv
  w += dw

  return v, w
end