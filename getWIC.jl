function getWIC(CP,rho,U)
  x = CP.x - 1/2*CP.chord
  y = CP.y
  z = CP.z

  # rename for convenience
  N = length(y)
  phi = CP.dihedral
  Lambda = CP.sweep
  chord = CP.chord
  tc = CP.tc
  ds = CP.ds

  # compute R matrix
  PhiM = ((cos(phi)*cos(phi)') + (sin(phi)*sin(phi)'))
  Rx = (sin(Lambda)*x' - (x.*sin(Lambda))*ones(1,N)).*PhiM
  Ry = cos(Lambda)*(y.*cos(phi))' - (y.*cos(Lambda))*cos(phi)'
  Rz = cos(Lambda)*(z.*sin(phi))' - (z.*cos(Lambda))*sin(phi)'
  R = Rx + Ry + Rz
  R = R - tril(R) # R(i,j) = 0 for j <= i

  # integrate along structural span
  ds_str = CP.ds'./cos(CP.sweep)

  BMM = rho*U*R.*(ones(N,1)*ds')
  WIC = 2*(ds_str./(tc'.*chord'))*BMM
  return WIC, BMM
end

