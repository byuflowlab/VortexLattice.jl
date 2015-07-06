function getWIC(CPfull,rho,U)
  CP_half = CP_def([0],[0],[0],[0],[0],[0],[0],[0],[0])
  # -------- take only half of wing -----------
  N = length(CPfull.y)
  if N % 2 == 1
    N_half = int(N/2)
  else
    N_half = int(N/2) + 1
  end
  CP_half.x = CPfull.x[N_half:end]
  CP_half.y = CPfull.y[N_half:end]
  CP_half.z = CPfull.z[N_half:end]
  CP_half.chord = CPfull.chord[N_half:end]
  CP_half.tc = CPfull.tc[N_half:end]
  CP_half.dihedral = CPfull.dihedral[N_half:end]
  CP_half.sweep = CPfull.sweep[N_half:end]
  CP_half.ds = CPfull.ds[N_half:end]
  # ------------------------------------------

  # assuming spar located at quarter chord
  x = CP_half.x - 1/2*CP_half.chord
  y = CP_half.y
  z = CP_half.z

  # rename for convenience
  N = length(y)
  phi = CP_half.dihedral
  Lambda = CP_half.sweep
  chord = CP_half.chord
  tc = CP_half.tc
  ds = CP_half.ds

  # compute R matrix
  PhiM = (cos(phi)'*cos(phi) + sin(phi)'*sin(phi))
  Rx = (sin(Lambda).*x' - (x.*sin(Lambda)).*ones(1,N)).*PhiM
  Ry = cos(Lambda)*(y'.*cos(phi)') - (y.*cos(Lambda))*cos(phi)'
  Rz = cos(Lambda)*(z.*sin(phi))' - (z.*cos(Lambda))*sin(phi)'
  R = Rx + Ry + Rz
  R = R - tril(R)   # R(i,j) = 0 for j <= i

  # integrate along structural span
  ds_str = CP_half.ds./cos(CP_half.sweep)

  BMM = rho*U*R.*(ones(1,N)*ds)   # Bending moment matrix * gamma = bending moment dist
  WIC = (ds_str./(tc.*chord))'*BMM

  # add other half
  WIC = [fliplr(WIC) WIC]
  BMM = [flipud(fliplr(BMM)) zeros(size(BMM))
         zeros(size(BMM)) BMM]
  return WIC, BMM
end

