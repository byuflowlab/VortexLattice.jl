function getViscous(cd1,cd2,CP,rho,U)
  D1 = cd1*rho*U*CP.ds
  D2 = cd2*2*rho./CP.chord.*CP.ds
  return [D1 D2]
end