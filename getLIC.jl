function getLIC(CP,rho,U)
  LIC = 2*rho*U*cos(CP.dihedral).*CP.ds
  return LIC
end

