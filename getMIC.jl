function getMIC(CP,rho,U,cg)

  xQC = CP.x - 1/2*CP.chord # quarter chord locations
  MIC = -(xQC - cg)*2*rho*U.*cos(CP.dihedral)'.*CP.ds
  return MIC
end

