function getAlpha2(L,LIC,AIC_inv,Vn,CP,U)

  qq = LIC*AIC_inv
  a = -U*qq*sin(CP.twist)'
  b = -U*qq*(cos(CP.twist).*cos(CP.dihedral))'
  c = L + qq*(Vn)'

  alpha = asin(c/b - a/b*(a*c + b*sqrt(a^2+b^2-c^2))/(a^2 + b^2))
  return alpha
end
function getclmax(CLmaxS,CP,LIC,AIC,gamma,rho,U)

  # temporary variables
  q = 1/2*rho*U^2
  tmpV = AIC\cos(CP.dihedral)'

  # add on the basic cl distribution
  CLIC = rho*U*(CLmaxS - LIC*gamma/q)*tmpV./(LIC*tmpV*CP.chord')

  # cl distribution at CLmax
  clmax = 2/U*gamma./CP.chord' + CLIC
  return clmax
end
function getROLLIC(CP,rho,U,cg)
  if (nargin == 3) # assume that cg is at quarter chord if not specified
    N = length(CP.y)
    M = N/2+1 # index for root panel of right half of wing
    cg.y = CP.y[M] - 1/2*CP.ds[M]*cos(CP.dihedral[M])
    cg.z = CP.z[M] - 1/2*CP.ds[M]*sin(CP.dihedral[M])
  end
  ROLLIC = -rho*U*CP.ds.*((CP.y-cg.y).*cos(CP.dihedral) + (CP.z-cg.z).*sin(CP.dihedral))
  return ROLLIC
end
function getYIC(CP,rho,U)
  YIC = -rho*U*sin(CP.dihedral).*CP.ds
  return YIC
end
