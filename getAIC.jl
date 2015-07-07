function getAIC(QC,TE,CP)

  v, w =  velocity3d(QC,TE,CP)

  m = length(CP.y)
  n = length(QC.y)-1

  AIC = zeros(m,n)

  for j=1:n
      AIC[:,j] = -v[:,j].*sin(CP.dihedral)+w[:,j].*cos(CP.dihedral)
  end
  return AIC
end

