function vortex(xA, yA, zA, xB, yB, zB, xC, yC, zC)
  m = length(yC)
  n = length(yA)
  v = zeros(m, n)
  w = zeros(m, n)

  for i = 1:m
    for j = 1:n
      x1 = xC[i] - xA[j]
      x2 = xC[i] - xB[j]
      x21 = xB[j] - xA[j]
      y1 = yC[i] - yA[j]
      y2 = yC[i] - yB[j]
      y21 = yB[j] - yA[j]
      z1 = zC[i] - zA[j]
      z2 = zC[i] - zB[j]
      z21 = zB[j] - zA[j]

      denom1 = (y1*z2-y2*z1)^2 + (x1*z2-x2*z1)^2 + (x1*y2-x2*y1)^2
      j1 = -x1*z2 + x2*z1
      k1 = x1*y2 - x2*y1
      frac1_1 = (x21*x1 + y21*y1 + z21*z1)/sqrt(x1^2 + y1^2 + z1^2)
      frac1_2 = (x21*x2 + y21*y2 + z21*z2)/sqrt(x2^2 + y2^2 + z2^2)
      frac1 = frac1_1 - frac1_2
      v[i,j] = 1/(4*pi)*j1/denom1*frac1
      w[i,j] = 1/(4*pi)*k1/denom1*frac1

    end
  end
  return v, w
end
