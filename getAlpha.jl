function getAlpha(Ltarget,CP,LIC,AIC,U)
  # solve quadratic equation for alpha
  qq = LIC/AIC;
  a = -U*qq*sin(CP.twist)';
  b = -U*qq*(cos(CP.twist).*cos(CP.dihedral))';
  c = Ltarget;
  # alpha = acos((a*c + b*sqrt(a^2+b^2-c^2))/(a^2+b^2));
  alpha = asin(c/b - a/b*(a*c + b*sqrt(a^2+b^2-c^2))/(a^2 + b^2));
  return alpha
end