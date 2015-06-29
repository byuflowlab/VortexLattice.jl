function Cdrag(CL,Lambda,tc,mach,supercrit)

  cosL = cos(Lambda)
  clp = CL/cosL^2
  tcp = tc/cosL

  # compute Mcc
  Mcc = 0.954-0.235*clp+0.0259*clp^2
  Mcc = Mcc - (1.963-1.078*clp+0.350*clp^2)*tcp
  Mcc = Mcc + (2.969-2.738*clp+1.469*clp^2)*tcp.^2
  Mcc = Mcc + supercrit*.06
  Mcc = Mcc/cosL

  # compute Cdc
  rm = mach / Mcc
  dm = rm-1.

  if rm < .5
      cdc = 0
  elseif (rm >= .5 && rm < .8)
      cdc = 1.3889e-4+5.5556e-4*dm+5.5556e-4*dm*dm
  elseif (rm >= .8 && rm < .95)
      cdc = 7.093e-4+6.733e-3*dm+.01956*dm*dm+.01185*dm*dm*dm
  elseif (rm >= .95 && rm < 1.0)
      cdc = .001000+.02727*dm+.4920*dm*dm+3.573*dm*dm*dm
  elseif (rm >= 1.0 && rm < 1.075)
      cdc = .001000+.02727*dm-.1952*dm*dm+19.09*dm*dm*dm
  else
      cdc = 0.01 + 0.33477*(dm-0.075) # linear extension
  end

  cdc = cdc*cosL^3
  return cdc
end

