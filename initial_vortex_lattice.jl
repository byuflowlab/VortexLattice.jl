type VLM
  surfaces
  panels
  nTotal

  gamma
  Sref
  cref
  # q = 0.5*Surface.rho*Surface.U*Surface.U

  function VLM(s, Sref, cref):
    surfaces = new Array
    surfaces.add(s)
    panels = new Array
    panels.add(s.getNPanels())
    nTotal
  end
end
