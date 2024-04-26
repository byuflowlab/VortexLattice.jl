"""
    `read_degengeom(filename::String)`

Read all geometry components from a DegenGeom file written out by OpenVSP

**Arguments**
- `filename::String`: DegenGeom filename

**Returns**
- `comp`: Vector of `vsp.VSPComponent` objects
"""
function read_degengeom(filename::String)
    comp = VSPGeom.readDegenGeom(filename)
    return comp
end

"""
    `import_vsp(comp::vsp.VSPComponent; geomType::String="", optargs...)

Imports properties from OpenVSP component to VortexLattice objects. Importing prop and duct geometries are under development.

**Arguments**
- `comp::VSPGeom.VSPComponent`: Single `VSPGeom.VSPComponent` object
- `geomType::String` : Geometry type may be one of - `wing`, `fuselage`, `prop`, `duct`
- `optargs` : Optional arguments that are passed into grid_to_surface_panels() called inside

**Returns**
- `grid`: Array with dimensions (3, i, j) containing the panel corners
- `surface`: Array with dimensions (i, j) containing generated panels
"""
function import_vsp(comp::VSPGeom.VSPComponent; geomType::String="", optargs...)

    # Infer type of geometry from VSPComponent if not specified
    if geomType == ""
        if lowercase(comp.type) == "lifting_surface"
            geomType = "wing"
        elseif lowercase(comp.type) == "body"
            geomType = "fuselage"
        elseif lowercase(comp.type) == "duct"
            geomType = "duct"
        end
    end

    if lowercase(geomType) == "wing"

        nXsecs, npts = VSPGeom.degenGeomSize(comp.plate)

        imax = npts
        jmax = nXsecs

        # 2 endpoints for nXsecs, each with 3 coordinates
        xyz = Array{Float64, 3}(undef, 3, imax, jmax)

        # x = comp.plate.x .+ comp.plate.zCamber .* comp.plate.nCamberx
        xyz[1, :, :] .= reshape(comp.plate.zCamber, imax, jmax)
        xyz[1, :, :] .*= reshape(comp.plate.nCamberx, imax, jmax)
        xyz[1, :, :] .+= reshape(comp.plate.x, imax, jmax)
        reverse!(view(xyz, 1, :, :), dims=1)

        # y = comp.plate.y .+ comp.plate.zCamber .* comp.plate.nCambery
        xyz[2, :, :] .= reshape(comp.plate.zCamber, imax, jmax)
        xyz[2, :, :] .*= reshape(comp.plate.nCambery, imax, jmax)
        xyz[2, :, :] .+= reshape(comp.plate.y, imax, jmax)
        reverse!(view(xyz, 2, :, :), dims=1)

        # z = comp.plate.z .+ comp.plate.zCamber .* comp.plate.nCamberz
        xyz[3, :, :] .= reshape(comp.plate.zCamber, imax, jmax)
        xyz[3, :, :] .*= reshape(comp.plate.nCamberz, imax, jmax)
        xyz[3, :, :] .+= reshape(comp.plate.z, imax, jmax)
        reverse!(view(xyz, 3, :, :), dims=1)

        return grid_to_surface_panels(xyz; optargs...)

    else
        error("The geomType \"$geomType\" is invalid")
        return Nothing
    end
end
