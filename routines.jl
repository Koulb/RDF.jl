using LinearAlgebra
using ProgressMeter
pi = Base.MathConstants.pi

"""
This function calculates the pair correlation function g(2)(r)
for a dense set of particles in a rectangular box. The particle
coordinates are supplied as [[x, y, z],...] lists.
All particles are used as central particles.
The spherical shell might extend beyond the known region.
An empty surrounding is assumed.
"""

function SphereCutVol(Rs::Float64, A::Float64, B::Float64)
    """
    Computes the volume of the wedge left when two planes
    at positions A and B intersect a sphere of radius Rs.

    Parameters
    ----------
    Rs: float, positive
        Sphere radius
    A, B: float, positive
        Distance of each box face from the origin

    Returns
    -------
    Vcut: float, positive
        Volume of the spherical cut
    """
    Root = sqrt(Rs^2 - A^2 - B^2)
    Vcut = 1/6 * Rs^3 * (pi - 2 * atan(A * B / (Rs * Root)))
    Vcut += 1/2 * (atan(A / Root) - pi / 2) * (Rs^2 * B - 1/3 * B^3)
    Vcut += 1/2 * (atan(B / Root) - pi / 2) * (Rs^2 * A - 1/3 * A^3)
    Vcut += 1/3 * A * B * Root
    return Vcut
end


function OctVolume(Rs::Float64, xb::Float64, yb::Float64, zb::Float64)
    """
    Compute the intersection volume between a sphere octant
    and a box corner.

    Parameters
    ----------
    Rs: float, positive
        Sphere radius
    xb, yb, zb: float, positive
        Distance of the box faces from the origin

    Returns
    -------
    VOctant: float, positive
        Intersection volume between the octant
        and the box corner
    """

    # if all boundaries are fully in the octant
    if xb^2 + yb^2 + zb^2 < Rs^2
        return xb * yb * zb
    end

    # if no boundary intersects we start with
    VOctant = 1/8 * 4/3 * pi * Rs^3

    # remove the spherical caps
    for B in [xb, yb, zb]
        if B < Rs
            VOctant -= pi/4 * (2/3 * Rs^3 - B * Rs^2 + 1/3 * B^3)
        end
    end

    # add the intersections of the caps
    for (a, b) in [(xb, yb), (xb, zb), (yb, zb)]
        if a^2 + b^2 < Rs^2
            VOctant += SphereCutVol(Rs, a, b)
        end
    end

    return VOctant
end


function SphereVolume(Rs::Float64, BoxBounds::Vector{Float64})
    """
    Computes the intersection volume of a sphere with a box.

    Parameters
    ----------
    Rs: float, positive
        Sphere radius
    BoxBounds: vector of 6 positive floats
        Distances of the box faces from the origin

    Returns
    -------
    VSphere: float, positive or zero
        Intersection volume of a sphere with a box
    """

    (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax) = BoxBounds

    VSphere = 0
    for xb in [Xmin, Xmax]
        for yb in [Ymin, Ymax]
            for zb in [Zmin, Zmax]
                # abs() mirrors the boundaries into the first octant
                VSphere += OctVolume(Rs, abs(xb), abs(yb), abs(zb))
            end
        end
    end

    return VSphere
end


function ShellVolume(Rmin::Float64, Rmax::Float64, BoxBounds::Vector{Float64})
    """
    Compute the intersection volume of a spherical shell and a box.

    Parameters
    ----------
    Rmin, Rmax: float, positive
        Inner and outer radius of the spherical shell
    BoxBounds: vector of 6 positive floats
        Distances of the box faces from the origin

    Returns
    -------
    Volume: float, positive or zero
        Intersection volume of a spherical shell and a box
    """

    # check for negative Rmin values
    Rmin = max(Rmin, 0)
    InnerShell = SphereVolume(Rmin, BoxBounds)
    OuterShell = SphereVolume(Rmax, BoxBounds)
    Volume = OuterShell - InnerShell

    return Volume
end

function RDF_AnalyticNorm(Particles::Matrix{Float64}, r, dr)
    """
    Computes g(r) from the particle set assuming that 
    the particles are bound by a rectangular box. The
    intersection volume between the radial bins and the 
    box are used to normalize g(r) correctly.

    Parameters
    ----------
    Particles: array of Float64
        Array with the individual particle coordinates
    r: array of Float64
        Center positions of the radial bins
    dr: Float64
        Width / thickness of each radial bin

    Returns
    -------
    Global_Gr: array of Float64
        g(r) values at the corresponding bin positions
    """

    # Gr averaged over all particles
    Global_Gr = zeros(length(r))

    # Keep track of the usefull shell volumes
    NonEmptyShells = zeros(length(r))

    # maximal radial distance
    MaxDist = r[end] + dr / 2

    # Box boundaries, could be supplied to the function
    XList = Particles[:,1]
    YList = Particles[:,2]
    ZList = Particles[:,3]

    # tight box around the particles
    BoxBounds = [minimum(XList), maximum(XList),
                 minimum(YList), maximum(YList),
                 minimum(ZList), maximum(ZList)]

    # box size
    Lx = BoxBounds[2] - BoxBounds[1]
    Ly = BoxBounds[4] - BoxBounds[3]
    Lz = BoxBounds[6] - BoxBounds[5]

    println("Lx = $Lx")
    println("Ly = $Ly")
    println("Lz = $Lz")

    MeanDensity = length(Particles[:,1]) / (Lx * Ly * Lz)

    # use every particle as the center once
    @showprogress 1 "Computing RDF..." for CentralP in 1:length(Particles[:,1])

        # local Gr around the current particle
        Local_Gr = zeros(length(r))

        # look at every other particle at most MaxDist away:
        for Neighbour in 1:length(Particles[:,1])

            if CentralP != Neighbour

                # calc the distance to the neighbour
                dx = Particles[CentralP,1] - Particles[Neighbour,1]
                dy = Particles[CentralP,2] - Particles[Neighbour,2]
                dz = Particles[CentralP,3] - Particles[Neighbour,3]

                d = sqrt(dx^2 + dy^2 + dz^2)

                # what bins is the particle in?
                IdxList = [k for k in 1:length(r) if abs(r[k] - d) <= dr / 2]

                # add one to every bin the particle is in
                for Pos in IdxList

                    # count the particle
                    Local_Gr[Pos] += 1

                end
            end
        end

        # shift the center of box cosy
        LocalBox = [BoxBounds[1] - Particles[CentralP,1],
                    BoxBounds[2] - Particles[CentralP,1],
                    BoxBounds[3] - Particles[CentralP,2],
                    BoxBounds[4] - Particles[CentralP,2],
                    BoxBounds[5] - Particles[CentralP,3],
                    BoxBounds[6] - Particles[CentralP,3]]

        # normalize with the shell volume
        for RIdx in 1:length(r)
            SVolume = ShellVolume(r[RIdx] - dr / 2, r[RIdx] + dr / 2, LocalBox)
            # check for safety
            if SVolume > 0.0
                Local_Gr[RIdx] /= SVolume
                NonEmptyShells[RIdx] += 1
            end
        end

        # normalize by the mean particle density
        Local_Gr .= Local_Gr / MeanDensity
        
        # save in the global g(r) for the average over particles
        Global_Gr .= Global_Gr + Local_Gr

        # println("Finished Particle $CentralP of $(length(Particles[:,1]))")
    end

    # final normalization considering the non empty shell volumes
    for k = 1:length(Global_Gr)
        if NonEmptyShells[k] != 0
            Global_Gr[k] /= NonEmptyShells[k]
        else
            println("All Shells at R = ", r[k], " are Empty!")
        end
    end    
        
    return  Global_Gr   
end             