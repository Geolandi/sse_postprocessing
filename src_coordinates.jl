# function save_variable(filename, varname, var)
#     FileIO.save(filename, varname, var)
#     return nothing
# end


function polyconic(Lat, Diff_long, Lat_Orig)
	
    #	Polyconic Projection
    #	Input:
    #		Lat 		= Latitude (decimal seconds)
    #		Diff_long 	= Differential Longitude (decimal seconds)
    #					relative to Central Meridian
    #		Lat_Orig	= Latitude of Origin (decimal seconds)
    #	Output: x =	Distance from Central Meridian
    #		y = 	Distance from Origin to Latitude
    
    p1 = Lat_Orig;
    p2 = Lat;
    il = Diff_long;

    arcone = 4.8481368e-6;
    esq = 6.7686580e-3;
    la = 6378206.4;
    a0 = 6367399.7;
    a2 = 32433.888;
    a4 = 34.4187;
    a6 = .0454;
    a8 = 6.0e-5;

    ip = p2 - p1;
    sinp2 = sin(p2 * arcone);
    cosp2 = cos(p2 * arcone);
    theta = il * sinp2;
    a = sqrt(1.0 - (esq * (2. * sinp2))) / (la * arcone);
    cot = cosp2 / sinp2;
    x = (cot * sin(theta * arcone)) / (a * arcone);
    ipr = ip * arcone;
    pr = ((p2 + p1) / 2.)*arcone;
    y = ((((a0*ipr) - ((a2*cos(2.0*pr))*sin(ipr))) + 
        ((a4*cos(4.0*pr))*sin(2.0*ipr))) - ((a6*cos(6.0*pr))*sin(3.0*ipr))) +
        ((a8*cos(8.0*pr))*sin(4.0*ipr));
    xy = [x, y];

    return xy
end


function llh2localxy(llh,ll_org)
    #LLH2LOCALXY   Convert lat-lon-height to local xy via a given origin
    #   [XY] = LLH2LOCALXY(LLH,LL_ORG) converts LLH matrices
    #   [Lat;Lon;Height] to locat coordinates XY = [X,Y] based on the
    #   origin vector LL_ORG = [LAT_ORIGIN, LONG_ORIGIN].
    #
    #   See also local2llh, PCAIM_driver.
    
    rows, nsta = size(llh);
    # change from decimal degrees to decimal seconds
    lat = 3600.0 * llh[1,:];
    lon = 3600.0 * llh[2,:];
    
    Lat_Orig = 3600.0 * ll_org[1];
    Diff_long  = 3600.0 * ll_org[2]*ones(size(lon)) - lon;

    xy = zeros(nsta,2);
    for i=1:nsta
        xy[i,:] = polyconic(lat[i], Diff_long[i], Lat_Orig);
    end
    # convert units from meter into kilometer and flip x-axis
    xy[:,1] = -xy[:,1] / 1000.0;
    xy[:,2] =  xy[:,2] / 1000.0;
    return xy
end


function llh2localxyz(lat,lon,height,origin)

    # llh2localxyz.m converts the geographic coordinates of a list of points
    # into local coordinates.
    # -------------------------------------------------------------------------
    # INPUT
    # lat: latitude (degrees)
    # lon: longitude (degrees)
    # height: altitude (km)
    # origin: [longitude, latitude] (degrees) of the origin of the new local
    #         coordinate reference system
    # -------------------------------------------------------------------------
    # OUTPUT
    # x: x in the local reference system (km), i.e. along the East direction
    # y: y in the local reference system (km), i.e. along the North direction
    # z: z in the local reference system (km), i.e. along the Vertical
    #    direction
    # -------------------------------------------------------------------------
    #
    # Adriano Gualandi - 19 Aug 2016
    # California Institute of Technology
    # Geological and Planetary Science Division
    

    xy = (llh2localxy(stack([lat, lon, height], dims=1),
        reverse(origin, dims = 2)));
    x = xy[:,1];
    y = xy[:,2];
    z = height;
    return x,y,z
end


function llh2localxyz_geodetic(geodetic_input,origin)

    # llh2localxyz_geodetic.m add to the geodetic structure the xE, yN, and zV
    # fields containing the catalog coordinate in the local reference system.
    # -------------------------------------------------------------------------
    # INPUT
    # geodetic_input: structure that must have the following fields:
    #   lat: latitude of the geodetic data (degrees)
    #   lon: longitude of the geodetic data (degrees)
    #   height: height of the geodetic data (km)
    # origin: [longitude, latitude] (degrees) of the origin of the new local
    #         coordinate reference system
    # -------------------------------------------------------------------------
    # OUTPUT
    # geodetic_output: structure equal to geodetic_input but with the
    # additional fields:
    #   xE: East location of the geodetic data in the local reference system
    #       (km)
    #   yN: North location of the geodetic data in the local reference system
    #       (km)
    #   zV: Vertical location of the geodetic data in the local reference
    #       system (km)
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi - 26 Aug 2016
    # California Institute of Technology
    # Geological and Planetary Science Division
    
    geodetic_output = geodetic_input;

    xE,yN,zV = llh2localxyz(geodetic_input["lat"],
                            geodetic_input["lon"],
                            geodetic_input["height"],
                            origin
                            );

    geodetic_output["xE"] = xE
    geodetic_output["yN"] = yN
    geodetic_output["zV"] = zV
    return geodetic_output
end

function llh2localxyz_seismicity(seismicity_input,origin)
    
    # llh2localxyz_seismicity.m add to the seismicity structure the x, y, and z
    # fields containing the catalog coordinate in the local reference system.
    # -------------------------------------------------------------------------
    # INPUT
    # seismicity_input: structure that must have the following fields:
    #   lat: latitude of the seismic catalog (degrees)
    #   lon: longitude of the seismic catalog (degrees)
    #   depth: depth of the seismic catalog (km)
    # origin: [longitude, latitude] (degrees) of the origin of the new local
    #         coordinate reference system
    # -------------------------------------------------------------------------
    # OUTPUT
    # seismicity_output: structure equal to seismicity_input but with the
    # additional fields:
    #   xE: East location of the seismic catalog in the local reference system
    #       (km)
    #   yN: North location of the seismic catalog in the local reference system
    #       (km)
    #   zV: Vertical location of the seismic catalog in the local reference
    #       system (km)
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi - 26 Aug 2016
    # California Institute of Technology
    # Geological and Planetary Science Division

    seismicity_output = seismicity_input;
    seismicity_output["xE"],seismicity_output["yN"],seismicity_output["zV"] = 
        llh2localxyz(seismicity_input["lat"],seismicity_input["lon"],
                    -seismicity_input["depth"], origin);
    return seismicity_output
end

function fault2local(vectors_fault,strike,dip)
    # fault2local takes as input n_vectors vectors in the Okada reference
    # frame, where x is along strike, y is along dip (updip), and z is
    # orthogonal to the fault plane, and transforms them in the local reference
    # frame where x is along East, y along North, and z is orthogonal to the
    # ground surface pointing outward the center of the Earth.
    # -------------------------------------------------------------------------
    # INPUT
    # vectors_fault: strike, updip, and tensile components of the vectors to be
    #                transformed; size: 3 x n_vectors.
    # strike: strike angle of the fault, defined starting from the North
    #         direction, increasing clockwise, until the strike direction is
    #         reached; size: scalar.
    # dip: dip angle of the fault, defined as the angle between the plane
    #      parallel to the ground and the updip direction of the fault; size:
    #      scalar.
    # -------------------------------------------------------------------------
    # OUTPUT
    # vectors_local: East, North, and Vertical components of the transformed
    #                vectors; size: 3 x n_vectors.
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi - 25 Aug 2016
    # California Institute of Technology
    # Geological and Planetary Science Division
    
    # Let us call [x;y;z] the system from which we start, i.e. the fault
    # reference system where x is along the strike direction, y is along the
    # updip direction, and z is along the tensile direction.
    # Let us call [x'';y'';z''] the system to which we end, i.e. the local
    # geographic reference system where x'' is along the East direction, y'' is
    # along the North direction, and z'' is along the Vertical outward
    # direction.
    # The transformation consists of two rotations. After the first rotation we
    # will be in the reference system defined by the intermediate coordinates
    # [x';y';z'].
    
    # The x and x'' directions, i.e. the strike and East directions, lie on a
    # plane parallel to the ground surface, i.e. parallel to the plane z''=0.
    # For this reason as first rotation we align the z axis with the z'' axis,
    # so that then we can rotate around z'(=z'') to align x' with x''. To do
    # that we have to perform the rotation around the x axis, so that x'=x. The
    # dip angle is defined as that angle that is formed between the plane
    # parallel to the ground and the fault plane. It follows that an analogue
    # definition can be the angle that goes from the vertical direction to the
    # tensile direction. Thus, we have to perform a rotation of -dip to align
    # the tensile axis z with the vertical axis z''. This rotation is performed
    # around the strike direction x.
    angle1 = -dip;
    
    # After the first rotation we now have to align the x'=x strike direction
    # with the East direction x''. We already know that now z has been aligned
    # with the Vertical direction z'', so we do not have to modify z'. Thus, we
    # can rotate around the axis z'=z''. The rotation must be of an angle equal
    # to the angle that separates the strike and East directions. The angle
    # between the strike direction and the East is given by (strike-90).
    # Indeed, the strike angle is increasing clockwise from the North. If we
    # subtract 90 degrees we get the angle between the East direction and the
    # strike direction clockwise increasing, i.e. moving the East axis towards
    # the strike direction. The rotation must be from the strike direction and
    # the East direction of the same angle, that will now increase
    # counterclockwise.
    angle2 = (strike-90);
    
    # The transformation matrices are calculated using the
    # rotation_reference_system function.
    # Rotation around the strike direction (x axis) of an angle equal to angle1
    R1 = rotation_reference_system("x",angle1);
    # Rotation around the Vertical direction (z'' axis) of an angle equal to
    # angle2
    R2 = rotation_reference_system("z",angle2);
    
    # The overall transformation is:
    # vectors_local = R2*R1*vectors_fault
    # where vectors_local and vectors_fault are the vectors in the local and
    # fault reference systems.
    vectors_local = R2*R1*vectors_fault;

    return vectors_local
end

function rotation_reference_system(axis,angle)
    # rotation_reference_system.m gives the linear operator that allows to
    # perform a rotation of the reference system around a given axis of a given
    # angle.
    # -------------------------------------------------------------------------
    # INPUT
    # axis: string that can assume the value 'x', 'y', or 'z'
    # angle: rotation angle in degrees, increasing counterclockwise in a right
    # hand system
    # -------------------------------------------------------------------------
    # OUTPUT
    # R: linear operator that can be applied to a 3-d vector to change the
    # reference system in the rotated one
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi
    # California Institute of Technology
    # Geological and Planetary Science Division
    
    if axis == "x"
        R = [ 1.0   0.0           0.0;
              0.0   cosd(angle)   sind(angle);
              0.0  -sind(angle)   cosd(angle)];
    elseif axis == "y"
        R = [cosd(angle)   0.0  -sind(angle);
                0.0        1.0   0.0;
             sind(angle)   0.0   cosd(angle)];
    elseif axis == "z"
        R = [ cosd(angle)   sind(angle)   0.0;
             -sind(angle)   cosd(angle)   0.0;
                0.0             0.0       1.0];
    else
        error("Please, select a string between 'x', 'y', and 'z' for the axis variable.");
    end
    
    return R
end

function plane_3points(A,B,C)
    # plane_3points.m calculates the coefficients of the general equation of a
    # plane passing through three non-aligned points A, B, and C. The equation
    # of the plane is;
    # 
    # a*x + b*y + c*z + d = 0                                           (eq. 1)
    # 
    # The coefficients are calculated using the following condition:
    # 
    #                                                                   (eq. 2)
    #     | x-xA,  y-yA,  z-zA|
    # det |xB-xA, yB-yA, zB-zA| = 0
    #     |xC-xA, yC-yA, zC-zA|
    # 
    # The coefficients a, b, and c represents the vector normal to the plane.
    # For this reason the output is given normalizing the coefficients a, b, c,
    # and d by the norm of the vector [a,b,c]'.
    # -------------------------------------------------------------------------
    # INPUT
    # A: x, y, and z coordinates of the first point
    # B: x, y, and z coordinates of the second point
    # C: x, y, and z coordinates of the third point
    # -------------------------------------------------------------------------
    # OUTPUT
    # a: First coefficient of (eq. 1)
    # b: Second coefficient of (eq. 1)
    # c: Third coefficient of (eq. 1)
    # d: Fourth coefficient of (eq. 1)
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi - 29 Oct 2016
    # California Institute of Technology
    # Geological and Planetary Sciences Division
    
    # Solve (eq. 1)
    a = (B[2]-A[2])*(C[3]-A[3]) - (B[3]-A[3])*(C[2]-A[2]);
    b = (B[3]-A[3])*(C[1]-A[1]) - (B[1]-A[1])*(C[3]-A[3]);
    c = (B[1]-A[1])*(C[2]-A[2]) - (B[2]-A[2])*(C[1]-A[1]);
    d = -a*A[1]-b*A[2]-c*A[3];
    # Norm of the normal vector
    n = norm([a,b,c]);
    # Normalize the coefficients
    a = a/n;
    b = b/n;
    c = c/n;
    d = d/n;

    return a, b, c, d
end