using MAT
using GeoGreensFunctions
using Statistics
using DelimitedFiles
import ProgressMeter
include("src_coordinates.jl")

function load_fault(dirs, options)
    ind_last_sep = findlast("/", options["fault"]["fault_file"])[end]
    fault_file = options["fault"]["fault_file"][ind_last_sep+1:end]
    # load fault matrix
    fault_matrix = readdlm(dirs["dir_fault"] * fault_file)
    # total number of patches
    n_patches, n_vertices = size(fault_matrix);

    # each row corresponds to a patch; for each patch there are 3 vertices,
    # and thus 3 values of lon, 3 values of lat, and 3 values of height
    fault = Dict()
    inds_lon = [1,4,7]
    inds_lat = [2,5,8]
    inds_height = [3,6,9]
    fault["lon"]    = fault_matrix[:,inds_lon]
    fault["lat"]    = fault_matrix[:,inds_lat]
    fault["height"] = fault_matrix[:,inds_height]

    fault["llh"] = transpose(stack([mean(fault["lon"], dims=2)[:,1],
        mean(fault["lat"], dims=2)[:,1],
        mean(fault["height"], dims=2)[:,1]],dims=1))
    fault_matrix = nothing

    # locate the origin at the center of the mesh
    origin = transpose([mean(mean(fault["lon"])), mean(mean(fault["lat"]))]);
    fault["origin"] = origin;

    # transform geografic (degrees) coordinates into local (km) coordinates
    fault["xE"] = zeros(n_patches,3);
    fault["yN"] = zeros(n_patches,3);
    fault["zV"] = zeros(n_patches,3);

    for i=1:3
        fault_i = Dict()
        fault_i["lon"] = fault["lon"][:,i];
        fault_i["lat"] = fault["lat"][:,i];
        fault_i["height"] = fault["height"][:,i];
        fault_i = llh2localxyz_geodetic(fault_i, origin);
        fault["xE"][:,i] = fault_i["xE"];
        fault["yN"][:,i] = fault_i["yN"];
        fault["zV"][:,i] = fault_i["zV"];
    end

    fault["xyz"] = transpose(stack([mean(fault["xE"], dims=2)[:,1],
                                    mean(fault["yN"], dims=2)[:,1],
                                    mean(fault["zV"], dims=2)[:,1]],dims=1));


    # find strike and rake of patches
    fault["strike"] = zeros(n_patches,1);
    fault["dip"]    = zeros(n_patches,1);
    fault["area"]   = zeros(n_patches,1);
    for i=1:n_patches
        # A: x, y, and z coordinates of the first point
        # B: x, y, and z coordinates of the second point
        # C: x, y, and z coordinates of the third point
        A = [fault["xE"][i,1], fault["yN"][i,1], fault["zV"][i,1]];
        B = [fault["xE"][i,2], fault["yN"][i,2], fault["zV"][i,2]];
        C = [fault["xE"][i,3], fault["yN"][i,3], fault["zV"][i,3]];
        a,b,c,d = plane_3points(A,B,C);
        fault["strike"][i] = 90 - atand(-a, b);
        # atand returns values between -180 and 180, the strike is thus
        # between -90 and 270
        # if -a/c<0, atand must be between 0 and 180, so the strike must be
        # between -90 and 90
        if -a/c < 0
            if atand(-a, b) < 0
                fault["strike"][i] = fault["strike"][i] + 180.0;
            end
        else
            if atand(-a, b) > 0
                fault["strike"][i] = fault["strike"][i] - 180.0;
            end
        end
        # let us keep the strike between 0 and 360
        fault["strike"][i] = mod(fault["strike"][i], 360);
        fault["dip"][i] = 90.0 - atand(c, sqrt(a^2 + b^2));
        if fault["dip"][i]>90.0
            fault["dip"][i] = 180.0 - fault["dip"][i];
        end
        AB = norm(B-A);
        AC = norm(C-A);
        theta = acosd(dot(B-A, C-A) / (AB * AC));
        fault["area"][i] = 0.5*AB*AC*sind(theta);
    end
    return fault
end

function _find_data_type(X)
    data_type = nothing
    if string(X["name"][1][end]) == "u"
        data_type = "1d";
    elseif string(X["name"][1][end]) == "e"
        if string(X["name"][3][end]) == "e"
            data_type = "2d";
        elseif string(X["name"][3][end]) == "u"
            data_type = "3d";
        end
    end
    return data_type
end

function _get_llh_stations(X, data_type)
    llh_stations = X["llh"][1:end,:];
    if data_type == "2d"
        llh_stations = X["llh"][1:2:end,:];
    elseif data_type == "3d"
        llh_stations = X["llh"][1:3:end,:];
    end
    return llh_stations
end

function create_greens_function(X, fault, options)
    nu = options["fault"]["nu"]
    origin = fault["origin"]

    km2m = 1e3

    n_patches = size(fault["lon"],1)

    v1 = Dict()
    v2 = Dict()
    v3 = Dict()
    
    v1["lon"] = zeros(n_patches,1);
    v1["lat"] = zeros(n_patches,1);
    v1["height"] = zeros(n_patches,1);
    v2["lon"] = zeros(n_patches,1);
    v2["lat"] = zeros(n_patches,1);
    v2["height"] = zeros(n_patches,1);
    v3["lon"] = zeros(n_patches,1);
    v3["lat"] = zeros(n_patches,1);
    v3["height"] = zeros(n_patches,1);
    ind_reorganize = zeros(Int, n_patches,3);

    for j=1:n_patches
        ind_lon = sortperm(fault["lon"][j,:]);
        ind1 = ind_lon[1];
        ind23 = [1,2,3];
        deleteat!(ind23, ind1)
        ind_lat = sortperm(fault["lat"][j,ind23]);
        ind_reorganize[j,1] = ind1;
        ind_reorganize[j,2:3] = ind23[ind_lat];
        v1["lon"][j,:]    .= fault["lon"][j,ind_reorganize[j,1]];
        v1["lat"][j,:]    .= fault["lat"][j,ind_reorganize[j,1]];
        v1["height"][j,:] .= fault["height"][j,ind_reorganize[j,1]];
        v2["lon"][j,:]    .= fault["lon"][j,ind_reorganize[j,2]];
        v2["lat"][j,:]    .= fault["lat"][j,ind_reorganize[j,2]];
        v2["height"][j,:] .= fault["height"][j,ind_reorganize[j,2]];
        v3["lon"][j,:]    .= fault["lon"][j,ind_reorganize[j,3]];
        v3["lat"][j,:]    .= fault["lat"][j,ind_reorganize[j,3]];
        v3["height"][j,:] .= fault["height"][j,ind_reorganize[j,3]];
    end
    v1 = llh2localxyz_geodetic(v1, origin);
    v2 = llh2localxyz_geodetic(v2, origin);
    v3 = llh2localxyz_geodetic(v3, origin);


    data_type = _find_data_type(X)
    llh_stations = _get_llh_stations(X, data_type)
    
    stations = Dict()
    stations["lon"] = llh_stations[:,1]
    stations["lat"] = llh_stations[:,2]
    stations["height"] = llh_stations[:,3]
    stations["origin"] = origin;
    stations = llh2localxyz_geodetic(stations, origin);

    n_stations = size(stations["lon"], 1);

    println("Creating Greens' function:")
    G = zeros(n_stations*3, n_patches*2);
    if options["flags"]["flag_parallel"] > 1

        # Construct a ParforProgressbar object
        progress = ProgressMeter.Progress(
            n_patches; desc = "Calculating Green's functions", enabled = true
            )
        # For every patch...
        Gstrike = zeros(n_stations*3, n_patches);
        Gdip = zeros(n_stations*3, n_patches);
        Threads.@threads for j in 1:n_patches
            Gj = zeros(n_stations*3, 2);
            
            x = [v1["xE"][j], v2["xE"][j], v3["xE"][j]]*km2m;
            y = [v1["yN"][j], v2["yN"][j], v3["yN"][j]]*km2m;
            z = [v1["zV"][j], v2["zV"][j], v3["zV"][j]]*km2m;

            for i = 1:n_stations
                sx = stations["xE"][i]*km2m;
                sy = stations["yN"][i]*km2m;
                sz = 0.0; #stations.zV(i)*km2m;
                ue_ss,un_ss,uv_ss = GeoGreensFunctions.disp_tri3_hs(sx,sy,sz,
                    [x[1], y[1], z[1]],
                    [x[2], y[2], z[2]],
                    [x[3], y[3], z[3]],1.0,0.0,0.0,nu);
                ue_ds,un_ds,uv_ds = GeoGreensFunctions.disp_tri3_hs(sx,sy,sz,
                    [x[1], y[1], z[1]],
                    [x[2], y[2], z[2]],
                    [x[3], y[3], z[3]],0.0,1.0,0.0,nu);
                
                Gj[i*3-2, 1] = ue_ss;
                Gj[i*3-1, 1] = un_ss;
                Gj[i*3, 1]   = uv_ss;
    
                Gj[i*3-2, 2] = ue_ds;
                Gj[i*3-1, 2] = un_ds;
                Gj[i*3, 2]   = uv_ds;
            end
            Gstrike[:,j] = Gj[:,1];
            Gdip[:,j]    = Gj[:,2];
            ProgressMeter.next!(progress)
        end
        G[:,1:n_patches] = Gstrike;
        G[:,n_patches+1:2*n_patches] = Gdip;
    else
        # For every patch...
        Gstrike = zeros(n_stations*3, n_patches);
        Gdip = zeros(n_stations*3, n_patches);
        for j in 1:n_patches
            Gj = zeros(n_stations*3, 2);
            
            x = [v1["xE"][j], v2["xE"][j], v3["xE"][j]]*km2m;
            y = [v1["yN"][j], v2["yN"][j], v3["yN"][j]]*km2m;
            z = [v1["zV"][j], v2["zV"][j], v3["zV"][j]]*km2m;

            for i = 1:n_stations
                sx = stations["xE"][i]*km2m;
                sy = stations["yN"][i]*km2m;
                sz = 0.0; #stations.zV(i)*km2m;
                ue_ss,un_ss,uv_ss = GeoGreensFunctions.disp_tri3_hs(sx,sy,sz,
                    [x[1], y[1], z[1]],
                    [x[2], y[2], z[2]],
                    [x[3], y[3], z[3]],1.0,0.0,0.0,nu);
                ue_ds,un_ds,uv_ds = GeoGreensFunctions.disp_tri3_hs(sx,sy,sz,
                    [x[1], y[1], z[1]],
                    [x[2], y[2], z[2]],
                    [x[3], y[3], z[3]],0.0,1.0,0.0,nu);
                
                Gj[i*3-2, 1] = ue_ss;
                Gj[i*3-1, 1] = un_ss;
                Gj[i*3, 1]   = uv_ss;
    
                Gj[i*3-2, 2] = ue_ds;
                Gj[i*3-1, 2] = un_ds;
                Gj[i*3, 2]   = uv_ds;
            end
            Gstrike[:,j] = Gj[:,1];
            Gdip[:,j]    = Gj[:,2];
            ProgressMeter.next!(progress)
        end
        G[:,1:n_patches] = Gstrike;
        G[:,n_patches+1:2*n_patches] = Gdip;
    end
    println("Done")
    return G
end

function downsample_V(decomp, Δt)
    decomp_downsampled = decomp;
    if Δt != 1
        n_samples, n_features = size(decomp["V"])
        
        n_t = Int(round(n_samples/Δt))

        t = zeros(1,n_t)
        V = zeros(n_t,n_features)
        var_V = zeros(n_t,n_features)
        Δt_mid = (Δt-1)/2
        for i=1:n_t
            j0 = Int((i-1)*Δt + 1)
            j1 = Int((i-1)*Δt + Δt)
            j_mid = Int((i-1)*Δt + Δt_mid)
            if i<n_t
                t[1,i] = decomp["timeline"][j_mid]
                V[i,:] = mean(decomp["V"][j0:j1, :], dims=1)
                var_V[i,:] = var(decomp["V"][j0:j1, :], dims=1)
            else
                t[1,i] = mean(decomp["timeline"][j0:end])
                V[i,:] = mean(decomp["V"][j0:end, :], dims=1)
                var_V[i,:] = var(decomp["V"][j0:end, :], dims=1)
            end
        end
        decomp_downsampled["timeline"] = t;
        decomp_downsampled["V"] = V;
        decomp_downsampled["var_V"] = var_V;
    end
    return decomp_downsampled
end