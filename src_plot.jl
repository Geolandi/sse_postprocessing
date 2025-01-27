using GLMakie
using GeometryBasics
using DSP
using GMT
using Dates
using StatsBase
include("src_selections.jl")

function plot_intro_map_gmt(X, fault, options)
    figsize = options["figsize"]
    limits = options["limits"]
    # proj = options["proj"]
    title = options["title"]
    ll_volcanoes = options["volcanoes"]
    boundaries_file = options["tect_boundaries"]
    stns2plot = options["stns2plot"]
    dir_output = options["dir_output"]
    output = options["output"]
    
    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    n_stns2plot = length(stns2plot)
    llh_stns = zeros(n_stns2plot,3)
    for i=1:n_stns2plot
        name = stns2plot[i]
        inds = [n[1:4] for n in X["name"]] .== name
        llh_stns[i,:] = X["llh"][inds[:],:][1,:]
    end

    vertices = get_vertices_gmt(fault)
    D = gmtread(boundaries_file);
    
    GMT.basemap(region=limits,
                proj=:Mercator,
                axes=:WSne,
                xaxis=(annot=2,ticks=1),
                yaxis=(annot=2,ticks=1),
                par=(:MAP_FRAME_TYPE,"fancy+"),
                show=false
                )
    
    
    # lc: line color; lw: line width
    GMT.imshow(D, region=limits, proj=:Mercator, axes=:WSne,
    lc="193/43/43", lw=2, show=false)

    # depth = -fault["llh"][:,3]
    # C = makecpt(cmap=:hot,
    #             inverse=true,
    #             range=(minimum(depth),
    #                     maximum(depth)));
    
    # GMT.colorbar!(pos=(inside=true, anchor=:BL, offset=(1,1), size=(4,0.4), horizontal=true),
    #             frame=(xlabel="Depth (km)",),
    #             xaxis=(annot=:auto, ticks=:auto),
    #             par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
    #                     (FONT_LABEL="16p,Helvetica,black"),),
    #             cmap=C,
    #             )
    # GMT.plot!(vertices,
    #             region=limits,
    #             level=(data=depth),
    #             pen="0.5p,gray",
    #             show=false
    #             )
    
    GMT.plot!(vertices,
                region=limits,
                pen="0.5p,gray",
                show=false
                )

    cont = copy(fault["llh"][:,:])
    cont[:,3] = -cont[:,3]
    GMT.contour!(cont,
                proj=:Mercator,
                limits=limits,
                pen="1p,30/150/65",
                annot=[i for i=15:15:80],
                show=false
                )
    
    GMT.coast!(shorelines="1/1p,black",
                resolution=:i,
                par=(:MAP_FRAME_TYPE,"fancy+"),
                show=false
                )
                
    GMT.arrows!([-128 48 20 1.4], 
        pen=(1,:black),
        fill="193/43/43",
        arrow4=(align=:middle,
                head=(arrowwidth="4p", headlength="18p", headwidth="7.5p"))
        )

    GMT.text!("40-45 mm/yr", x=-127.7, y=47.5,
                par=((FONT_ANNOT_PRIMARY="12p,Helvetica,193/43/43"),
                (FONT_LABEL="12p,Helvetica,193/43/43"),),
                justify=:CB, show=false)


    GMT.arrows!([-126 41 5 0.87], 
        pen=(1,:black),
        fill="193/43/43",
        arrow4=(align=:middle,
                head=(arrowwidth="4p", headlength="18p", headwidth="7.5p"))
        )
    GMT.text!("25 mm/yr", x=-126.2, y=41.3,
                par=((FONT_ANNOT_PRIMARY="12p,Helvetica,193/43/43"),
                (FONT_LABEL="12p,Helvetica,193/43/43"),),
                justify=:CB, show=false)

    GMT.arrows!([-123.95 39.9 -58 0.8], pen=2,
        arrow=(len=0.5, half=:left, angle=30), fill="193/43/43")
    GMT.arrows!([-123.95 39.1 122 0.8], pen=2,
        arrow=(len=0.5, half=:left, angle=30), fill="193/43/43")
    
    GMT.text!("North", x=limits[2]-4.4, y=limits[4]-1.3,
                par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                (FONT_LABEL="16p,Helvetica,black"),),
                justify=:CB, show=false)
    GMT.text!("America", x=limits[2]-4.4, y=limits[4]-1.6,
                par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                (FONT_LABEL="16p,Helvetica,black"),),
                justify=:CB, show=false)
    GMT.text!("plate", x=limits[2]-4.4, y=limits[4]-1.9,
                par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                (FONT_LABEL="16p,Helvetica,black"),),
                justify=:CB, show=false)
    GMT.text!("Juan de Fuca", x=limits[1]+1.7, y=limits[3]+7.0,
                par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                (FONT_LABEL="16p,Helvetica,black"),),
                justify=:CB, show=false)
    GMT.text!("plate", x=limits[1]+1.7, y=limits[3]+6.7,
                par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                (FONT_LABEL="16p,Helvetica,black"),),
                justify=:CB, show=false)

    
    GMT.scatter!(X["llh"][1:3:end,1], X["llh"][2:3:end,2],
            #color="43/43/193",
            color="0/84/147",
            markersize=0.3,
            markeredgecolor=:black,
            show=false
            )

    GMT.scatter!(ll_volcanoes[:,1],ll_volcanoes[:,2],
            color="255/165/0",
            #color="30/150/65",
            marker=:triangle,
            markersize=0.3,
            markeredgecolor=:black,
            show=false,
            )

    GMT.scatter!(llh_stns[:,1], llh_stns[:,2],
            marker=:diamond,
            color="232/56/137",
            markersize=0.5,
            markeredgecolor=:black,
            show=false
            )

    GMT.basemap!(
            inset = (; anchor = :TR, size=(4,2), offset = (-0.1,-0.1),
                        save="xx000.dat"),
            box = (; fill = :white, pen = 0, clearance = 0.1, shaded = true),
            show=false,
            #savefig=output*".png"
            )
    t = gmtread("xx000.dat").data
    GMT.coast!(region="-180/-20/0/70",
            proj=:Mercator,
            xaxis=(annot=30,grid=10),
            yaxis=(annot=20,grid=10),
            par=(:MAP_FRAME_TYPE,"fancy+"),
            land=:lightgray,
            shore=:thinnest,
            figsize=t[3],
            xshift=t[1], yshift=t[2],
            show=false,
            )
    rect = [limits[1] limits[3];
            limits[1] limits[4];
            limits[2] limits[4];
            limits[2] limits[3];
            limits[1] limits[3]]
    GMT.plot!(rect,
            region="-180/-20/0/70", lw=2, lc=:blue,
            show=false,
            savefig=dir_output * output * ".png")
    rm("xx000.dat")
    if isfile("./gmt.history")
        rm("./gmt.history")
    end
    return nothing
end


function plot_ts_gmt(X, options)
    figsize = options["figsize"]
    name = options["name"]
    dir_output = options["dir_output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    timeline = X["timeline"][1,:]
    t0 = timeline[1]
    t1 = timeline[end]
    n_samples = length(timeline)

    inds = [n[1:4] for n in X["name"]] .== name
    enu = X["ts"][inds[:],:]
    var_enu = X["var_ts"][inds[:],:]
    e = enu[1,:]
    n = enu[2,:]
    u = enu[3,:]
    var_e = var_enu[1,:]
    var_n = var_enu[2,:]
    var_u = var_enu[3,:]
    
    Ae = maximum(e) - minimum(e)
    ye_min = minimum(e) - 0.1*Ae
    ye_max = maximum(e) + 0.1*Ae
    An = maximum(n) - minimum(n)
    yn_min = minimum(n) - 0.1*An
    yn_max = maximum(n) + 0.1*An
    Au = maximum(u) - minimum(u)
    yu_min = minimum(u) - 0.1*Au
    yu_max = maximum(u) + 0.1*Au

    e_vars = [timeline, e, zeros(n_samples), sqrt.(var_e)]
    n_vars = [timeline, n, zeros(n_samples), sqrt.(var_n)]
    u_vars = [timeline, u, zeros(n_samples), sqrt.(var_u)]

    subplot(grid=(3,1), dims=(size=(15,15), frac=((1,), (1,1,1)) ),
            title=name, row_axes=(left=false), col_axes=(bott=false),
            margins=0, autolabel=true,
            savefig=dir_output * name * ".png")
        GMT.basemap(region=(t0,t1,ye_min,ye_max),
                    xaxis=(annot=5,ticks=1),
                    yaxis=(annot=5,ticks=1),
                    par=(:MAP_FRAME_TYPE,"fancy+"),
                    axes=:WsNe,
                    ylabel="East (mm)",
                    panel=(1,1)
                    )
        GMT.plot(e_vars,
                error_bars=(x=:x, pen=1),
                )
        GMT.scatter(timeline, e,
                color=:red,
                markersize=0.02
                )

        GMT.basemap(region=(t0,t1,yn_min,yn_max),
                xaxis=(annot=5,ticks=1),
                yaxis=(annot=5,ticks=1),
                par=(:MAP_FRAME_TYPE,"fancy+"),
                axes=:Wsne,
                ylabel="North (mm)",
                panel=(2,1)
                )
        GMT.plot(n_vars,
                error_bars=(x=:x, pen=1),
                )
        GMT.scatter(timeline, n,
                color=:red,
                markersize=0.02
                )

        GMT.basemap(region=(t0,t1,yu_min,yu_max),
                xaxis=(annot=5,ticks=1),
                yaxis=(annot=5,ticks=1),
                par=(:MAP_FRAME_TYPE,"fancy+"),
                axes=:WSne,
                xlabel="Time (yr)",
                ylabel="Vertical (mm)",
                panel=(3,1)
                )
        GMT.plot(u_vars,
                error_bars=(x=:x, pen=1),
                )
        GMT.scatter(timeline, u,
                color=:red,
                markersize=0.02
                )
    subplot()
    
    if isfile("./gmt.history")
        rm("./gmt.history")
    end

    return nothing
end


function plot_comp_gmt(decomp, options)
    figsize = options["figsize"]
    limits = options["limits"]
    xy_legend = options["ll_legend"]
    proj = options["proj"]
    title = options["title"]
    show_progress = options["show_progress"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    x = decomp["llh"][1:3:end,1]
    y = decomp["llh"][1:3:end,2]
    append!(x, xy_legend[1])
    append!(y, xy_legend[2])
    
    
    timeline = decomp["timeline"][1,:]
    t0 = timeline[1]
    t1 = timeline[end]
    U = decomp["U"]
    S = decomp["S"]
    V = decomp["V"]
    V_var = decomp["var_V"]
    min_V = minimum(V, dims=1)
    max_V = maximum(V, dims=1)

    norm_facts = max_V .- min_V
    temp_funcs = (V .- min_V) ./ norm_facts
    temp_funcs_err = sqrt.(decomp["var_V"]) ./ norm_facts
    
    spat_maps = (U * S) .* norm_facts
    n_samples = length(timeline)
    n_comps = size(decomp["S"])[1]
    
    progress = ProgressMeter.Progress(
        n_comps; desc = "Plotting components: ", enabled = show_progress
    )

    for i=1:n_comps
        
        time_func = temp_funcs[:,i]
        time_func_err = temp_funcs_err[:,i]

        time_plot_vars = [timeline, time_func, zeros(n_samples), time_func_err]
        
        spat_map = spat_maps[:,i]

        u = spat_map[1:3:end]
        v = spat_map[2:3:end]
        w = spat_map[3:3:end]
        u_legend = Integer(round(maximum(sqrt.(u.^2 + v.^2)) / 2))
        append!(u, u_legend)
        append!(v, 0)
        
        n_stn = length(w)
        err_u = zeros(n_stn+1)
        err_v = zeros(n_stn+1)
        
        subplot(grid=(2,1), dims=(size=(7.1,15), frac=((1,), (1,3)) ),
            title=title*string(i), row_axes=(left=false), col_axes=(bott=false),
            margins=0, autolabel=true,
            savefig=dir_output*output*string(i)*".png",
            par=((FONT_HEADING="24p,Helvetica,black"),)
            )
            GMT.basemap(region=(t0,t1,-0.1,1.1),
                    xaxis=(annot=5,ticks=1),
                    yaxis=(annot=0.5,ticks=0.1),
                    par=(:MAP_FRAME_TYPE,"fancy+"),
                    axes=:WsNe,
                    ylabel="Time function",
                    panel=(1,1)
                    )
            GMT.plot(time_plot_vars,
                    error_bars=(x=:x, pen=1),
                    #xlabel="Time (yr)",
                    )
            GMT.scatter(timeline, time_func,
                    color=:red,
                    markersize=0.02
                    #xlabel="Time (yr)",
                    )
            GMT.basemap(region=limits,
                    proj=:Mercator,
                    axes=:WSne,
                    xaxis=(annot=2,ticks=1),
                    yaxis=(annot=2,ticks=1),
                    par=(:MAP_FRAME_TYPE,"fancy+"),
                    panel=(2,1)
                    )
            GMT.coast(shorelines="1/0.5p,black",
                    resolution=:i,
                    )
            C = makecpt(cmap=:polar,
                        range=(-maximum(abs.(spat_map)),
                                maximum(abs.(spat_map))));
            GMT.scatter(x[1:end-1],y[1:end-1],
                zcolor=w,              # Assign color to each symbol
                size=0.25,                  # The symbl sizes
                markeredgecolor=:black,
                cmap=C
                #panel=(2,1)
                )
            
            GMT.colorbar!(pos=(justify=:BL, offset=(-6.7,-8.5), size=(4,0.4)),
                        frame=(xlabel="Vertical (mm)",),
                        xaxis=(annot=:auto, ticks=:auto),
                        par=((FONT_ANNOT_PRIMARY="16p,Helvetica,black"),
                                (FONT_LABEL="16p,Helvetica,black"),),
                        cmap=C,
                        )
                        
            velo_vecs = hcat(x, y, u, v, err_u, err_v, zeros(n_stn+1))
            GMT.velo(mat2ds(velo_vecs, repeat([" "], n_stn+1)),
                    pen=(0.6,:black),
                    Se="0.2/0.39/18",
                    #arrow="0.3c+p1p+e+gblack+n",
                    arrow="0.3c+p1p+e+gblack",
                    region=limits,
                    )
            GMT.text([x[end]+0.75 y[end]+0.7], text="Horizontal",
                par=((FONT_ANNOT_PRIMARY="8p,Helvetica,black"),
                (FONT_LABEL="8p,Helvetica,black"),), justify=:CB)
            GMT.text([x[end]+0.75 y[end]+0.3], text=string(u_legend)*" mm",
                par=((FONT_ANNOT_PRIMARY="8p,Helvetica,black"),
                (FONT_LABEL="8p,Helvetica,black"),), justify=:CB)
        subplot()
        ProgressMeter.next!(progress)
    end

    if isfile("./gmt.history")
        rm("./gmt.history")
    end

    return nothing
end

function plot_select_comps_and_smoothing(f, cs, misfit_comps, options)

    sigmas = options["sigma"]
    ind_comps = options["ind_comps"]
    ind_comps_complex = options["ind_comps_complex"]
    ind_comps_complex_removed = options["ind_comps_complex_removed"]
    ind_comps_common = options["ind_comps_common"]
    ind_comps_common_removed = options["ind_comps_common_removed"]
    
    titlea = options["titlea"]
    titleb = options["titleb"]
    title = options["title"]
    dir_output = options["dir_output"]
    output = options["output"]

    n_comp = size(cs)[1]
    n_comps_inverted = size(misfit_comps)[2]
    x_min_axb = minimum(sigmas)
    x_max_axb = maximum(sigmas)
    y_min_axb = minimum(misfit_comps)
    y_max_axb = maximum(misfit_comps)

    x_min_axa = minimum(f)
    x_max_axa = maximum(f)
    y_min_axa = minimum(cs)
    y_max_axa = maximum(cs)

    subplot(grid=(1,2), dims=(size=(30,15),
        frac=(ntuple(x -> 1, 2), ntuple(x -> 1, 1)) ),
            title=title, row_axes=(left=false), col_axes=(bott=false),
            margins=1, autolabel=true, savefig=dir_output*output*".png")

        GMT.basemap(region=(x_min_axa, x_max_axa, 0.9*y_min_axa, y_max_axa*1.1),
                axes=:WSne,
                xlabel="Frequency (1/yr)",
                ylabel="Cumulative PSD",
                # xaxis=(annot=0.5),
                # yaxis=(annot=0.2),
                title=titlea,
                panel=(1,1)
                )
        
        for i=1:n_comp
            GMT.plot!(f, cs[i,:], lw=0.5, lc=:black)
        end
        n_ind_comps = length(ind_comps)
        # for i=1:n_ind_comps
        #     GMT.plot!(f, cs[ind_comps[i],:], pen="2p,193/43/43")
        # end
        n_ind_comps_complex = length(ind_comps_complex)
        for i=1:n_ind_comps_complex
            GMT.plot!(f, cs[ind_comps_complex[i],:],
                    pen="2p,193/43/43")
        end
        n_ind_comps_complex_removed = length(ind_comps_complex_removed)
        for i=1:n_ind_comps_complex_removed
            GMT.plot!(f, cs[ind_comps_complex_removed[i],:],
                    pen="2p,255/165/0,dashed")
        end
        n_ind_comps_common = length(ind_comps_common)
        for i=1:n_ind_comps_common
            GMT.plot!(f, cs[ind_comps_common[i],:],
                    pen="2p,0/84/147,dashed")
        end

        GMT.basemap(region=(x_min_axb, x_max_axb, 0.9*y_min_axb, y_max_axb*1.1),
                proj=:logxy,
                axes=:WSne,
                xlabel="@~s@~@-m0@-",
                ylabel="Misfit spatial distribution (non dimensional)",
                xaxis=(annot=1, scale=:pow),
                yaxis=(annot=1, scale=:pow),
                title=titleb,
                panel=(1,2)
                )
        for i=1:n_comps_inverted
            GMT.plot!(sigmas, misfit_comps[:,i],lw=2,lc="0/84/147")
        end
        for i=1:n_comps_inverted
            ind_min = argmin(misfit_comps[:,i])[1]
            if i<n_comps_inverted    
                GMT.scatter!(sigmas[ind_min], misfit_comps[ind_min,i],
                    fill="255/165/0",
                    marker=:circle,
                    markersize=0.3,
                    markeredgecolor=:black)
            else
                GMT.scatter!(sigmas[ind_min], misfit_comps[ind_min,i],
                    fill="255/165/0",
                    marker=:circle,
                    markersize=0.3,
                    markeredgecolor=:black,
                    )
            end
        end
    subplot()
    
    if isfile("./gmt.history")
        rm("./gmt.history")
    end

    return nothing
end

function plot_psd_gmt(f, psds, options)

    n_comps = size(psds)[2]

    ind_comps = options["ind_comps"]
    ind_comps_complex = options["ind_comps_complex"]
    ind_comps_complex_removed = options["ind_comps_complex_removed"]
    ind_comps_common = options["ind_comps_common"]
    ind_comps_common_removed = options["ind_comps_common_removed"]

    color_list = options["color_list"]
    label = options["label"]
    title = options["title"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    n_subplots_x, n_subplots_y = options["subplot_grid"]
    n_subplots = n_subplots_x * n_subplots_y
    n_per_subplot = Integer(floor(n_comps / n_subplots));
    n_extra = Integer(mod(n_comps, n_subplots));

    f0 = f[1]
    f1 = f[end]
    y_min = 0
    i = 0
    k = 0
    x_legend = f[end]*0.9
    dx = -0.1*f[end]
    subplot(grid=(n_subplots_x,n_subplots_y), dims=(size=(15,15),
        frac=(ntuple(x -> 1, n_subplots_y), ntuple(x -> 1, n_subplots_x)) ),
            title=title, row_axes=(left=false), col_axes=(bott=false),
            margins=0, autolabel=true, savefig=dir_output*output*".png")
            
        for ix=1:n_subplots_x
            for iy=1:n_subplots_y
                i = i+1;
                if i<=n_extra
                    n_per_subplot_i = copy(n_per_subplot+1);
                else
                    n_per_subplot_i = copy(n_per_subplot);
                end
                
                y_max = maximum(psds[:,1+k:k+n_per_subplot_i])
                y_legend = y_max*0.95
                dy = -0.15*y_max
                if ix == n_subplots_x && iy == n_subplots_y
                    ann = Integer(floor((1.1*y_max - y_min) / 3))
                    GMT.basemap(region=(f0,f1,y_min,y_max*1.1),
                                xaxis=(annot=ann,ticks=ann),
                                yaxis=(annot=ann,ticks=ann),
                                par=(:MAP_FRAME_TYPE,"fancy+"),
                                axes=:wSnE,
                                xlabel="Frequency (1/yr)",
                                panel=(ix,iy)
                                )
                elseif ix < n_subplots_x && iy == n_subplots_y
                    ann = Integer(floor((1.1*y_max - y_min) / 3))
                    GMT.basemap(region=(f0,f1,y_min,y_max*1.1),
                                xaxis=(annot=ann,ticks=ann),
                                yaxis=(annot=ann,ticks=ann),
                                par=(:MAP_FRAME_TYPE,"fancy+"),
                                axes=:wsnE,
                                panel=(ix,iy)
                                )
                elseif ix < n_subplots_x && iy < n_subplots_y
                    ann = Integer(floor((1.1*y_max - y_min) / 3))
                    GMT.basemap(region=(f0,f1,y_min,y_max*1.1),
                                xaxis=(annot=ann,ticks=ann),
                                yaxis=(annot=ann,ticks=ann),
                                par=(:MAP_FRAME_TYPE,"fancy+"),
                                axes=:Wsne,
                                ylabel="PSD",
                                panel=(ix,iy)
                                )
                elseif ix == n_subplots_x && iy < n_subplots_y
                    ann = Integer(floor((1.1*y_max - y_min) / 3))
                    GMT.basemap(region=(f0,f1,y_min,y_max*1.1),
                                xaxis=(annot=ann,ticks=ann),
                                yaxis=(annot=ann,ticks=ann),
                                par=(:MAP_FRAME_TYPE,"fancy+"),
                                axes=:WSne,
                                xlabel="Frequency (1/yr)",
                                ylabel="PSD",
                                panel=(ix,iy)
                                )
                end
                for j=1:n_per_subplot_i
                    k = k + 1;
                    GMT.plot!(f,psds[:,k],lw=2,lc=color_list[j])
                    GMT.text!(label*string(k), x=x_legend, y=y_legend+(j-1)*dy,
                        par=((FONT_ANNOT_PRIMARY="12p,Helvetica,"*color_list[j]),
                        (FONT_LABEL="12p,Helvetica,"*color_list[j]),),
                        justify=:CM)
                    # if !isempty(intersect(k, ind_comps))
                    #     GMT.scatter!(x_legend+dx, y_legend+(j-1)*dy,
                    #         fill="193/43/43",
                    #         markersize=0.3,
                    #         marker=:star,
                    #         markeredgecolor=:black,
                    #         )
                    # end
                    if !isempty(intersect(k, ind_comps_complex))
                        GMT.scatter!(x_legend+dx, y_legend+(j-1)*dy,
                            fill="193/43/43",
                            markersize=0.3,
                            marker=:star,
                            markeredgecolor=:black,
                            )
                    end
                    if !isempty(intersect(k, ind_comps_complex_removed))
                        GMT.scatter!(x_legend+dx, y_legend+(j-1)*dy,
                            fill="255/165/0",
                            markersize=0.3,
                            marker=:star,
                            markeredgecolor=:black,
                            )
                    end
                    if !isempty(intersect(k, ind_comps_common))
                        GMT.scatter!(x_legend+1.5*dx, y_legend+(j-1)*dy,
                            fill="0/84/147",
                            markersize=0.3,
                            marker=:star,
                            markeredgecolor=:black,
                            )
                    end
                end
            end
        end
    subplot()

    if isfile("./gmt.history")
        rm("./gmt.history")
    end

    return nothing
end


function plot_comp(decomp, options)

    figsize = options["figsize"]
    limits = options["limits"]
    xy_legend = options["ll_legend"]
    proj = options["proj"]
    title = options["title"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    x = decomp["llh"][1:3:end,1]
    y = decomp["llh"][1:3:end,2]
    append!(x, xy_legend[1])
    append!(y, xy_legend[2])
    
    
    timeline = decomp["timeline"][1,:]
    U = decomp["U"]
    S = decomp["S"]
    V = decomp["V"]
    V_var = decomp["var_V"]
    min_V = minimum(V, dims=1)
    max_V = maximum(V, dims=1)

    norm_facts = max_V .- min_V
    temp_funcs = (V .- min_V) ./ norm_facts
    temp_funcs_err = sqrt.(decomp["var_V"]) ./ norm_facts
    
    spat_maps = (U * S) .* norm_facts

    n_comps = size(ICA["S"])[1]
    for i=1:n_comps
        time_func = temp_funcs[:,i]
        time_func_err = temp_funcs_err[:,i]
        
        spat_map = spat_maps[:,i]

        fig = Figure(size=figsize)
        gt = fig[0, 1] = GridLayout()
        gab = fig[1:2, 1] = GridLayout()
        ga = gab[1, 1] = GridLayout()
        gb = gab[2, 1] = GridLayout()
        
        #############
        ### TITLE ###
        #############
        titlelayout = GridLayout(gt[1, 1],
                                halign = :center,
                                tellwidth = false
                                )
        Label(titlelayout[1, 1],
            title*string(i),
            halign = :center,
            fontsize=36,
            font = "Times New Roman"
            )
        
        ###################
        ### TIME SERIES ###
        ###################
        sum_tr = sum(tremors["R"], dims=1)[1,:]
        axa = Axis(ga[1,1],
                xlabel=L"Time ($yr$)",
                ylabel=L"Normalized\ntime evolution",
                xlabelsize=24,
                ylabelsize=24,
                xticklabelsize=24,
                yticklabelsize=24,
                xaxisposition=:top,
                )
        # V_plot = Makie.lines!(
        #         axa,
        #         timeline, time_func,
        #         color=:black,
        #         )
        V_plot = Makie.errorbars!(
                axa,
                timeline, time_func,
                time_func_err,
                whiskerwidth = 3,
                color=:black,
                )
        ###########
        ### MAP ###
        ###########
        lon_min = limits[1]
        lon_max = limits[2]
        lat_min = limits[3]
        lat_max = limits[4]
        xticks = collect(lon_min-1:2:lon_max+1)
        xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
        yticks = collect(lat_min-1:2:lat_max+1)
        yticks = [round(Int, ytick) for ytick in yticks] # Ensure ticks are integers

        
        axb = GeoAxis(
            gb[1,1],
            dest=proj,
            limits=limits,
            xgridwidth=0.1,
            ygridwidth=0.1,
            xticklabelsize=24,
            yticklabelsize=24,
            xticks=xticks,
            yticks=yticks,
            xgridstyle=:solid,
            )
        stations_ver = Makie.scatter!(
            axb,
            decomp["llh"][:,1], decomp["llh"][:,2],
            marker=:circle,
            color=spat_map[3:3:end],
            markersize=20,
            alpha=1.0,
            colormap=:bwr,
            )

        u = spat_map[1:3:end]
        v = spat_map[2:3:end]
        append!(u, 10)
        append!(v, 0)

        max_vec_norm = maximum(sqrt.(u.^2 .+ v.^2))
        lengthscale = 1 / max_vec_norm
        stations_hor = Makie.arrows!(
                axb,
                x, y,
                u, v,
                lengthscale = lengthscale,
                linewidth = 1.2,
                #arrowsize = 1.2*max_vec_norm,
                )
        
        # Coastline
        # ne_10m_coastline = GeoMakie.coastlines(10)
        # save_object("ne_10m_coastline.jld2", ne_10m_coastline)
        coastline_file = dirs["dir_coastlines"] * "ne_10m_coastline.jld2"
        ne_10m_coastline = load_object(coastline_file)
        cl_obj = GeoMakie.lines!(axb, ne_10m_coastline, color=:black)
        translate!(cl_obj, 0, 0, 200)

        max_abs_value = maximum(abs.(spat_map[3:3:end]))
        tick_step = ceil(Int, max_abs_value / 4)
        ticks = collect(-4*tick_step:tick_step:4*tick_step)
        ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
        Colorbar(fig[3, 1],
                limits=(-maximum(abs.(spat_map[3:3:end])), maximum(abs.(spat_map[3:3:end]))),
                colormap=:bwr,
                flipaxis=false,
                label="Vertical (mm)",
                vertical=false,
                labelsize=24,
                ticks=ticks,
                ticklabelsize=24
                )
                
        #################################
        ### GENERAL LAYOUT ADJUSTMENT ###
        #################################
        rowgap!(gab, -50)
        rowsize!(gab, 1, Auto(0.3))

        save(output_dir*output*string(i)*".png", fig)
    end
    return nothing
end

# get fault vertices
function get_vertices(fault, proj="2d")

    vertices = Polygon[]
    n_features = size(fault["llh"])[1];

    if proj=="2d"
        for i = 1:n_features
            p1 = Point2f(fault["lon"][i,1], fault["lat"][i,1])
            p2 = Point2f(fault["lon"][i,2], fault["lat"][i,2])
            p3 = Point2f(fault["lon"][i,3], fault["lat"][i,3])
            append!(vertices,[Polygon(vec([p1,p2,p3]))])
        end
    else
        for i = 1:n_features
            p1 = Point2f(fault["lon"][i,1], fault["lat"][i,1], fault["height"][i,1])
            p2 = Point2f(fault["lon"][i,2], fault["lat"][i,2], fault["height"][i,2])
            p3 = Point2f(fault["lon"][i,3], fault["lat"][i,3], fault["height"][i,3])
            append!(vertices,[Polygon(vec([p1,p2,p3]))])
        end
    end
    return vertices
end


# get fault vertices for GeoMakie plot
function get_vertices_GeoMakie(fault)

    n_features = size(fault["llh"])[1];
    vertices = NaN .* zeros(n_features*4, 2)

    j = 0
    for i = 1:4:4*n_features
        j = j+1
        vertices[i+1,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
        vertices[i+2,:] = [fault["lon"][j, 2], fault["lat"][j, 2]]
        vertices[i+3,:] = [fault["lon"][j, 3], fault["lat"][j, 3]]
    end
    return vertices
end

function get_vertices_gmt(fault)
    n_features = size(fault["lon"])[1]
    vertices = NaN .* zeros(n_features*5, 2)

    j = 0
    for i = 1:5:5*n_features
        j = j+1
        vertices[i+1,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
        vertices[i+2,:] = [fault["lon"][j, 2], fault["lat"][j, 2]]
        vertices[i+3,:] = [fault["lon"][j, 3], fault["lat"][j, 3]]
        vertices[i+4,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
    end

    return vertices
end


function make_figure_map_lat_time(obs_var, fault, options)

    figsize = options["figsize"]
    title = options["title"]
    n_lats_edges = options["n_lats_edges"]
    color_thresh  = options["color_thresh_perc"]
    n_future_days = options["n_future_days"]
    output = options["output"]
    dir_output = options["dir_output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    timeline = obs_var["timeline"][1,:]
    obs = obs_var["obs"]

    tremors_timeline = tremors["timeline"]
    tremors_timeline_R = tremors["timeline_R"][1,:]

    n_samples = length(timeline)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    color_max = maximum(abs.(color))
    color_threshold = color_thresh*color_max
    color[abs.(color) .< color_threshold] .= 0
    colorrange = (-ceil(maximum(abs.(color))), ceil(maximum(abs.(color))))

    n_lats = n_lats_edges - 1;
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    lats_edges = range(lat_min, stop=lat_max, length=n_lats_edges);
    lats = lats_edges[1:end-1] + (lats_edges[2:end]-lats_edges[1:end-1])/2;
    color_lats_min = zeros(n_lats,n_samples);
    color_lats_max = zeros(n_lats,n_samples);

    for i=1:n_lats
        ind_lats = fault["llh"][:,2] .> lats_edges[i] .&&
                   fault["llh"][:,2] .<= lats_edges[i+1];
        color_lats_min[i,:] = minimum(color[ind_lats,:], dims=1);
        color_lats_max[i,:] = maximum(color[ind_lats,:], dims=1);
    end
    color_lats_minmax = color_lats_max;
    inds = color_lats_max .< abs.(color_lats_min);
    color_lats_minmax[inds] = color_lats_min[inds];

    fig = Figure(size=figsize)

    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(fig[0, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    xticks = collect(timeline[1]-1:2:timeline[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    ax = Axis(fig[1, 1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Latitude ($\degree$)",
               xaxisposition=:top,
               xlabelsize=24,
               ylabelsize=24,
               xgridwidth=0.1,
               ygridwidth=0.1,
               xticks=xticks,
               xticklabelsize=24,
               yticklabelsize=24
               )
    
    hm = heatmap!(ax, timeline, lats, color_lats_minmax',
            colormap=:bwr, colorrange=colorrange)
    translate!(hm, 0, 0, 100) # move above surface plot
    tr_map = Makie.scatter!(
        ax,
        tremors_timeline, tremors["lat"],
        color=:black,
        markersize=1,
        alpha=1.0,
        )
    translate!(tr_map, 0, 0, 200) # move above surface plot

    xlims!(ax, timeline[1],
                DateFormats.yeardecimal(today() + Day(n_future_days)))

    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[2, 1],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )
    save(dir_output*output*"_"*title*".png", fig)
    return fig
end

function make_figure_map_lat_time_ts(obs_var, fault, options)

    figsize = options["figsize"]
    title = options["title"]
    n_lats_edges = options["n_lats_edges"]
    color_thresh  = options["color_thresh_perc"]
    n_future_days = options["n_future_days"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end
    
    timeline = obs_var["timeline"][1,:]
    obs = copy(obs_var["obs"])

    tremors_timeline = tremors["timeline"]
    tremors_timeline_R = tremors["timeline_R"][1,:]

    #obs[fault["xyz"][:,3].<-30 .&& fault["xyz"][:,3].>0, :] .= 0
    n_samples = length(timeline)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    color_max = maximum(abs.(color))
    color_threshold = color_thresh*color_max
    color[abs.(color) .< color_threshold] .= 0
    colorrange = (-ceil(maximum(abs.(color))), ceil(maximum(abs.(color))))

    n_lats = n_lats_edges - 1;
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    lats_edges = range(lat_min, stop=lat_max, length=n_lats_edges);
    lats = lats_edges[1:end-1] + (lats_edges[2:end]-lats_edges[1:end-1])/2;
    color_lats_min = zeros(n_lats,n_samples);
    color_lats_max = zeros(n_lats,n_samples);

    for i=1:n_lats
        ind_lats = fault["llh"][:,2] .> lats_edges[i] .&&
                   fault["llh"][:,2] .<= lats_edges[i+1];
        color_lats_min[i,:] = minimum(color[ind_lats,:], dims=1);
        color_lats_max[i,:] = maximum(color[ind_lats,:], dims=1);
    end
    color_lats_minmax = color_lats_max;
    inds = color_lats_max .< abs.(color_lats_min);
    color_lats_minmax[inds] = color_lats_min[inds];

    fig = Figure(size=figsize)
    gt = fig[0, 1] = GridLayout()
    gab = fig[1:2, 1] = GridLayout()
    ga = gab[1, 1] = GridLayout()
    gb = gab[2, 1] = GridLayout()
    
    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(gt[1, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ###################
    ### TIME SERIES ###
    ###################
    sum_tr = sum(tremors["R"], dims=1)[1,:]
    max_obs = maximum(color_lats_minmax, dims=1)[1,:]
    min_obs = minimum(color_lats_minmax, dims=1)[1,:]
    xticks = collect(timeline[1]-1:2:timeline[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    axa = Axis(ga[1,1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xticks=xticks,
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    sum_tr_line = Makie.lines!(
        axa,
        tremors_timeline_R, sum_tr ./ maximum(sum_tr),
        color=:black,
        )
    max_obs_norm = max_obs ./ maximum(max_obs)
    max_obs_line = Makie.lines!(
        axa,
        timeline, max_obs_norm,
        color=:red,
        )

    # min_obs_norm = (abs.(min_obs) ./ maximum(abs.(min_obs)))
    # min_obs_line = Makie.lines!(
    #     axa,
    #     timeline, min_obs_norm,
    #     color=:blue,
    #     )
    #xlims!(axa, timeline[1], timeline[end])
    xlims!(axa, timeline[1],
                DateFormats.yeardecimal(today() + Day(n_future_days)))
    ylims!(axa, minimum([0, minimum(max_obs_norm)]), 1)
    # ylims!(axa, minimum([0, minimum(max_obs_norm), minimum(min_obs_norm)]), 1)
    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    xticks = collect(timeline[1]-1:2:timeline[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    axb = Axis(gb[1, 1],
               ylabel=L"Latitude ($\degree$)",
               xlabelsize=24,
               ylabelsize=24,
               xgridwidth=0.1,
               ygridwidth=0.1,
               xticks=xticks,
               xticklabelsize=24,
               yticklabelsize=24
               )
    hidexdecorations!(axb, grid=false)
    hm = heatmap!(axb, timeline, lats, color_lats_minmax',
            colormap=:bwr, colorrange=colorrange)
    translate!(hm, 0, 0, 100) # move above surface plot
    tr_map = Makie.scatter!(
        axb,
        tremors_timeline, tremors["lat"],
        color=:black,
        markersize=1,
        alpha=1.0,
        )
    translate!(tr_map, 0, 0, 200) # move above surface plot
    #xlims!(axb, timeline[1], timeline[end])
    xlims!(axb, timeline[1], 
                DateFormats.yeardecimal(today() + Day(n_future_days)))
    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[3, 1],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )

    #################################
    ### GENERAL LAYOUT ADJUSTMENT ###
    #################################
    rowgap!(gab, 0)
    rowsize!(gab, 1, Auto(0.3))
    
    save(dir_output*output*"_ts_"*title*".png", fig)
    return fig, r
end

function get_color(obs, color_thresh)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    color_max = maximum(abs.(color))
    color_threshold = color_thresh*color_max
    color[abs.(color) .< color_threshold] .= 0
    colorrange = (-ceil(color_max), ceil(color_max))
    return color, colorrange, color_exp
end

function get_color_lontime_map(color, fault, n_samples, n_lats_edges)
    n_lats = n_lats_edges - 1;
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    lats_edges = range(lat_min, stop=lat_max, length=n_lats_edges);
    lats = lats_edges[1:end-1] + (lats_edges[2:end]-lats_edges[1:end-1])/2;
    color_lats_min = zeros(n_lats,n_samples);
    color_lats_max = zeros(n_lats,n_samples);

    for i=1:n_lats
        ind_lats = fault["llh"][:,2] .> lats_edges[i] .&&
                   fault["llh"][:,2] .<= lats_edges[i+1];
        color_lats_min[i,:] = minimum(color[ind_lats,:], dims=1);
        color_lats_max[i,:] = maximum(color[ind_lats,:], dims=1);
    end
    color_lats_minmax = color_lats_max;
    inds = color_lats_max .< abs.(color_lats_min);
    color_lats_minmax[inds] = color_lats_min[inds];
    return color_lats_minmax, lats
end

function plot_map_lattimets(obs_var, tremors, fault, options)

    figsize = options["figsize"]
    title = options["title"]
    n_lats_edges  = options["n_lats_edges"]
    t0            = options["t0"] 
    t1            = options["t1"] 
    color_thresh  = options["color_thresh_perc"]
    Δt = options["Δt"]
    t_max_crosscorr = options["t_max_crosscorr"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    t0_decyear, t0_date = get_time_decyear_and_date(t0)
    t1_decyear, t1_date = get_time_decyear_and_date(t1)

    n_patches = size(fault["llh"])[1]

    ##############################
    ### FIND VARIABLES TO PLOT ###
    ##############################
    # Observable from GNSS solution (e.g., slip potency rate)
    obs_ind_t0 = findfirst(obs_var["timeline"][1,:] .- t0_decyear .>= 0)
    obs_ind_t1 = findlast(t1_decyear .- obs_var["timeline"][1,:] .>= 0)
    obs_timeline = obs_var["timeline"][1,obs_ind_t0:obs_ind_t1]
    obs_dates = Date.(DateFormats.yeardecimal.(obs_timeline))
    obs = obs_var["obs"][:,obs_ind_t0:obs_ind_t1]
    n_samples_obs = length(obs_timeline)
    # Find color and maxmin values along latitude for the lat-time map
    obs_on_fault, colorrange, color_exp = get_color(obs, color_thresh)
    obs_timelats, lats = get_color_lontime_map(obs_on_fault,
                            fault, n_samples_obs, n_lats_edges)
    max_obs_timelats = maximum(obs_timelats, dims=1)[1,:]
    n_lats = n_lats_edges - 1;
    # Tremors on original timeline (non-uniform)
    tremors_ind_t0 = findfirst(tremors["timeline"] .- t0_decyear .>= 0)
    tremors_ind_t1 = findlast(t1_decyear .- tremors["timeline"] .>= 0)
    tremors_timeline = tremors["timeline"][tremors_ind_t0:tremors_ind_t1]
    tremors_dates = Date.(DateFormats.yeardecimal.(tremors_timeline))
    lat_tremors = tremors["lat"][tremors_ind_t0:tremors_ind_t1]
    # Lon-lat of all tremors in considered time range
    tremors_ll = Point2f.(tremors["lon"][tremors_ind_t0:tremors_ind_t1],
                          tremors["lat"][tremors_ind_t0:tremors_ind_t1])
    # Tremors on daily timeline
    tremors_R_ind_t0 = findfirst(tremors["timeline_R"] .- t0_decyear .>= 0)
    tremors_R_ind_t1 = findlast(t1_decyear .- tremors["timeline_R"] .>= 0)
    tremors_timeline_R =tremors["timeline_R"][tremors_R_ind_t0:tremors_R_ind_t1]
    tremors_dates_R = Date.(DateFormats.yeardecimal.(tremors_timeline_R))
    sumR = sum(tremors["R"], dims=1)[1,tremors_R_ind_t0:tremors_R_ind_t1]
    
    ##################################################
    ### CREATE TIMELINE OF TIMES AND DATES TO PLOT ###
    ##################################################
    timeline2plot, dates2plot = create_timeline(t0_date, t1_date)
    n_samples2plot = length(dates2plot)

    # Initialize to NaN variables to plot
    obs_on_fault2plot = zeros(n_patches, n_samples2plot) .* NaN
    obs_timelats2plot = zeros(n_lats, n_samples2plot) .* NaN
    max_obs_timelats2plot = zeros(n_samples2plot) .* NaN
    sumR2plot = zeros(n_samples2plot) .* NaN
    lat_tremors2plot = zeros(n_samples2plot) .* NaN

    # Find indices of dates shared with the dates to plot
    ind_obs = findall(x -> x in dates2plot, obs_dates)
    ind_tremors = findall(x -> x in dates2plot, tremors_dates)
    ind_tremors_R = findall(x -> x in dates2plot, tremors_dates_R)

    # Fill with values the variables to plot
    obs_on_fault2plot[:,ind_obs] = obs_on_fault[:,:]
    obs_timelats2plot[:,ind_obs] = obs_timelats[:,:]
    max_obs_timelats2plot[ind_obs] = max_obs_timelats[:]
    sumR2plot[ind_tremors_R] = sumR[:]
    
    ind_notnan_max_obs_timelats2plot = findall(x -> !isnan(x),
                                            max_obs_timelats2plot)
    norm_max_obs_timelats2plot = max_obs_timelats2plot ./ 
                maximum(max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot])
    ind_notnan_sumR2plot = findall(x -> !isnan(x), sumR2plot)
    norm_sumR2plot = sumR2plot ./ maximum(sumR2plot[ind_notnan_sumR2plot])
    
    min_ylims_axa = minimum([0,
        minimum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        minimum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    max_ylims_axa = maximum([1,
        maximum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        maximum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    
    # Find ticks for x ans y axes
    xticks = collect(timeline2plot[1]-1:2:timeline2plot[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    yticks_axb = collect(lat_min-1:2:lat_max+1)
    yticks_axb = [round(Int, ytick) for ytick in yticks_axb] # Ensure ticks are integers

    
    ind_nan_norm_max_obs_timelats2plot = findall(!isnan, 
                                            norm_max_obs_timelats2plot)
    ind_nan_norm_sumR2plot = findall(!isnan, norm_sumR2plot)
    inds = intersect(ind_nan_norm_max_obs_timelats2plot,ind_nan_norm_sumR2plot)
    x = norm_max_obs_timelats2plot[inds]
    y = norm_sumR2plot[inds]

    ## Explanation of crosscor
    ## (see https://discourse.julialang.org/t/cross-correlation-with-statsbase-crosscor-answer-does-not-make-sense/44374)
    # Julia 1.5
    # using LinearAlgebra
    # numone, numtwo = [2,3,4,5], [1,2,3,4]
    # x, y = numone[1:end-1], numtwo[2:end]
    # ra = dot(x,y)/sqrt(dot(numone,numone)*dot(numtwo,numtwo)) # what crosscor does
    # rb = dot(x,y)/sqrt(dot(x,x)*dot(y,y)) # what you expect
    # [ra rb]

    inds_crosscorr = range(-t_max_crosscorr,t_max_crosscorr,step=1)
    r = crosscor(x,y,inds_crosscorr; demean=false)
    r_max = maximum(r)
    ind_r_max = argmax(r)
    n_days_max_crosscorr = inds_crosscorr[ind_r_max]

    ##############
    ### FIGURE ###
    ##############
    fig = Figure(size=figsize)
    gt = fig[0, 1] = GridLayout()
    gab = fig[1:2, 1] = GridLayout()
    ga = gab[1, 1] = GridLayout()
    gb = gab[2, 1] = GridLayout()
    
    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(gt[1, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ###################
    ### TIME SERIES ###
    ###################
    axa = Axis(ga[1,1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xlabelsize=24,
               ylabelsize=24,
               xticks=xticks,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    xlims!(axa, timeline2plot[1], timeline2plot[end] + Δt)
    ylims!(axa, min_ylims_axa, max_ylims_axa)

    # tremors_ts_scatter = Makie.scatter!(
    #     axa,
    #     timeline2plot, norm_sumR2plot,
    #     color=:black,
    #     markersize=5,
    #     alpha=1,
    #     )
    tremors_ts_lines = Makie.lines!(
        axa,
        timeline2plot, norm_sumR2plot,
        color=:black,
        )

    # obs_ts_lines = Makie.scatter!(
    #     axa,
    #     timeline2plot, norm_max_obs_timelats2plot,
    #     color=:red,
    #     markersize=5,
    #     alpha=1,
    #     )
    obs_ts_lines = Makie.lines!(
        axa,
        timeline2plot, norm_max_obs_timelats2plot,
        color=:red,
        )
    
    dx = (timeline2plot[end] + Δt - timeline2plot[1])*0.02
    dy = (max_ylims_axa - min_ylims_axa)*0.15
    Makie.text!(axa, timeline2plot[1]+dx, max_ylims_axa-dy;
            text="corr = "*string(round(r_max, sigdigits=2)), fontsize=16)
    if abs(n_days_max_crosscorr) == 1
        str_n_days_max_crosscorr = "Δt = "*string(n_days_max_crosscorr)*" day"
    else
        str_n_days_max_crosscorr = "Δt = "*string(n_days_max_crosscorr)*" days"
    end
    Makie.text!(axa, timeline2plot[1]+dx, max_ylims_axa-2*dy;
            text=str_n_days_max_crosscorr, fontsize=16)
    

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axb = Axis(gb[1, 1],
               ylabel=L"Latitude ($\degree$)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks_axb
               )
    hidexdecorations!(axb, grid = false)

    hm = heatmap!(axb,
                  timeline2plot,
                  lats,
                  obs_timelats2plot',
                  colormap=:bwr,
                  colorrange=colorrange
                  )
    
    tr_map = Makie.scatter!(
        axb,
        tremors_timeline,
        lat_tremors,
        color=:black,
        markersize=1,
        alpha=1.0,
        )
    translate!(tr_map, 0, 0, 100) # move above surface plot
    
    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[3, 1],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )
    
    #################################
    ### GENERAL LAYOUT ADJUSTMENT ###
    #################################
    rowgap!(gab, 0)
    rowsize!(gab, 1, Auto(0.3))
    
    save(dir_output*output*"_ts_"*title*".png", fig)
    
    return r
end


function plot_map_lattimets_notremors(obs_var, fault, options)

    figsize = options["figsize"]
    title = options["title"]
    n_lats_edges  = options["n_lats_edges"]
    t0            = options["t0"] 
    t1            = options["t1"] 
    color_thresh  = options["color_thresh_perc"]
    Δt = options["Δt"]
    t_max_crosscorr = options["t_max_crosscorr"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    t0_decyear, t0_date = get_time_decyear_and_date(t0)
    t1_decyear, t1_date = get_time_decyear_and_date(t1)

    n_patches = size(fault["llh"])[1]
    
    ##############################
    ### FIND VARIABLES TO PLOT ###
    ##############################
    # Observable from GNSS solution (e.g., slip potency rate)
    obs_ind_t0 = findfirst(obs_var["timeline"][1,:] .- t0_decyear .>= 0)
    obs_ind_t1 = findlast(t1_decyear .- obs_var["timeline"][1,:] .>= 0)
    obs_timeline = obs_var["timeline"][1,obs_ind_t0:obs_ind_t1]
    obs_dates = Date.(DateFormats.yeardecimal.(obs_timeline))
    obs = obs_var["obs"][:,obs_ind_t0:obs_ind_t1]
    n_samples_obs = length(obs_timeline)
    # Find color and maxmin values along latitude for the lat-time map
    obs_on_fault, colorrange, color_exp = get_color(obs, color_thresh)
    obs_timelats, lats = get_color_lontime_map(obs_on_fault,
                            fault, n_samples_obs, n_lats_edges)
    max_obs_timelats = maximum(obs_timelats, dims=1)[1,:]
    n_lats = n_lats_edges - 1;
    
    ##################################################
    ### CREATE TIMELINE OF TIMES AND DATES TO PLOT ###
    ##################################################
    timeline2plot, dates2plot = create_timeline(t0_date, t1_date)
    n_samples2plot = length(dates2plot)

    # Initialize to NaN variables to plot
    obs_on_fault2plot = zeros(n_patches, n_samples2plot) .* NaN
    obs_timelats2plot = zeros(n_lats, n_samples2plot) .* NaN
    max_obs_timelats2plot = zeros(n_samples2plot) .* NaN
    
    # Find indices of dates shared with the dates to plot
    ind_obs = findall(x -> x in dates2plot, obs_dates)
    
    # Fill with values the variables to plot
    obs_on_fault2plot[:,ind_obs] = obs_on_fault[:,:]
    obs_timelats2plot[:,ind_obs] = obs_timelats[:,:]
    max_obs_timelats2plot[ind_obs] = max_obs_timelats[:]
    
    ind_notnan_max_obs_timelats2plot = findall(x -> !isnan(x),
                                            max_obs_timelats2plot)
    norm_max_obs_timelats2plot = max_obs_timelats2plot ./ 
                maximum(max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot])
    
    min_ylims_axa = minimum([0,
        minimum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot])])
    max_ylims_axa = maximum([1,
        maximum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot])])
    
    # Find ticks for x ans y axes
    xticks = collect(timeline2plot[1]-1:2:timeline2plot[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    yticks_axb = collect(lat_min-1:2:lat_max+1)
    yticks_axb = [round(Int, ytick) for ytick in yticks_axb] # Ensure ticks are integers

    ##############
    ### FIGURE ###
    ##############
    fig = Figure(size=figsize)
    gt = fig[0, 1] = GridLayout()
    gab = fig[1:2, 1] = GridLayout()
    ga = gab[1, 1] = GridLayout()
    gb = gab[2, 1] = GridLayout()
    
    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(gt[1, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ###################
    ### TIME SERIES ###
    ###################
    axa = Axis(ga[1,1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xlabelsize=24,
               ylabelsize=24,
               xticks=xticks,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    xlims!(axa, timeline2plot[1], timeline2plot[end] + Δt)
    ylims!(axa, min_ylims_axa, max_ylims_axa)

    # obs_ts_scatter = Makie.scatter!(
    #     axa,
    #     timeline2plot, norm_max_obs_timelats2plot,
    #     color=:red,
    #     markersize=5,
    #     alpha=1,
    #     )
    obs_ts = Makie.lines!(
        axa,
        timeline2plot, norm_max_obs_timelats2plot,
        color=:red,
        )
    
    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axb = Axis(gb[1, 1],
               ylabel=L"Latitude ($\degree$)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks_axb
               )
    hidexdecorations!(axb, grid = false)

    hm = heatmap!(axb,
                  timeline2plot,
                  lats,
                  obs_timelats2plot',
                  colormap=:bwr,
                  colorrange=colorrange
                  )
    
    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[3, 1],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )
    
    #################################
    ### GENERAL LAYOUT ADJUSTMENT ###
    #################################
    rowgap!(gab, 0)
    rowsize!(gab, 1, Auto(0.3))
    
    save(dir_output*output*"_ts_notremors_"*title*".png", fig)
    
    return nothing
end


function plot_map_lattimets_northsouth(obs_var, tremors, fault, options)

    figsize         = options["figsize"]
    title           = options["title"]
    n_lats_edges    = options["n_lats_edges"]
    t0              = options["t0"] 
    t1              = options["t1"] 
    color_thresh    = options["color_thresh_perc"]
    Δt              = options["Δt"]
    t_max_crosscorr = options["t_max_crosscorr"]
    lat_split       = options["lat_split"] 
    dir_output      = options["dir_output"]
    output          = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    t0_decyear, t0_date = get_time_decyear_and_date(t0)
    t1_decyear, t1_date = get_time_decyear_and_date(t1)

    n_patches = size(fault["llh"])[1]

    ##############################
    ### FIND VARIABLES TO PLOT ###
    ##############################
    # Observable from GNSS solution (e.g., slip potency rate)
    obs_ind_t0 = findfirst(obs_var["timeline"][1,:] .- t0_decyear .>= 0)
    obs_ind_t1 = findlast(t1_decyear .- obs_var["timeline"][1,:] .>= 0)
    obs_timeline = obs_var["timeline"][1,obs_ind_t0:obs_ind_t1]
    obs_dates = Date.(DateFormats.yeardecimal.(obs_timeline))
    obs = obs_var["obs"][:,obs_ind_t0:obs_ind_t1]
    n_samples_obs = length(obs_timeline)
    # Find color and maxmin values along latitude for the lat-time map
    obs_on_fault, colorrange, color_exp = get_color(obs, color_thresh)
    obs_timelats, lats = get_color_lontime_map(obs_on_fault,
                            fault, n_samples_obs, n_lats_edges)
    # Split North and South observations
    ind_lats_north = lats .>= lat_split
    ind_lats_south = lats .< lat_split
    obs_timelats_north = obs_timelats[ind_lats_north,:]
    obs_timelats_south = obs_timelats[ind_lats_south,:]
    max_obs_timelats_north = maximum(obs_timelats_north, dims=1)[1,:]
    max_obs_timelats_south = maximum(obs_timelats_south, dims=1)[1,:]
    n_lats = n_lats_edges - 1;
    # Tremors on original timeline (non-uniform)
    tremors_ind_t0 = findfirst(tremors["timeline"] .- t0_decyear .>= 0)
    tremors_ind_t1 = findlast(t1_decyear .- tremors["timeline"] .>= 0)
    tremors_timeline = tremors["timeline"][tremors_ind_t0:tremors_ind_t1]
    tremors_dates = Date.(DateFormats.yeardecimal.(tremors_timeline))
    lat_tremors = tremors["lat"][tremors_ind_t0:tremors_ind_t1]
    # Lon-lat of all tremors in considered time range
    tremors_ll = Point2f.(tremors["lon"][tremors_ind_t0:tremors_ind_t1],
                          tremors["lat"][tremors_ind_t0:tremors_ind_t1])
    # Tremors on daily timeline
    tremors_R_ind_t0 = findfirst(tremors["timeline_R"] .- t0_decyear .>= 0)
    tremors_R_ind_t1 = findlast(t1_decyear .- tremors["timeline_R"] .>= 0)
    tremors_timeline_R =tremors["timeline_R"][tremors_R_ind_t0:tremors_R_ind_t1]
    tremors_dates_R = Date.(DateFormats.yeardecimal.(tremors_timeline_R))
    
    ind_lats_fault_north = fault["llh"][:,2] .>= lat_split
    ind_lats_fault_south = fault["llh"][:,2] .< lat_split
    tremorsR_north = tremors["R"][ind_lats_fault_north,:]
    tremorsR_south = tremors["R"][ind_lats_fault_south,:]
    sumR_north=sum(tremorsR_north, dims=1)[1,tremors_R_ind_t0:tremors_R_ind_t1]
    sumR_south=sum(tremorsR_south, dims=1)[1,tremors_R_ind_t0:tremors_R_ind_t1]
    
    ##################################################
    ### CREATE TIMELINE OF TIMES AND DATES TO PLOT ###
    ##################################################
    timeline2plot, dates2plot = create_timeline(t0_date, t1_date)
    n_samples2plot = length(dates2plot)

    # Initialize to NaN variables to plot
    obs_on_fault2plot = zeros(n_patches, n_samples2plot) .* NaN
    obs_timelats2plot = zeros(n_lats, n_samples2plot) .* NaN
    max_obs_timelats2plot_north = zeros(n_samples2plot) .* NaN
    max_obs_timelats2plot_south = zeros(n_samples2plot) .* NaN
    sumR2plot_north = zeros(n_samples2plot) .* NaN
    sumR2plot_south = zeros(n_samples2plot) .* NaN
    lat_tremors2plot = zeros(n_samples2plot) .* NaN

    # Find indices of dates shared with the dates to plot
    ind_obs = findall(x -> x in dates2plot, obs_dates)
    ind_tremors = findall(x -> x in dates2plot, tremors_dates)
    ind_tremors_R = findall(x -> x in dates2plot, tremors_dates_R)

    # Fill with values the variables to plot
    obs_on_fault2plot[:,ind_obs] = obs_on_fault[:,:]
    obs_timelats2plot[:,ind_obs] = obs_timelats[:,:]
    max_obs_timelats2plot_north[ind_obs] = max_obs_timelats_north[:]
    max_obs_timelats2plot_south[ind_obs] = max_obs_timelats_south[:]
    sumR2plot_north[ind_tremors_R] = sumR_north[:]
    sumR2plot_south[ind_tremors_R] = sumR_south[:]
    
    ind_notnan_max_obs_timelats2plot_north = findall(x -> !isnan(x),
                                            max_obs_timelats2plot_north)
    norm_max_obs_timelats2plot_north = max_obs_timelats2plot_north ./ 
        maximum(
            max_obs_timelats2plot_north[ind_notnan_max_obs_timelats2plot_north])
    
    ind_notnan_max_obs_timelats2plot_south = findall(x -> !isnan(x),
                                            max_obs_timelats2plot_south)
    norm_max_obs_timelats2plot_south = max_obs_timelats2plot_south ./ 
        maximum(
            max_obs_timelats2plot_south[ind_notnan_max_obs_timelats2plot_south])
    
    ind_notnan_sumR2plot_north = findall(x -> !isnan(x), sumR2plot_north)
    norm_sumR2plot_north = sumR2plot_north ./
                    maximum(sumR2plot_north[ind_notnan_sumR2plot_north])

    ind_notnan_sumR2plot_south = findall(x -> !isnan(x), sumR2plot_south)
    norm_sumR2plot_south = sumR2plot_south ./
                    maximum(sumR2plot_south[ind_notnan_sumR2plot_south])
    
    min_ylims_axa = minimum([0,
        minimum(norm_max_obs_timelats2plot_north[
            ind_notnan_max_obs_timelats2plot_north]),
        minimum(norm_sumR2plot_north[ind_notnan_sumR2plot_north])]
        )
    max_ylims_axa = maximum([1,
        maximum(norm_max_obs_timelats2plot_north[
            ind_notnan_max_obs_timelats2plot_north]),
        maximum(norm_sumR2plot_north[ind_notnan_sumR2plot_north])]
        )

    min_ylims_axc = minimum([0,
        minimum(norm_max_obs_timelats2plot_south[
            ind_notnan_max_obs_timelats2plot_south]),
        minimum(norm_sumR2plot_south[ind_notnan_sumR2plot_south])]
        )
    max_ylims_axc = maximum([1,
        maximum(norm_max_obs_timelats2plot_south[
            ind_notnan_max_obs_timelats2plot_south]),
        maximum(norm_sumR2plot_south[ind_notnan_sumR2plot_south])]
        )
    
    # Find ticks for x ans y axes
    xticks = collect(timeline2plot[1]-1:2:timeline2plot[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    yticks_axb = collect(lat_min-1:2:lat_max+1)
    yticks_axb = [round(Int, ytick) for ytick in yticks_axb] # Ensure ticks are integers

    
    ind_nan_norm_max_obs_timelats2plot_north = findall(!isnan, 
                                            norm_max_obs_timelats2plot_north)
    ind_nan_norm_max_obs_timelats2plot_south = findall(!isnan, 
                                            norm_max_obs_timelats2plot_south)
    ind_nan_norm_sumR2plot_north = findall(!isnan, norm_sumR2plot_north)
    ind_nan_norm_sumR2plot_south = findall(!isnan, norm_sumR2plot_south)
    inds_north = intersect(ind_nan_norm_max_obs_timelats2plot_north,
                            ind_nan_norm_sumR2plot_north)
    inds_south = intersect(ind_nan_norm_max_obs_timelats2plot_south,
                            ind_nan_norm_sumR2plot_south)
    x_north = norm_max_obs_timelats2plot_north[inds_north]
    x_south = norm_max_obs_timelats2plot_south[inds_south]
    y_north = norm_sumR2plot_north[inds_north]
    y_south = norm_sumR2plot_south[inds_south]

    ## Explanation of crosscor
    ## (see https://discourse.julialang.org/t/cross-correlation-with-statsbase-crosscor-answer-does-not-make-sense/44374)
    # Julia 1.5
    # using LinearAlgebra
    # numone, numtwo = [2,3,4,5], [1,2,3,4]
    # x, y = numone[1:end-1], numtwo[2:end]
    # ra = dot(x,y)/sqrt(dot(numone,numone)*dot(numtwo,numtwo)) # what crosscor does
    # rb = dot(x,y)/sqrt(dot(x,x)*dot(y,y)) # what you expect
    # [ra rb]

    inds_crosscorr = range(-t_max_crosscorr,t_max_crosscorr,step=1)
    r_north = crosscor(x_north,y_north,inds_crosscorr; demean=false)
    r_south = crosscor(x_south,y_south,inds_crosscorr; demean=false)
    r_north_max = maximum(r_north)
    r_south_max = maximum(r_south)
    ind_r_north_max = argmax(r_north)
    ind_r_south_max = argmax(r_south)
    n_days_max_crosscorr_north = inds_crosscorr[ind_r_north_max]
    n_days_max_crosscorr_south = inds_crosscorr[ind_r_south_max]

    ##############
    ### FIGURE ###
    ##############
    fig = Figure(size=figsize)
    gt = fig[0, 1] = GridLayout()
    gabc = fig[1:3, 1] = GridLayout()
    ga = gabc[1, 1] = GridLayout()
    gb = gabc[2, 1] = GridLayout()
    gc = gabc[3, 1] = GridLayout()
    
    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(gt[1, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ###################
    ### TIME SERIES ###
    ###################
    axa = Axis(ga[1,1],
               xlabel=L"Time (yr)$$",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xlabelsize=24,
               ylabelsize=24,
               xticks=xticks,
               yticks=[0,0.5,1],
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    xlims!(axa, timeline2plot[1], timeline2plot[end] + Δt)
    ylims!(axa, min_ylims_axa, max_ylims_axa)

    # tremors_ts_scatter = Makie.scatter!(
    #     axa,
    #     timeline2plot, norm_sumR2plot_north,
    #     color=:black,
    #     markersize=5,
    #     alpha=1,
    #     )
    tremors_ts_lines = Makie.lines!(
        axa,
        timeline2plot, norm_sumR2plot_north,
        color=:black,
        )

    # obs_ts_lines = Makie.scatter!(
    #     axa,
    #     timeline2plot, norm_max_obs_timelats2plot_north,
    #     color=:red,
    #     markersize=5,
    #     alpha=1,
    #     )
    obs_ts_lines = Makie.lines!(
        axa,
        timeline2plot, norm_max_obs_timelats2plot_north,
        color=:red,
        )
    
    dx = (timeline2plot[end] + Δt - timeline2plot[1])*0.02
    dy = (max_ylims_axa - min_ylims_axa)*0.17
    Makie.text!(axa, timeline2plot[1]+dx, max_ylims_axa-dy;
            text="corr = "*string(round(r_north_max, sigdigits=2)), fontsize=16)
    if abs(n_days_max_crosscorr_north) == 1
        str_n_days_max_crosscorr_north = "Δt = " * 
                        string(n_days_max_crosscorr_north)*" day"
    else
        str_n_days_max_crosscorr_north = "Δt = " * 
                            string(n_days_max_crosscorr_north)*" days"
    end
    Makie.text!(axa, timeline2plot[1]+dx, max_ylims_axa-2*dy;
            text=str_n_days_max_crosscorr_north, fontsize=16)
    

    axc = Axis(gc[1,1],
               xlabel=L"Time (yr)$$",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xlabelsize=24,
               ylabelsize=24,
               xticks=xticks,
               yticks=[0,0.5,1],
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:bottom,
               )
    xlims!(axc, timeline2plot[1], timeline2plot[end] + Δt)
    ylims!(axc, min_ylims_axc, max_ylims_axc)

    # tremors_ts_scatter = Makie.scatter!(
    #     axc,
    #     timeline2plot, norm_sumR2plot_south,
    #     color=:black,
    #     markersize=5,
    #     alpha=1,
    #     )
    tremors_ts_lines = Makie.lines!(
        axc,
        timeline2plot, norm_sumR2plot_south,
        color=:black,
        )

    # obs_ts_lines = Makie.scatter!(
    #     axc,
    #     timeline2plot, norm_max_obs_timelats2plot_south,
    #     color=:red,
    #     markersize=5,
    #     alpha=1,
    #     )
    obs_ts_lines = Makie.lines!(
        axc,
        timeline2plot, norm_max_obs_timelats2plot_south,
        color=:red,
        )
    
    dy = (max_ylims_axc - min_ylims_axc)*0.17
    Makie.text!(axc, timeline2plot[1]+dx, max_ylims_axc-dy;
            text="corr = "*string(round(r_south_max, sigdigits=2)), fontsize=16)
    if abs(n_days_max_crosscorr_south) == 1
        str_n_days_max_crosscorr_south = "Δt = " * 
                        string(n_days_max_crosscorr_south)*" day"
    else
        str_n_days_max_crosscorr_south = "Δt = " * 
                            string(n_days_max_crosscorr_south)*" days"
    end
    Makie.text!(axc, timeline2plot[1]+dx, max_ylims_axc-2*dy;
            text=str_n_days_max_crosscorr_south, fontsize=16)

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axb = Axis(gb[1, 1],
               ylabel=L"Latitude ($\degree$)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks_axb
               )
    hidexdecorations!(axb, grid = false)

    hm = heatmap!(axb,
                  timeline2plot,
                  lats,
                  obs_timelats2plot',
                  colormap=:bwr,
                  colorrange=colorrange
                  )
    
    tr_map = Makie.scatter!(
        axb,
        tremors_timeline,
        lat_tremors,
        color=:black,
        markersize=1,
        alpha=1.0,
        )
    translate!(tr_map, 0, 0, 100) # move above surface plot
    
    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[2, 2],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=true,
             label=L"%$(label_text) (m$^3$/yr) $\times 10^{%$(color_exp)}$",
             vertical=true,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )
    
    #################################
    ### GENERAL LAYOUT ADJUSTMENT ###
    #################################
    rowgap!(gabc, 0)
    rowsize!(gabc, 1, Auto(0.2))
    rowsize!(gabc, 3, Auto(0.2))
    
    save(dir_output*output*"_ts_"*title*".png", fig)
    
    return r_north, r_south
end