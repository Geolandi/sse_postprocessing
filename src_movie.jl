using GLMakie
using GeoMakie
using JLD2

function make_video_map(obs_var, tremors, fault, options)

    limits        = options["limits"]
    figsize       = options["figsize"]
    proj          = options["proj"]
    title         = options["title"]
    output        = options["output"]
    show_progress = options["show_progress"]
    framerate     = options["framerate"]
    n_lats_edges  = options["n_lats_edges"]
    Δt            = options["Δt"] 
    t0            = options["t0"] 
    t0_time       = options["t0_time"]
    color_thresh  = options["color_thresh_perc"]

    timeline = obs_var["timeline"][1,:]
    obs = obs_var["obs"]

    if t0_time == false
        ind_t0 = t0
    else
        ind_t0 = findfirst(timeline .- t0 .> 0)
    end
    n_samples = length(timeline)
    timestamps = range(ind_t0, n_samples, step=1)


    vertices = get_vertices(fault)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    color_max = maximum(abs.(color))
    color_threshold = color_thresh*color_max
    color[abs.(color) .< color_threshold] .= 0
    #color[fault["llh"][:,3] .< -50, :] .= 0
    colorrange = (-ceil(color_max), ceil(color_max))

    dates = Date.(DateFormats.yeardecimal.(timeline))  
    dates_tremors = Date.(DateFormats.yeardecimal.(tremors["timeline"]))

    ind_t_obs = Observable(ind_t0);
    date_obs_str = @lift(string(dates[$ind_t_obs]));
    color_obs = @lift(color[:,$ind_t_obs]);
    t1 = @lift(timeline[$ind_t_obs])
    ind_tr_obs = Observable(dates_tremors .== dates[1])
    points_all = Point2f.(tremors["lon"], tremors["lat"])
    points = @lift(points_all[$ind_tr_obs])
    
    

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
    gt = fig[0, 1:2] = GridLayout()
    ga = fig[1:2, 1] = GridLayout()
    gbc = fig[1:2, 2] = GridLayout()
    gb = gbc[1, 1] = GridLayout()
    gc = gbc[2, 1] = GridLayout()

    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(gt[1,1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1:2],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ######################
    ### SUBDUCTION MAP ###
    ######################
    lon_min = minimum(minimum(fault["lon"]));
    lon_max = maximum(maximum(fault["lon"]));
    xticks = collect(lon_min-1:2:lon_max+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    yticks = collect(lat_min-1:2:lat_max+1)
    yticks = [round(Int, ytick) for ytick in yticks] # Ensure ticks are integers
    axa = GeoAxis(
        ga[1,1],
        dest=proj,
        limits=limits,
        xgridwidth=0.1,
        ygridwidth=0.1,
        xticklabelsize=24,
        yticklabelsize=24,
        xticks=xticks,
        yticks=yticks,
        xgridstyle=:solid,
        title=date_obs_str,
        titlesize=36,
        titlefont="Times New Roman",
        titlegap=50,
        )
    # # Date
    # Makie.text!(ax1, lon_min+0.2, lat_min+0.2; text=date_obs_str, fontsize=24)
    # Observable on fault
    trg = poly!(
        axa,
        vertices,
        color=color_obs,
        strokecolor=:black,
        strokewidth=0.0,
        colormap=:bwr,
        colorrange=colorrange,
        )
    translate!(trg, 0, 0, 100) # move above surface plot
    # Tremors
    tr = GeoMakie.scatter!(
            axa,
            points,
            color=:black,
            markersize=5,
            )
    translate!(tr, 0, 0, 250) # move above surface plot
    # Coastline
    # ne_10m_coastline = GeoMakie.coastlines(10)
    # save_object("ne_10m_coastline.jld2", ne_10m_coastline)
    coastline_file = dirs["dir_coastlines"] * "ne_10m_coastline.jld2"
    ne_10m_coastline = load_object(coastline_file)
    cl_obj = GeoMakie.lines!(axa, ne_10m_coastline, color=:black)
    translate!(cl_obj, 0, 0, 200)

    ###################
    ### TIME SERIES ###
    ###################
    sum_tr = sum(tremors["R"], dims=1)[1,:]
    axb = Axis(gb[1,1],
               xlabel=L"Time ($yr$)",
               ylabel="tbd",
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    sum_tr_sc = Makie.scatter!(
        axb,
        timeline, sum_tr,
        color=:black,
        markersize=5,
        alpha=1,
        )
    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axc = Axis(gc[1, 1],
               ylabel="Latitude (°)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks
               )
    hidexdecorations!(axc, grid = false)

    hm = heatmap!(axc, timeline, lats, color_lats_minmax',
            colormap=:bwr, colorrange=colorrange)
    tr_map = Makie.scatter!(
        axc,
        tremors["timeline"], tremors["lat"],
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
    Colorbar(fig[3, 1:2],
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
    rowsize!(fig.layout, 1, Auto(0.5))
    colsize!(ga, 1, Auto(0.5))
    rowsize!(gbc, 1, Auto(0.5))
    rowgap!(gbc, 0)

    #############
    ### MOVIE ###
    #############
    progress = ProgressMeter.Progress(
        n_samples-ind_t0; desc = "Movie time step: ", enabled = show_progress
    )
    record(fig, output*".mp4", timestamps; framerate = framerate) do time_idx
        ind_t_obs[] = time_idx
        ind_tr_obs[] = dates_tremors .== dates[time_idx]
        xlims!(axb, t1.val-Δt, t1.val)
        xlims!(axc, t1.val-Δt, t1.val)
        ProgressMeter.next!(progress)
    end

end


function make_video_map_dθ(obs_var, tremors, fault, d, θ, options)

    limits        = options["limits"]
    figsize       = options["figsize"]
    proj          = options["proj"]
    title         = options["title"]
    output        = options["output"]
    show_progress = options["show_progress"]
    framerate     = options["framerate"]
    n_lats_edges  = options["n_lats_edges"]
    Δt            = options["Δt"] 
    t0            = options["t0"] 
    t0_time       = options["t0_time"]
    n_tail        = options["n_tail"]

    timeline = obs_var["timeline"][1,:]
    obs = obs_var["obs"]

    if t0_time == false
        ind_t0 = t0
    else
        ind_t0 = findfirst(timeline .- t0 .> 0)
    end
    n_samples = length(timeline)
    timestamps = range(ind_t0, n_samples, step=1)


    vertices = get_vertices(fault)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    #color[fault["llh"][:,3] .< -50, :] .= 0
    colorrange = (-ceil(maximum(abs.(color))), ceil(maximum(abs.(color))))

    dates = Date.(DateFormats.yeardecimal.(timeline))  
    dates_tremors = Date.(DateFormats.yeardecimal.(tremors["timeline"]))

    ind_t_obs = Observable(ind_t0);
    date_obs_str = @lift(string(dates[$ind_t_obs]));
    color_obs = @lift(color[:,$ind_t_obs]);
    t1 = @lift(timeline[$ind_t_obs])
    ind_tr_obs = Observable(dates_tremors .== dates[1])
    points_all = Point2f.(tremors["lon"], tremors["lat"])
    points = @lift(points_all[$ind_tr_obs])
    
    

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
    gt = fig[0, 1:2] = GridLayout()
    gab = fig[1, 1:2] = GridLayout()
    ga = gab[1, 1] = GridLayout()
    gb = gab[1, 2] = GridLayout()
    gcd = fig[2, 1:2] = GridLayout()
    gc = gcd[1, 1] = GridLayout()
    gd = gcd[1, 2] = GridLayout()

    #############
    ### TITLE ###
    #############
    if title=="date"
        title = date_obs_str
    end
    titlelayout = GridLayout(gt[1,1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1:2],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    ################
    ### d-θ plot ###
    ################
    color_dθ_all = maximum(obs, dims=1)[1,:] ./ 10^color_exp
    points_dθ_all = Point2f.(d, θ)
    point_dθ = @lift(points_dθ_all[$ind_t_obs])
    
    ax = Axis(ga[1,1],
               xlabel=L"$d_1$",
               ylabel=L"$\theta$",
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               )
    dθ_sc_all = Makie.scatter!(
        ax,
        d, θ,
        color=color_dθ_all,
        markersize=5,
        alpha=0.3,
        )
    dθ_sc = Makie.scatter!(
        ax,
        point_dθ,
        color=:red,
        markersize=10,
        alpha=1.0,
        )
    translate!(dθ_sc, 0, 0, 100) # move above surface plot

    label_text = "Max Slip Potency Rate"
    max_abs_value = maximum(color_dθ_all)
    tick_step = ceil(Int, max_abs_value / 4)
    ticks = collect(0:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(ga[0, 1],
             limits=(0, max_abs_value),
             colormap=:viridis,
             flipaxis=true,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=18,
             ticks=ticks,
             ticklabelsize=18
             )
    
    ###################
    ### TIME SERIES ###
    ###################
    axb = Axis(gb[1,1],
               xlabel=L"Time ($yr$)",
               ylabel="tbd",
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               )
    sc = Makie.scatter!(
        axb,
        timeline, color_dθ_all  ./ maximum(color_dθ_all),
        color=:black,
        markersize=5,
        alpha=1,
        )
    sc_d = Makie.scatter!(
        axb,
        timeline, d ./ maximum(d),
        color=:red,
        markersize=5,
        alpha=1,
        )
    sc_θ = Makie.scatter!(
        axb,
        timeline, θ ./ maximum(θ),
        color=:blue,
        markersize=5,
        alpha=1,
        )

    ######################
    ### SUBDUCTION MAP ###
    ######################
    lon_min = minimum(minimum(fault["lon"]));
    lon_max = maximum(maximum(fault["lon"]));
    xticks = collect(lon_min-1:2:lon_max+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    yticks = collect(lat_min-1:2:lat_max+1)
    yticks = [round(Int, ytick) for ytick in yticks] # Ensure ticks are integers
    ax1 = GeoAxis(
        gc[1, 1],
        dest=proj,
        limits=limits,
        xgridwidth=0.1,
        ygridwidth=0.1,
        xticklabelsize=24,
        yticklabelsize=24,
        xticks=xticks,
        yticks=yticks,
        # title=date_obs_str,
        # titlesize=36,
        # titlefont="Times New Roman",
        # titlegap=50,
        )
    # # Date
    # Makie.text!(ax1, lon_min+0.2, lat_min+0.2; text=date_obs_str, fontsize=24)
    # Observable on fault
    trg = poly!(
        ax1,
        vertices,
        color=color_obs,
        strokecolor=:black,
        strokewidth=0.0,
        colormap=:bwr,
        colorrange=colorrange,
        )
    translate!(trg, 0, 0, 100) # move above surface plot
    # Tremors
    tr = GeoMakie.scatter!(
            ax1,
            points,
            color=:black,
            markersize=5,
            )
    translate!(tr, 0, 0, 250) # move above surface plot
    # Coastline
    # ne_10m_coastline = GeoMakie.coastlines(10)
    # save_object("ne_10m_coastline.jld2", ne_10m_coastline)
    ne_10m_coastline = load_object("./coastlines/ne_10m_coastline.jld2")
    cl_obj = GeoMakie.lines!(ax1, ne_10m_coastline, color=:black)
    translate!(cl_obj, 0, 0, 200)

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    ax2 = Axis(gd[1, 1],
               xlabel="Time (yr)",
               ylabel="Latitude (°)",
               xaxisposition=:top,
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               yticks=yticks
               )
    
    hm = heatmap!(ax2, timeline, lats, color_lats_minmax',
            colormap=:bwr, colorrange=colorrange)
    tr_map = Makie.scatter!(
        ax2,
        tremors["timeline"], tremors["lat"],
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
    Colorbar(fig[3, 1:2],
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
    rowsize!(fig.layout, 1, Auto(0.5))
    colsize!(gab, 1, Auto(0.5))
    colsize!(gcd, 1, Auto(0.5))
    rowgap!(fig.layout, 1)

    #############
    ### MOVIE ###
    #############
    progress = ProgressMeter.Progress(
        n_samples-ind_t0; desc = "Movie time step: ", enabled = show_progress
    )
    record(fig, output*".mp4", timestamps; framerate = framerate) do time_idx
        ind_t_obs[] = time_idx
        ind_tr_obs[] = dates_tremors .== dates[time_idx]
        xlims!(axb, t1.val-Δt, t1.val)
        xlims!(ax2, t1.val-Δt, t1.val)
        ProgressMeter.next!(progress)
    end

end