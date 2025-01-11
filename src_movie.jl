using GLMakie
using GeoMakie
using JLD2
using DateFormats
include("src_time.jl")


function make_video_map_2(obs_var, tremors, fault, options)

    limits        = options["limits"]
    figsize       = options["figsize"]
    proj          = options["proj"]
    title         = options["title"]
    show_progress = options["show_progress"]
    framerate     = options["framerate"]
    n_lats_edges  = options["n_lats_edges"]
    Δt            = options["Δt"] 
    t0            = options["t0"] 
    t1            = options["t1"] 
    color_thresh  = options["color_thresh_perc"]
    n_future_days = options["n_future_days"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    t0_decyear, t0_date = get_time_decyear_and_date(t0)
    t1_decyear, t1_date = get_time_decyear_and_date(t1)

    n_patches = size(fault["llh"])[1]
    vertices = get_vertices(fault)
    lon_min = minimum(minimum(fault["lon"]));
    lon_max = maximum(maximum(fault["lon"]));
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    xticks_axa = collect(lon_min-1:2:lon_max+1)
    xticks_axa = [round(Int, xtick) for xtick in xticks_axa] # Ensure ticks are integers
    yticks_axa = collect(lat_min-1:2:lat_max+1)
    yticks_axa = [round(Int, ytick) for ytick in yticks_axa] # Ensure ticks are integers
    yticks_axc = copy(yticks_axa)

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

    # println(size(dates2plot))
    # println(size(tremors_dates_R))
    # println(ind_tremors_R[end])
    # println(n_samples2plot)
    # println(size(sumR))
    # stop
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
    
    min_ylims_axb = minimum([0,
        minimum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        minimum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    max_ylims_axb = maximum([1,
        maximum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        maximum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    timestamps = range(1, n_samples2plot, step=1)

    # println(max_obs_timelats2plot)
    # println(maximum(max_obs_timelats2plot[~isnan(max_obs_timelats2plot)]))
    

    # dates = Date.(DateFormats.yeardecimal.(obs_timeline))  
    # dates_tremors = Date.(DateFormats.yeardecimal.(tremors_timeline))
    # dates_tremors_R = Date.(DateFormats.yeardecimal.(tremors_timeline_R))  

    # Start from first index of the timeline
    ind_t_now = Observable(1)
    date_now_str = @lift(string(dates2plot[$ind_t_now]));
    obs_on_fault_now = @lift(obs_on_fault2plot[:,$ind_t_now]);

    # println(size(timeline2plot))
    # println(size(norm_sumR2plot))
    timelinenormsumR = Point2f.(timeline2plot,norm_sumR2plot)
    timelinenormsumR_1now = @lift(timelinenormsumR[1:$ind_t_now])

    timelinenormmaxobstimelats = Point2f.(timeline2plot,
                                            norm_max_obs_timelats2plot)
    timelinenormmaxobstimelats_1now = @lift(
                                timelinenormmaxobstimelats[1:$ind_t_now])
    
    obs_timelats2plot_1now_tr = @lift(obs_timelats2plot[:,1:$ind_t_now]')
    timeline2plot_1now = @lift(timeline2plot[1:$ind_t_now])
    # obs_timelats2plot_1now = @lift(obs_timelats2plot[:,1:$ind_t_now])
    # max_obs_timelats2plot_1now = @lift(max_obs_timelats2plot[1:$ind_t_now])
    # sumR2plot_1now = @lift(sumR2plot[1:$ind_t_now])

    ind_tremors_now = Observable(tremors_dates .== dates2plot[1])
    tremors_ll_now = @lift(tremors_ll[$ind_tremors_now])

    ind_tremors_1now = Observable(tremors_dates .<= dates2plot[1])
    lat_tremors2plot = Point2f.(tremors_timeline, lat_tremors)
    lat_tremors2plot_1now = @lift(lat_tremors2plot[$ind_tremors_1now])

    #########################
    ### CREATE THE FIGURE ###
    #########################
    fig = Figure(size=figsize)
    gt = fig[0, 1:2]  = GridLayout()
    ga = fig[1:2, 1]  = GridLayout()
    gatop = ga[1, 1]  = GridLayout()
    gabot = ga[2, 1]  = GridLayout()
    gbc = fig[1:2, 2] = GridLayout()
    gb = gbc[1, 1]    = GridLayout()
    gc = gbc[2, 1]    = GridLayout()

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

    ############
    ### DATE ###
    ############
    axatop = Axis(gatop[1,1])
    hidedecorations!(axatop) # hides ticks, grid and lables
    hidespines!(axatop) # hide the frame
    Makie.text!(axatop, 0, 0;
            text=date_now_str,
            fontsize=36,
            font="Times New Roman",
            align=(:center,:center),
            )

    ######################
    ### SUBDUCTION MAP ###
    ######################
    
    axa = GeoAxis(
        gabot[1,1],
        dest=proj,
        limits=limits,
        xgridwidth=0.1,
        ygridwidth=0.1,
        xticklabelsize=24,
        yticklabelsize=24,
        xticks=xticks_axa,
        yticks=yticks_axa,
        xgridstyle=:solid,
        # title=date_now_str,
        # titlesize=36,
        # titlefont="Times New Roman",
        # titlegap=50,
        )
    hidexdecorations!(axa, grid = false)
    # Observable on fault
    trg = poly!(
        axa,
        vertices,
        color=obs_on_fault_now,
        strokecolor=:black,
        strokewidth=0.0,
        colormap=:bwr,
        colorrange=colorrange,
        )
    translate!(trg, 0, 0, 100) # move above surface plot
    # Tremors
    tr = GeoMakie.scatter!(
            axa,
            tremors_ll_now,
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
    axb = Axis(gb[1,1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    ylims!(axb, min_ylims_axb, max_ylims_axb)
    
    normsumR_scatter = Makie.scatter!(
        axb,
        timelinenormsumR_1now,
        color=:black,
        markersize=5,
        alpha=1,
        )
    normsumR_lines = Makie.lines!(
        axb,
        timelinenormsumR_1now,
        color=:black,
        )

    normmaxobstimelats_scatter = Makie.scatter!(
        axb,
        timelinenormmaxobstimelats_1now,
        color=:red,
        markersize=5,
        alpha=1,
        )
    
    normmaxobstimelats_lines = Makie.lines!(
        axb,
        timelinenormmaxobstimelats_1now,
        color=:red,
        )

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axc = Axis(gc[1, 1],
               ylabel=L"Latitude ($\degree$)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks_axc
               )
    hidexdecorations!(axc, grid = false)

    hm = heatmap!(axc,
                  timeline2plot_1now,
                  lats,
                  obs_timelats2plot_1now_tr,
                  colormap=:bwr,
                  colorrange=colorrange
                  )
    
    tr_map = Makie.scatter!(
        axc,
        lat_tremors2plot_1now,
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
    #rowsize!(fig.layout, 1, Auto(0.5))
    #colsize!(ga, 1, Auto(0.5))
    rowsize!(ga, 1, Auto(0.17))
    rowsize!(gbc, 1, Auto(0.3))
    rowgap!(ga, 0)
    rowgap!(gbc, 0)
    colgap!(fig.layout, 0)

    #############
    ### MOVIE ###
    #############
    progress = ProgressMeter.Progress(
        n_samples2plot; desc = "Movie time step: ", enabled = show_progress
    )
    record(fig, dir_output*output*".mp4", timestamps;
                        framerate = framerate) do time_idx
        ind_t_now[] = time_idx
        ind_tremors_now[] = tremors_dates .== dates2plot[time_idx]
        ind_tremors_1now[] = tremors_dates .<= dates2plot[time_idx]
        # ind_tr_obs[] = dates_tremors .<= dates[time_idx]
        
        x1 = DateFormats.yeardecimal(dates2plot[time_idx] + Day(n_future_days))
        xlims!(axb, timeline2plot[time_idx]-Δt, x1)
        xlims!(axc, timeline2plot[time_idx]-Δt, x1)
        ProgressMeter.next!(progress)
    end
    # save(dir_output*output*"_ts_"*title*".png", fig)
end








function make_video_map(obs_var, tremors, fault, options)

    limits        = options["limits"]
    figsize       = options["figsize"]
    proj          = options["proj"]
    title         = options["title"]
    show_progress = options["show_progress"]
    framerate     = options["framerate"]
    n_lats_edges  = options["n_lats_edges"]
    Δt            = options["Δt"] 
    t0            = options["t0"] 
    t0_time       = options["t0_time"]
    color_thresh  = options["color_thresh_perc"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

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
    record(fig, dir_output*output*".mp4", timestamps;
                        framerate = framerate) do time_idx
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
    show_progress = options["show_progress"]
    framerate     = options["framerate"]
    n_lats_edges  = options["n_lats_edges"]
    Δt            = options["Δt"] 
    t0            = options["t0"] 
    t0_time       = options["t0_time"]
    n_future_days = options["n_future_days"]
    n_tail        = options["n_tail"]
    dir_output = options["dir_output"]
    output = options["output"]

    if ~isdir(dir_output)
        mkdir(dir_output)
    end

    timeline = obs_var["timeline"][1,:]
    obs = obs_var["obs"]

    timeline_tremors = tremors["timeline"]
    timeline_tremors_R = tremors["timeline_R"][1,:]

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
    dates_tremors = Date.(DateFormats.yeardecimal.(tremors["timeline_R"]))

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
        timeline_tremors, tremors["lat"],
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
    record(fig, dir_output*output*".mp4", timestamps;
                    framerate = framerate) do time_idx
        ind_t_obs[] = time_idx
        ind_tr_obs[] = dates_tremors .== dates[time_idx]

        x1 = DateFormats.yeardecimal(dates[time_idx] + Day(n_future_days))
        xlims!(axb, t1.val-Δt, x1)
        xlims!(axc, t1.val-Δt, x1)
        # xlims!(axb, t1.val-Δt, t1.val)
        # xlims!(ax2, t1.val-Δt, t1.val)
        ProgressMeter.next!(progress)
    end

end