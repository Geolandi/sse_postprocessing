using GLMakie
using GeoMakie
using JLD2
using DateFormats
include("src_time.jl")


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
    t1            = options["t1"] 
    color_thresh  = options["color_thresh_perc"]
    n_future_days = options["n_future_days"]
    dir_coastlines = options["dir_coastlines"]
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
    obs_timelats, lats = get_color_lattime_map(obs_on_fault,
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
    ind_nan_sumR2plot = findall(x -> isnan(x), sumR2plot)
    norm_sumR2plot = sumR2plot ./ maximum(sumR2plot[ind_notnan_sumR2plot])
    norm_sumR2plot[ind_nan_sumR2plot] .= 0.0

    inds_tremors_nan = findall((timeline2plot .- tremors["timeline"][1] .< 0) .|
                               (timeline2plot .- tremors["timeline"][end] .> 0))
    norm_sumR2plot[inds_tremors_nan] .= NaN
    
    min_ylims_axb = minimum([0,
        minimum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        minimum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    max_ylims_axb = maximum([1,
        maximum(norm_max_obs_timelats2plot[ind_notnan_max_obs_timelats2plot]),
        maximum(norm_sumR2plot[ind_notnan_sumR2plot])]
        )
    timestamps = range(1, n_samples2plot, step=1)

    # Start from first index of the timeline
    ind_t_now = Observable(1)
    date_now_str = @lift(string(dates2plot[$ind_t_now]));
    obs_on_fault_now = @lift(obs_on_fault2plot[:,$ind_t_now]);

    timelinenormsumR = Point2f.(timeline2plot,norm_sumR2plot)
    timelinenormsumR_1now = @lift(timelinenormsumR[1:$ind_t_now])

    timelinenormmaxobstimelats = Point2f.(timeline2plot,
                                            norm_max_obs_timelats2plot)
    timelinenormmaxobstimelats_1now = @lift(
                                timelinenormmaxobstimelats[1:$ind_t_now])
    
    obs_timelats2plot_1now_tr = @lift(obs_timelats2plot[:,1:$ind_t_now]')
    timeline2plot_1now = @lift(timeline2plot[1:$ind_t_now])
    
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
    coastline_file = dir_coastlines * "ne_10m_coastline.jld2"
    ne_10m_coastline = load_object(coastline_file)
    cl_obj = GeoMakie.lines!(axa, ne_10m_coastline, color=:black)
    translate!(cl_obj, 0, 0, 200)

    ###################
    ### TIME SERIES ###
    ###################
    xticks = collect(round(timeline2plot[1])-1:0.5:round(timeline2plot[end]+1))
    #xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    axb = Axis(gb[1,1],
               xlabel=L"Time ($yr$)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xticks=xticks,
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


function make_video_map_northsouth(obs_var, tremors, fault, options)

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
    lat_split     = options["lat_split"] 
    dir_coastlines = options["dir_coastlines"]
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
    obs_timelats, lats = get_color_lattime_map(obs_on_fault,
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
    ind_nan_sumR2plot_north = findall(x -> isnan(x), sumR2plot_north)
    norm_sumR2plot_north = sumR2plot_north ./
                    maximum(sumR2plot_north[ind_notnan_sumR2plot_north])
    norm_sumR2plot_north[ind_nan_sumR2plot_north] .= 0.0

    ind_notnan_sumR2plot_south = findall(x -> !isnan(x), sumR2plot_south)
    ind_nan_sumR2plot_south = findall(x -> isnan(x), sumR2plot_south)
    norm_sumR2plot_south = sumR2plot_south ./
                    maximum(sumR2plot_south[ind_notnan_sumR2plot_south])
    norm_sumR2plot_south[ind_nan_sumR2plot_south] .= 0.0

    inds_tremors_nan = findall((timeline2plot .- tremors["timeline"][1] .< 0) .|
                               (timeline2plot .- tremors["timeline"][end] .> 0))
    norm_sumR2plot_north[inds_tremors_nan] .= NaN
    norm_sumR2plot_south[inds_tremors_nan] .= NaN

    min_ylims_axb = minimum([0,
        minimum(norm_max_obs_timelats2plot_north[
            ind_notnan_max_obs_timelats2plot_north]),
        minimum(norm_sumR2plot_north[ind_notnan_sumR2plot_north])]
        )
    max_ylims_axb = maximum([1,
        maximum(norm_max_obs_timelats2plot_north[
            ind_notnan_max_obs_timelats2plot_north]),
        maximum(norm_sumR2plot_north[ind_notnan_sumR2plot_north])]
        )

    min_ylims_axd = minimum([0,
        minimum(norm_max_obs_timelats2plot_south[
            ind_notnan_max_obs_timelats2plot_south]),
        minimum(norm_sumR2plot_south[ind_notnan_sumR2plot_south])]
        )
    max_ylims_axd = maximum([1,
        maximum(norm_max_obs_timelats2plot_south[
            ind_notnan_max_obs_timelats2plot_south]),
        maximum(norm_sumR2plot_south[ind_notnan_sumR2plot_south])]
        )

    timestamps = range(1, n_samples2plot, step=1)

    # Start from first index of the timeline
    ind_t_now = Observable(1)
    date_now_str = @lift(string(dates2plot[$ind_t_now]));
    obs_on_fault_now = @lift(obs_on_fault2plot[:,$ind_t_now]);

    timelinenormsumR_north = Point2f.(timeline2plot,norm_sumR2plot_north)
    timelinenormsumR_north_1now = @lift(timelinenormsumR_north[1:$ind_t_now])

    timelinenormmaxobstimelats_north = Point2f.(timeline2plot,
                                            norm_max_obs_timelats2plot_north)
    timelinenormmaxobstimelats_north_1now = @lift(
                                timelinenormmaxobstimelats_north[1:$ind_t_now])
    
    obs_timelats2plot_1now_tr = @lift(obs_timelats2plot[:,1:$ind_t_now]')
    timeline2plot_1now = @lift(timeline2plot[1:$ind_t_now])
    
    ind_tremors_now = Observable(tremors_dates .== dates2plot[1])
    tremors_ll_now = @lift(tremors_ll[$ind_tremors_now])

    ind_tremors_1now = Observable(tremors_dates .<= dates2plot[1])
    lat_tremors2plot = Point2f.(tremors_timeline, lat_tremors)
    lat_tremors2plot_1now = @lift(lat_tremors2plot[$ind_tremors_1now])

    timelinenormsumR_south = Point2f.(timeline2plot,norm_sumR2plot_south)
    timelinenormsumR_south_1now = @lift(timelinenormsumR_south[1:$ind_t_now])

    timelinenormmaxobstimelats_south = Point2f.(timeline2plot,
                                            norm_max_obs_timelats2plot_south)
    timelinenormmaxobstimelats_south_1now = @lift(
                                timelinenormmaxobstimelats_south[1:$ind_t_now])

    #########################
    ### CREATE THE FIGURE ###
    #########################
    fig = Figure(size=figsize)
    gt = fig[0, 1:2]  = GridLayout()
    ga = fig[1:3, 1]  = GridLayout()
    gatop = ga[1, 1]  = GridLayout()
    gamid = ga[2, 1]  = GridLayout()
    gabot = ga[3, 1]  = GridLayout()
    gbcd = fig[1:3, 2:3] = GridLayout()
    gb = gbcd[1, 1]    = GridLayout()
    gc = gbcd[2, 1]    = GridLayout()
    gd = gbcd[3, 1]    = GridLayout()
    gcb = gbcd[2, 2]   = GridLayout()

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
        gamid[1,1],
        dest=proj,
        limits=limits,
        xgridwidth=0.1,
        ygridwidth=0.1,
        xticklabelsize=24,
        yticklabelsize=24,
        xticks=xticks_axa,
        yticks=yticks_axa,
        xgridstyle=:solid,
        ygridstyle=:solid,
        )
    # hidexdecorations!(axa, grid = false)
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
    coastline_file = dir_coastlines * "ne_10m_coastline.jld2"
    ne_10m_coastline = load_object(coastline_file)
    cl_obj = GeoMakie.lines!(axa, ne_10m_coastline, color=:black)
    translate!(cl_obj, 0, 0, 200)

    ###################
    ### TIME SERIES ###
    ###################
    xticks = collect(round(timeline2plot[1])-1:0.5:round(timeline2plot[end]+1))
    axb = Axis(gb[1,1],
               xlabel="Time (yr)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xticks=xticks,
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:top,
               )
    ylims!(axb, min_ylims_axb, max_ylims_axb)
    
    normsumR_scatter = Makie.scatter!(
        axb,
        timelinenormsumR_north_1now,
        color=:black,
        markersize=5,
        alpha=1,
        )
    normsumR_lines = Makie.lines!(
        axb,
        timelinenormsumR_north_1now,
        color=:black,
        )

    normmaxobstimelats_scatter = Makie.scatter!(
        axb,
        timelinenormmaxobstimelats_north_1now,
        color=:red,
        markersize=5,
        alpha=1,
        )
    
    normmaxobstimelats_lines = Makie.lines!(
        axb,
        timelinenormmaxobstimelats_north_1now,
        color=:red,
        )

    
    axd = Axis(gd[1,1],
               xlabel="Time (yr)",
               ylabel=L"Norm. $R$ and $\dot{p}$",
               xticks=xticks,
               xlabelsize=24,
               ylabelsize=24,
               xticklabelsize=24,
               yticklabelsize=24,
               xaxisposition=:bottom,
               )
    ylims!(axd, min_ylims_axd, max_ylims_axd)
    
    normsumR_scatter = Makie.scatter!(
        axd,
        timelinenormsumR_south_1now,
        color=:black,
        markersize=5,
        alpha=1,
        )
    normsumR_lines = Makie.lines!(
        axd,
        timelinenormsumR_south_1now,
        color=:black,
        )

    normmaxobstimelats_scatter = Makie.scatter!(
        axd,
        timelinenormmaxobstimelats_south_1now,
        color=:red,
        markersize=5,
        alpha=1,
        )
    
    normmaxobstimelats_lines = Makie.lines!(
        axd,
        timelinenormmaxobstimelats_south_1now,
        color=:red,
        )

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    axc = Axis(gc[1, 1],
               #ylabel=L"Latitude ($\degree$)",
               ylabelsize=24,
               yticklabelsize=24,
               yticks=yticks_axc,
               xticks=xticks
               )
    axc.ytickformat = yticks_axc -> string.(Int.(yticks_axc)) .* "°"
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
    Colorbar(gcb[1,1],
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
    rowsize!(ga, 1, Auto(0.04))
    rowsize!(ga, 3, Auto(0.08))
    rowsize!(gbcd, 1, Auto(0.2))
    rowsize!(gbcd, 3, Auto(0.2))
    rowgap!(ga, 0)
    rowgap!(gbcd, 0)
    colgap!(gbcd, 10)

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
        xlims!(axd, timeline2plot[time_idx]-Δt, x1)
        ProgressMeter.next!(progress)
    end
    # save(dir_output*output*"_ts_"*title*".png", fig)
end