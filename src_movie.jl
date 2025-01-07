using GLMakie
using GeoMakie
using JLD2

function make_video_map(obs_var, tremors, fault, options)

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

    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(fig[0, 1:2],
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
    ax1 = GeoAxis(
        fig[1,1],
        dest=proj,
        limits=limits,
        xgridwidth=0.1,
        ygridwidth=0.1,
        xticklabelsize=24,
        yticklabelsize=24,
        xticks=xticks,
        yticks=yticks,
        title=date_obs_str,
        titlesize=24,
        titlefont="Times New Roman",
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
    ax2 = Axis(fig[1, 2],
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
    Colorbar(fig[2, 1:2],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )

    #############
    ### MOVIE ###
    #############
    progress = ProgressMeter.Progress(
        n_samples-ind_t0; desc = "Movie time step: ", enabled = show_progress
    )
    record(fig, output*".mp4", timestamps; framerate = framerate) do time_idx
        ind_t_obs[] = time_idx
        ind_tr_obs[] = dates_tremors .== dates[time_idx]
        xlims!(ax2, t1.val-Δt, t1.val)
        ProgressMeter.next!(progress)
    end

end


# options["plot"]["map_lat_time"] = Dict()

# options["plot"]["map_lat_time"]["n_lats_edges"] = 101;
# options["plot"]["map_lat_time"]["colormap"] = :bwr
# options["plot"]["map_lat_time"]["FontSize"] = 16;
# options["plot"]["map_lat_time"]["fig_pos"] = (10,10);
# options["plot"]["map_lat_time"]["fig_size"] = (1800, 800);
# options["plot"]["map_lat_time"]["output_name"] = "Fig_slipratemap";

# n_lats_edges   = options["plot"]["map_lat_time"]["n_lats_edges"];
# colormap_color = options["plot"]["map_lat_time"]["colormap"];
# FontSize       = options["plot"]["map_lat_time"]["FontSize"];
# fig_pos        = options["plot"]["map_lat_time"]["fig_pos"];
# fig_size       = options["plot"]["map_lat_time"]["fig_size"];
# output_name    = options["plot"]["map_lat_time"]["output_name"];

# n_lats = n_lats_edges - 1;
# lat_min = minimum(minimum(fault["lat"]));
# lat_max = maximum(maximum(fault["lat"]));
# lats_edges = linspace(lat_min, lat_max, n_lats_edges);
# lats = lats_edges[1:end-1] + (lats_edges[2:end]-lats_edges[1:end-1])/2;
# color_lats_min = zeros(n_lats,n_samples);
# color_lats_max = zeros(n_lats,n_samples);

# for i=1:n_lats
#     ind_lats = fault["llh"][:,2] .> lats_edges[i] .&& fault["llh"][:,2] .<= lats_edges[i+1];
#     color_lats_min[i,:] = minimum(color[ind_lats,:], dims=1);
#     color_lats_max[i,:] = maximum(color[ind_lats,:], dims=1);
# end

# color_lats_minmax = color_lats_max;
# inds = color_lats_max .< abs.(color_lats_min);
# color_lats_minmax[inds] = color_lats_min[inds];


# fig = Figure(size=(1200,800))
# ax = Axis(fig[1, 1], xlabel="Time (yr)", ylabel="Latitude (°)")
# heatmap!(ax, timeline_dot[1,:], lats, color_lats_minmax',
#         colormap=:bwr, colorrange=colorrange)
# label_text = "Slip Potency Rate"
# Colorbar(fig[2, 1], limits = colorrange, colormap=:bwr, flipaxis=false,
#     label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$", vertical=false,)
# display(fig)