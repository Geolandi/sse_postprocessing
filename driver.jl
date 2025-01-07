using MAT
using DynamicalSystems
using CairoMakie
using GLMakie
#using DateFormats
using GeoMakie
using NaturalEarth
using JLD2

include("src_load.jl")
include("src_inversion.jl")
include("src_rate.jl")
include("src_localindices.jl")
include("src_plot.jl")
include("src_movie.jl")
include("src_tremors.jl")
# # Add ffmpeg directory to PATH
# ENV["PATH"] = ENV["PATH"] * ":/opt/homebrew/bin"


#solution_date = "2024-10-27/"
# dirs["dir_results2load"] = 
#     "../../Matlab/vbICA_Cascadia/Scenarios/Cascadia/SSE/ES/results/" * 
#     solution_date

solution_date = "2025-01-03/"
dirs = Dict()
dirs["dir_data"] = "/Users/ag2347/Work/Data/"
dirs["dir_case"] = "Cascadia/"
dirs["dir_results2load"] = "./matfiles/" * dirs["dir_case"] * solution_date
dirs["dir_fault"]   = dirs["dir_data"] * "Faults/" * dirs["dir_case"]
dirs["dir_tremors"] = dirs["dir_data"] * "Tremors/" * dirs["dir_case"] * "PNSN/"

# # load indices of components to invert
# ind_comps    = matread(dirs["dir_results2load"]*"ind_comps.mat")["ind_comps"]
# ind_comps = round.(Int, ind_comps)[:,1]
# # load misfit of inverted components
# misfit_comps = matread(
#     dirs["dir_results2load"]*"misfit_comps.mat")["misfit_comps"]
# # load options
# options      = matread(dirs["dir_results2load"]*"options.mat")["options"]
# # load position time series used as input for vbICA
# X            = matread(dirs["dir_results2load"]*"X.mat")["Xd"]
# # load vbICA results
# ICA          = matread(dirs["dir_results2load"]*"ICA.mat")["ICA_essential"]
# ICA["type"]     = X["type"];
# ICA["llh"]      = X["llh"];
# ICA["timeline"] = X["timeline"];
# ICA["decmode"]  = X["decmode"];

# # load fault
# fault = load_fault(dirs, options)

# # calculate Greens' functions
# G = create_greens_function(X, fault, options)

# # find smoothing parameters for inversion
# min_misfit_comps, ind_sigma0_comps_Cartesian = findmin(misfit_comps, dims=1)
# ind_sigma0_comps = getindex.(ind_sigma0_comps_Cartesian, 1)

# # invert components
# m, Cm = invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, options)
# println("Done")

# # calculate co-variance matrix for slip potency components
# n_ICs2invert = length(ind_comps)
# n_patches = length(fault["area"])
# mp = repeat(fault["area"],2,1) .* m;
# Cmp = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
# for i=1:n_ICs2invert
#     Cmp[i] = repeat(fault["area"],2,1) .* Cm[i];
# end

# options["slip_rate"] = Dict()
# options["slip_rate"]["rake"] = 90;
# options["slip_rate"]["windowsize"] = 7;


# options["smooth"] = Dict()
# # options["smooth"]["name"] = "loess"
# #options["smooth"]["name"] = "rollmean"
# options["smooth"]["name"] = "lowpass"

# if options["smooth"]["name"] == "loess"
#     options["smooth"]["span"] = 0.1
# elseif options["smooth"]["name"] == "rollmean"
#     options["smooth"]["windowsize"] = 3
#     options["smooth"]["centered_output"] = false
# elseif options["smooth"]["name"] == "lowpass"
#     options["smooth"]["filter_cutofffreq"] = 1/28
#     options["smooth"]["fs"] = 1
#     options["smooth"]["filter_type"] = "hanning"
#     options["smooth"]["filter_window"] = 365
#     options["smooth"]["scale"] = true
#     options["smooth"]["causal"]  = false
# end

# V = ICA["V"]
# V_smooth, timeline_smooth = create_xsmooth(
#     V', ICA["timeline"], options["smooth"])

# V_dot, timeline_dot = calc_derivative(V_smooth', timeline_smooth,
#     options["slip_rate"]["windowsize"], true)


# # update and load tremors
# update_tremors(dirs, timeline_dot[end]);
# tremors = load_tremors(dirs, fault["origin"]);
# tremors = select_tremors(tremors, timeline_dot, fault)

# options["inversion"]["rake_pos"] = 90;

# ICA_smooth = ICA;
# ICA_smooth["V"] = V_smooth;
# n_samples_smooth = length(timeline_smooth)
# ICA_smooth["var_V"] = ICA["var_V"][end-n_samples_smooth+1:end,:];
# ICA_smooth["timeline"] = timeline_smooth;
# slip_smooth = create_model(m, Cm, ICA_smooth, fault, options, ind_comps)

# ICA_dot = ICA;
# ICA_dot["V"] = V_dot;
# n_samples_dot = length(timeline_dot)
# ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
# ICA_dot["timeline"] = timeline_dot;
# slip_rate = create_model(m, Cm, ICA_dot, fault, options, ind_comps)

# #########################
# ### SLIP POTENCY RATE ###
# #########################
# # calculate slip potency rate
# slip_potency_rate = create_model(mp, Cmp, ICA_dot, fault, options, ind_comps)



# First, make a surface plot
options["plot"]["video"] = nothing

options["plot"]["video"]["map_video"]                  = Dict()
options["plot"]["video"]["map_video"]["limits"]        = (-128.2, # lon_min
                                                          -121.0, # lon_max
                                                          39.0,   # lat_min
                                                          51.0)   # lat_max
options["plot"]["video"]["map_video"]["figsize"]       = (1000, 1000)
options["plot"]["video"]["map_video"]["proj"]          = "+proj=merc"
options["plot"]["video"]["map_video"]["title"]         = "Cascadia"
options["plot"]["video"]["map_video"]["output"]        =  "./movies/" * 
                                                        dirs["dir_case"] * 
                                                        "slip_potency_rate_map"
options["plot"]["video"]["map_video"]["show_progress"] = true
options["plot"]["video"]["map_video"]["framerate"]     = 20
options["plot"]["video"]["map_video"]["n_lats_edges"]  = 101
options["plot"]["video"]["map_video"]["Δt"]            = 1.0
options["plot"]["video"]["map_video"]["t0"]            = 2024.5
options["plot"]["video"]["map_video"]["t0_time"]       = true

# fig = make_figure_map_lat_time(slip_potency_rate, fault, (1510, 900), title)

make_video_map(slip_potency_rate,
               tremors,
               fault,
               options["plot"]["video"]["map_video"]
               )

X = StateSpaceSet(slip_potency_rate["obs"]');
p = 0.99
est = :exp
d, θ = extremevaltheory_dims_persistences(X, Exceedances(p, est))
make_video_dtheta(options["plot"]["movie"])

# fig = Figure(size=figsize)
#     ax = GeoAxis(
#         fig[1,1],
#         dest=proj,
#         limits=limits,
#         xgridwidth=0,
#         ygridwidth=0,
#         title=title,
#         titlefont="Times New Roman",
#         titlesize=36,
#         xlabel = "Longitude (°)",
#         )
# GeoMakie.text!(-125, 38, text="test")
# display(fig)


# # Plot the triangle
# trg = poly!(
#     ax,
#     vertices,
#     color=color[:,20],
#     strokecolor=:black,
#     strokewidth=0.0,
#     colormap=:bwr,
#     colorrange=colorrange,
#     )
# translate!(trg, 0, 0, 100) # move above surface plot


# tr = GeoMakie.scatter!(
#     ax,
#     tremors["lon"][ind_tr],
#     tremors["lat"][ind_tr],
#     color=:black,
#     markersize=5,
#     )
# translate!(tr, 0, 0, 250) # move above surface plot

# display(fig)
# stop




# tr = GeoMakie.scatter!(
#     ax,
#     tremors_obs[],
#     color=:black,
#     markersize=5,
#     )


# translate!(tr, 0, 0, 250) # move above surface plot

# tr = GeoMakie.scatter!(
#     ax,
#     tremors["lon"],
#     tremors["lat"],
#     color=:black,
#     markersize=5,
#     )
# translate!(tr, 0, 0, 250) # move above surface plot



# progress = ProgressMeter.Progress(
#         n_samples; desc = "Movie time step: ", enabled = show_progress
#     )
# if isnothing(output) == true
#     for i = 1:n_samples
#         #poly!(v, color = ṗᵢ_def_sse_filt[i,:], colormap=colormap, colorrange=(-max_ṗᵢ_sse_filt, max_ṗᵢ_sse_filt))
#         time_obs[] = timeline_dot[i];
#         color_obs[] = color[:,i];
#         push!(traj[], Point2f(xa[i], ya[i]))
#             traj[] = traj[] # <- important! Updating in-place the value of an
#                             # `Observable` does not trigger an update!
            
#         sleep(1/framerate);
#         ProgressMeter.next!(progress)
#     end
# else
#     record(fig, "movie_test.mp4", dates; framerate = framerate) do t
#         time[] = t
#     end

#     record(f, string(output,".mp4"), range(1,n_samples);
#             framerate = framerate) do i
#         time[] = t[i];
#         if isnothing(xb) == true | isnothing(yb) == true
#             color[] = colorb[i,:];
#         else
#             color[] = colorb[i,:,:]';
#         end
#         pointa[] = Point2f(xa[i], ya[i]);
#         push!(traj[], Point2f(xa[i], ya[i]))
#         traj[] = traj[] # <- important! Updating in-place the value of an
#                         # `Observable` does not trigger an update!
#         ProgressMeter.next!(progress)
#     end
# end

# region = [-128.2, -121, 39, 51]
# # Create the map
# fig = coast(region=region,
#             proj=:Mercator,
#             shorelines="1/0.5p,black",
#             resolution=:i,
#             xaxis=(annot=2,ticks=1),
#             yaxis=(annot=2,ticks=1),
#             scale=1,
#             par=(:MAP_FRAME_TYPE,"fancy+"),
#             title="Cascadia"
#             )

# # Define the file path
# vertices_file = "./vertices.txt"

# # Open the file in write mode and write some content
# open(vertices_file, "w") do file
#     for i=1:n_patches
#         write(file, "> -Z"*string.(slip_potency_rate["obs"][i,1])*"\n")
#         write(file, string(fault["lon"][i, 1])*" "*string(fault["lat"][i, 1])*"\n")
#         write(file, string(fault["lon"][i, 2])*" "*string(fault["lat"][i, 2])*"\n")
#         write(file, string(fault["lon"][i, 3])*" "*string(fault["lat"][i, 3])*"\n")
#         write(file, string(fault["lon"][i, 1])*" "*string(fault["lat"][i, 1])*"\n")
#     end
# end



# function create_frame(color, t, date)
#     fig = basemap(region=region,
#             proj=:Mercator,
#             xaxis=(annot=2,ticks=1),
#             yaxis=(annot=2,ticks=1),
#             scale=1,
#             par=(:MAP_FRAME_TYPE,"fancy+"),
#             title="Cascadia"
#             )
#     GMT.plot!(vertices_file,
#             region=region,
#             level=(data=color[:,t]),
#             pen="0p,black",
#             C=cpt
#             )

#     GMT.text!(string(date), x=region[1] + 0.5, y=region[3] + 0.5, font=18, justify=:LB)

#     GMT.coast!(region=region,
#             proj=:Mercator,
#             shorelines="1/0.5p,black",
#             resolution=:i,
#             scale=1,
#             par=(:MAP_FRAME_TYPE,"fancy+")
#             )
#     # Add a colorbar next to the plot
#     GMT.colorbar!(
#         pos=(anchor=:BC,length=(8,0.4), horizontal=true, offset=(0,1.0)),
#         frame=(annot=3, ticks=1, xlabel="Slip Potency Rate (m@+3@+/yr) x 10@+"*string(color_exp)*"@+",),
#         par=(FONT_ANNOT_PRIMARY=18,),
#         savefig="./figures/frame_$(lpad(t, 4, '0')).png"
#         )
#     return fig
# end

# # dates_dot = decyear2date(timeline_dot)
# dates_dot = Date.(DateFormats.yeardecimal.(timeline_dot))[1,:]

# dates_tremors = Date.(DateFormats.yeardecimal.(tremors["timeline"]))

# n_samples = length(timeline_dot)
# tremors["timeline_daily"] = timeline_dot
# tremors["N"] = zeros(n_samples)
# for i=1:n_samples
#     dates_tremors .== dates_dot[i]
#     if any(dates_tremors .== dates_dot[i])
#         tremors["N"][i] = sum(dates_tremors .== dates_dot[i])
#     end
# end

# # Define the file path
# file4movie = "./file4movie.txt"
# # Open the file in write mode and write some content
# open(file4movie, "w") do file
#     for i=1:10
#         write(file, string.(dates_dot[i]))
#         write(file, "\n")
#         # write(file, " ")
#         # for j=1:n_patches
#         #     write(file, string.(color[j,i]))
#         #     if j<n_patches
#         #         write(file, " ")
#         #     else
#         #         write(file, "\n")
#         #     end
#         # end
#     end
# end

# # 
# cmap_range_file = "./gmt_files/cmap_range.txt"
# open(cmap_range_file, "w") do file
#     write(file, string.(colorrange[1])*" "*
#                 string.(colorrange[2])*" "*
#                 string.(colorstep)*"\n"
#                 )
# end

# # 
# cmap_range_file = "./gmt_files/region.txt"
# open(cmap_range_file, "w") do file
#     write(file, string.(region[1])*" "*
#                 string.(region[2])*" "*
#                 string.(region[3])*" "*
#                 string.(region[4])*"\n"
#                 )
# end

# function create_frame()
#     fig = GMT.basemap(region=region,
#             proj=:Mercator,
#             xaxis=(annot=2,ticks=1),
#             yaxis=(annot=2,ticks=1),
#             scale=1,
#             par=(:MAP_FRAME_TYPE,"fancy+"),
#             title="Cascadia"
#             )
    
#     cpt = makecpt(cmap=:polar,
#             range=(colorrange[1], colorrange[2], 1), continuous=true
#             )
    
#     # #gmt makecpt -Cpolar -T-$dataverlim/$dataverlim/$dataverstep -Z > color.cpt
    
#     GMT.plot!(vertices_file,
#             region=region,
#             level=(data=color[:,1]),
#             pen="0p,black",
#             C=cpt
#             )
    
#     GMT.coast!(region=region,
#             proj=:Mercator,
#             shorelines="1/0.5p,black",
#             resolution=:i,
#             scale=1,
#             par=(:MAP_FRAME_TYPE,"fancy+")
#             )
#     # Add a colorbar next to the plot
#     # GMT.colorbar!(
#     #     pos=(anchor=:BC,length=(8,0.4), horizontal=true, offset=(0,1.0)),
#     #     frame=(annot=3, ticks=1, xlabel="Slip Potency Rate (m@+3@+/yr) x 10@+"*string(color_exp)*"@+",),
#     #     par=(FONT_ANNOT_PRIMARY=18,),
#     #     savefig="./figures/frame_$(lpad(t, 4, '0')).png"
#     #     )

#     #GMT.text!(x=region[1] + 0.5, y=region[3] + 0.5, font=18, justify=:LB, text="MOVIE_FRAME")
#     GMT.text!(offset=(shift=(0.5,0.5)), region=region, proj=:merc, figscale=1, font=18,
#         region_justify=:LB, text="MOVIE_COL0", axes=:noannot)
#     #return fig
# end

# #run(`./main_script.sh`)

# GMT.movie(create_frame,
#           canvas="10cx20cx30",
#           name=:cascadia_movie,
#           frames=file4movie,
#           format=:mp4,
#           frame_rate=4,
#           clean=true
#           )



######################
# dates_dot = decyear2date(timeline_dot)
# n_samples = length(color[1,:])
# for t=n_samples-3:n_samples
#     create_frame(color, t, dates_dot[t])
# end

# #ffmpeg -n -framerate 10 -pattern_type glob -i 'frame_*.png' test.mp4

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



# options["plot"]["map_lat_time"]["n_lons_edges"] = 41;
# n_lons_edges   = options["plot"]["map_lat_time"]["n_lons_edges"];
# n_lons = n_lons_edges - 1;
# lon_min = minimum(minimum(fault["lon"]));
# lon_max = maximum(maximum(fault["lon"]));
# lons_edges = linspace(lon_min, lon_max, n_lons_edges);
# lons = lons_edges[1:end-1] + (lons_edges[2:end]-lons_edges[1:end-1])/2;
# color_lons_min = zeros(n_lons,n_samples);
# color_lons_max = zeros(n_lons,n_samples);

# for i=1:n_lons
#     ind_lons = fault["llh"][:,1] .> lons_edges[i] .&& fault["llh"][:,1] .<= lons_edges[i+1];
#     color_lons_min[i,:] = minimum(color[ind_lons,:], dims=1);
#     color_lons_max[i,:] = maximum(color[ind_lons,:], dims=1);
# end

# color_lons_minmax = color_lons_max;
# inds = color_lons_max .< abs.(color_lons_min);
# color_lons_minmax[inds] = color_lons_min[inds];

# fig = Figure(size=(600, 1200))
# ax = Axis(fig[1, 1], xlabel="Longitude (°)", ylabel="Time (yr)")
# heatmap!(ax, lons, timeline_dot[1,:], color_lons_minmax,
#         colormap=:bwr, colorrange=colorrange)
# label_text = "Slip Potency Rate"
# Colorbar(fig[2, 1], limits = colorrange, colormap=:bwr, flipaxis=false,
#     label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$", vertical=false,)
# display(fig)