include("src_load.jl")
include("src_inversion.jl")
include("src_rate.jl")
include("src_localindices.jl")
include("src_plot.jl")
include("src_movie.jl")
include("src_tremors.jl")

# Add ffmpeg directory to PATH to be able to plot with GMT
ENV["PATH"] = ENV["PATH"] * ":/opt/homebrew/bin"

#solution_date = "2024-08-18/"
solution_date = "2025-01-06/"
dirs = Dict()
dirs["dir_data"] = "/Users/ag2347/Work/Data/"
dirs["dir_results"] = "/Users/ag2347/Work/Results/"
dirs["dir_case"] = "Cascadia/"
dirs["dir_results2load"] = dirs["dir_results"] * "Slowquakes/real-time/" *
                           dirs["dir_case"] * "matfiles/" * solution_date
dirs["dir_fault"]   = dirs["dir_data"] * "Faults/" * dirs["dir_case"]
dirs["dir_tremors"] = dirs["dir_data"] * "Tremors/" * dirs["dir_case"] * "PNSN/"
dirs["dir_coastlines"] = dirs["dir_data"] * "NaturalEarth/coastlines_jl/"

# load indices of components to invert
ind_comps    = matread(dirs["dir_results2load"]*"ind_comps.mat")["ind_comps"]
ind_comps = round.(Int, ind_comps)[:,1]
# load misfit of inverted components
misfit_comps = matread(
    dirs["dir_results2load"]*"misfit_comps.mat")["misfit_comps"]
# load options
options      = matread(dirs["dir_results2load"]*"options.mat")["options"]
# load position time series used as input for vbICA
X            = matread(dirs["dir_results2load"]*"X.mat")["Xd"]
# load vbICA results
ICA          = matread(dirs["dir_results2load"]*"ICA.mat")["ICA_essential"]
ICA["type"]     = X["type"];
ICA["llh"]      = X["llh"];
ICA["timeline"] = X["timeline"];
ICA["decmode"]  = X["decmode"];

# load fault
fault = load_fault(dirs, options)

# calculate Greens' functions
G = create_greens_function(X, fault, options)

# find smoothing parameters for inversion
min_misfit_comps, ind_sigma0_comps_Cartesian = findmin(misfit_comps, dims=1)
ind_sigma0_comps = getindex.(ind_sigma0_comps_Cartesian, 1)

# invert components
m, Cm = invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, options)
println("Done")

# calculate co-variance matrix for slip potency components
n_ICs2invert = length(ind_comps)
n_patches = length(fault["area"])
mp = repeat(fault["area"],2,1) .* m;
Cmp = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
for i=1:n_ICs2invert
    Cmp[i] = repeat(fault["area"],2,1) .* Cm[i];
end

options["slip_rate"] = Dict()
options["slip_rate"]["rake"] = 90;
options["slip_rate"]["windowsize"] = 7;


options["smooth"] = Dict()
# options["smooth"]["name"] = "loess"
# options["smooth"]["name"] = "rollmean"
options["smooth"]["name"] = "lowpass"

if options["smooth"]["name"] == "loess"
    options["smooth"]["span"] = 0.1
elseif options["smooth"]["name"] == "rollmean"
    options["smooth"]["windowsize"] = 3
    options["smooth"]["centered_output"] = false
elseif options["smooth"]["name"] == "lowpass"
    options["smooth"]["filter_cutofffreq"] = 1/28
    options["smooth"]["fs"] = 1
    options["smooth"]["filter_type"] = "hanning"
    options["smooth"]["filter_window"] = 365
    options["smooth"]["scale"] = true
    options["smooth"]["causal"]  = false
end

V = ICA["V"]
V_smooth, timeline_smooth = create_xsmooth(
    V', ICA["timeline"], options["smooth"])

V_dot, timeline_dot = calc_derivative(V_smooth', timeline_smooth,
    options["slip_rate"]["windowsize"], true)


# update and load tremors
update_tremors(dirs, ICA["timeline"][end]);
tremors = load_tremors(dirs, fault["origin"]);
tremors = select_tremors(tremors, timeline_dot, fault)

options["inversion"]["rake_pos"] = 90;

ICA_smooth = copy(ICA);
ICA_smooth["V"] = V_smooth;
n_samples_smooth = length(timeline_smooth)
ICA_smooth["var_V"] = ICA["var_V"][end-n_samples_smooth+1:end,:];
ICA_smooth["timeline"] = timeline_smooth;
slip_smooth = create_model(m, Cm, ICA_smooth, fault, options, ind_comps)

ICA_dot = copy(ICA);
ICA_dot["V"] = V_dot;
n_samples_dot = length(timeline_dot)
ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
ICA_dot["timeline"] = timeline_dot;
slip_rate = create_model(m, Cm, ICA_dot, fault, options, ind_comps)

#########################
### SLIP POTENCY RATE ###
#########################
# calculate slip potency rate
slip_potency_rate = create_model(mp, Cmp, ICA_dot, fault, options, ind_comps)


##############
### MOVIES ###
##############
# First, make a surface plot
options["plot"]["video"] = Dict()

options["plot"]["video"]["map_video"]                  = Dict()
options["plot"]["video"]["map_video"]["limits"]        = (-128.2, # lon_min
                                                          -121.0, # lon_max
                                                          39.0,   # lat_min
                                                          51.0)   # lat_max
options["plot"]["video"]["map_video"]["figsize"]       = (1000, 1000)
options["plot"]["video"]["map_video"]["proj"]          = "+proj=merc"
options["plot"]["video"]["map_video"]["title"]         = "Cascadia"
options["plot"]["video"]["map_video"]["show_progress"] = true
options["plot"]["video"]["map_video"]["framerate"]     = 20
options["plot"]["video"]["map_video"]["n_lats_edges"]  = 101
options["plot"]["video"]["map_video"]["Δt"]            = 1.0
options["plot"]["video"]["map_video"]["t0"]            = 2019.5
options["plot"]["video"]["map_video"]["t0_time"]       = true
options["plot"]["video"]["map_video"]["dir_output"]    =  dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] * 
                                                        solution_date
options["plot"]["video"]["map_video"]["output"] = "slip_potency_rate_map"

# if abs.(color) .< color_thresh_perc, set color to 0
options["plot"]["video"]["map_video"]["color_thresh_perc"] = 0 #0.15

make_video_map(slip_potency_rate,
               tremors,
               fault,
               options["plot"]["video"]["map_video"]
               )

###############
### FIGURES ###
###############
options["plot"]["figures"] = Dict()

# Intro figure
options["plot"]["figures"]["intro"] = Dict()
options["plot"]["figures"]["intro"]["figsize"] = (1000, 1000)
options["plot"]["figures"]["intro"]["limits"]  = (-129.0, # lon_min
                                                  -117.5, # lon_max
                                                   39.0,   # lat_min
                                                   51.5)   # lat_max
# options["plot"]["figures"]["intro"]["proj"]  = :Mercator
options["plot"]["figures"]["intro"]["title"] = "Cascadia"
options["plot"]["figures"]["intro"]["volcanoes"] = 
                                options["scen"]["select"]["origin"][2:end,:]
options["plot"]["figures"]["intro"]["stns2plot"] = ["ALBH", "P161"]
# Tectonic plate boundaries in json format downloaded from
# https://github.com/fraxen/tectonicplates/tree/master?tab=readme-ov-file
# (original files from http://peterbird.name/oldFTP/PB2002/)
options["plot"]["figures"]["intro"]["tect_boundaries"] = dirs["dir_data"] *
                "Plates/tectonicplates-master/GeoJSON/PB2002_boundaries.json"
options["plot"]["figures"]["intro"]["dir_output"] = dirs["dir_results"] *
                                            "Slowquakes/real-time/"*
                                            dirs["dir_case"] * solution_date
options["plot"]["figures"]["intro"]["output"] = "intro"
fig = plot_intro_map_gmt(X, fault, options["plot"]["figures"]["intro"])


# Independent Components
options["plot"]["figures"]["comps"] = Dict()
options["plot"]["figures"]["comps"]["figsize"] = (500, 1000)
options["plot"]["figures"]["comps"]["ll_legend"] = (-128.2, 40.0)
options["plot"]["figures"]["comps"]["limits"]        = (-129.0, # lon_min
                                                        -117.5, # lon_max
                                                        39.0,   # lat_min
                                                        51.5)   # lat_max
options["plot"]["figures"]["comps"]["proj"]  = "+proj=merc"
options["plot"]["figures"]["comps"]["title"] = "IC "
options["plot"]["figures"]["comps"]["show_progress"] = true
options["plot"]["figures"]["comps"]["dir_output"] = dirs["dir_results"] *
                                                "Slowquakes/real-time/"*
                                                dirs["dir_case"] * solution_date
options["plot"]["figures"]["comps"]["output"] = "IC"
# fig = plot_comp(ICA, options["plot"]["figures"]["comps"])
fig = plot_comp_gmt(ICA, options["plot"]["figures"]["comps"])


# Power spectrum of the ICs
# options["PSD"] = Dict()
#  # sample rate of the signal
# options["PSD"]["fs"] = 365.25; # number of samples per year
# options["PSD"]["subplot_grid"] = (4,2)
# options["PSD"]["output"] = dirs["dir_results"] * "Slowquakes/real-time/"* 
#                             dirs["dir_case"] * "PSD"
# options["PSD"]["color_list"] = ["193/43/43", "0/84/147", "30/150/65", "232/56/137"]
# # fig = plot_psd_gmt(ICA, options["PSD"])


# Map latitude-time and latitude-time with rate time series
options["plot"]["figures"]["map_ts"] = Dict()
options["plot"]["figures"]["map_ts"]["figsize"] = (1510, 900)
options["plot"]["figures"]["map_ts"]["title"]   = "Cascadia"
options["plot"]["figures"]["map_ts"]["n_lats_edges"]  = 101
options["plot"]["figures"]["map_ts"]["color_thresh_perc"] = 0.0
options["plot"]["figures"]["map_ts"]["dir_output"]   = dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    solution_date
options["plot"]["figures"]["map_ts"]["output"] = "slip_potency_rate"
fig = make_figure_map_lat_time(slip_potency_rate, fault,
    options["plot"]["figures"]["map_ts"])
fig = make_figure_map_lat_time_ts(slip_potency_rate, fault,
    options["plot"]["figures"]["map_ts"])


# Selected time series
options["plot"]["figures"]["ts"] = Dict()
options["plot"]["figures"]["ts"]["figsize"] = (1000, 1000)
options["plot"]["figures"]["ts"]["dir_output"] = dirs["dir_results"] *
                                                "Slowquakes/real-time/"*
                                                dirs["dir_case"] * solution_date
options["plot"]["figures"]["ts"]["name"] = "ALBH"
fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
options["plot"]["figures"]["ts"]["name"] = "P161"
fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])




###############################
### DYNAMICAL LOCAL INDICES ###
###############################
# options["local_indices"] = Dict()
# options["local_indices"]["p"]   = 0.99
# options["local_indices"]["est"] = :exp

# d, θ = calc_d_θ(slip_potency_rate, options["local_indices"])

# make_video_map_dθ(slip_potency_rate, tremors, fault, d, θ,
#                     options["plot"]["video"]["map_video"])

# make_video_dtheta(options["plot"]["movie"])