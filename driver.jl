include("src_load.jl")
include("src_inversion.jl")
include("src_rate.jl")
include("src_options.jl")
include("src_plot.jl")
include("src_movie.jl")
include("src_tremors.jl")
include("src_time.jl")

# Add ffmpeg directory to PATH to be able to plot with GMT
ENV["PATH"] = ENV["PATH"] * ":/opt/homebrew/bin"

solution_date = "2025-01-26/"

println(" ")
println(solution_date)

mm2m = 1e-3
km2m = 1e3

#######################
### SET DIRECTORIES ###
#######################
dirs = Dict()
dirs["dir_data"] = "/Users/ag2347/Work/Data/"
dirs["dir_results"] = "/Users/ag2347/Work/Results/"
dirs["dir_case"] = "Cascadia/"
dirs["dir_results2load"] = dirs["dir_results"] * "Slowquakes/real-time/" *
                        dirs["dir_case"] * "matfiles/" * solution_date
dirs["dir_fault"]   = dirs["dir_data"] * "Faults/" * dirs["dir_case"]
dirs["dir_tremors"] = dirs["dir_data"]*"Tremors/"*dirs["dir_case"]*"PNSN/"
dirs["dir_coastlines"] = dirs["dir_data"] * "NaturalEarth/coastlines_jl/"

######################
### LOAD VARIABLES ###
######################
# load options
options = matread(dirs["dir_results2load"]*"options.mat")["options"]
options["solution_date"] = solution_date
# load position time series used as input for vbICA
X = matread(dirs["dir_results2load"]*"X.mat")["Xd"]
# load vbICA results
ICA = matread(dirs["dir_results2load"]*"ICA.mat")["ICA_essential"]
ICA["type"]     = X["type"];
ICA["llh"]      = X["llh"];
ICA["timeline"] = X["timeline"];
ICA["decmode"]  = X["decmode"];
if options["scen"]["unit_output"] == "mm"
    # convert S in m, so that later on the slip potency is in m^3 (SI units)
    ICA["S"] = ICA["S"] .* mm2m
end

##################################################
### LOAD FAULT AND CALCULATE GREENS' FUNCTIONS ###
##################################################
# load fault
fault = load_fault(dirs, options)
# calculate Greens' functions
G = create_greens_function(X, fault, options)

#########################
### SELECT COMPONENTS ###
#########################
# refine selection of components to use in SSEs reconstruction
options = read_options_select_comps(options)
comps_selected = select_comps(dirs, ICA, options)

#########################
### INVERT COMPONENTS ###
#########################
# invert selected components
ind_comps = comps_selected["ind_comps"]
misfit_comps = comps_selected["misfit_comps"]
ind_sigma0_comps = comps_selected["ind_sigma0_comps"]
m, Cm = invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, options)
println("Done")

###############################
### SLIP POTENCY COMPONENTS ###
###############################
n_ICs2invert = length(ind_comps)
n_patches = length(fault["area"])
# slip potency components
mp = repeat(fault["area"] .* km2m^2,2,1) .* m;
# co-variance matrix for slip potency components
Cmp = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
for i=1:n_ICs2invert
    Cmp[i] = repeat(fault["area"] .* km2m^2,2,1) .* Cm[i];
end

#########################
### SMOOTH COMPONENTS ###
#########################
options = read_options_smoothing(options)
V_smooth, timeline_smooth = create_xsmooth(
    ICA["V"]', ICA["timeline"], options["smooth"])

####################################
### CALCULATE RATE OF COMPONENTS ###
####################################
options = read_options_sliprate(options)
V_dot, timeline_dot = calc_derivative(V_smooth', timeline_smooth,
    options["slip_rate"]["windowsize"], true)


########################################
### CREATE SLIP AND SLIP RATE MODELS ###
########################################
# # create slip model
# slip = create_model(m, Cm, ICA, fault, options, ind_comps)

# option to specify what rake direction is positive
options = read_options_inversion(options)
# set smooth ICA dictionary
ICA_smooth = copy(ICA);
ICA_smooth["V"] = V_smooth;
n_samples_smooth = length(timeline_smooth)
ICA_smooth["var_V"] = ICA["var_V"][end-n_samples_smooth+1:end,:];
ICA_smooth["timeline"] = timeline_smooth;
# create slip model with smooth components
# slip_smooth = create_model(m, Cm, ICA_smooth, fault, options, ind_comps)

# set time derivative ICA dictionary
ICA_dot = copy(ICA);
ICA_dot["V"] = V_dot;
n_samples_dot = length(timeline_dot)
ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
ICA_dot["timeline"] = timeline_dot;
# create slip rate model with rate of smoothed components
# slip_rate = create_model(m, Cm, ICA_dot, fault, options, ind_comps)

#########################
### SLIP POTENCY RATE ###
#########################
# calculate slip potency rate
slip_potency_rate = create_model(mp,Cmp,ICA_dot,fault,options,ind_comps)

###############
### TREMORS ###
###############
# update and load tremors
update_tremors(dirs, "yesterday");
options = read_options_tremors(options, fault)
tremors = load_tremors(dirs, options["tremors"]);
t0_decyear, t0_date = get_time_decyear_and_date(
                                options["scen"]["first_epoch"])
t1_decyear, t1_date = get_time_decyear_and_date("today")
timeline2plot, dates2plot = create_timeline(t0_date, t1_date)
tremors = select_tremors(tremors, timeline2plot, fault)

####################################
### FLAGS FOR MOVIES AND FIGURES ###
####################################
flag_plot_video_map            = false
flag_plot_video_map_northsouth = false
flag_plot_fig_intro            = false
flag_plot_fig_ts               = false
flag_plot_fig_2ts              = true
flag_plot_fig_ts_misfit        = false
flag_plot_fig_ICs              = false
flag_plot_fig_ICs_fit          = false
flag_plot_fig_PSD              = false
flag_plot_fig_selectparams     = false
flag_plot_fig_map_lat_time     = false
flag_plot_fig_map_lat_time_ns  = false

##############
### MOVIES ###
##############
if flag_plot_video_map == true
    options = read_options_movie_map(options, dirs)
    make_video_map(slip_potency_rate,
                tremors,
                fault,
                options["plot"]["video"]["map_video"]
                )
end
if flag_plot_video_map_northsouth == true
    options = read_options_movie_map_northsouth(options, dirs)
    make_video_map_northsouth(slip_potency_rate,
                tremors,
                fault,
                options["plot"]["video"]["map_video"]
                )
end
            
###############
### FIGURES ###
###############
options["plot"]["figures"] = Dict()

# intro figure
if flag_plot_fig_intro == true
    options = read_options_fig_intro(options, dirs)
    fig = plot_intro_map_gmt(X, fault, options["plot"]["figures"]["intro"])
end

# selected time series
if flag_plot_fig_ts == true
    options = read_options_fig_ts(options, dirs)
    options["plot"]["figures"]["ts"]["name"] = "ALBH"
    fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
    options["plot"]["figures"]["ts"]["name"] = "P160"
    fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
    # options["plot"]["figures"]["ts"]["name"] = "PGC5"
    # fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
end
if flag_plot_fig_2ts == true
    options = read_options_fig_2ts(options, dirs)
    fig = plot_2ts_gmt(X, options["plot"]["figures"]["2ts"])
end

# independent components
if flag_plot_fig_ICs == true
    options = read_options_fig_comps(options, dirs)
    fig = plot_comp_gmt(ICA, options["plot"]["figures"]["comps"])
end

# Power Spectral Density
if flag_plot_fig_PSD == true
    options = read_options_fig_psd(options, dirs, comps_selected)
    fig = plot_psd_gmt(comps_selected["PSD"]["f"],
                    comps_selected["PSD"]["psds"], options["plot"]["PSD"])
end

# selection of ICs and smoothing parameters
if flag_plot_fig_selectparams == true
    options = read_options_selection_params(options, dirs, comps_selected)
    fig = plot_select_comps_and_smoothing(
                        comps_selected["PSD"]["fselected"],
                        comps_selected["PSD"]["cumulative_psds"],
                        misfit_comps, options["plot"]["selection_params"])
end

# fit inverted independent components
if flag_plot_fig_ICs_fit == true
    options = read_options_fig_comps_fit(options, dirs, comps_selected)
    model = Dict()
    model["U"] = G*m
    n_ts = size(G)[1]
    # co-variance matrix for modelled displacement
    model["var_U"] = [zeros(n_ts,n_ts) for i=1:n_ICs2invert]
    for i=1:n_ICs2invert
        model["var_U"][i] = G*Cm[i]*G'
    end
    fig = plot_comp_fit_gmt(ICA, model, options["plot"]["figures"]["comps_fit"])
end

# selected time series: SSEs ICs, model and misfit
if flag_plot_fig_ts_misfit == true
    X_ICs_ind_comps = copy(X)
    X̂ = Dict()
    if options["scen"]["unit_output"]=="mm"
        S = ICA["S"] ./ mm2m
    else
        S = ICA["S"]
    end
    X_ICs_ind_comps["ts"] = ICA["U"][:,ind_comps]*
            S[ind_comps,ind_comps]*ICA["V"][:,ind_comps]'
    X̂["ts"] = G*m*S[ind_comps,ind_comps]*ICA["V"][:,ind_comps]'
    
    
    X_ICs_ind_comps["var_ts"] = zeros(size(X["ts"]))
    X̂["var_ts"] = zeros(size(X̂["ts"]))
    for (i,j) in enumerate(ind_comps)
        X_ICs_ind_comps["var_ts"] = X_ICs_ind_comps["var_ts"] .+
                S[j,j].^2 * (
                    ICA["var_U"][:,j]*ICA["var_V"][:,j]' .+
                    ICA["var_U"][:,j].*(ICA["V"][:,j].^2)' .+
                    (ICA["U"][:,j].^2).*ICA["var_V"][:,j]')
        
        û = G*m[:,i]
        var_û = diag(G*Cm[i])
        X̂["var_ts"] = X̂["var_ts"] .+
                S[j,j].^2 * (
                    var_û*ICA["var_V"][:,j]' .+
                    var_û.*(ICA["V"][:,j].^2)' .+
                    (û.^2).*ICA["var_V"][:,j]')
    end
    
    options = read_options_fig_ts_fit(options, dirs, comps_selected)
    options["plot"]["figures"]["ts_fit"]["name"] = "ALBH"
    fig = plot_ts_fit_gmt(X_ICs_ind_comps, X̂,
        options["plot"]["figures"]["ts_fit"])
    options["plot"]["figures"]["ts_fit"]["name"] = "P160"
    fig = plot_ts_fit_gmt(X_ICs_ind_comps, X̂,
        options["plot"]["figures"]["ts_fit"])
end

# slip potency rate lat vs. time map
# map latitude-time and latitude-time with rate time series
if flag_plot_fig_map_lat_time == true
    options = read_options_map_ts(options, dirs)
    r = plot_map_lattimets(slip_potency_rate, tremors, fault,
        options["plot"]["figures"]["map_ts"])
    fig = plot_map_lattimets_notremors(slip_potency_rate, fault,
        options["plot"]["figures"]["map_ts"])

    options["plot"]["figures"]["map_ts"]["figsize"] = (755, 900)
    options["plot"]["figures"]["map_ts"]["xticks_dist"] = 0.5
    options["plot"]["figures"]["map_ts"]["t0"] = Date(2024,1,27)
    options["plot"]["figures"]["map_ts"]["t1"] = Date(2025,4,25) #Date(2025,1,26)
    options["plot"]["figures"]["map_ts"]["Δt"] = 0/365.25
    options["plot"]["figures"]["map_ts"]["t_max_crosscorr"] = 0
    options["plot"]["figures"]["map_ts"]["output"] = "slip_potency_rate_zoom"
    r = plot_map_lattimets(slip_potency_rate, tremors, fault,
        options["plot"]["figures"]["map_ts"])
    fig = plot_map_lattimets_notremors(slip_potency_rate, fault,
        options["plot"]["figures"]["map_ts"])
end

# slip potency rate lat vs. time map
# map latitude-time and latitude-time with rate time series
if flag_plot_fig_map_lat_time_ns == true
    options = read_options_map_ts_northsouth(options, dirs)
    r_north, r_south = plot_map_lattimets_northsouth(slip_potency_rate,
        tremors, fault,
        options["plot"]["figures"]["map_ts"])
    options = read_options_map_ts_northsouth_zoom(options, dirs)
    r_north_zoom, r_south_zoom = plot_map_lattimets_northsouth(
        slip_potency_rate, tremors, fault,
        options["plot"]["figures"]["map_ts_zoom"])
end