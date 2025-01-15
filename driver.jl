include("src_load.jl")
include("src_inversion.jl")
include("src_rate.jl")
include("src_localindices.jl")
include("src_options.jl")
include("src_plot.jl")
include("src_movie.jl")
include("src_tremors.jl")
include("src_time.jl")

# Add ffmpeg directory to PATH to be able to plot with GMT
ENV["PATH"] = ENV["PATH"] * ":/opt/homebrew/bin"

solutions_dates = ["2024-08-01/",
                   "2024-08-02/",
                   "2024-08-03/",
                   "2024-08-04/",
                   "2024-08-06/",
                   "2024-08-07/",
                   "2024-08-08/",
                   "2024-08-09/",
                   "2024-08-10/",
                   "2024-08-11/",
                   "2024-08-12/",
                   "2024-08-13/",
                   "2024-08-14/",
                   "2024-08-16/",
                   "2024-08-17/",
                   "2024-08-18/",
                   "2024-08-19/",
                   "2024-08-20/",
                   "2024-08-21/",
                   "2024-08-22/",
                   "2024-08-23/",
                   "2024-08-24/",
                   "2024-08-26/",
                   "2024-08-28/",
                   "2024-08-29/",
                   "2024-08-30/",
                   "2024-08-31/",
                   ]
solution_dates = ["2024-08-24/"]
# solution_date = "2024-08-24/"
for solution_date in solutions_dates

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
    options      = matread(dirs["dir_results2load"]*"options.mat")["options"]
    options["solution_date"] = solution_date
    # load position time series used as input for vbICA
    X            = matread(dirs["dir_results2load"]*"X.mat")["Xd"]
    # load vbICA results
    ICA          = matread(dirs["dir_results2load"]*"ICA.mat")["ICA_essential"]
    ICA["type"]     = X["type"];
    ICA["llh"]      = X["llh"];
    ICA["timeline"] = X["timeline"];
    ICA["decmode"]  = X["decmode"];

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
    mp = repeat(fault["area"],2,1) .* m;
    # co-variance matrix for slip potency components
    Cmp = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
    for i=1:n_ICs2invert
        Cmp[i] = repeat(fault["area"],2,1) .* Cm[i];
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
    # option to specify what rake direction is positive
    options = read_options_inversion(options)
    # set smooth ICA dictionary
    ICA_smooth = copy(ICA);
    ICA_smooth["V"] = V_smooth;
    n_samples_smooth = length(timeline_smooth)
    ICA_smooth["var_V"] = ICA["var_V"][end-n_samples_smooth+1:end,:];
    ICA_smooth["timeline"] = timeline_smooth;
    # create slip model with smooth components
    slip_smooth = create_model(m, Cm, ICA_smooth, fault, options, ind_comps)

    # set time derivative ICA dictionary
    ICA_dot = copy(ICA);
    ICA_dot["V"] = V_dot;
    n_samples_dot = length(timeline_dot)
    ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
    ICA_dot["timeline"] = timeline_dot;
    # create slip rate model with rate of smoothed components
    slip_rate = create_model(m, Cm, ICA_dot, fault, options, ind_comps)

    #########################
    ### SLIP POTENCY RATE ###
    #########################
    # calculate slip potency rate
    slip_potency_rate = create_model(mp,Cmp,ICA_dot,fault,options,ind_comps)

    ###############
    ### TREMORS ###
    ###############
    # update and load tremors
    #update_tremors(dirs, "today");
    tremors = load_tremors(dirs, fault["origin"]);
    t0_decyear, t0_date = get_time_decyear_and_date(timeline_dot[1])
    t1_decyear, t1_date = get_time_decyear_and_date("today")
    timeline2plot, dates2plot = create_timeline(t0_date, t1_date)
    tremors = select_tremors(tremors, timeline2plot, fault)
    # tremors = select_tremors(tremors, timeline_dot, fault)


    ####################################
    ### FLAGS FOR MOVIES AND FIGURES ###
    ####################################
    flag_plot_video_map        = true
    flag_plot_fig_intro        = true
    flag_plot_fig_ts           = true
    flag_plot_fig_ICs          = true
    flag_plot_fig_PSD          = true
    flag_plot_fig_selectparams = true
    flag_plot_fig_map_lat_time = true

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
                
    ###############
    ### FIGURES ###
    ###############
    options["plot"]["figures"] = Dict()

    # intro figure
    if flag_plot_fig_intro == true
        options = read_options_fig_intro(options, dirs)
        fig = plot_intro_map_gmt(X, fault, options["plot"]["figures"]["intro"])
    end

    # selecter time series
    if flag_plot_fig_ts == true
        options = read_options_fig_ts(options, dirs)
        options["plot"]["figures"]["ts"]["name"] = "ALBH"
        fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
        options["plot"]["figures"]["ts"]["name"] = "P161"
        fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
        options["plot"]["figures"]["ts"]["name"] = "PGC5"
        fig = plot_ts_gmt(X, options["plot"]["figures"]["ts"])
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


    # slip potency rate lat vs. time map
    # map latitude-time and latitude-time with rate time series
    if flag_plot_fig_map_lat_time == true
        options = read_options_map_ts(options, dirs)
        r = plot_map_lattimets(slip_potency_rate, tremors, fault,
            options["plot"]["figures"]["map_ts"])
        fig = plot_map_lattimets_notremors(slip_potency_rate, fault,
            options["plot"]["figures"]["map_ts"])
    end
end



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

# using GLMakie
# me = zeros(3,size(ICA["U"])[2])
# st = zeros(3,size(ICA["U"])[2])
# sk = zeros(3,size(ICA["U"])[2])
# for i=1:size(ICA["U"])[2]
#     for j=1:3
#         fig_hist = Makie.Figure()
#         ax_hist = Makie.Axis(fig_hist[1,1])
#         Makie.hist!(ICA["U"][j:3:end,i])
#         #save("./tmp/hist_U"*string(i)*"_"*string(j)*".png", fig_hist)
#         me[j,i] = mean(ICA["U"][j:3:end,i])
#         st[j,i] = std(ICA["U"][j:3:end,i])
#         sk[j,i] = skewness(ICA["U"][j:3:end,i])
#     end
# end

# fig = Makie.Figure()
# ax = Makie.Axis(fig[1,1])
# Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), st[1,:])
# Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), st[2,:])
# Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), st[3,:])
# # Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), sk[1,:])
# # Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), sk[2,:])
# # Makie.lines!(ax, range(1,size(ICA["U"])[2],step=1), sk[3,:])
# display(fig)