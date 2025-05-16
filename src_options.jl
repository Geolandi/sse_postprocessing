using Dates

#########################
### SELECT COMPONENTS ###
#########################
function read_options_select_comps(options)
    # options["inversion"]["select_comps"]["frequency_analysis"]["f_threshold"] =
    #                                                             3.5
    # options["inversion"]["select_comps"]["frequency_analysis"]["f_skip"] = 
    #                                                             [1.0 2.0 3.0]
    # options["inversion"]["select_comps"]["frequency_analysis"]["sigmaf_skip"] =
    #                                                         [0.05 0.05 0.05]
    # options["inversion"]["select_comps"]["frequency_analysis"]["cs_psd1"] = 0.5
    options["inversion"]["select_comps"]["frequency_analysis"]["f_psd_st"] = 0.1
    options["inversion"]["select_comps"]["frequency_analysis"]["cs_psd_st"]= 0.9
    
    return options
end

###########
### PSD ###
###########
function read_options_psd(options)
    options["PSD"] = Dict()
    # sample rate of the signal
    options["PSD"]["fs"] = 365.25; # number of samples per year
    options["PSD"]["n_freq"] =
        options["inversion"]["select_comps"]["frequency_analysis"]["n_freq"];
    options["PSD"]["n_sample"] =
        options["inversion"]["select_comps"]["frequency_analysis"]["n_sample"];
    options["PSD"]["f_last"] = 3.655466864042958 # to reproduce Matlab results

    return options
end

####################
### COMMON MODES ###
####################
function read_options_common_modes(options)
    # options["inversion"]["select_comps"]["common_mode_perc"] = 0.95
    options["inversion"]["select_comps"]["common_mode_stddist"] = 3
    return options
end

#################
### INVERSION ###
#################
function read_options_inversion(options)
    # rake_pos determines the rake for which we need to consider the slip as
    # positive; opposite directions will be negative;
    # 90 = thrust; -90 = normal
    options["inversion"]["rake_pos"] = 90;

    return options
end

#################
### SMOOTHING ###
#################
function read_options_smoothing(options)
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
        options["smooth"]["filter_cutofffreq"] = 1/7 #1/28
        options["smooth"]["fs"] = 1
        options["smooth"]["filter_type"] = "hamming" #"hanning"
        options["smooth"]["filter_window"] = 365
        options["smooth"]["scale"] = true
        options["smooth"]["causal"]  = false
    end
    return options
end

#################
### SLIP RATE ###
#################
function read_options_sliprate(options)
    options["slip_rate"] = Dict()
    options["slip_rate"]["rake"] = 90;
    options["slip_rate"]["windowsize"] = 7;

    return options
end

###############
### TREMORS ###
###############
function read_options_tremors(options, fault)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

    options["tremors"] = Dict()
    options["tremors"]["t0"] = options["scen"]["first_epoch"]
    options["tremors"]["t1"] = Date(yyyy,mm,dd) + Day(1)
    options["tremors"]["origin"] = fault["origin"]

    return options
end

##############
### MOVIES ###
##############
function read_options_movie_map(options, dirs)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

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
    options["plot"]["video"]["map_video"]["t0"] = options["scen"]["first_epoch"]
    options["plot"]["video"]["map_video"]["t1"] = Date(yyyy,mm,dd) + Day(1)
    
    options["plot"]["video"]["map_video"]["n_future_days"]  = 60
    options["plot"]["video"]["map_video"]["dir_coastlines"] = 
                                                    dirs["dir_coastlines"]
    # if abs.(color) .< color_thresh_perc, set color to 0
    options["plot"]["video"]["map_video"]["color_thresh_perc"] = 0 #0.15
    options["plot"]["video"]["map_video"]["dir_output"] =  dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] * 
                                                        options["solution_date"]
    options["plot"]["video"]["map_video"]["output"] = "slip_potency_rate_map"

    return options
end

function read_options_movie_map_northsouth(options, dirs)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

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
    options["plot"]["video"]["map_video"]["t0"] = options["scen"]["first_epoch"]
    options["plot"]["video"]["map_video"]["t1"] = Date(yyyy,mm,dd) + Day(30)
    
    options["plot"]["video"]["map_video"]["n_future_days"] = 60
    options["plot"]["video"]["map_video"]["lat_split"] = 44
    options["plot"]["video"]["map_video"]["dir_coastlines"] = 
                                                    dirs["dir_coastlines"]
    # if abs.(color) .< color_thresh_perc, set color to 0
    options["plot"]["video"]["map_video"]["color_thresh_perc"] = 0 #0.15
    options["plot"]["video"]["map_video"]["dir_output"] =  dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] * 
                                                        options["solution_date"]
    options["plot"]["video"]["map_video"]["output"] = "slip_potency_rate_map_ns"

    return options
end

###############
### FIGURES ###
###############
# INTRO
function read_options_fig_intro(options, dirs)
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
    options["plot"]["figures"]["intro"]["stn_legend"] = [-128.0 40.0]
    options["plot"]["figures"]["intro"]["stn_example_legend"] = [-128.0 39.6]
    options["plot"]["figures"]["intro"]["ll_volcanoes_legend"] = [-128.0 39.2]
    options["plot"]["figures"]["intro"]["stns2plot"] = ["ALBH", "P160"]
    options["plot"]["figures"]["intro"]["dx_stns2plot"] = [-1.0 -0.8]
    options["plot"]["figures"]["intro"]["dy_stns2plot"] = [0.15 0.2]
    # Tectonic plate boundaries in json format downloaded from
    # https://github.com/fraxen/tectonicplates/tree/master?tab=readme-ov-file
    # (original files from http://peterbird.name/oldFTP/PB2002/)
    options["plot"]["figures"]["intro"]["tect_boundaries"] = dirs["dir_data"] *
                    "Plates/tectonicplates-master/GeoJSON/PB2002_boundaries.json"
    options["plot"]["figures"]["intro"]["dir_output"] = dirs["dir_results"] *
                                                "Slowquakes/real-time/"*
                                                dirs["dir_case"] *
                                                options["solution_date"]
    options["plot"]["figures"]["intro"]["output"] = "intro"
    return options
end

# TIME SERIES
function read_options_fig_ts(options, dirs)
    options["plot"]["figures"]["ts"] = Dict()
    options["plot"]["figures"]["ts"]["figsize"] = (1000, 1000)
    options["plot"]["figures"]["ts"]["annot_y_e"] = 5
    options["plot"]["figures"]["ts"]["annot_y_n"] = 5
    options["plot"]["figures"]["ts"]["annot_y_u"] = 10
    options["plot"]["figures"]["ts"]["ticks_y_e"] = 1
    options["plot"]["figures"]["ts"]["ticks_y_n"] = 1
    options["plot"]["figures"]["ts"]["ticks_y_u"] = 2
    options["plot"]["figures"]["ts"]["ms"] = 0.02
    options["plot"]["figures"]["ts"]["ml"] = 0.01
    options["plot"]["figures"]["ts"]["p"] = 1
    options["plot"]["figures"]["ts"]["alpha"] = 0.5
    options["plot"]["figures"]["ts"]["dir_output"] = dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    options["solution_date"]
    return options
end

# TIME SERIES
function read_options_fig_2ts(options, dirs)
    options["plot"]["figures"]["2ts"] = Dict()
    options["plot"]["figures"]["2ts"]["figsize"] = (1000, 1000)
    options["plot"]["figures"]["2ts"]["names"] = ["ALBH" "P160"]
    options["plot"]["figures"]["2ts"]["annot_y_e"] = [5 5]
    options["plot"]["figures"]["2ts"]["annot_y_n"] = [5 5]
    options["plot"]["figures"]["2ts"]["annot_y_u"] = [10 10]
    options["plot"]["figures"]["2ts"]["ticks_y_e"] = [1 1]
    options["plot"]["figures"]["2ts"]["ticks_y_n"] = [1 1]
    options["plot"]["figures"]["2ts"]["ticks_y_u"] = [2 2]
    options["plot"]["figures"]["2ts"]["ms"] = 0.02
    options["plot"]["figures"]["2ts"]["ml"] = 0.01
    options["plot"]["figures"]["2ts"]["p"] = 1
    options["plot"]["figures"]["2ts"]["alpha"] = 0.5
    options["plot"]["figures"]["2ts"]["dir_output"] = dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    options["solution_date"]
    return options
end

# TIME SERIES
function read_options_fig_ts_fit(options, dirs, comps_selected)
    options["plot"]["figures"]["ts_fit"] = Dict()
    options["plot"]["figures"]["ts_fit"]["figsize"] = (1000, 1000)
    options["plot"]["figures"]["ts_fit"]["annot_y_e"] = 2
    options["plot"]["figures"]["ts_fit"]["annot_y_n"] = 1
    options["plot"]["figures"]["ts_fit"]["annot_y_u"] = 5
    options["plot"]["figures"]["ts_fit"]["ticks_y_e"] = 1
    options["plot"]["figures"]["ts_fit"]["ticks_y_n"] = 0.5
    options["plot"]["figures"]["ts_fit"]["ticks_y_u"] = 1
    options["plot"]["figures"]["ts_fit"]["ms"] = 0.02
    options["plot"]["figures"]["ts_fit"]["ml"] = 0.01
    options["plot"]["figures"]["ts_fit"]["p"] = 1
    options["plot"]["figures"]["ts_fit"]["alpha"] = 0.5
    options["plot"]["figures"]["ts_fit"]["comps_selected"] = comps_selected
    options["plot"]["figures"]["ts_fit"]["dir_output"] = dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    options["solution_date"]
    return options
end

# COMPONENTS
function read_options_fig_comps(options, dirs)
    options["plot"]["figures"]["comps"] = Dict()
    options["plot"]["figures"]["comps"]["figsize"] = (500, 1000)
    options["plot"]["figures"]["comps"]["ll_legend"] = (-128.2, 40.0)
    options["plot"]["figures"]["comps"]["limits"]        = (-129.0, # lon_min
                                                            -117.5, # lon_max
                                                            39.0,   # lat_min
                                                            51.5)   # lat_max
    options["plot"]["figures"]["comps"]["proj"]  = "+proj=merc"
    options["plot"]["figures"]["comps"]["title"] = "IC "
    options["plot"]["figures"]["comps"]["unit_output"] = 
                                                options["scen"]["unit_output"]
    options["plot"]["figures"]["comps"]["ms"] = 0.02
    options["plot"]["figures"]["comps"]["ml"] = 0.01
    options["plot"]["figures"]["comps"]["p"] = 1
    options["plot"]["figures"]["comps"]["alpha"] = 0.5
    options["plot"]["figures"]["comps"]["show_progress"] = true
    options["plot"]["figures"]["comps"]["dir_output"] = dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    options["solution_date"]
    options["plot"]["figures"]["comps"]["output"] = "IC"

    return options
end

# COMPONENTS FIT
function read_options_fig_comps_fit(options, dirs, comps_selected)
    options["plot"]["figures"]["comps_fit"] = Dict()
    options["plot"]["figures"]["comps_fit"]["figsize"] = (500, 1000)
    options["plot"]["figures"]["comps_fit"]["ll_legend"] = (-127.7, 40.6)
    options["plot"]["figures"]["comps_fit"]["ll_legend_model"] = (-127.7, 39.8)
    options["plot"]["figures"]["comps_fit"]["limits"]    = (-129.0, # lon_min
                                                            -117.5, # lon_max
                                                            39.0,   # lat_min
                                                            51.5)   # lat_max
    options["plot"]["figures"]["comps_fit"]["proj"]  = "+proj=merc"
    options["plot"]["figures"]["comps_fit"]["title"] = "Fit IC "
    options["plot"]["figures"]["comps_fit"]["unit_output"] = 
                                                options["scen"]["unit_output"]
    options["plot"]["figures"]["comps_fit"]["show_progress"] = true
    options["plot"]["figures"]["comps_fit"]["comps_selected"] = comps_selected
    options["plot"]["figures"]["comps_fit"]["dir_output"] = 
                                                    dirs["dir_results"] *
                                                    "Slowquakes/real-time/"*
                                                    dirs["dir_case"] *
                                                    options["solution_date"]
    options["plot"]["figures"]["comps_fit"]["output"] = "fit_IC"

    return options
end

# PSD
function read_options_fig_psd(options, dirs, comps_selected)
    options["plot"]["PSD"] = copy(options["PSD"])
    options["plot"]["PSD"]["subplot_grid"] = (4,2)
    options["plot"]["PSD"]["ind_comps"] = comps_selected["ind_comps"]
    options["plot"]["PSD"]["ind_comps_complex"] = 
                                    comps_selected["PSD"]["ind_comps_complex"]
    options["plot"]["PSD"]["ind_comps_complex_removed"] = 
                            comps_selected["PSD"]["ind_comps_complex_removed"]
    options["plot"]["PSD"]["ind_comps_common"] =
                            comps_selected["common_modes"]["ind_comps_common"]
    options["plot"]["PSD"]["ind_comps_common_removed"] =
                    comps_selected["common_modes"]["ind_comps_common_removed"]

    options["plot"]["PSD"]["label"] = "IC "
    options["plot"]["PSD"]["title"] = "Power Spectral Density"
    options["plot"]["PSD"]["dir_output"] = dirs["dir_results"] *
                                            "Slowquakes/real-time/"* 
                                            dirs["dir_case"] *
                                            options["solution_date"]
    options["plot"]["PSD"]["output"] = "PSD"
    options["plot"]["PSD"]["color_list"] = ["0/84/147", "255/165/0",
                                            "30/150/65", "193/43/43",
                                            "232/56/137"]
                                            
    return options
end

# PARAMETERS SELECTION
function read_options_selection_params(options, dirs, comps_selected)
    options["plot"]["selection_params"] = Dict()
    # options["plot"]["selection_params"]["PSD"] = options["PSD"]
    options["plot"]["selection_params"]["sigma"] = 
                                        options["inversion"]["sigma"][1,:]
    options["plot"]["selection_params"]["ind_comps"] = 
                                        comps_selected["ind_comps"]
    options["plot"]["selection_params"]["ind_comps_complex"] = 
                            comps_selected["PSD"]["ind_comps_complex"]
    options["plot"]["selection_params"]["ind_comps_complex_removed"] = 
                            comps_selected["PSD"]["ind_comps_complex_removed"]
    options["plot"]["selection_params"]["ind_comps_common"] = 
                            comps_selected["common_modes"]["ind_comps_common"]
    options["plot"]["selection_params"]["ind_comps_common_removed"] = 
                    comps_selected["common_modes"]["ind_comps_common_removed"]
    # options["plot"]["selection_params"]["ind_comps_removed"] = ind_comps_removed
    options["plot"]["selection_params"]["titlea"] = "ICs selection"
    options["plot"]["selection_params"]["titleb"] = "Smoothing parameter selection"
    options["plot"]["selection_params"]["title"] = "Automatic selections"
    options["plot"]["selection_params"]["dir_output"] = dirs["dir_results"] *
                                                        "Slowquakes/real-time/"* 
                                                        dirs["dir_case"] *
                                                        options["solution_date"]
    options["plot"]["selection_params"]["output"] = "automatic_selection"
    return options
end

# MAP LAT-TIME
function read_options_map_ts(options, dirs)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

    options["plot"]["figures"]["map_ts"] = Dict()
    options["plot"]["figures"]["map_ts"]["figsize"] = (1510, 900)
    options["plot"]["figures"]["map_ts"]["title"]   = "Cascadia"
    options["plot"]["figures"]["map_ts"]["n_lats_edges"]  = 101
    options["plot"]["figures"]["map_ts"]["color_thresh_perc"] = 0.0
    options["plot"]["figures"]["map_ts"]["Δt"] = 60/365.25
    options["plot"]["figures"]["map_ts"]["t_max_crosscorr"] = 14
    options["plot"]["figures"]["map_ts"]["t0"] = options["scen"]["first_epoch"]
        options["plot"]["figures"]["map_ts"]["t1"] = Date(yyyy,mm,dd) + Day(1)
    options["plot"]["figures"]["map_ts"]["xticks_dist"] = 2
    options["plot"]["figures"]["map_ts"]["dir_output"] = dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] *
                                                        options["solution_date"]
    options["plot"]["figures"]["map_ts"]["output"] = "slip_potency_rate"
    return options
end

# MAP LAT-TIME
function read_options_map_ts_northsouth(options, dirs)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

    options["plot"]["figures"]["map_ts"] = Dict()
    options["plot"]["figures"]["map_ts"]["figsize"] = (1510, 900)
    options["plot"]["figures"]["map_ts"]["title"]   = "Cascadia"
    options["plot"]["figures"]["map_ts"]["n_lats_edges"]  = 101
    options["plot"]["figures"]["map_ts"]["color_thresh_perc"] = 0.0
    options["plot"]["figures"]["map_ts"]["Δt"] = 60/365.25
    options["plot"]["figures"]["map_ts"]["t_max_crosscorr"] = 14
    options["plot"]["figures"]["map_ts"]["t0"] = options["scen"]["first_epoch"]
    options["plot"]["figures"]["map_ts"]["t1"] = Date(yyyy,mm,dd) + Day(1)
    options["plot"]["figures"]["map_ts"]["lat_split"] = 44
    options["plot"]["figures"]["map_ts"]["xticks_dist"] = 2
    options["plot"]["figures"]["map_ts"]["dir_output"] = dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] *
                                                        options["solution_date"]
    options["plot"]["figures"]["map_ts"]["output"] = "slip_potency_rate_ns"
    return options
end

# MAP LAT-TIME ZOOM
function read_options_map_ts_northsouth_zoom(options, dirs)
    yyyy = parse(Int, options["solution_date"][1:4])
    mm = parse(Int, options["solution_date"][6:7])
    dd = parse(Int, options["solution_date"][9:10])

    options["plot"]["figures"]["map_ts_zoom"] = Dict()
    options["plot"]["figures"]["map_ts_zoom"]["figsize"] = (1510, 900)
    options["plot"]["figures"]["map_ts_zoom"]["title"]   = "Cascadia"
    options["plot"]["figures"]["map_ts_zoom"]["n_lats_edges"]  = 101
    options["plot"]["figures"]["map_ts_zoom"]["color_thresh_perc"] = 0.0
    options["plot"]["figures"]["map_ts_zoom"]["Δt"] = 60/365.25
    options["plot"]["figures"]["map_ts_zoom"]["t_max_crosscorr"] = 0
    options["plot"]["figures"]["map_ts_zoom"]["t0"] = 2023.0
    options["plot"]["figures"]["map_ts_zoom"]["t1"] = Date(yyyy,mm,dd) + Day(1)
    options["plot"]["figures"]["map_ts_zoom"]["lat_split"] = 44
    options["plot"]["figures"]["map_ts_zoom"]["xticks_dist"] = 1
    options["plot"]["figures"]["map_ts_zoom"]["dir_output"] = 
                                                        dirs["dir_results"] *
                                                        "Slowquakes/real-time/"*
                                                        dirs["dir_case"] *
                                                        options["solution_date"]
    options["plot"]["figures"]["map_ts_zoom"]["output"] = 
                                                "slip_potency_rate_ns_zoom"
    return options
end