include("src_load.jl")
include("src_inversion.jl")
include("src_options.jl")
include("src_selections.jl")
include("src_time.jl")
include("src_rate.jl")

solution_date = "2025-01-26/"
# solution_date = "2025-04-25/"

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
# option to specify what rake direction is positive
options = read_options_inversion(options)
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

#################################
### CREATE SLIP POTENCY MODEL ###
#################################
# create slip potency model
slip_potency = create_model(mp, Cmp, ICA, fault, options, ind_comps)


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

######################################
### CREATE SLIP POTENCY RATE MODEL ###
######################################
# option to specify what rake direction is positive
options = read_options_inversion(options)
# set time derivative ICA dictionary
ICA_dot = copy(ICA);
ICA_dot["V"] = V_dot;
n_samples_dot = length(timeline_dot)
ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
ICA_dot["timeline"] = timeline_dot;

# calculate slip potency rate
slip_potency_rate = create_model(mp,Cmp,ICA_dot,fault,options,ind_comps)