include("src_load.jl")
include("src_inversion.jl")
include("src_rate.jl")
include("src_options.jl")
include("src_time.jl")
include("src_selections.jl")

using DataFrames, CSV


solution_date = "2024-08-01/"#"2025-03-31/"#"2024-08-08/"#

mm2m = 1e-3
km2m = 1e3

println(" ")
println(solution_date)
#######################
### SET DIRECTORIES ###
#######################
dirs = Dict()
dirs["dir_results"] = "./results/"
dirs["dir_results2load"] = dirs["dir_results"]  * solution_date
dirs["dir_fault"]   = "./fault/"

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
mp = repeat(fault["area"] .* km2m^2,2,1) .* m .*mm2m;
# co-variance matrix for slip potency components
Cmp = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
for i=1:n_ICs2invert
    Cmp[i] = repeat(fault["area"] .* km2m^2,2,1) .* Cm[i] .* mm2m^2;
end

#########################
### SMOOTH COMPONENTS ###
#########################
# options = read_options_smoothing(options)
# V_smooth, timeline_smooth = create_xsmooth(
#     ICA["V"]', ICA["timeline"], options["smooth"])

####################################
### CALCULATE RATE OF COMPONENTS ###
####################################
# options = read_options_sliprate(options)
# V_dot, timeline_dot = calc_derivative(V_smooth', timeline_smooth,
#     options["slip_rate"]["windowsize"], true)


########################################
### CREATE SLIP AND SLIP RATE MODELS ###
########################################
# # create slip model
# slip = create_model(m, Cm, ICA, fault, options, ind_comps)

# option to specify what rake direction is positive
options = read_options_inversion(options)
# # set smooth ICA dictionary
# ICA_smooth = copy(ICA);
# ICA_smooth["V"] = V_smooth;
# n_samples_smooth = length(timeline_smooth)
# ICA_smooth["var_V"] = ICA["var_V"][end-n_samples_smooth+1:end,:];
# ICA_smooth["timeline"] = timeline_smooth;
# create slip model with smooth components
# slip_smooth = create_model(m, Cm, ICA_smooth, fault, options, ind_comps)

# set time derivative ICA dictionary
# ICA_dot = copy(ICA);
# ICA_dot["V"] = V_dot;
# n_samples_dot = length(timeline_dot)
# ICA_dot["var_V"] = ICA["var_V"][end-n_samples_dot+1:end,:];
# ICA_dot["timeline"] = timeline_dot;
# create slip rate model with rate of smoothed components
# slip_rate = create_model(m, Cm, ICA_dot, fault, options, ind_comps)

#########################
### SLIP POTENCY RATE ###
#########################
# calculate slip potency rate
# slip_potency_rate = create_model(mp,Cmp,ICA_dot,fault,options,ind_comps)
slip_potency = create_model(mp,Cmp,ICA,fault,options,ind_comps)
# println("Done")

#########################
print("Saving results...")
# Extract the slip strike and dip components
# slip_strike_r = slip_potency_rate["obs_strike"]
# slip_dip_r = slip_potency_rate["obs_dip"]
# Consider using the slip_potency
slip_strike = slip_potency["obs_strike"]
slip_dip = slip_potency["obs_dip"]
timeline = slip_potency["timeline"]


df_slip_strike =  DataFrame( transpose(slip_strike), :auto)
df_slip_dip =  DataFrame( transpose(slip_dip), :auto)

# df_slip_strike_r =  DataFrame( transpose(slip_strike_r), :auto)
# df_slip_dip_r =  DataFrame( transpose(slip_dip_r), :auto)

df_timeline =  DataFrame( timeline, :auto)
lambda0 = options["inversion"]["lambda0"]
df_lambda0 = DataFrame(Value = [lambda0])

CSV.write(dirs["dir_results"]  * solution_date *"slip_strike.csv", df_slip_strike)
CSV.write(dirs["dir_results"]  * solution_date *"slip_dip.csv", df_slip_dip)

# CSV.write(dirs["dir_results"]  * solution_date *"slip_strike_r.csv", df_slip_strike_r)
# CSV.write(dirs["dir_results"]  * solution_date *"slip_dip_r.csv", df_slip_dip_r)
CSV.write(dirs["dir_results"]  * solution_date *"timeline.csv", df_timeline)
CSV.write(dirs["dir_results"]  * solution_date *"lambda0.csv", df_lambda0)