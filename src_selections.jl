using DSP

function calc_psd(decomp, options)
    
    fs = options["fs"]
    n_sample = options["n_sample"]
    f_last = options["f_last"]
    
    timeline = decomp["timeline"][1,:]
    n_samples = length(timeline)
    V = decomp["V"]
    n_comps = size(V)[2]
    
    p = periodogram(V[:,1]; nfft=Integer(round(n_samples*n_sample)), fs=fs)
    f = DSP.freq(p)[:]
    # f_last = 1/n_samples*(round(365.25*n_samples/n_freq))

    f0 = f[1]
    ind_last = findlast(f .< f_last)
    f1 = f[ind_last]
    f = f[1:ind_last]
    n_freqs = length(f)
    psds = zeros(n_freqs,n_comps)
    
    progress = ProgressMeter.Progress(
        n_comps; desc = "Calculating variance", enabled = true
        )
    for i=1:n_comps
        p = periodogram(V[:,i]; nfft=Integer(round(n_samples*n_sample)), fs=fs)
        psds[:,i] = DSP.power(p)[1:ind_last]
        ProgressMeter.next!(progress)
    end

    return f, psds
end

function select_complex_psds(f, psds, options)

    # load options
    # f_skip are the frequencies to skip in the calculation of the
    # cumulative PSD
    f_skip      = options["f_skip"]
    # sigmaf_skip is the uncertainty around f_skip
    # (I do not use all frequencies in the range 
    # _skip - sigmaf_skip:f_skip + sigmaf_skip)
    sigmaf_skip = options["sigmaf_skip"]
    # 
    f_threshold = options["f_threshold"]
    cs_psd1     = options["cs_psd1"]
    cs_psd2     = options["cs_psd2"]

    # number of components
    n_comp = size(psds)[2]
    # find boundaries of frequency regions to avoid in the calculation of the
    # cumulative PSD
    f_low  = f_skip - sigmaf_skip
    f_high = f_skip + sigmaf_skip    
    ind_f2keep = fill(true,size(f))
    for i=1:size(f_low)[2]
        ind_f2rm_tmp = f .>= f_low[i] .&& f .<= f_high[i]
        ind_f2keep[ind_f2rm_tmp] .= false
    end
    # keep only the frequencies in the desired regions
    fselected = f[ind_f2keep]
    # find how many frequencies are left
    n_freq_f2keep = sum(ind_f2keep)
    # initialize the cumulative PSD variable
    cpsd = zeros(n_comp,n_freq_f2keep)
    # for every component...
    for i=1:n_comp
        # normalize the psd and retain only the values at desired frequencies
        psd_norm = psds[ind_f2keep,i] ./ sum(psds[:,i])
        # calculate the cumulative sum
        cpsd[i,:] = cumsum(psd_norm)
    end
    # 
    f_thr = minimum(abs.(f[ind_f2keep] .- f_threshold))
    ind_f_thr = argmin(abs.(f[ind_f2keep] .- f_threshold))

    ind_comps_cpsd = findall(cpsd[:,ind_f_thr] .>= cs_psd1 .&&
                           cpsd[:,ind_f_thr] .<= cs_psd2)

    ind_comps2rm   = []
    ind_comps2keep = []
    n_ind_comps_cpsd = length(ind_comps_cpsd)
    n_f_skip = length(f_skip)
    for i=1:n_ind_comps_cpsd
        j = ind_comps_cpsd[i]
        psd_max = maximum(psds[:,j])
        ind_f_max = argmax(psds[:,j])
        f_max = f[ind_f_max]
        for k=1:n_f_skip
            if (f_max .>= f_skip[k] .- sigmaf_skip[k]) .&& 
                    (f_max .<= f_skip[k] .+ sigmaf_skip[k])
                append!(ind_comps2rm, i)
            else
                append!(ind_comps2keep, i)
            end
        end
    end

    ind_comps_complex = copy(ind_comps_cpsd)
    ind_comps_removed = ind_comps_cpsd[ind_comps2rm]
    deleteat!(ind_comps_complex, ind_comps2rm)

    return fselected, cpsd, ind_comps_complex, ind_comps_removed
end

function select_common_modes(decomp, options)
    common_mode_perc = options["common_mode_perc"]
    U = decomp["U"]
    n_ts, n_comp = size(U)
    n_stn = n_ts/3
    ind_comps_common   = []
    for i=1:n_comp
        flag_u_common = 0
        for j=1:3
            u = U[j:3:end,i]
            u_plus = sum(u .> 0)
            u_minus = sum(u .< 0)
            u_common = maximum([u_plus, u_minus]) ./ n_stn
            if u_common > common_mode_perc
                flag_u_common = flag_u_common + 1
            end
        end
        if flag_u_common > 0
            append!(ind_comps_common, i)
        end
    end
    return ind_comps_common
end

function select_comps2invert(decomp, f, psds, options)
    
    fselected, cpsd, ind_comps_complex, ind_comps_complex_removed = 
            select_complex_psds(f, psds, options["frequency_analysis"])
    ind_comps_common = select_common_modes(decomp, options)

    ind_comps_common_removed = intersect(ind_comps_complex, ind_comps_common)
    
    if !isempty(ind_comps_common_removed)
        ind_comps =
                ind_comps_complex[ind_comps_common_removed.!=ind_comps_complex]
    else
        ind_comps = ind_comps_complex
    end

    return fselected, cpsd, ind_comps, ind_comps_complex,
            ind_comps_complex_removed, ind_comps_common,
            ind_comps_common_removed
end


function select_comps(dirs, decomp, options)

    # load misfit of inverted components
    misfit_comps_calculated = matread(
        dirs["dir_results2load"]*"misfit_comps.mat")["misfit_comps"]
    
    # find smoothing parameters for inversion
    min_misfit_comps_calculated, ind_sigma0_comps_Cartesian_calculated = 
                                        findmin(misfit_comps_calculated, dims=1)
    ind_sigma0_comps_calculated = getindex.(
                        ind_sigma0_comps_Cartesian_calculated, 1)

    # load indices of inverted components
    ind_comps_calculated = matread(
                    dirs["dir_results2load"]*"ind_comps.mat")["ind_comps"]
    ind_comps_calculated = round.(Int, ind_comps_calculated)[:,1]
    
    # further the selection removing possible common mode components
    options = read_options_psd(options)

    freqs, psds = calc_psd(decomp, options["PSD"])

    # freqs_selected, cumulative_psds, ind_comps_complex,
    #     ind_comps_complex_removed = 
    #         select_complex_psds(
    #             freqs,
    #             psds,
    #             options["inversion"]["select_comps"]["frequency_analysis"]
    #             )

    options = read_options_common_modes(options)
    
    freqs_selected, cumulative_psds, ind_comps, ind_comps_complex,
        ind_comps_complex_removed, ind_comps_common,
        ind_comps_common_removed = select_comps2invert(
                                        decomp, freqs, psds,
                                        options["inversion"]["select_comps"]
                                        )

    intersected_comps = intersect(ind_comps_calculated,ind_comps)
    ind_in_ind_comps_calculated = findall(x -> x in intersected_comps, 
                                        ind_comps_calculated)
    
    misfit_comps = misfit_comps_calculated[:,ind_in_ind_comps_calculated]
    min_misfit_comps = min_misfit_comps_calculated[ind_in_ind_comps_calculated]
    ind_sigma0_comps = ind_sigma0_comps_calculated[ind_in_ind_comps_calculated]

    results = Dict()
    results["misfit_comps"] = misfit_comps
    results["ind_sigma0_comps"] = ind_sigma0_comps
    results["ind_comps"] = ind_comps
    results["PSD"] = Dict()
    results["PSD"]["f"] = freqs
    results["PSD"]["psds"] = psds
    results["PSD"]["fselected"] = freqs_selected
    results["PSD"]["cumulative_psds"] = cumulative_psds
    results["PSD"]["ind_comps_complex"] = ind_comps_complex
    results["PSD"]["ind_comps_complex_removed"] = ind_comps_complex_removed
    results["common_modes"] = Dict()
    results["common_modes"]["ind_comps_common"] = ind_comps_common
    results["common_modes"]["ind_comps_common_removed"] = 
                                                ind_comps_common_removed

    return results
end