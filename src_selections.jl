using DSP

function calc_psd(decomp, options)
    
    fs = options["fs"]
    n_sample = options["n_sample"]
    f_last = options["f_last"]
    
    timeline = decomp["timeline"][1,:]
    n_samples = length(timeline)
    V = decomp["V"]
    
    p = periodogram(V[:,1]; nfft=Integer(round(n_samples*n_sample)), fs=fs)
    f = DSP.freq(p)[:]
    f1 = 1/n_samples*(0:1:round(365.25*n_samples/n_freq));

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