import ProgressMeter
using GLM
using RollingFunctions
using DataFrames
using Loess
using DSP

function create_xsmooth(X, timeline, method)
    n_features, n_samples = size(X)
    method_name = method["name"]
    if method_name == "rollmean"
        windowsize = method["windowsize"]
        centered = method["centered_output"]
        Xsmooth = zeros(n_features, n_samples-windowsize+1)
        for i = 1:n_features
            Xsmooth[i,:] = rollmean(X[i,:], windowsize)
        end

        timelinesmooth = zeros(1, n_samples-windowsize+1)
        if centered==true
            half_windowsize = Int((windowsize-1)/2)
            timelinesmooth[1,:] = timeline[1,1+half_windowsize:end-half_windowsize]
        else
            timelinesmooth = timeline[1,windowsize+1:end]
        end
    elseif method_name == "loess"
        span = method["span"]
        Xsmooth = zeros(n_features, n_samples)
        for i=1:n_features
            model = loess(timeline[1,:], X[i,:], span=span)
            Xsmooth[i,:] = predict(model, timeline[1,:])
        end
        timelinesmooth = timeline
    elseif method_name == "lowpass"
        fs = method["fs"]
        Wn = method["filter_cutofffreq"]
        filter_win = method["filter_window"]
        responsetype = Lowpass(Wn; fs)
        if method["filter_type"] == "hanning"
            designmethod = FIRWindow(hanning(filter_win; zerophase=false))
        end
        Xsmooth = zeros(n_features, n_samples)
        for i=1:n_features
            if method["causal"] == true
                Xsmooth[i,:] = filt(digitalfilter(responsetype, designmethod),
                            X[i,:])
            else
                x_tr = X[i,end]
                Xsmooth[i,:] = x_tr .+ filtfilt(digitalfilter(
                            responsetype, designmethod), X[i,:] .- x_tr)
            end
        end
        timelinesmooth = timeline
    end

    return Xsmooth', timelinesmooth
end

function calc_derivative(X, timeline, windowsize=7, centered=true)
    n_features, n_samples = size(X)
    Vdot = zeros(n_features, n_samples-windowsize+1)
    for i = 1:n_features
        for t = 1:n_samples-windowsize+1
            df = DataFrame(x = timeline[1,t:t+windowsize-1] .- timeline[1,t])
            df.y = X[i,t:t+windowsize-1]
            model = lm(@formula(y ~ 1 + x), df)
            v = coef(model)[2]
            Vdot[i,t] = v
            df = nothing
        end
    end
    timelinedot = zeros(1, n_samples-windowsize+1)
    if centered==true
        half_windowsize = Int((windowsize-1)/2)
        timelinedot[1,:] = timeline[1,1+half_windowsize:end-half_windowsize]
    else
        timelinedot = timeline
    end
    return Vdot', timelinedot
end