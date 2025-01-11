using Dates
using DateFormats

function get_time_decyear_and_date(t)
    if typeof(t) == Date
        t_date = t
        t_decyear = DateFormats.yeardecimal.(t_date)
    elseif typeof(t) == String
        if t == "today"
            t_date = today()
            t_decyear = DateFormats.yeardecimal.(t_date)
        else
            println("`t` not recognized")
        end
    else
        t_decyear = t
        t_date = Date(DateFormats.yeardecimal.(t_decyear));
    end
    return t_decyear, t_date
end

function create_timeline(t0,t_end)

    if typeof(t0) == Date
        date0 = t0
    elseif typeof(t0) == String
        println("`t0` must be either a Date or a Float (decimal year).")
    else
        date0 = Date(DateFormats.yeardecimal.(t0));
    end

    if typeof(t_end) == Date
        date_end = t_end
    elseif typeof(t_end) == String
        if t_end == "today"
            date_end = today()
        else
            println("`t_end` not recognized")
        end
    else
        date_end = Date(DateFormats.yeardecimal.(t_end));
    end

    dates = collect.(range.(date0, date_end));
    n_dates = length(dates)
    
    timeline = DateFormats.yeardecimal.(dates)
    dates = Date.(DateFormats.yeardecimal.(timeline) .+ Hour.(ones(n_dates)*12))

    return timeline, dates
end