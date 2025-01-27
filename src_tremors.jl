using HTTP
using Dates
using GeoDataFrames
using GeoInterface
using DateFormats
using PolygonOps
using StaticArrays
include("src_plot.jl")
include("src_time.jl")

function update_tremors(dirs, t_end)

    files_tremors = readdir(dirs["dir_tremors"]);
    if length(files_tremors[end]) == 21
        date_last_tremor_file = Date(files_tremors[end][8:17], "yyyy-mm-dd");
        date_first_tremor_file2download = date_last_tremor_file + Day(1)
    elseif length(files_tremors[end]) == 18
        date_last_tremor_file = Date(files_tremors[end][8:14], "yyyy-mm");
        date_first_tremor_file2download = date_last_tremor_file + Month(1)
    end
    if typeof(t_end) == String
        if t_end == "today"
            date_last_epoch2download = today()
        elseif t_end == "yesterday"
            date_last_epoch2download = today() - Day(1)
        else
            println("`t_end` not recognized")
        end
    else
        date_last_epoch2download = Date(DateFormats.yeardecimal.(t_end));
    end
    
    n_days = Dates.value(date_last_epoch2download -
                            date_first_tremor_file2download)
    for i=0:ceil(n_days)
        date2download = string(date_first_tremor_file2download + Day(i));
        # date2download = string(date_last_tremor_file + Day(i));
        download_tremors_day(date2download, dirs)
    end
    
    return nothing
end

function download_tremors_day(date2download, dirs)

    print(date2download)
    tremor_url = 
            "https://tremorapi.pnsn.org/api/v3.0/events?starttime="*
            date2download*"T00:00:00&endtime="*date2download*"T23:59:59";
    try
        HTTP.download(tremor_url,
                    dirs["dir_tremors"]*"tremor_"*date2download*".shp";
                    update_period=Inf);
        println(" ")
    catch
        println(" No data")
    end
    
end



function download_tremors_past(year_in, year_end, dirs)

    for i=year_in:year_end
        if isleapyear(i)
            days_in_month = [31,29,31,30,31,30,31,31,30,31,30,31];
        else
            days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31];
        end
    
        for j=1:12
            j_str = string(lpad(j,2,"0"))
            println(string(i) * "-" *j_str)
            tremor_url = 
                    "https://tremorapi.pnsn.org/api/v3.0/events?" *
                    "starttime=" * string(i) * "-" * j_str * "-01T00:00:00&" *
                    "endtime=" * string(i) * "-" * j_str * "-" *
                    string(days_in_month[j]) * "T23:59:59"
            try
                HTTP.download(tremor_url,
                    dirs["dir_tremors"]*"tremor_"*string(i)*'-'*j_str*".shp";
                    update_period=Inf);
                println(" ")
            catch
                println(" No data")
            end
        end
    end
    return nothing
end


function load_tremors(dirs, options)
    
    t0 = options["t0"]
    t1 = options["t1"]
    origin = options["origin"]

    t0_decyear, t0_date = get_time_decyear_and_date(t0)
    t1_decyear, t1_date = get_time_decyear_and_date(t1)

    files = readdir(dirs["dir_tremors"]);
    tremors = Dict();

    tremors["timeline"] = [];
    tremors["lon"]      = [];
    tremors["lat"]      = [];
    tremors["depth"]    = [];
    tremors["mag"]      = [];
    n_files = length(files);
    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        n_files; desc = "Loading tremors", enabled = true
        )
    for i=1:n_files
        try
            if length(files[i]) == 21
                date_last_tremor_file = Date(files[i][8:17], "yyyy-mm-dd");
            elseif length(files[i]) == 18
                date_last_tremor_file = Date(files[i][8:14], "yyyy-mm");
            end
            t_tremor_file_decyear, t_tremor_file_date = 
                get_time_decyear_and_date(date_last_tremor_file)
            if t_tremor_file_decyear >= t0_decyear &&
                    t_tremor_file_decyear <= t1_decyear
                tremors_i = load_tremors_shp(dirs["dir_tremors"]*files[i]);
                tremors["timeline"] = [tremors["timeline"];
                                        tremors_i["timeline"]];
                tremors["lon"]      = [tremors["lon"]; tremors_i["lon"]];
                tremors["lat"]      = [tremors["lat"]; tremors_i["lat"]];
                tremors["depth"]    = [tremors["depth"]; tremors_i["depth"]];
                tremors["mag"]      = [tremors["mag"]; tremors_i["mag"]];
            end
        catch
            #println("Skipped.")
        end
        ProgressMeter.next!(progress)
    end

    if ~isnothing(origin)
        tremors = llh2localxyz_seismicity(tremors, origin);
    end

    println("Done")

    return tremors
end


function load_tremors_shp(file)
    #print("Loading tremors from "*file*"... ")
    tremors = Dict();
    try
        data = GeoDataFrames.read(file);
        
        n_tremors = length(data.time);
        tremors["timeline"] = zeros(1, n_tremors);
        tremors["lon"] = zeros(n_tremors);
        tremors["lat"] = zeros(n_tremors);
        tremors["depth"] = zeros(n_tremors);
        tremors["mag"] = NaN*zeros(n_tremors);

        a = [d[6:25] for d in data.time]
        tremors["timeline"] = DateFormats.yeardecimal.(
                                    DateTime.(a, "dd u yyyy HH:MM:SS"))
        ll = GeoInterface.coordinates.(data.geometry)
        tremors["lon"] = [l[1] for l in ll]
        tremors["lat"] = [l[2] for l in ll]
        tremors["depth"] = data.depth
        try 
            tremors["mag"] = data.magnitude
        catch
            tremors["mag"] = NaN*zeros(n_tremors)
        end
        #println(" Done.")
    catch
        #print("File not loaded: ")
    end
    return tremors
end

function select_tremors(tremors_input, timeline, fault)
    tremors_output = tremors_input

    tremors_output["timeline_R"] = timeline
    dates = Date.(DateFormats.yeardecimal.(timeline))
    
    dates_tremors = Date.(DateFormats.yeardecimal.(tremors_input["timeline"]))
    
    polygon = SVector.(hcat(fault["xE"], fault["xE"][:,1]),
                       hcat(fault["yN"], fault["yN"][:,1])
                       )
    
    n_dates = length(dates)
    n_patches = size(fault["xE"])[1]
    tremors_output["R"] = zeros(n_patches, n_dates) .* NaN
    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        n_dates; desc = "Finding tremors in patches", enabled = true
        )
    Threads.@threads for i=1:n_dates
        ind_dates_tremors = dates_tremors .== dates[i]

        points = vec(SVector.(tremors_input["xE"][ind_dates_tremors],
                              tremors_input["yN"][ind_dates_tremors])
                    )
        if ~isempty(points)
            for j in 1:n_patches
                tremors_output["R"][j,i] = sum([inpolygon(p, polygon[j,:];
                                in=true, on=false, out=false) for p in points])
                
            end
        end
        ProgressMeter.next!(progress)
    end
    ind_not_selected = tremors["timeline_R"] .< tremors["timeline"][1] .||
                            tremors["timeline_R"] .> tremors["timeline"][end]
    tremors["R"][:,ind_not_selected] .= NaN
    
    return tremors_output
end

# function select_tremors_space(tremors, fault)
#     n_patches = size(fault["xE"])[1]
#     n_tremors = length(tremors["timeline"])
#     polygon = SVector.(hcat(fault["xE"], fault["xE"][:,1]),
#                        hcat(fault["yN"], fault["yN"][:,1]))
#     points = vec(SVector.(tremors["xE"],tremors["yN"]))
#     ind_ll_tremors = zeros(n_patches, n_tremors)


#     # Construct a ParforProgressbar object
#     progress = ProgressMeter.Progress(
#         n_patches; desc = "Finding tremors in patches", enabled = true
#         )
#     Threads.@threads for j in 1:n_patches
#         ind_ll_tremors[j,:] = [inpolygon(p, polygon[j,:];
#                         in=true, on=false, out=false) for p in points]
#         ProgressMeter.next!(progress)
#     end
#     return ind_ll_tremors
# end