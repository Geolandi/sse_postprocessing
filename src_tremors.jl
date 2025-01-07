using HTTP
using Dates
using GeoDataFrames
using GeoInterface
using DateFormats
using PolygonOps
include("src_plot.jl")

function update_tremors(dirs, t_end)

    files_tremors = readdir(dirs["dir_tremors"]);
    if length(files_tremors[end]) == 21
        date_last_tremor_file = Date(files_tremors[end][8:17], "yyyy-mm-dd");
    elseif length(files_tremors[end]) == 18
        date_last_tremor_file = Date(files_tremors[end][8:17], "yyyy-mm-dd");
    end
    date_last_epoch2download = Date(DateFormats.yeardecimal.(t_end));
    
    n_days = Dates.value(date_last_epoch2download - date_last_tremor_file);
    # n_days = daysdif(date_last_tremor_file,date_last_epoch2download);
    println(date_last_tremor_file)
    for i=0:ceil(n_days)
        date2download = string(date_last_tremor_file + Day(i));
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



# function download_tremors_past(year_in, year_end, dirs)

#     for i=year_in:year_end
#         if isleap(i)
#             days_in_month = [31,29,31,30,31,30,31,31,30,31,30,31];
#         else
#             days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31];
#         end
    
#         for j=1:12
#             if j<10
#                 j_str = string(sprintf('0%d',j));
#             else
#                 j_str = string(sprintf('%d',j));
#             end
#             fprintf('%d-%s',i,j_str)
#             tremor_url = append(...
#                     "https://tremorapi.pnsn.org/api/v3.0/events?",...
#                     "starttime=",string(i),"-",j_str,"-01T00:00:00&",...
#                     "endtime=",string(i),"-",j_str,"-",...
#                     string(days_in_month(j)),"T23:59:59");
#             try
#                 websave(append(dirs.dir_data, dirs.dir_case,'PNSN/tremor_',...
#                     string(i),'-',j_str,'.shp'), tremor_url);
#                 fprintf("\n")
#             catch
#                 fprintf(' No data\n')
#             end
#         end
#     end


function load_tremors(dirs, origin=nothing)
    
    files = readdir(dirs["dir_tremors"]);
    tremors = Dict();

    tremors["timeline"] = [];
    tremors["lon"]      = [];
    tremors["lat"]      = [];
    tremors["depth"]    = [];
    tremors["mag"]      = [];
    n_files = length(files);
    for i=1:n_files
        try
            tremors_i = load_tremors_shp(dirs["dir_tremors"]*files[i]);
            tremors["timeline"] = [tremors["timeline"]; tremors_i["timeline"]];
            tremors["lon"]      = [tremors["lon"]; tremors_i["lon"]];
            tremors["lat"]      = [tremors["lat"]; tremors_i["lat"]];
            tremors["depth"]    = [tremors["depth"]; tremors_i["depth"]];
            tremors["mag"]      = [tremors["mag"]; tremors_i["mag"]];
        catch
            println("Skipped.")
        end
    end
    
    if ~isnothing(origin)
        tremors = llh2localxyz_seismicity(tremors, origin);
    end
    return tremors
end


function load_tremors_shp(file)
    print("Loading tremors from "*file*"... ")
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
        println(" Done.")
    catch
        print("File not loaded: ")
    end
    return tremors
end

function select_tremors(tremors_input, timeline, fault)
    tremors_output = tremors_input

    dates = Date.(DateFormats.yeardecimal.(timeline))[1,:]
    
    dates_tremors = Date.(DateFormats.yeardecimal.(tremors_input["timeline"]))
    
    polygon = SVector.(hcat(fault["xE"], fault["xE"][:,1]),
                       hcat(fault["yN"], fault["yN"][:,1])
                       )
    
    n_dates = length(dates)
    n_patches = size(fault["xE"])[1]
    tremors_output["N"] = zeros(n_patches, n_dates)
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
                tremors_output["N"][j,i] = sum([inpolygon(p, polygon[j,:];
                                in=true, on=false, out=false) for p in points])
                
            end
        end
        ProgressMeter.next!(progress)
    end
    return tremors_output
end

function select_tremors_space(tremors, fault)
    n_patches = size(fault["xE"])[1]
    n_tremors = length(tremors["timeline"])
    polygon = SVector.(hcat(fault["xE"], fault["xE"][:,1]),
                       hcat(fault["yN"], fault["yN"][:,1]))
    points = vec(SVector.(tremors["xE"],tremors["yN"]))
    ind_ll_tremors = zeros(n_patches, n_tremors)


    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        n_patches; desc = "Finding tremors in patches", enabled = true
        )
    Threads.@threads for j in 1:n_patches
        ind_ll_tremors[j,:] = [inpolygon(p, polygon[j,:];
                        in=true, on=false, out=false) for p in points]
        ProgressMeter.next!(progress)
    end
    return ind_ll_tremors
end