using CairoMakie
using GLMakie
using GeometryBasics

# get fault vertices
function get_vertices(fault, proj="2d")

    vertices = Polygon[]
    n_features = size(fault["llh"])[1];

    if proj=="2d"
        for i = 1:n_features
            p1 = Point2f(fault["lon"][i,1], fault["lat"][i,1])
            p2 = Point2f(fault["lon"][i,2], fault["lat"][i,2])
            p3 = Point2f(fault["lon"][i,3], fault["lat"][i,3])
            append!(vertices,[Polygon(vec([p1,p2,p3]))])
        end
    else
        for i = 1:n_features
            p1 = Point2f(fault["lon"][i,1], fault["lat"][i,1], fault["height"][i,1])
            p2 = Point2f(fault["lon"][i,2], fault["lat"][i,2], fault["height"][i,2])
            p3 = Point2f(fault["lon"][i,3], fault["lat"][i,3], fault["height"][i,3])
            append!(vertices,[Polygon(vec([p1,p2,p3]))])
        end
    end
    return vertices
end


# get fault vertices for GeoMakie plot
function get_vertices_GeoMakie(fault)

    n_features = size(fault["llh"])[1];
    vertices = NaN .* zeros(n_features*4, 2)

    j = 0
    for i = 1:4:4*n_features
        j = j+1
        vertices[i+1,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
        vertices[i+2,:] = [fault["lon"][j, 2], fault["lat"][j, 2]]
        vertices[i+3,:] = [fault["lon"][j, 3], fault["lat"][j, 3]]
    end
    return vertices
end

function get_vertices_gmt(fault)
    n_features = size(fault["lon"])[1]
    vertices = NaN .* zeros(n_features*5, 2)

    j = 0
    for i = 1:5:5*n_features
        j = j+1
        vertices[i+1,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
        vertices[i+2,:] = [fault["lon"][j, 2], fault["lat"][j, 2]]
        vertices[i+3,:] = [fault["lon"][j, 3], fault["lat"][j, 3]]
        vertices[i+4,:] = [fault["lon"][j, 1], fault["lat"][j, 1]]
    end

    return string.(vertices)
end


function make_figure_map_lat_time(obs_var, fault, figsize, title)

    timeline = obs_var["timeline"][1,:]
    obs = obs_var["obs"]
    n_samples = length(timeline)
    color_exp = Int(minimum([ceil(log10(abs(minimum(obs)))),
                        ceil(log10(abs(maximum(obs))))]) - 1)
    color = obs ./ 10^color_exp
    colorrange = (-ceil(maximum(abs.(color))), ceil(maximum(abs.(color))))

    n_lats = n_lats_edges - 1;
    lat_min = minimum(minimum(fault["lat"]));
    lat_max = maximum(maximum(fault["lat"]));
    lats_edges = range(lat_min, stop=lat_max, length=n_lats_edges);
    lats = lats_edges[1:end-1] + (lats_edges[2:end]-lats_edges[1:end-1])/2;
    color_lats_min = zeros(n_lats,n_samples);
    color_lats_max = zeros(n_lats,n_samples);

    for i=1:n_lats
        ind_lats = fault["llh"][:,2] .> lats_edges[i] .&&
                   fault["llh"][:,2] .<= lats_edges[i+1];
        color_lats_min[i,:] = minimum(color[ind_lats,:], dims=1);
        color_lats_max[i,:] = maximum(color[ind_lats,:], dims=1);
    end
    color_lats_minmax = color_lats_max;
    inds = color_lats_max .< abs.(color_lats_min);
    color_lats_minmax[inds] = color_lats_min[inds];

    fig = Figure(size=figsize)

    #############
    ### TITLE ###
    #############
    titlelayout = GridLayout(fig[0, 1],
                             halign = :center,
                             tellwidth = false
                             )
    Label(titlelayout[1, 1],
          title,
          halign = :center,
          fontsize=36,
          font = "Times New Roman"
          )

    #########################
    ### MAP LATITUDE-TIME ###
    #########################
    xticks = collect(timeline[1]-1:2:timeline[end]+1)
    xticks = [round(Int, xtick) for xtick in xticks] # Ensure ticks are integers
    ax = Axis(fig[1, 1],
               xlabel="Time (yr)",
               ylabel="Latitude (°)",
               xaxisposition=:top,
               xlabelsize=24,
               ylabelsize=24,
               xgridwidth=0.1,
               ygridwidth=0.1,
               xticks=xticks,
               xticklabelsize=24,
               yticklabelsize=24
               )
    
    hm = heatmap!(ax, timeline, lats, color_lats_minmax',
            colormap=:bwr, colorrange=colorrange)
    translate!(hm, 0, 0, 100) # move above surface plot
    tr_map = Makie.scatter!(
        ax,
        tremors["timeline"], tremors["lat"],
        color=:black,
        markersize=1,
        alpha=1.0,
        )
    translate!(tr_map, 0, 0, 200) # move above surface plot
    
    ################
    ### COLORBAR ###
    ################
    label_text = "Slip Potency Rate"
    max_abs_value = maximum(abs.(colorrange))
    tick_step = ceil(Int, max_abs_value / 3)
    ticks = collect(-3*tick_step:tick_step:3*tick_step)
    ticks = [round(Int, tick) for tick in ticks]  # Ensure ticks are integers
    Colorbar(fig[2, 1],
             limits=colorrange,
             colormap=:bwr,
             flipaxis=false,
             label=L"%$(label_text) $(m^3/yr)$ $\times 10^{%$(color_exp)}$",
             vertical=false,
             labelsize=24,
             ticks=ticks,
             ticklabelsize=24
             )

    return fig
end