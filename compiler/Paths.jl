module Paths

using LightGraphs

export top_path, bottom_path


function top_path(graph, qi1, qi2, distmx=weights(graph))
    info = dijkstra_shortest_paths(graph, qi2, distmx, allpaths=true)
    path = [qi1]
    next = info.predecessors[path[end]]
    while length(next) > 0
        push!(path, minimum(next))
        next = info.predecessors[path[end]]
    end
    path[end] == qi2 ? (path, info.dists[qi1]) : (nothing, nothing)
end
function bottom_path(graph, qi1, qi2, distmx=weights(graph))
    info = dijkstra_shortest_paths(graph, qi2, distmx, allpaths=true)
    path = [qi1]
    next = info.predecessors[path[end]]
    while length(next) > 0
        push!(path, maximum(next))
        next = info.predecessors[path[end]]
    end
    path[end] == qi2 ? (path, info.dists[qi1]) : (nothing, nothing)
end


end

