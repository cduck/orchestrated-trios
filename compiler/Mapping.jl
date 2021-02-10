module Mapping

import Base: show, getindex, size
using LightGraphs
using MetaGraphs

using Tweaks.BiMaps, Tweaks.GraphUtil
#using DrawSvg  # required by draw_mapping

using ..Circuits, ..Connectivity

export make_interaction_cost_matrix, map_to_hardware, draw_mapping


struct AncDistMatrix{GraphT} <: AbstractArray{Float64,2}
    dev_graph::GraphT
    hw_ancilla_ids::Set{Int}
end
AncDistMatrix(dev_con::AncAssistDeviceConnectivity) = (
    AncDistMatrix(dev_con.graph, dev_con.hw_ancilla_ids))
getindex(mat::AncDistMatrix, i, j)::Float64 = (
    i in mat.hw_ancilla_ids || j in mat.hw_ancilla_ids ? 2/3 : 1)
size(mat::AncDistMatrix) = (l=length(vertices(mat.dev_graph)); (l, l))

function make_interaction_cost_matrix(cost_func::Function,
                                      dev_con::DeviceConnectivity)
    if dev_con isa AncAssistDeviceConnectivity
        spaths = floyd_warshall_shortest_paths(dev_con.graph,
                                               AncDistMatrix(dev_con))
    else
        spaths = floyd_warshall_shortest_paths(dev_con.graph)
    end
    cost_func.(spaths.dists)
end

function largest_edges(graph)
    max_weight = maximum(get_weight(graph, e)
                         for e in edges(graph))
    [e for e in edges(graph)
       if get_weight(graph, e) == max_weight]
end

function largest_edge(graph)
    @assert length(edges(graph)) >= 1
    argmin(largest_edges(graph)) do e
        e.src, e.dst
    end
end

function largest_vertices(graph)
    v_weights = [
        (qi,
         length(neighbors(graph, qi)) > 0
         ? sum(get_weight(graph, qi, qi_n, 0) for qi_n in neighbors(graph, qi))
         : 0)
        for qi in vertices(graph)
    ]
    max_weight = maximum(w for (qi, w) in v_weights)
    [qi for (qi, w) in v_weights
        if w == max_weight]
end

function largest_vertex(graph)
    @assert length(edges(graph)) >= 1
    minimum(largest_vertices(graph))
end

function map_to_hardware(dev_con::DeviceConnectivity,
                         interaction_graph::AbstractMetaGraph,
                         interaction_cost_matrix::Matrix,
                         qids,
                         ; max_qubit_first=false)
    if max_qubit_first
        first_qubit_ids = [largest_vertex(interaction_graph)]
    else
        e = largest_edge(interaction_graph)
        first_qubit_ids = [e.src, e.dst]
    end
    if dev_con isa AncAssistDeviceConnectivity
        ignore_hw_ids = dev_con.hw_ancilla_ids
    else
        ignore_hw_ids = Set{Int}()
    end
    map_to_hardware(dev_con, ignore_hw_ids, interaction_graph,
                    interaction_cost_matrix, qids, first_qubit_ids)
end
function map_to_hardware(dev_con::DeviceConnectivity,
                         ignore_hw_ids::Set{Int},
                         interaction_graph::AbstractMetaGraph,
                         interaction_cost_matrix::Matrix,
                         qids,
                         first_qubit_ids)
    @assert(length(qids) <= length(dev_con.qubit_map) - length(ignore_hw_ids),
            "Too many qubit to map to device: "
            * "$(length(qids)) > "
            * "$(length(dev_con.qubit_map) - length(ignore_hw_ids))")
    frontier_qhi_set = Set{Int}()
    circuit_to_device_map = BiMap{Int, Int}()
    weights_to_placed = Dict{Int, Float64}()
    function place_qubit_id(qi)
        @assert(!haskey(circuit_to_device_map, qi),
                "Qubit already placed: $qi")
        # Pick a hardware qubit id to place the circuit qubit id
        if length(circuit_to_device_map) <= 0
            qh = device_center_qubit(dev_con, ignore=ignore_hw_ids)
            qhi = dev_con.qubit_map[qh]
            @assert haskey(dev_con.qubit_map, qh)
            @assert !hasval(circuit_to_device_map, qhi)
        else
            qhi = argmin(frontier_qhi_set) do qhi_f
                cost = sum(begin
                        dist_cost = interaction_cost_matrix[qhi_f, qhi_p]
                        w = get_weight(interaction_graph, qi, qi_p, 0)
                        w * dist_cost
                    end
                    for (qi_p, qhi_p) in pairs(circuit_to_device_map)
                )
            end
        end
        @assert !(qhi in ignore_hw_ids) "$(rev(dev_con.qubit_map)[qhi])"
        circuit_to_device_map[qi] = qhi
        # Update frontier
        delete!(frontier_qhi_set, qhi)
        potential_frontier = Set{Int}(neighbors(dev_con.graph, qhi))
        while length(potential_frontier) > 0
            qhi_n = pop!(potential_frontier)
            hasval(circuit_to_device_map, qhi_n) && continue
            if qhi_n in ignore_hw_ids
                union!(potential_frontier, neighbors(dev_con.graph, qhi_n))
                continue
            end
            push!(frontier_qhi_set, qhi_n)
        end
        # Update weights_to_placed
        delete!(weights_to_placed, qi)
        for qi2 in neighbors(interaction_graph, qi)
            haskey(circuit_to_device_map, qi2) && continue
            w = get_weight(interaction_graph, qi, qi2, 0)
            weights_to_placed[qi2] = get(weights_to_placed, qi2, 0) + w
        end
    end

    # Place starting qubits
    for qi in first_qubit_ids
        place_qubit_id(qi)
    end
    # Place the rest of the qubits
    for _ in length(first_qubit_ids)+1:length(qids)
        if length(weights_to_placed) > 0
            max_w = maximum(w for (qi, w) in weights_to_placed
                              if !haskey(circuit_to_device_map, qi))
            next_qi = minimum(qi for (qi, w) in weights_to_placed
                                 if (w >= max_w
                                     && !haskey(circuit_to_device_map, qi)))
        else
            next_qi = minimum(qi for qi in qids
                                 if !haskey(circuit_to_device_map, qi))
        end
        place_qubit_id(next_qi)
    end
    @assert length(circuit_to_device_map) == length(qids)

    circuit_to_device_map
end

function draw_mapping(dev_con::DeviceConnectivity,
                      interaction_graph::AbstractMetaGraph,
                      circuit_to_device_map::BiMap;
                      render_w=800, render_h=300)
    if dev_con isa AncAssistDeviceConnectivity
        hw_ancilla_ids = dev_con.hw_ancilla_ids
    else
        hw_ancilla_ids = Set{Int}()
    end

    x_min, x_max, y_min, y_max = device_bounds(dev_con)
    w, h = x_max-x_min, y_max-y_min
    d = draw.Drawing(w+1, h+1, origin=(x_min-0.5, -y_max-0.5))

    for (src, dst) in edges(dev_con.graph)
        q1 = rev(dev_con.qubit_map)[src]
        q2 = rev(dev_con.qubit_map)[dst]
        d.append(draw.Line(q1.col, -q1.row, q2.col, -q2.row,
                           stroke_width=0.05, stroke="#ccc"))
    end
    for (q, qhi) in pairs(dev_con.qubit_map)
        fill = qhi in hw_ancilla_ids ? "#ccc" : "#afa"
        d.append(draw.Circle(q.col, -q.row, 0.3, fill=fill))
        d.append(draw.Text([string(q)], 0.25, q.col, -q.row, center=true,
                           fill="#777"))
    end

    max_weight = maximum(get_weight(interaction_graph, e, 0)
                         for e in edges(interaction_graph))
    for (qi1, qi2) in edges(interaction_graph)
        qhi1 = circuit_to_device_map[qi1]
        qhi2 = circuit_to_device_map[qi2]
        qh1 = rev(dev_con.qubit_map)[qhi1]
        qh2 = rev(dev_con.qubit_map)[qhi2]
        line_len = hypot(qh1.col-qh2.col, qh1.row-qh2.row)
        width = 0.1/max_weight * get_weight(interaction_graph, qi1, qi2)
        stroke = "black"
        if qhi1 in hw_ancilla_ids || qhi2 in hw_ancilla_ids
            stroke = "red"
        end
        d.append(draw.Path(stroke_width=width, stroke=stroke, opacity=0.4,
                           fill="none"
            ).M(qh1.col, -qh1.row
            ).Q((qh1.col+qh2.col)/2+0.1*line_len,
                -(qh1.row+qh2.row)/2+0.2*line_len,
                qh2.col, -qh2.row))
    end

    d.setRenderSize(w=render_w, h=render_h)
end


end
