module Scheduling

import Base: show, iterate, length, size, getindex
using LightGraphs
using PyCall

using Tweaks.BiMaps, Tweaks.GraphUtil

using ..Circuits, ..Connectivity, ..Paths
import ..Circuits.to_cirq

export HardwareGate, ScheduledCircuit, to_cirq, schedule_for_static_mapping

function __init__()
    _init_cirq()
end


include("scheduling_objs.jl")
include("scheduling_util.jl")
include("scheduling_gate2.jl")
include("scheduling_gate3.jl")


"""
schedule_for_static_mapping(c_inter, dev_con, circuit_to_device_map, mode)

Return a `ScheduledCircuit`.

Modes: `:swap`, swap_one_way, anc_swap, anc_swap_one_way, anc_assisted,
    const_every2
"""
function schedule_for_static_mapping(
        c_inter::CircuitInteractions,
        dev_con::DeviceConnectivity,
        circuit_to_device_map::BiMap{Int, Int},
        mode::Symbol,
        alt_mode::Symbol,
        dist_threshold::Int=0,
        use_anc::Bool=true)
    orig_qubit_map = copy(dev_con.qubit_map)
    if dev_con isa AncAssistDeviceConnectivity
        orig_hw_ancilla_ids = copy(dev_con.hw_ancilla_ids)
    end
    gate_list = HardwareGate[]
    real_mode = mode
    if dev_con isa AncAssistDeviceConnectivity && use_anc
        distmat = AncChainWeights(dev_con.graph, dev_con.hw_ancilla_ids)
    else
        distmat = weights(dev_con.graph)
    end
    for qi_list in c_inter
        qhi_list = [circuit_to_device_map[qi] for qi in qi_list]
        length(qhi_list) >= 2 || continue
        if length(qhi_list) == 2
            # Schedule CNOT, CZ
            qhi1, qhi2 = qhi_list
            path, _ = top_path(dev_con.graph, qhi1, qhi2, distmat)
            dist = length(path) - 1
            real_mode = dist < dist_threshold ? alt_mode : mode
            hardware_gates_for_interaction!(gate_list, dev_con, path, real_mode)
        elseif length(qhi_list) == 3
            # Schedule Toffoli, CCZ
            qhi1, qhi2, qhi3 = qhi_list
            best_cost = (typemax(Int),)
            best_paths = nothing
            for (qhi_a, qhi_b, qhi_c) in [
                    (qhi1, qhi2, qhi3), (qhi2, qhi3, qhi1), (qhi3, qhi1, qhi2),
                    (qhi3, qhi2, qhi1), (qhi1, qhi3, qhi2), (qhi2, qhi1, qhi3)]
                path1, dist1 = top_path(dev_con.graph, qhi_a, qhi_b, distmat)
                mat_removed = RemovedNodeDistMatrix(
                    distmat, Set(@view path1[1:end-1]), 1<<24)
                path2, dist2 = top_path(dev_con.graph, qhi_b, qhi_c,
                                        mat_removed)
                if path2 !== nothing
                    cost = (dist1 + dist2,
                            length(path1) + length(path2),
                            abs(length(path1) - length(path2)),
                            qhi_a, qhi_b, qhi_c)
                    if cost < best_cost
                        best_paths = path1, path2
                        best_cost = cost
                    end
                end
            end
            @assert(best_paths !== nothing, "no double path found for $qi_list")
            path1, path2 = best_paths
            dist1 = length(path1) - 1
            dist2 = length(path2) - 1
            real_mode = (dist1 < dist_threshold && dist2 < dist_threshold
                         ? alt_mode : mode)
            hardware_gates_for_ccz!(gate_list, dev_con, path1, path2, real_mode)
        else
            @assert(2 <= length(qhi_list) <= 3,
                    "gates on more than 3 qubits not implemented")
        end
    end
    copy!(dev_con.qubit_map, orig_qubit_map)  # Restore initial map
    if dev_con isa AncAssistDeviceConnectivity
        # Restore initial ancilla
        copy!(dev_con.hw_ancilla_ids, orig_hw_ancilla_ids)
    end
    ScheduledCircuit(gate_list, values(circuit_to_device_map),
                     dev_con.qubit_map)
end


end
