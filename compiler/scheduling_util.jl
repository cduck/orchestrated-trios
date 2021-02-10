hardware_gates_for_swaps(path::AbstractVector{Int}, inv::Bool=false,
                         mutate_ancilla::AbstractSet{Int}=Set{Int}(),
                         mutate_map::Union{Nothing, BiMap}=nothing,
                         cancel_last_cnot::Bool=false) = (
    hardware_gates_for_swaps!(HardwareGate[], path, inv,
                              mutate_ancilla, mutate_map, cancel_last_cnot))
"""
    hardware_gates_for_swaps(gates, path[, inv[, mutate_ancilla[, mutate_map]]])

Append CNOT gates to `gates` that swap `path[1]` to `path[end]` (or in reverse
if inv=true).

Uses the 2-CNOT version if swapping with an ancilla listed in
`[mutate_]ancilla`.  If `mutate_map` is given, update it and `mutate_ancilla`
with the result of the swaps.
"""
function hardware_gates_for_swaps!(gates::AbstractVector{HardwareGate},
                                   path::AbstractVector{Int}, inv::Bool,
                                   mutate_ancilla::AbstractSet{Int},
                                   mutate_map::BiMap,
                                   cancel_last_cnot::Bool=false)
    inv && (path = reverse(path))
    itr = 1:length(path)-1
    did_append = false
    for i in itr
        qi1, qi2 = path[i], path[i+1]
        if qi2 in mutate_ancilla
            append!(gates, [
                HardwareGate(:cnot, qi1, qi2, is_comm=true)
                HardwareGate(:cnot, qi2, qi1, is_comm=true)
            ][cancel_last_cnot && !did_append ? (2:end) : (:)])
        else
            c1 = HardwareGate(:cnot, qi2, qi1, is_comm=true)
            c2 = HardwareGate(:cnot, qi1, qi2, is_comm=true)
            seq = inv ? [c2, c1, c2] : [c1, c2, c1]
            append!(gates, seq[cancel_last_cnot && !did_append ? (2:end) : (:)])
        end
        did_append = true
        # Mutate the qubit map
        q1 = pop!(rev(mutate_map), qi1)
        q2 = pop!(rev(mutate_map), qi2)
        mutate_map[q2] = qi1
        mutate_map[q1] = qi2
        # Mutate the ancilla set
        anc_is_1 = qi1 in mutate_ancilla
        anc_is_2 = qi2 in mutate_ancilla
        delete!(mutate_ancilla, qi1)
        delete!(mutate_ancilla, qi2)
        anc_is_1 && push!(mutate_ancilla, qi2)
        anc_is_2 && push!(mutate_ancilla, qi1)
    end
    if cancel_last_cnot && did_append && !inv
        pop!(gates)
    end
    gates
end
function hardware_gates_for_swaps!(gates::AbstractVector{HardwareGate},
                                   path::AbstractVector{Int}, inv::Bool=false,
                                   ancilla::AbstractSet{Int}=Set{Int}(),
                                   mutate_map::Nothing=nothing,
                                   cancel_last_cnot::Bool=false)
    itr = 1:length(path)-1
    inv && Iterators.reverse(itr)
    did_append = false
    for i in itr
        qi1, qi2 = path[i], path[i+1]
        if qi2 in ancilla
            append!(gates, [
                HardwareGate(:cnot, qi1, qi2, is_comm=true)
                HardwareGate(:cnot, qi2, qi1, is_comm=true)
            ][cancel_last_cnot && !did_append ? (2:end) : (:)])
        else
            append!(gates, [
                HardwareGate(:cnot, qi2, qi1, is_comm=true)
                HardwareGate(:cnot, qi1, qi2, is_comm=true)
                HardwareGate(:cnot, qi2, qi1, is_comm=true)
            ][cancel_last_cnot && !did_append ? (2:end) : (:)])
        end
        did_append = true
    end
    if cancel_last_cnot && did_append && !inv
        pop!(gates)
    end
    gates
end

function linear_ccz_gates!(gates::AbstractVector{HardwareGate},
                    qhi1::Int, qhi2::Int, qhi3::Int; is_comm::Bool=false)
    # TODO: Single qubit gates
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    gates
end
function complete_ccz_gates!(gates::AbstractVector{HardwareGate},
                             qhi1::Int, qhi2::Int, qhi3::Int;
                             is_comm::Bool=false)
    # TODO: Single qubit gates
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi2, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi3, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    push!(gates, HardwareGate(:cnot, qhi1, qhi2, is_comm=is_comm))
    gates
end
function either_ccz_gates!(dev_con::DeviceConnectivity,
                           gates::AbstractVector{HardwareGate},
                           qhi1::Int, qhi2::Int, qhi3::Int; is_comm::Bool=false)
    if has_edge(dev_con.graph, qhi1, qhi3)
        complete_ccz_gates!(gates, qhi1, qhi2, qhi3)
    else
        linear_ccz_gates!(gates, qhi1, qhi2, qhi3)
    end
end
