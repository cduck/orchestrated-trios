module Circuits

import Base: show, print
using Random
using PyCall
using LightGraphs
using MetaGraphs

using Tweaks.BiMaps, Tweaks.GraphUtil

export GridQubit, CircuitInteractions, to_cirq, make_simple_interaction_graph,
    shuffle_circuit_qubits

global cirq
function __init__()
    global cirq = pyimport("cirq")
end


struct GridQubit{T}
    row::T
    col::T
end

Base.print(io::IO, q::GridQubit) = print(io, "($(q.row),$(q.col))")


"""
Represents a sequence of two-qubit interactions.
"""
struct CircuitInteractions
    qubit_map::BiMap{Any, Int}
    interactions::Vector{Tuple{Vararg{Int}}}
end

"""
Create from the two-qubits gates in a Cirq circuit.
"""
function CircuitInteractions(circuit::PyCall.PyObject;
                             allow_toffoli::Bool=false,
                             decompose_toffoli::Bool=false,
                             include_single::Bool=false)
    qubit_map = BiMap(Dict{Any, Int}(
        q=>i for (i, q) in enumerate(sort!(collect(circuit.all_qubits())))
       ))
    qubit_i = 0
    interactions = Tuple{Vararg{Int}}[]
    for op in circuit.all_operations()
        qi_list = ((qubit_map[q] for q in op.qubits)...,)
        _push_interaction!(interactions, qi_list, allow_toffoli,
                          decompose_toffoli, include_single)
    end
    CircuitInteractions(qubit_map, interactions)
end

"""
Create from a circuit interactions `*.cint` file.
"""
function CircuitInteractions(fname::AbstractString;
                             allow_toffoli::Bool=false,
                             decompose_toffoli::Bool=false,
                             include_single::Bool=false)
    open(fname) do f
        for line in eachline(f)
            # Skip pre-header before "===" line
            line == "===" && break
        end
        n_given = nothing
        depth = nothing
        gates = nothing
        # Read header info
        for line in eachline(f)
            line == "---" && break
            startswith(line, "n=") && (n_given = parse(Int, line[3:end]))
            startswith(line, "depth=") && (depth = parse(Int, line[7:end]))
            startswith(line, "gates=") && (gates = parse(Int, line[7:end]))
        end
        # Read content
        n = 0
        interactions = Tuple{Vararg{Int}}[]
        gates !== nothing && sizehint!(interactions, gates)
        for line in eachline(f)
            (line == "---" || line == "===") && break
            length(line) == 0 && continue
            line[1] in "0123456789" || continue
            qs_list = split(line)
            qi_list = ((parse(Int, qs)+1 for qs in qs_list)...,)
            _push_interaction!(interactions, qi_list, allow_toffoli,
                              decompose_toffoli, include_single)
            n = max(n, maximum(qi_list))
        end
        if n > n_given
            throw(ErrorException(
                "File parse consistency error: Given n is smaller than actual n"
                * " ($_given < $n_g)"))
        end
        qubit_map = BiMap(Dict{Any, Int}(
            i=>i+1 for i in 0:n-1
        ))
        CircuitInteractions(qubit_map, interactions)
    end
end

function _push_interaction!(interactions, qi_list, allow_toffoli,
                            decompose_toffoli, include_single)
    length(qi_list) >= 2 - include_single || return interactions
    length(qi_list) <= 2 + (allow_toffoli || decompose_toffoli) || (
        throw(ArgumentError(
            "unsupported gate on more than two (or three) qubits")))
    if decompose_toffoli && length(qi_list) == 3
        qi1, qi2, qi3 = qi_list
        # TODO: Single qubit gates
        push!(interactions, (qi2, qi3))
        push!(interactions, (qi1, qi3))
        push!(interactions, (qi2, qi3))
        push!(interactions, (qi1, qi3))
        push!(interactions, (qi1, qi2))
        push!(interactions, (qi1, qi2))
    else
        push!(interactions, qi_list)
    end
    interactions
end

function Base.iterate(c_inter::CircuitInteractions)
    length(c_inter.interactions) >= 1 || return nothing
    c_inter.interactions[1], 2
end
function Base.iterate(c_inter::CircuitInteractions, i::Int)
    length(c_inter.interactions) >= i || return nothing
    c_inter.interactions[i], i+1
end
Base.length(c_inter::CircuitInteractions) = length(c_inter.interactions)

to_cirq(list::AbstractVector) = [to_cirq(x) for x in list]
to_cirq(list::Tuple) = ((to_cirq(x) for x in list)...,)
to_cirq(qubit::Real) = cirq.LineQubit(qubit)
to_cirq(qubit::Tuple{Real, Real}) = cirq.GridQubit(row=qubit[1], col=qubit[2])
to_cirq(qubit::GridQubit) = cirq.GridQubit(row=qubit.row, col=qubit.col)
function to_cirq(c_inter::CircuitInteractions)
    c = cirq.Circuit()
    for qi_list in c_inter
        length(qi_list) >= 1 || continue
        q_list = [rev(c_inter.qubit_map)[qi] for qi in qi_list]
        cq_list = [to_cirq(q) for q in q_list]
        op = cirq.Z(cq_list[end])
        for cq in cq_list[end-1:-1:1]
            op = op.controlled_by(cq)
        end
        c.append(op)
    end
    c
end


function make_simple_interaction_graph(c_inter::CircuitInteractions)
    interaction_graph = MetaGraph{Int, Float64}(length(c_inter.qubit_map))
    defaultweight!(interaction_graph, 0.0)
    for qi_list in c_inter
        length(qi_list) >= 2 || continue
        @assert(length(qi_list) <= 3,
                "interaction over >=3 qubits not implemented")
        w = length(qi_list) == 2 ? 1 : 2
        for (qi1, qi2) in Iterators.product(qi_list, qi_list)
            qi1 < qi2 || continue
            add_edge!(interaction_graph, Edge(qi1, qi2))
            set_prop!(interaction_graph, qi1, qi2, :weight,
                      get_weight(interaction_graph, Edge(qi1, qi2)) + w)
        end
    end
    interaction_graph
end


function shuffle_circuit_qubits(c::CircuitInteractions, seed)
    shuffle_circuit_qubits(c, MersenneTwister(seed))
end
function shuffle_circuit_qubits(c::CircuitInteractions,
                                rng::Union{Nothing, Random.AbstractRNG})
    perm = 1:length(c.qubit_map)
    rng !== nothing && (perm = shuffle(rng, 1:length(c.qubit_map)))
    map = BiMap(
        k => perm[v] for (k, v) in pairs(c.qubit_map)
    )
    interactions = [
        (rng === nothing ? ((perm[qi] for qi in qi_list)...,)
            : (shuffle!(rng, [perm[qi] for qi in qi_list])...,))
        for qi_list in c.interactions
    ]
    CircuitInteractions(map, interactions)
end


end
