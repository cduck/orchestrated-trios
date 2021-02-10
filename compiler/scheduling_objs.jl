global cirq, feedback_gate, feedback_xor_gate
function _init_cirq()
    global cirq = pyimport("cirq")
    # Define some cirq gates
    py"""
    import cirq
    class Feedback(cirq.Gate):
        def __init__(self, size):
            self._size = size
        def _num_qubits_(self):
            return self._size
        def _circuit_diagram_info_(self, args):
            return ('F',) * args.known_qubit_count
    class FeedbackXor(cirq.Gate):
        def __init__(self, size):
            self._size = size
        def _num_qubits_(self):
            return self._size
        def _circuit_diagram_info_(self, args):
            return ('F^',) * args.known_qubit_count
    class Barrier(cirq.Gate):
        def __init__(self, size):
            self._size = size
        def _num_qubits_(self):
            return self._size
        def _circuit_diagram_info_(self, args):
            return ('|B',) * args.known_qubit_count
    def feedback_gate(*qubits):
        return Feedback(len(qubits))(*qubits)
    def feedback_xor_gate(*qubits):
        return FeedbackXor(len(qubits))(*qubits)
    def barrier_gate(*qubits):
        return Barrier(len(qubits))(*qubits)
    """
    global feedback_gate = py"feedback_gate"
    global feedback_xor_gate = py"feedback_xor_gate"
    global barrier_gate = py"barrier_gate"
end


struct HardwareGate{N}
    kind::Symbol
    qubits::Tuple{Vararg{Int, N}}
    is_comm::Bool
end
function HardwareGate(kind::Symbol, qubits::Int...; is_comm)
    # Expensive check
    #@assert(length(Set(qubits)) == length(qubits),
    #        "Duplicate qubits ($qubits)")
    HardwareGate{length(qubits)}(kind, qubits, is_comm)
end


struct ScheduledCircuit{QubitT}
    moments_rev::Vector{Vector{HardwareGate}}
    qubit_map::BiMap{QubitT, Int}
end
ScheduledCircuit(gate_list::Vector{HardwareGate},
                 used_hardware_qids,
                 qubit_map::BiMap{QubitT, Int} where QubitT) = (
    ScheduledCircuit([gate_list], used_hardware_qids, qubit_map))
function ScheduledCircuit(gate_lists::Vector{Vector{HardwareGate}},
                          used_hardware_qids,
                          qubit_map::BiMap{QubitT, Int} where QubitT)
    used = sizehint!(Set{Int}(used_hardware_qids),
                     length(used_hardware_qids))
    moments_rev = Vector{HardwareGate}[]
    for gate_list in Iterators.reverse(gate_lists)
        for gate in Iterators.reverse(gate_list)
            if any(qhi in used for qhi in gate.qubits)
                empty!(used)
                push!(moments_rev, HardwareGate[])
            end
            push!(moments_rev[end], gate)
            union!(used, gate.qubits)
        end
    end
    ScheduledCircuit(moments_rev, qubit_map)
end

Base.length(c_sched::ScheduledCircuit) = length(c_sched.moments_rev)
Base.iterate(c_sched::ScheduledCircuit) = (
    iterate(Iterators.reverse(c_sched.moments_rev)))
Base.iterate(c_sched::ScheduledCircuit, state) = (
    iterate(Iterators.reverse(c_sched.moments_rev), state))

function to_cirq(c_sched::ScheduledCircuit)
    map = rev(c_sched.qubit_map)
    c = cirq.Circuit()
    for moment in c_sched
        cirq_gates = sizehint!([], length(moment))
        for gate in moment
            qubits = [to_cirq(get(map, qhi, qhi))
                      for qhi in gate.qubits]
            g = if gate.kind == :measure
                cirq.measure(qubits...)
            elseif gate.kind == :cnot
                cirq.CNOT(qubits...)
            elseif gate.kind in (:gate2, :cz)
                cirq.CZ(qubits...)
            elseif gate.kind == :gate1
                cirq.I(qubits...)
            elseif gate.kind == :feedback
                feedback_gate(qubits...)
            elseif gate.kind == :feedback_xor
                feedback_xor_gate(qubits...)
            elseif gate.kind == :barrier
                barrier_gate(qubits...)
            else
                throw(ArgumentError("unknown gate kind: $(gate.kind)"))
            end::PyObject
            push!(cirq_gates, g)
        end
        c.append(cirq.Moment(cirq_gates))
    end
    c
end


struct AncChainWeights <: AbstractMatrix{Int}
    graph::AbstractGraph{Int}
    ancilla_ids::Set{Int}
end
Base.size(w::AncChainWeights) = (l=length(vertices(w.graph)); (l, l))
function Base.getindex(w::AncChainWeights, i, j)
    (i in w.ancilla_ids || j in w.ancilla_ids) && return 1
    return 1<<12
end


struct RemovedNodeDistMatrix{WeightT, WeightsT, NodeT} <: AbstractArray{WeightT,2}
    weights::WeightsT
    removed_nodes::Set{NodeT}
    missing_weight::WeightT
end
RemovedNodeDistMatrix(weights::WeightsT, removed_nodes::Set{NodeT},
                      missing_weight) where {WeightsT, NodeT} = (
    RemovedNodeDistMatrix{eltype(WeightsT), WeightsT, NodeT}(
        weights, removed_nodes, missing_weight))
Base.getindex(mat::RemovedNodeDistMatrix, i, j) = (
    (i in mat.removed_nodes || j in mat.removed_nodes
        ? mat.missing_weight
        : mat.weights[i, j]))
Base.size(mat::RemovedNodeDistMatrix) = size(mat.weights)
