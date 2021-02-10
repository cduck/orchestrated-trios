module Evaluation

import Base: show, +, //
using LightGraphs

using Tweaks.BiMaps, Tweaks.GraphUtil

using ..Circuits, ..Connectivity, ..Scheduling

export CircuitCost, ErrorModel, schedule_cost, compact_schedule_cost,
    total_error_probability, total_success_probability


struct CircuitCost
    num_gates1::Int
    num_gates2::Int
    num_comm_gates2::Int
    num_meas::Int
    num_feedback::Int
    num_feedback_xor::Int
    total_active_time::Float64
    time::Float64
    function CircuitCost(; num_gates1=0, num_gates2=0, num_comm_gates2=0,
                         num_meas=0, num_feedback=0, num_feedback_xor=0,
                         total_active_time=NaN, time=NaN)
        new(num_gates1, num_gates2, num_comm_gates2, num_meas, num_feedback,
            num_feedback_xor, total_active_time, time)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", cost::CircuitCost)
    print(io, """CircuitCost(
        num_gates1=$(cost.num_gates1),
        num_gates2=$(cost.num_gates2),
        num_comm_gates2=$(cost.num_comm_gates2),
        num_meas=$(cost.num_meas),
        num_feedback=$(cost.num_feedback),
        num_feedback_xor=$(cost.num_feedback_xor),
        total_active_time=$(cost.total_active_time),
        time=$(cost.time),
    )""")
end

function +(cost1::CircuitCost, cost2::CircuitCost)
    CircuitCost(
        num_gates1=cost1.num_gates1 + cost2.num_gates1,
        num_gates2=cost1.num_gates2 + cost2.num_gates2,
        num_comm_gates2=cost1.num_comm_gates2 + cost2.num_comm_gates2,
        num_meas=cost1.num_meas + cost2.num_meas,
        num_feedback=cost1.num_feedback + cost2.num_feedback,
        num_feedback_xor=cost1.num_feedback_xor + cost2.num_feedback_xor,
        total_active_time=cost1.total_active_time + cost2.total_active_time,
        time=cost1.time + cost2.time,
    )
end

function //(cost1::CircuitCost, cost2::CircuitCost)
    CircuitCost(
        num_gates1=cost1.num_gates1 + cost2.num_gates1,
        num_gates2=cost1.num_gates2 + cost2.num_gates2,
        num_comm_gates2=cost1.num_comm_gates2 + cost2.num_comm_gates2,
        num_meas=cost1.num_meas + cost2.num_meas,
        num_feedback=cost1.num_feedback + cost2.num_feedback,
        num_feedback_xor=cost1.num_feedback_xor + cost2.num_feedback_xor,
        total_active_time=cost1.total_active_time + cost2.total_active_time,
        time=max(cost1.time, cost2.time),
    )
end


struct ErrorModel
    gate1_err::Float64
    gate2_err::Float64
    meas_err::Float64
    dur2_idle_err::Float64

    gate1_time::Float64
    gate2_time::Float64
    meas_time::Float64
    feedback_time::Float64
    feedback_xor_time::Float64

    extra_idle_qubits::Int

    function ErrorModel(; gate1_err=0, gate2_err=0, meas_err=0, dur2_idle_err=0,
                        gate1_time=0, gate2_time=0, meas_time=0,
                        feedback_time=0, feedback_xor_time=0,
                        extra_idle_qubits=0)
        new(gate1_err, gate2_err, meas_err, dur2_idle_err, gate1_time,
            gate2_time, meas_time, feedback_time, feedback_xor_time,
            extra_idle_qubits)
    end
    function ErrorModel(old::ErrorModel; gate1_err=old.gate1_err,
                        gate2_err=old.gate2_err, meas_err=old.meas_err,
                        dur2_idle_err=old.dur2_idle_err,
                        gate1_time=old.gate1_time, gate2_time=old.gate2_time,
                        meas_time=old.meas_time,
                        feedback_time=old.feedback_time,
                        feedback_xor_time=old.feedback_xor_time,
                        extra_idle_qubits=old.extra_idle_qubits)
        new(gate1_err, gate2_err, meas_err, dur2_idle_err, gate1_time,
            gate2_time, meas_time, feedback_time, feedback_xor_time,
            extra_idle_qubits)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", model::ErrorModel)
    print(io, """ErrorModel(
        gate1_err=$(model.gate1_err),
        gate2_err=$(model.gate2_err),
        meas_err=$(model.meas_err),
        dur2_idle_err=$(model.dur2_idle_err),
        gate1_time=$(model.gate1_time),
        gate2_time=$(model.gate2_time),
        meas_time=$(model.meas_time),
        feedback_time=$(model.feedback_time),
        feedback_xor_time=$(model.feedback_xor_time),
        extra_idle_qubits=$(model.extra_idle_qubits),
    )""")
end


struct _GateErrorInfo
    gate_err::Float64
    time::Float64
    counter_i::Int
    starts_activity::Bool
    stops_activity::Bool
    count_is_per_qubit::Bool
end

const _gate1_i, _gate2_i, _meas_i, _feedback_i, _feedback_xor_i, _comm_i,
    _barrier_i = 1:7

function _gate_info_table(model::ErrorModel)
    Dict(
        :measure => _GateErrorInfo(model.meas_err, model.meas_time, _meas_i,
                            true, true, true),
        :cnot => _GateErrorInfo(model.gate2_err, model.gate2_time, _gate2_i,
                            true, false, false),
        :gate2 => _GateErrorInfo(model.gate2_err, model.gate2_time, _gate2_i,
                            true, false, false),
        :cz => _GateErrorInfo(model.gate2_err, model.gate2_time, _gate2_i,
                            true, false, false),
        :gate1 => _GateErrorInfo(model.gate1_err, model.gate1_time, _gate1_i,
                            true, false, false),
        :feedback => _GateErrorInfo(0, model.feedback_time, _feedback_i,
                            false, false, false),
        :feedback_xor => _GateErrorInfo(0, model.feedback_xor_time, _feedback_xor_i,
                            false, false, false),
        :barrier => _GateErrorInfo(0, 0, _barrier_i,
                            false, false, false)
    )
end

function schedule_cost(c_sched::ScheduledCircuit, model::ErrorModel)
    info_table = _gate_info_table(model)
    counters = repeat([0], 7)
    circuit_time = 0
    # Qid -> (Currently active, total time, last touch time)
    qubit_activity = Dict{Int, Tuple{Bool, Float64, Float64}}()
    for moment in c_sched
        moment_time = 0
        for gate in moment
            @assert length(gate.qubits) > 0
            @assert length(Set(gate.qubits)) == length(gate.qubits)
            info = info_table[gate.kind]
            # Update gate type counters
            counters[info.counter_i] += (info.count_is_per_qubit
                                         ? length(gate.qubits) : 1)
            if info.counter_i == _gate2_i && gate.is_comm
                counters[_comm_i] += (info.count_is_per_qubit
                                      ? length(gate.qubits) : 1)
            end
            # Calculate moment duration
            moment_time = max(moment_time, info.time)
            # Update qubit active time
            for qid in gate.qubits
                is_active, q_total, q_last = get(
                        qubit_activity, qid, (false, 0.0, 0.0))
                is_active && (q_total += circuit_time-q_last)
                qubit_activity[qid] = ((is_active || info.starts_activity)
                                            && !info.stops_activity,
                                       q_total + info.time,
                                       circuit_time + info.time)
            end
        end
        circuit_time += moment_time
    end
    total_active_time = if length(qubit_activity) > 0
        sum(q_total for (_, q_total, _) in values(qubit_activity))
    else
        0.0
    end::Float64
    CircuitCost(
        num_gates1=counters[_gate1_i],
        num_gates2=counters[_gate2_i],
        num_comm_gates2=counters[_comm_i],
        num_meas=counters[_meas_i],
        num_feedback=counters[_feedback_i],
        num_feedback_xor=counters[_feedback_xor_i],
        total_active_time=total_active_time,
        time = circuit_time,
    )
end


function compact_schedule_cost(c_sched::ScheduledCircuit, model::ErrorModel)
    info_table = _gate_info_table(model)
    counters = repeat([0], 7)
    circuit_time = 0
    # Qid -> (Currently active, total time, last touch time)
    qubit_activity = Dict{Int, Tuple{Float64, Float64, Float64}}()
    for moment in c_sched.moments_rev  # Reverse iteration
        for gate in moment
            (length(gate.qubits) == 0 && gate.kind == :feedback_xor) && continue
            @assert length(gate.qubits) > 0
            @assert length(Set(gate.qubits)) == length(gate.qubits)
            info = info_table[gate.kind]
            # Update gate type counters
            counters[info.counter_i] += (info.count_is_per_qubit
                                         ? length(gate.qubits) : 1)
            if info.counter_i == _gate2_i && gate.is_comm
                counters[_comm_i] += (info.count_is_per_qubit
                                      ? length(gate.qubits) : 1)
            end
            # Actually latest because it's flipped
            earliest_time = maximum(
                get(qubit_activity, q, (NaN, 0.0, 0.0))[3]
                for q in gate.qubits
            )
            # Update qubit active time
            for qid in gate.qubits
                last_active, q_total, q_last = get(
                        qubit_activity, qid, (NaN, 0.0, 0.0))
                if (isfinite(last_active) && info.starts_activity
                        && !info.stops_activity)
                    q_total += earliest_time-last_active
                end
                info.starts_activity && (last_active = earliest_time+info.time)
                qubit_activity[qid] = (last_active,
                                       q_total + info.time,
                                       earliest_time + info.time)
            end
        end
    end
    if length(c_sched) > 0
        circuit_time = maximum(end_time
                               for (_, _, end_time) in values(qubit_activity))
        total_active_time = sum(q_total
                                for (_, q_total, _) in values(qubit_activity))
    else
        circuit_time = 0.0
        total_active_time = 0.0
    end
    CircuitCost(
        num_gates1=counters[_gate1_i],
        num_gates2=counters[_gate2_i],
        num_comm_gates2=counters[_comm_i],
        num_meas=counters[_meas_i],
        num_feedback=counters[_feedback_i],
        num_feedback_xor=counters[_feedback_xor_i],
        total_active_time=total_active_time,
        time = circuit_time,
    )
end


total_error_probability(model::ErrorModel, cost::CircuitCost) = (
    1 - total_success_probability(model, cost))
function total_success_probability(model::ErrorModel, cost::CircuitCost)
    return prod([
        (1-model.gate1_err) ^ cost.num_gates1,
        (1-model.gate2_err) ^ cost.num_gates2,
        (1-model.meas_err) ^ cost.num_meas,
        (1-model.dur2_idle_err) ^ (
            cost.total_active_time != 0
            ? cost.total_active_time/model.gate2_time
            : 0),
        (1-model.dur2_idle_err) ^ (
            cost.time != 0
            ? model.extra_idle_qubits*cost.time/model.gate2_time
            : 0),
    ])
end

end
