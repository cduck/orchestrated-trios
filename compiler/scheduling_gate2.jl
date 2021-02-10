"""
hardware_gates_for_interaction(dev_con, path, mode)

Modes: `:swap`, ...
"""
function hardware_gates_for_interaction end
# Dispatch based on value of `mode`
function hardware_gates_for_interaction(dev_con::DeviceConnectivity,
                                        path::AbstractVector{Int}, mode::Symbol
                                       )::AbstractVector{HardwareGate}
    hardware_gates_for_interaction!(HardwareGate[], dev_con, path, mode)
end
"""
hardware_gates_for_interaction!(gates, dev_con, path, mode)

Modes: `:swap`, ...
"""
function hardware_gates_for_interaction!(gates::AbstractVector{HardwareGate},
                                         dev_con::DeviceConnectivity,
                                         path::AbstractVector{Int},
                                         mode::Symbol
                                        )::AbstractVector{HardwareGate}
    @assert(length(path) >= 2, "path is too short")
    hardware_gates_for_interaction!(Val(mode), gates, dev_con, path)
end

# Implementation for mode :swap
_swap_warned = false
function hardware_gates_for_interaction!(::Val{:swap},
                                         gates::AbstractVector{HardwareGate},
                                         dev_con::DeviceConnectivity,
                                         path::AbstractVector{Int})
    global _swap_warned
    if dev_con isa AncAssistDeviceConnectivity && !_swap_warned
        @warn("Wasted qubits!  Using an ancilla assited mapping with regular "
              * "swaps.")
        _swap_warned = true
    end
    hardware_gates_for_swaps!(gates, (@view path[1:end-1]), false)
    push!(gates, HardwareGate(:gate2, path[end-1], path[end], is_comm=false))
    hardware_gates_for_swaps!(gates, (@view path[1:end-1]), true)
    gates
end

# Implementation for mode :swap_one_way
_swap_one_way_warned = false
"""
Mutates dev_con.qubit_map in the process!
"""
function hardware_gates_for_interaction!(::Val{:swap_one_way},
                                         gates::AbstractVector{HardwareGate},
                                         dev_con::DeviceConnectivity,
                                         path::AbstractVector{Int})
    global _swap_one_way_warned
    if dev_con isa AncAssistDeviceConnectivity && !_swap_one_way_warned
        @warn("Wasted qubits!  Using an ancilla assited mapping with regular "
              * "swaps.")
        _swap_one_way_warned = true
    end
    hardware_gates_for_swaps!(gates, (@view path[1:end-1]), false,
                              Set{Int}(), dev_con.qubit_map)
    push!(gates, HardwareGate(:gate2, path[end-1], path[end], is_comm=false))
    gates
end
