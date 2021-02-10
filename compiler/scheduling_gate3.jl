function hardware_gates_for_ccz! end
# Dispatch based on value of `mode`
function hardware_gates_for_ccz!(gates::AbstractVector{HardwareGate},
                                 dev_con::DeviceConnectivity,
                                 path1::Vector{Int}, path2::Vector{Int},
                                 mode::Symbol)
    @assert(length(path1) >= 2, "path 1 is too short")
    @assert(length(path2) >= 2, "path 2 is too short")
    # Orient paths so path1 ends where path2 ends
    if path1[end] == path2[1]
        path2 = reverse(path2)
    elseif path1[1] == path2[end]
        path1 = reverse(path1)
    elseif path1[1] == path2[1]
        path1 = reverse(path1)
        path2 = reverse(path2)
    end
    # Ensure path2 is the same or longer than path1
    if length(path1) > length(path2)
        path1, path2 = path2, path1
    end
    @assert(path1[end] == path2[end],
            "path 1 and path 2 must have a common end")
    @assert(intersect!(Set(path1), path2) == Set([path1[end]]),
            "path 1 and path 2 must not intersect")
    hardware_gates_for_ccz!(Val(mode), gates, dev_con, path1, path2)
end

# Implementation for mode :swap
function hardware_gates_for_ccz!(::Val{:swap},
                                 gates::AbstractVector{HardwareGate},
                                 dev_con::DeviceConnectivity,
                                 path1::Vector{Int}, path2::Vector{Int})
    hardware_gates_for_swaps!(gates, (@view path1[1:end-1]), false)
    hardware_gates_for_swaps!(gates, (@view path2[1:end-1]), false)
    either_ccz_gates!(dev_con, gates, path1[end-1], path1[end], path2[end-1])
    hardware_gates_for_swaps!(gates, (@view path1[1:end-1]), true)
    hardware_gates_for_swaps!(gates, (@view path2[1:end-1]), true)
    gates
end

# Implementation for mode :swap_one_way
function hardware_gates_for_ccz!(::Val{:swap_one_way},
                                 gates::AbstractVector{HardwareGate},
                                 dev_con::DeviceConnectivity,
                                 path1::Vector{Int}, path2::Vector{Int})
    hardware_gates_for_swaps!(gates, (@view path1[1:end-1]), false,
                              Set{Int}(), dev_con.qubit_map)
    hardware_gates_for_swaps!(gates, (@view path2[1:end-1]), false,
                              Set{Int}(), dev_con.qubit_map)
    either_ccz_gates!(dev_con, gates, path1[end-1], path1[end], path2[end-1])
    gates
end

