module Connectivity

import Base: show
using LightGraphs
using MetaGraphs
using GraphPlot

using Tweaks.BiMaps, Tweaks.GraphUtil

using ..Circuits: GridQubit

export DeviceConnectivity, SimpleDeviceConnectivity,
    AncAssistDeviceConnectivity, device_bounds, device_center_qubit


abstract type DeviceConnectivity{QubitT} end

struct SimpleDeviceConnectivity{QubitT} <: DeviceConnectivity{QubitT}
    qubit_map::BiMap{QubitT, Int}
    graph::AbstractGraph{Int}
end

struct AncAssistDeviceConnectivity{QubitT} <: DeviceConnectivity{QubitT}
    qubit_map::BiMap{QubitT, Int}
    graph::AbstractGraph{Int}
    hw_ancilla_ids::Set{Int}
end

function Base.show(io, mime::MIME"image/svg+xml",
                   dev_con::DeviceConnectivity)
    qubits = [rev(dev_con.qubit_map)[i]
              for i in vertices(dev_con.graph)]
    nlabels = [string(q) for q in qubits]
    args = []
    if first(keys(dev_con.qubit_map)) isa GridQubit
        xs = [q.col for q in qubits]
        ys = [q.row for q in qubits]
        push!(args, xs, ys)
        nlabels = 0:length(qubits)-1
    end
    kwargs = Dict{Symbol, Any}()
    if dev_con isa AncAssistDeviceConnectivity
        nodefillc = [i in dev_con.hw_ancilla_ids ? "silver" : "turquoise"
                     for i in vertices(dev_con.graph)]
        kwargs[:nodefillc] = nodefillc
    end
    context = gplot(dev_con.graph, args...; nodelabel=nlabels, kwargs...)
    show(io, mime, context)
end

raw"""
# Example layout
```
o-o o
|X \|
o-o-o
 /|
o-o-o
```
"""
function DeviceConnectivity(layout::AbstractString; hw_ancilla=nothing,
                            ancilla_char='A')
    lines = split(layout, '\n')
    #line_arrs = [collect(line) for line in lines]
    line_width = maximum(length(line) for line in lines)
    char_grid = fill(' ', (length(lines), line_width))
    for (y, line) in enumerate(lines)
        char_grid[y, 1:length(line)] .= collect(line)
    end
    w, h = div(line_width+1, 2), div(length(lines)+1, 2)
    DeviceConnectivity(char_grid, hw_ancilla=hw_ancilla,
                       ancilla_char=ancilla_char)
end
function DeviceConnectivity(layout::Matrix{Char}; hw_ancilla=nothing,
                            ancilla_char='A')
    if ancilla_char in layout
        @assert(hw_ancilla === nothing,
                "Don't specify ancilla in both layout and hw_ancilla.")
        hw_ancilla = []
    end

    h, w = size(layout)
    h, w = div(h+1, 2), div(w+1, 2)
    num_qubits = count(c->!(c in " -|/\\X\n\r"), layout)
    map = BiMap{GridQubit, Int}()
    g = Graph{Int}(num_qubits)
    # Assign qubit ids
    for y in 1:h, x in 1:w
        c = layout[2y-1, 2x-1]
        c == ' ' && continue
        qh = GridQubit(y, x)
        map[GridQubit(y, x)] = length(map) + 1
        c == ancilla_char && push!(hw_ancilla, qh)
    end
    # Horizontal connections
    for y in 1:h, x in 1:w-1
        c = layout[2y-1, 2x]
        @assert(c in " -",
                "Unknown horizontal connection '$c' at ($(2y-1), $(2x))")
        c == ' ' && continue
        q1, q2 = GridQubit(y, x), GridQubit(y, x+1)
        add_edge!(g, Edge(map[q1], map[q2]))
    end
    # Vertical connections
    for y in 1:h-1, x in 1:w
        c = layout[2y, 2x-1]
        @assert(c in " |",
                "Unknown vertical connection '$c' at ($(2y), $(2x-1))")
        c == ' ' && continue
        q1, q2 = GridQubit(y, x), GridQubit(y+1, x)
        add_edge!(g, Edge(map[q1], map[q2]))
    end
    # Diagonal connections
    for y in 1:h-1, x in 1:w-1
        c = layout[2y, 2x]
        @assert(c in " X\\/",
                "Unknown diagonal connection '$c' at ($(2y), $(2x-1))")
        q1, q2 = GridQubit(y, x), GridQubit(y, x+1)
        q3, q4 = GridQubit(y+1, x), GridQubit(y+1, x+1)
        c in "\\X" && add_edge!(g, Edge(map[q1], map[q4]))
        c in "/X" && add_edge!(g, Edge(map[q2], map[q3]))
    end

    DeviceConnectivity(map, g, hw_ancilla=hw_ancilla)
end

function DeviceConnectivity(width::Int, height::Int; hw_ancilla=nothing,
                            periodic_ancilla::Union{Nothing, Int}=nothing,
                            anc_off_x::Int=0, anc_off_y::Int=0)
    layout = fill(' ', (height*2-1, width*2-1))
    layout[1:2:end, 1:2:end] .= 'o'
    layout[2:2:end, 1:2:end] .= '|'
    layout[1:2:end, 2:2:end] .= '-'
    if periodic_ancilla !== nothing
        p = periodic_ancilla
        layout[2p-1+2anc_off_y:2p:end, 2p-1+2anc_off_x:2p:end] .= 'A'
    end
    DeviceConnectivity(layout, hw_ancilla=hw_ancilla, ancilla_char='A')
end

function DeviceConnectivity(map::BiMap{T, Int} where T, g::AbstractGraph{Int};
                            hw_ancilla=nothing)
    if hw_ancilla !== nothing
        hw_ancilla_ids = Set{Int}(map[q] for q in hw_ancilla)
        AncAssistDeviceConnectivity(map, g, hw_ancilla_ids)
    else
        SimpleDeviceConnectivity(map, g)
    end
end

## Test create
#g = Graph{Int}(2)
#add_edge!(g, Edge(1, 2))
#dev_con = DeviceConnectivity(
#    BiMap(GridQubit(1, 2)=>1, GridQubit(2, 5)=>2),
#    g)
## Grid create
#dev_con = DeviceConnectivity(4, 3)
## Any layout
dev_con_demo = DeviceConnectivity(raw"""
o-o o
|X \|
o-A-o
 /|
o-o-o
""")


function device_bounds(dev_con::DeviceConnectivity{GridQubit})
    x_min, x_max = extrema(q.col for q in keys(dev_con.qubit_map))
    y_min, y_max = extrema(q.row for q in keys(dev_con.qubit_map))
    x_min, x_max, y_min, y_max
end

function device_center_qubit(dev_con::DeviceConnectivity{GridQubit};
                             ignore=Set{Int}())
    x_min, x_max, y_min, y_max = device_bounds(dev_con)
    center_trunc = GridQubit(div(y_min+y_max, 2), div(x_min+x_max, 2))
    if (mod(x_min+x_max, 2) == 0 && mod(y_min+y_max, 2) == 0
            && haskey(dev_con.qubit_map, center_trunc)
            && !(dev_con.qubit_map[center_trunc] in ignore))
        center_trunc
    else
        xc, yc = (x_min+x_max)/2-1e-10, (y_min+y_max)/2-1e-10
        argmin(q for (q, qhi) in pairs(dev_con.qubit_map)
                 if !(qhi in ignore)) do q
            (hypot(xc-q.col, yc-q.row), q.col, q.row)
        end
    end
end

# Test
@assert device_center_qubit(DeviceConnectivity(3, 9)) == GridQubit(5, 2)
@assert device_center_qubit(dev_con_demo) == GridQubit(2, 2)
@assert device_center_qubit(DeviceConnectivity(3, 10)) == GridQubit(5, 2)
@assert device_center_qubit(DeviceConnectivity(10, 2)) == GridQubit(1, 5)


function find_device_path(dev_con::DeviceConnectivity)
    dev_con.qubit
end

end
