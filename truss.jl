
using LinearAlgebra, SparseArrays
using MKLSparse
using Plots
mutable struct PoleArgs
    E::Float64
    A::Float64
end
abstract type Abstract_Pole end
abstract type Abstract_Node end

mutable struct Node <: Abstract_Node
    id::Int64
    poles::Vector{Abstract_Pole}
    external_force::Vector{Tuple{Float64, Float64}}
    position::Tuple{Float64, Float64}
    constraint::Tuple{Bool, Bool}
end

mutable struct Pole <: Abstract_Pole
    id::Int64
    vertex::Tuple{Abstract_Node, Abstract_Node}
    args::PoleArgs
end

mutable struct System
    node_list::Vector{Abstract_Node}
    pole_list::Vector{Abstract_Pole}
    ε::Vector{Float64}
    K::SparseMatrixCSC{Float64,Int64}
    P::Vector{Float64}
    δ::Vector{Float64}
    C::SparseMatrixCSC{Float64,Int64}
    reduced_rows::Vector{Int64}
end

function System()
    return System(
        Vector{Abstract_Node}(),
        Vector{Abstract_Pole}(),
        Vector{Float64}(),
        spzeros(0,0),
        Vector{Float64}(),
        Vector{Float64}(),
        spzeros(0,0),
        Vector{Int64}(),
    )
end

function Node()
    return Node(0, Vector{Pole}(), Vector{Tuple{Float64, Float64}}(), (0.0, 0.0) , (false, false))
end

function Pole(node_1::Node, node_2::Node, args::PoleArgs)
    return Pole(-1,(node_1, node_2), args)
end

function Kᵉ(p::Pole)
    A = p.args.A
    E = p.args.E

    x_1 = p.vertex[1].position[1]
    y_1 = p.vertex[1].position[2]
    x_2 = p.vertex[2].position[1]
    y_2 = p.vertex[2].position[2]
    
    Δx = x_2-x_1
    Δy = y_2-y_1

    l = Δx^2 + Δy^2 |> sqrt
    θ = Δx/l |> acos |> Float64

    k = A*E/l
    K̂ᵉ = zeros(Float64, 4, 4)
    K̂ᵉ[1,1] = K̂ᵉ[3,3] = k
    K̂ᵉ[1,3] = K̂ᵉ[3,1] = -k
    T_unit = [
        cos(θ)  -sin(θ) ;
        sin(θ)  cos(θ)  ;
    ]
    T = zeros(Float64, 4, 4)
    T[1:2,1:2] = T[3:4,3:4] = T_unit
    Kᵉ = T * K̂ᵉ * transpose(T)
    return Kᵉ
end

function fill_K_with(K::SparseMatrixCSC{Float64,Int64}, p::Pole)
    Nᵢ = p.vertex[1]
    Nⱼ = p.vertex[2]
    A = p.args.A
    E = p.args.E

    i = Nᵢ.id
    j = Nⱼ.id

    xᵢ = Nᵢ.position[1]
    yᵢ = Nᵢ.position[2]
    xⱼ = Nⱼ.position[1]
    yⱼ = Nⱼ.position[2]
    
    Δx = xⱼ-xᵢ
    Δy = yⱼ-yᵢ

    l = sqrt(Δx^2 + Δy^2)

    k = A*E/l
    α = Δx/l
    β = Δy/l

    Kₑ = k*[α^2 α*β; α*β β^2] 
    K[2i-1:2i, 2i-1:2i] += Kₑ
    K[2j-1:2j, 2j-1:2j] += Kₑ
    K[2i-1:2i, 2j-1:2j] -= Kₑ
    K[2j-1:2j, 2i-1:2i] -= Kₑ
end

function node(sys::System, positions::Vector{Tuple{Float64, Float64}})
    index = 0
    for p in positions
        index += 1
        n = Node(index, Vector{Pole}(), Vector{Tuple{Float64, Float64}}(), p, (false, false))
        push!(sys.node_list, n)
    end
    return sys
end

function constraint(sys::System, node_id::Int64, x::Bool, y::Bool)
    sys.node_list[node_id].constraint = (x,y)
    return sys
end

function external_force(sys::System, node_id::Int64, f::Tuple{Float64, Float64})
    push!(sys.node_list[node_id].external_force, f)
    return sys
end

function pole(sys::System, node_index::Tuple{Int64, Int64}, args::PoleArgs)
    (id_1, id_2) = node_index
    n1 = sys.node_list[id_1]
    n2 = sys.node_list[id_2]
    p = Pole(n1, n2, args)
    push!(n1.poles, p)
    push!(n2.poles, p)
    push!(sys.pole_list, p)
    p.id = length(sys.pole_list)
    return sys
end

function pole(sys::System, node_indexs_list::Vector{Tuple{Int64, Int64}}, args::PoleArgs)
    for nodes_index in node_indexs_list
        pole(sys, nodes_index, args)
    end
    return sys
end

function generate_K(sys::System)
    N = length(sys.node_list)
    K = spzeros(2N,2N)
    for p in sys.pole_list
        fill_K_with(K, p)
    end
    sys.K = K
    return sys
end

function generate_P(sys::System)
    N = sys.node_list |> length
    P = zeros(2N)
    for n in sys.node_list
        px = 0
        py = 0
        for p in n.external_force
            px += p[1]
            py += p[2]
        end
        i = n.id
        P[2i-1] = px
        P[2i] = py
    end
    sys.P = P;
    return sys
end

function generate_C(sys::System)
    i = 1
    reduced_rows = Vector{Int64}()
    N = length(sys.node_list)
    for n in sys.node_list
        if n.constraint[1]
            push!(reduced_rows,i)
        end
        if n.constraint[2]
            push!(reduced_rows,i+1)
        end
        i += 2
    end
    println(reduced_rows)
    s = length(reduced_rows)
    C = spzeros(2N-s,2N)

    r = 1
    i_r = 1
    for c in 1:2N
        if i_r <= s && reduced_rows[i_r] != c
            C[r,c] = 1
            r += 1
            println(C)
            println("r=$(r)")
        else
            i_r += 1
            println("ir=$(i_r)")
        end
    end

    sys.C = C
    sys.reduced_rows = reduced_rows
    return sys
end;

function solve_δ(sys::System)
    δ̂ = (sys.C*sys.K*transpose(sys.C))\(sys.C*sys.P)
    N = length(sys.node_list)
    s = length(sys.reduced_rows)
    δ = zeros(2N)
    i_r = 1
    i = 1
    for r in 1:2N
        if i_r<=s && r!=sys.reduced_rows[i_r]
            δ[r] = δ̂[i]
            i += 1
        else
            i_r += 1
            δ[r] = 0
        end
    end
    sys.δ = δ
    sys.P = sys.K*sys.δ
    return sys
end

function sys_solve(sys::System)
    sys |> generate_K |> generate_P |> solve_delta
    return sys
end

function caculate_ε(sys::System)
    ε = zeros(length(sys.pole_list))
    index = 1
    for p in sys.pole_list
        Nᵢ = p.vertex[1]
        Nⱼ = p.vertex[2]
        A = p.args.A
        E = p.args.E
    
        i = Nᵢ.id
        j = Nⱼ.id
    
        xᵢ = Nᵢ.position[1]
        yᵢ = Nᵢ.position[2]
        xⱼ = Nⱼ.position[1]
        yⱼ = Nⱼ.position[2]

        uᵢ = sys.δ[2i-1]
        vᵢ = sys.δ[2i]
        uⱼ = sys.δ[2j-1]
        vⱼ = sys.δ[2j]

        Δx = xⱼ-xᵢ
        Δy = yⱼ-yᵢ
    
        l = sqrt(Δx^2 + Δy^2)
    
        α = Δx/l
        β = Δy/l

        ε[index] = α*(uⱼ-uᵢ)/l + β*(vⱼ-vᵢ)/l
        index += 1
    end
    sys.ε = ε
    return sys
end
E = 210e9;
A = 0.005;
args = PoleArgs(E,A); # 杆件参数
nodes = [
    (0.0,   0.0),
    (5.0,   7.0),
    (5.0,   0.0),
    (10.0,  7.0),
    (10.0,  0.0),
    (15.0,  0.0),
]; # 节点

poles = [
    (1,2),
    (1,3),
    (2,3),
    (2,4),
    (2,5),
    (3,5),
    (4,5),
    (4,6),
    (5,6),
]; # 杆件

truss = System(); # 创建系统
node(truss, nodes); # 添加节点 
pole(truss, poles, args); # 添加杆件
constraint(truss, 1, true, true); # 添加约束
constraint(truss, 6, true, true); # 添加约束
external_force(truss, 2, (20e3, 0.00)); # 添加外力
truss |> generate_K |> generate_P |> generate_C |> solve_δ |> caculate_ε; # 生成, 求解

δ_mm = 1000* truss.δ;
σ = E*truss.ε;