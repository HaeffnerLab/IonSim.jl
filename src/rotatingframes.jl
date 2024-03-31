using QuantumOptics
using SparseArrays
export RotatingFrame
export rotatingphase
export unitary

mutable struct RotatingFrame
    chamber::Chamber
    unitary::Vector{Float64}
    
    """
    In general, a change of basis to the Hamiltonian does not affect the underlying physics. Sometimes, such
    a change of basis is beneficial computationally because if the change of basis itself is a function of
    time, the frequency of temporal oscillations in the hamiltonian can be reduced (sometimes to zero). This
    greatly simplifies the problem and improves computational speed.
    
    A rotating frame transformation may be thought of in terms of either operators or states acquiring time-
    dependent phases. For experimentalists, it is often easiest to think of the phase as being attached to some
    operator which appears in the Hamiltonian. For ion-trap experiments, those are most commonly a particular
    projection operator (onto an atomic eigenstate) or a creation operator which acts on a vibrational mode.
    
    IonSim's RotatingFrame object keeps track of such time dependent phases that the user may choose to apply
    in order to reduce the time dependence in their Hamiltonian. When constructing a Hamiltonian, a 
    RotatingFrame can be passed for this purpose.
    
    The RotatingFrame object should be interacted with through the use of operators, in agreement with the 
    philosophy outlined above, which are easily constructed natively by IonSim. The constructor then generates
    a list of the phases which are applied to each eigenstate of the bare Hamiltonian and stores it in the
    an attribute called unitary.

    The user constructs a RotatingFrame by passing a vector of tuples.
    Each element of the vector has the form (operator, phase) where operator is an Operator type and phase is 
    a Float64. The Operator should be diagonal (either a projection operator of an atomoic state or the 
    number operator for a vibrational mode). For example, if the user wishes to apply the unitary U given by
    
    U(t) = exp(i * |S⟩⟨S| * ϕ₁t) × exp(i * a†a * ϕ₂t)
    
    then they should call
    
    RotatingFrame(chamber, [(|S⟩⟨S|, ϕ₁),(a†a, ϕ₂)]).
    
    Each of the operators needs to have the same dimensions as the Hamiltonian. This should be natively handled
    properly by all of IonSim's operator construction functions.
    """
    
    function RotatingFrame(chamber::Chamber,
        op_phase_tuples::Any
    )
        @assert typeof(op_phase_tuples) <: Vector "Argument op_phase_tuples must be of type Vector."
        hdim = length(basis(chamber))
        u = zeros(hdim)
        
        for (j,tup) in enumerate(op_phase_tuples)
            @assert typeof(tup)<:Tuple{Operator,Float64} "Elements of op_phase_tuples must be of type Tuple{Operator,Float64}."
            op = tup[1]
            phase = tup[2]
            @assert √(length(op))==hdim "Operator $j has incorrect dimensions."
            for k in 1:hdim
                u[k] = u[k] + phase*op.data[k,k]
            end
        end
                
        new(chamber, u)
    end
end

# Consider: T = X₁ ⊗ X₂ ⊗ ... ⊗ X_n (Xᵢ ∈ ℝ{dims[i]×dims[i]}), and indices:
# indxs[1], indxs[2], ..., indsx[N] = (i1, j1), (i2, j2), ..., (iN, jN).
# This function returns (k, l) such that: T[k, l] = X₁[i1, j1] * X₂[i2, j2] *...* X_N[iN, jN]
function get_kron_indxs(indxs::Vector{Tuple{Int64, Int64}}, dims::Vector{Int64})
    L = length(indxs)
    rowcol = Int64[0, 0]
    @assert indxs[L][1] <= dims[L] "indxs[$L][1] > dims[$L]"
    @assert indxs[L][2] <= dims[L] "indxs[$L][2] > dims[$L]"
    for i in 1:(L-1)
        @assert indxs[i][1] <= dims[i] "indxs[$i][1] > dims[$i]"
        @assert indxs[i][2] <= dims[i] "indxs[$i][2] > dims[$i]"
        rowcol .+= prod(view(dims,(i+1):L)) .* (indxs[i] .- 1)
    end
    rowcol .+= indxs[L]
    return rowcol
end

# The inverse of get_kron_indxs. If T = X₁ ⊗ X₂ ⊗ X₃ and X₁, X₂, X₃ are M×M, N×N and L×L
# dimension matrices, then we should input dims=(M, N, L).
"""function inv_get_kron_indxs(indxs, dims)
    row, col = indxs
    N = length(dims)
    ret_rows = Array{Int64}(undef, N)
    ret_cols = Array{Int64}(undef, N)
    for i in 1:N
        tensor_N = prod(dims[i:N])
        M = tensor_N ÷ dims[i]
        rowflag = false
        colflag = false
        for j in 1:dims[i]
            jM = j * M
            if !rowflag && row <= jM
                @inbounds ret_rows[i] = j
                row -= jM - M
                rowflag = true
            end
            if !colflag && col <= jM
                @inbounds ret_cols[i] = j
                col -= jM - M
                colflag = true
            end
            rowflag && colflag && break
        end
    end
    return Tuple(ret_rows), Tuple(ret_cols)
end"""

function rotatingphase(rf::RotatingFrame,sublvlindices::Vector{Int64},modeindices::Vector{Int64})
    @assert length(sublvlindices)≡length(ions(rf.chamber)) "Number of sublevel indices must match number of ions."
    @assert length(modeindices)≡length(modes(rf.chamber)) "Number of mode occupation numbers must match number of modes."
    ch = rf.chamber
    indxs = Vector{Tuple{Int64,Int64}}()
    dims = Vector{Int64}()
    
    for (i,ion) in enumerate((ions(ch)))
        push!(dims,shape(ion)[1])
    end

    for (m,mode) in enumerate(modes(ch))
        push!(dims,shape(mode)[1])
    end

    for sublvlidx in (sublvlindices)
        push!(indxs, (sublvlidx,sublvlidx))
    end    

    for modeidx in modeindices
        push!(indxs, (modeidx+1,modeidx+1))
    end

    dims = reverse(dims)
    indxs = reverse(indxs)
    
    kron_idxs = get_kron_indxs(indxs,dims)

    return unitary(rf)[kron_idxs[1]]
end

function rotatingphase(rf::RotatingFrame, σ::Any)
    @assert typeof(σ) <: Operator "Must provide a transition operator."
    (nzrows, nzcols, nzvals) = findnz(σ.data)
    @assert length(nzrows) ≥ 1 "Must provide a nonzero operator."
    return unitary(rf)[nzrows[1]] - unitary(rf)[nzcols[1]]
end

function rotatingphase(rf::RotatingFrame, ionidx::Int64, transition::Tuple{Tuple{String, Rational},Tuple{String, Rational}}, modeidx::Int64, Δphonons::Int64)
    ch = rf.chamber
    allions = ions(ch)
    allmodes = modes(ch)
    sl1 = 0
    sl2 = 0
    for (i,sublevel) in enumerate(sublevels(allions[ionidx]))
        if transition[1]==sublevel
            sl1 = i
        end
        if transition[2]==sublevel
            sl2 = i
        end
    end
    sl1vector = [1 for i in 1:length(allions)]
    sl1vector[ionidx] = sl1
    sl2vector = [1 for i in 1:length(allions)]
    sl2vector[ionidx] = sl2
    
    mode1vector = [0 for i in 1:length(allmodes)]
    mode2vector = [0 for i in 1:length(allmodes)]
    mode2vector[modeidx] = Δphonons
    
    return rotatingphase(rf, sl2vector, mode2vector) - rotatingphase(rf, sl1vector, mode1vector)
end

function rotatingphase(rf::RotatingFrame, ionidx::Int64, transition::Tuple{Tuple{String, Rational},Tuple{String, Rational}})
    return rotatingphase(rf, ionidx, transition, 1, 0)
end

function unitary(rf::RotatingFrame)
    return rf.unitary
end