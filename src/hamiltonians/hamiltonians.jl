using LaTeXStrings

export chamber
export Hamiltonian
export displacement
export lambdickeorder
export timescale
export rwacutoff
export timedependenteta
export prettyprint_labframe_matrixelement

"""
This block in triple quotes is Joe's outline for what we might eventually want.
mutable struct Hamiltonian
    chamber::Chamber  # Pull user-defined params from here

    bare_hamiltonian
    interacting_hamiltonian
    rotatingframe
    rotatingframe_hamiltonian


    timescale::Real  # Scaling factor to be applied to all instances of time
    lambdickeorder::Union{Vector{Int}, Int}  # Only terms order η^lambdickeorder are retained
    rwacutoff::Real  # Time-dep terms with frequencies > rwacutoff are discarded
    displacement::String  # Approximation technique to be applied to displacement operators
    timedependenteta::Bool # Should η's reflect time-dependence of trap frequencies?
end

function Hamiltonian(chamber::Chamber, 
        bare_hamiltonian, interacting_hamiltonian, 
        rotating_frame, rotating_frame_hamiltonian
    )
    Hamiltonian(
        chamber,
        bare_hamiltonian, interacting_hamiltonian, 
        rotating_frame, rotating_frame_hamiltonian,
        1, 1, Inf, "truncated", false
    )
end
"""

mutable struct Hamiltonian
    chamber::Chamber
    timescale::Real
    lambdickeorder::Union{Vector{Int}, Int}
    rwacutoff::Real
    displacement::String
    timedependenteta::Bool
    labframe_nz_elements::Vector{Any}
    
    function Hamiltonian(chamber::Chamber; timescale::Real=1.0, lambdickeorder::Union{Vector{Int}, Int}=1,
            rwacutoff::Real=1e-6, displacement::String="truncated", timedependenteta::Bool=false)
        
        labframe_nz_elements = findmatrixelements(chamber, lambdickeorder)
        new(chamber, timescale, lambdickeorder, rwacutoff, displacement, timedependenteta, labframe_nz_elements)
    end
     
end

function chamber(h::Hamiltonian)
    return h.chamber
end

function labframe_nz_elements(h::Hamiltonian)
    return h.labframe_nz_elements
end

function timescale(h::Hamiltonian)
    return h.timescale
end

function rwacutoff(h::Hamiltonian)
    return h.rwacutoff
end

function timedependenteta(h::Hamiltonian)
    return h.timedependenteta
end

function displacement(h::Hamiltonian)
    return h.displacement
end

function lambdickeorder(h::Hamiltonian)
    return h.lambdickeorder
end

function prettyprint_labframe_matrixelement(h::Hamiltonian,(n,k,i,p,q))
    @assert n in range(1,length(ions(chamber(h)))) "ion index n = $n is out of bounds"
    @assert k in range(1,length(subleveltransitions(ions(chamber(h))[n]))) "transition index k = $k is out of bounds"
    @assert i in range(1,length(modes(chamber(h)))) "vibrational mode index i = $i is out of bounds"
    @assert p in range(0,modecutoff(modes(chamber(h))[i])) "initial vibrational state index p = $p is out of bounds"
    @assert q in range(0,modecutoff(modes(chamber(h))[i])) "final vibrational state index q = $q is out of bounds"
    s = findall(x->x[1]==(n,k,i,p,q),labframe_nz_elements(h))
    method = displacement(h)
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"

    if length(s)==0
        return latexstring(0)
    end
    texstr = ""
    #print(labframe_nz_elements(h)[s][1])
    for m in labframe_nz_elements(h)[s][1][2]
        newtexstr = L"\Omega_{%$n,%$k,%$m}(e^{i \nu_%$m t}+e^{-i \nu_%$m t})"
        Dstring = ""
        if method ≡ "analytic"
            Dstring = L"{\frac{%$p !}{%$q !}}(i\eta)^{%$(q-p)}e^{-\eta^2/2}L_%$p^{%$(q-p)}\left(\eta^2\right)"
        else
            for j in range(abs(q-p),step=2,lambdickeorder(h))
                if j == abs(q-p)
                    Dstring = latexstring("(i\\eta_{$n,$i,$p,$q,$m})^$j")
                else
                    Dstring = latexstring(Dstring,"+ (i\\eta_{$n,$i,$p,$q,$m})^$j")
                end
            end
        end
        newtexstr = latexstring(Dstring,newtexstr)
        if m != labframe_nz_elements(h)[s][1][2][end]
            newtexstr = latexstring(newtexstr,"+")
        end
        texstr = latexstring(texstr,newtexstr)
    end
    return texstr
    
end

function findmatrixelements(
        chamber,
        lambdickeorder
    )
    nonzero_indices = []
    mymodes = modes(chamber)
    myions = ions(chamber)
    mylasers = lasers(chamber)
    L = length(mymodes)
    N = length(myions)
    M = length(mylasers)
    νlist = Tuple([frequency(mode) for mode in mymodes])
    mode_dims = [modecutoff(mode) + 1 for mode in mymodes]
    lambdickeorder = _checklambdickeorder(lambdickeorder, L)
    for n in eachindex(myions), m in eachindex(mylasers)
        for k in eachindex(subleveltransitions(myions[n]))
            Ω0(chamber,1e-6,n,m,k)==0 && continue
            
            for (i,d) in enumerate(mode_dims)
                for q in 1:d, p in 1:q
                    q-p > lambdickeorder[i] && continue
                    this_element_index = findall(x->x[1]==(n,k,i,p,q),nonzero_indices)
                    if length(this_element_index)==0
                        push!(nonzero_indices,[(n,k,i,p,q),[m]])
                    else
                        push!(nonzero_indices[this_element_index][1][2],m)
                    end
                end
            end
        end
    end
    return nonzero_indices
end

function _checklambdickeorder(lambdickeorder, L)
    if typeof(lambdickeorder) <: Int
        return [lambdickeorder for _ in 1:L]
    else
        @assert(
            length(lambdickeorder) == L,
            "if typeof(lambdickeorder)<:Vector, then length of lambdickeorder must ",
            "equal number of modes"
        )
        reverse(lambdickeorder)
    end
end