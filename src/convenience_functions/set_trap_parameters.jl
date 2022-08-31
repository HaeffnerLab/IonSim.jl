export global_beam!

"""
    global_beam!(T::Trap, laser::Laser)
Set `laser` to shine with full intensity on all ions in `Trap`.
"""
function global_beam!(T::Trap, laser::Laser)
    for n in eachindex(T.configuration.ions)
        push!(laser.pointing, (n, 1.0))
    end
end
