function vorticity(u, v, Δx=1, Δy=1)
    g₁, g₂ = parent.(Kernel.bickley())
    ImageFiltering.mapwindow(StructArray(; u, v), (3, 3)) do w
        mapreduce(+, eachindex(w)) do i
            w[i].v * g₁[i] / Δx - w[i].u * g₂[i] / Δy
        end
    end
end
