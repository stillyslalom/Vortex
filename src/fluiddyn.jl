function vorticity(u, v, Δx=1, Δy=1)
    g₁, g₂ = parent.(Kernel.bickley())
    ImageFiltering.mapwindow(StructArray(; u, v), (3, 3)) do w
        mapreduce(+, eachindex(w)) do i
            w[i].v * g₁[i] / Δx - w[i].u * g₂[i] / Δy
        end
    end
end

function vortex_integrals(ω, r, z)
    dr = step(r)
    dz = step(z)
    dΓ = ω .* dr * dz
    Γ = sum(dΓ)
    R = sum(dΓ .* r) / Γ
    Z = sum(dΓ .* z') / Γ
    R₂² = sum(dΓ .* r.^2) / Γ
    R_diff = R₂² - R^2
    a = R_diff < 0u"m^2" ? R/2 : sqrt(R_diff)
    (; Γ, R, Z, R₂², a)
end