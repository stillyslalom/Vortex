using JLD2
"Distance between two points"
dist2d(x, y) = sqrt(sum((x .- y).^2))
dist2d(x) = dist2d(x...)

function add_registrations!(runlist)
    calibrations = collect_results(datadir("calibrations"), verbose=false)
    calibrations.path = basename.(calibrations.path)
    rename!(calibrations, :path => :registration_path)
    leftjoin!(runlist, calibrations, on=:registration_path, matchmissing=:equal)
end
