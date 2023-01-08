function rawdatadir(date, args...)
    basepath = raw"S:\users\shocktube\runs"
    date = split(string(date), '-')
    return joinpath(basepath, date..., args...)
end

TSIname(basename, i, AB) = i == "missing" ? missing : basename * lpad(i, 6, '0') * ".T000.D000.P000.H000.L$AB.TIF"

runname(runmeta) = string(runmeta.Date, '_', runmeta.ID)


"""
    validatepaths(runmeta)
Schema:
Date	ID	timings_path	ptrace_path	TSI_ID	TSI_idx	cine_ID	TSI_bg_path

Check that all files listed in the metadata actually exist and return true if so, false otherwise.
Check:
    - rawdatadir(runmeta.Date, runmeta.ID): base directory containing all data for this run
    - runmeta.timings_path: timings .tsv
    - runmeta.ptrace_path: pressure trace .lvm
    - runmeta.TSI_LA_path: first TSI image
    - runmeta.cine_ID*".cine": cine file
    - datadir("PIV", "bg", runmeta.TSI_bg_path): background image
"""
function validatepaths(runmeta)
    basedir = rawdatadir(runmeta.Date, runmeta.ID)
    isdir(basedir) || (@warn("Directory $basedir does not exist") && return false)
    for path in (runmeta.timings_path, 
                 runmeta.ptrace_path, 
                 runmeta.TSI_LA_path,
                 runmeta.cine_ID*".cine")
        (isnothing(path) || ismissing(path)) && continue
        if !isfile(joinpath(basedir, path))
            @warn("File $path does not exist for run $(runname(runmeta))")
            return false
        end
    end
    isfile(datadir("PIV", "bg", runmeta.TSI_bg_path)) || (@warn("Background image $(runmeta.TSI_bg_path) does not exist") && return false)
    return true
end