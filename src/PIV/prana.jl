using MAT, MATLAB, TOML, ImageIO, TiffImages, ImageCore, ImageTransformations, Statistics, Printf

"""
    config2toml(config::Dict, tomlpath)

Write a dictionary of prana config data to a TOML file.
"""
function config2toml(config, tomlpath)
    open(tomlpath, "w") do io
        TOML.print(io, config; sorted=true) do x
            x isa Matrix{Any} && return vec(string.(x))
            return x
        end
    end
end

"""
    toml2config(tomlpath)
"""
function toml2config(tomlpath)
    TOML.parsefile(tomlpath)
end

"""
    prana2config(pranapath)

Convert a prana config file saved in the .MAT format to a human-readable TOML file.
"""
function prana2config(pranapath, tomlpath)
    prana_cfg_m = MatFile(pranapath)
    prana_cfg = sortdict(get_variable(prana_m, "Data"))
    close(prana_cfg_m)

    prana_cfg
end

"""
    config2prana(config, pranapath)

Convert a dict of Prana configuration data to a prana config file saved in the .MAT format.
"""
function config2prana(config, pranapath)
    function reshapecfg!(cfg) # reshape all vectors to 1xn
        for (k, v) in cfg
            v isa Vector && (cfg[k] = reshape(v, 1, :))
            v isa Dict && reshapecfg!(v)
            v isa Number && (cfg[k] = string(v))
        end
    end
    reshapecfg!(config)

    prana_cfg_m = MatFile(pranapath, "w")
    put_variable(prana_cfg_m, "Data", config)
    close(prana_cfg_m)
end

# Need to assign several run-specific parameters to the prana config file 
# and allow some config parameters to be assigned as keyword arguments, but others
# can be left as default.
# For each run, we need to specify paths to the:
# - per-run data output folder
# - per-run processing folder (typically a temp folder)
# - Prana template file (PIV.toml)
#
# From specified inputs, need to do a few things:
# - background-subtract raw images and save in processing folder
# - save intensity-adjusted summary image as .png in data output folder
# - find region of interest from pixel size in calibration file
# - write config fields:
#   - exp_ROI
#   - exp_date
#   - imdirec
#   - imdirec2
#   - outdirec
#   - passes
#   - wrmag (pixel size in μm)
#   - wrsep (pulse separation in Hz)
#   - per-PIV-pass fields
# 
# After setting up PIV passes, copy PIV1 settings to PIV0, write config to TOML file in data output folder,
# and convert to .MAT file in processing folder. Return path to .MAT file.
"""
    setupprana(runmeta; settings=Dict{String,Any}();
        processingdir=mktempdir(),
        templatepath=srcdir("PIV", "PIV.toml"),
        )
"""
function setupprana(runmeta, settings=Dict{String,Any}();
    processingdir=mktempdir(),
    templatepath=srcdir("PIV", "PIV.toml"),
    )
    # Create output and processing directories if they don't exist
    outdir = runmeta.outdir
    !isdir(outdir) && mkpath(outdir)
    runmeta.outdir = outdir
    !isdir(processingdir) && mkpath(processingdir)

    # Load raw images and background-subtract
    LA_path = rawdatadir(runmeta.Date, runmeta.ID, runmeta.TSI_LA_path)
    LB_path = replace(LA_path, "LA" => "LB")
    BGA_path = datadir("PIV", "bg", runmeta.TSI_bg_path)
    BGB_path = replace(BGA_path, "LA" => "LB")
    LA, LB = load.((LA_path, LB_path))
    if isfile(BGA_path) && isfile(BGB_path)
        BGA, BGB = load.((BGA_path, BGB_path))
        clamp01!(Gray{Float32}.(LA) .- BGA)
        clamp01!(Gray{Float32}.(LB) .- BGB)
    else
        @warn("No background images found for run $(runname(runmeta)).")
    end

    # Save intensity-adjusted summary image
    summaryimg = overlapimages(imadjust(LA, qmax=0.999), imadjust(LB, qmax=0.999)) |> restrict |> rotl90
    save(joinpath(outdir, "particles.png"), summaryimg)

    # Save background-subtracted images to processing folder
    for (path, img) in zip(joinpath.(Ref(processingdir), ("LA000001.tif", "LA000002.tif")), (LA, LB))
        isfile(path) && rm(path)
        save(path, rotl90(Gray{N0f16}.(img)))
    end

    # calculate camera FOV size
    pixel_size = runmeta.dx
    fov_size_μm = size(LA) .* pixel_size
    fov_size_m = fov_size_μm ./ 1e6
    pulse_separation = runmeta.dt

    # Load template config file
    prana_cfg = toml2config(templatepath)
    prana_cfg["exp_date"] = string(runmeta.Date)
    prana_cfg["imdirec"] = processingdir
    prana_cfg["imdirec2"] = processingdir
    prana_cfg["outdirec"] = processingdir
    prana_cfg["wrmag"] = pixel_size
    prana_cfg["wrsep"] = pulse_separation
    prana_cfg["exp_ROI"] = @sprintf("%f,%f", fov_size_m...)

    mergewith!(prana_cfg, settings) do x, y
        x isa Dict && y isa Dict && return merge!(x, y)
        return y
    end
    prana_cfg["PIV0"] = prana_cfg["PIV1"]
    prana_cfg["passes"] = string(count(prana_cfg) do (k, v)
        k isa String && startswith(k, "PIV") && v isa Dict
    end - 1) # subtract 1 for PIV0
    
    # Save config file to TOML and .MAT formats
    prana_cfg_toml = joinpath(outdir, "PIV.toml")
    config2toml(prana_cfg, prana_cfg_toml)
    prana_cfg_mat = joinpath(processingdir, "PIV.mat")
    config2prana(prana_cfg, prana_cfg_mat)

    return prana_cfg_mat
end

function runprana(runmeta, prana_cfg_path, 
        prana_dir=raw"C:\Users\aames\Documents\MATLAB\prana-master")
    outdir = runmeta.outdir
    # Launch a MATLAB session and run Prana
    s = MSession()
    put_variable(s, :prana_cfg_path, prana_cfg_path)
    put_variable(s, :prana_dir, prana_dir)
    eval_string(s, """
        addpath(prana_dir)
        prana_cfg = load(prana_cfg_path)
        pranaPIVcode(prana_cfg.Data)
    """)
    close(s)

    # Move Prana output files to output directory
    prana_cfg_dir = dirname(prana_cfg_path)
    for f in readdir(prana_cfg_dir)
        if startswith(f, "PIVpass") && endswith(f, ".mat")
            mv(joinpath(prana_cfg_dir, f), joinpath(outdir, f), force=true)
        end
        if startswith(f, "ExpSummary")
            mv(joinpath(prana_cfg_dir, f), joinpath(outdir, f), force=true)
        end
    end
end