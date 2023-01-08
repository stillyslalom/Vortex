using ImageIO, TiffImages, MATLAB, ImageCore, Statistics

include(srcdir("pathutils.jl"))
include(srcdir("imageutils.jl"))

function load_targets(tgtdata)
    TSI = rotl90(load(rawdatadir(tgtdata.Date, tgtdata.TSI_tgt_path)))
    Phantom = Gray.(rotl90(load(rawdatadir(tgtdata.Date, tgtdata.Phantom_tgt_path))))
    return (; TSI, Phantom)
end

function run_cpselect(targets, wd=mktempdir())
    cd(wd) do
        TSIpath = joinpath(wd, "TSI.png")
        Phantompath = joinpath(wd, "Phantom.png")

        save(TSIpath, imadjust(targets.TSI))
        save(Phantompath, imadjust(targets.Phantom))

        mat"[m, f] = cpselect($TSIpath, $Phantompath, 'Wait', true)"
        mat"trans = fitgeotrans(m, f, 'projective')"
        mat"$tform = trans.T"
        tform
    end
end