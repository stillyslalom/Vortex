using Unitful, Dates, CSV, DataFramesMeta

"""
    parse_MA1(s)

Parse the MA1 field of the NOAA data.
Format: [sea level pressure, quality code, barometric pressure, quality code]
Quality code.
0 = Passed gross limits check
1 = Passed all quality control checks
2 = Suspect
3 = Erroneous
4 = Passed gross limits check, data originate from an NCEI data source
5 = Passed all quality control checks, data originate from an NCEI data source
6 = Suspect, data originate from an NCEI data source
7 = Erroneous, data originate from an NCEI data source
M = Manual change made to value based on information provided by NWS or FAA
9 = Passed gross limits check if element is present
"""
function parse_MA1(s)
    SLP, q_SLP, BP, q_BP = split(s, ',')
    return parse(Int, SLP)*10u"Pa", parse(Int, BP)*10u"Pa", q_SLP, q_BP
end

## Extract barometric pressure from NOAA data
function barometric_pressure_record(year)
    url = "https://www.ncei.noaa.gov/data/local-climatological-data/access/$year/72641014837.csv"
    path = joinpath(tempdir(), "$year.csv")
    download(url, path)

    wx = CSV.read(path, DataFrame)
    filter!(:HourlyStationPressure => p -> !ismissing(p) && p != "", wx)
    select!(wx, [:DATE, :HourlyStationPressure])
    wx.P_inchHg = tryparse.(Float64, wx.HourlyStationPressure)
    filter!(:P_inchHg => !isnothing, wx)
    wx.P = wx.P_inchHg*3386.39u"Pa"

    select!(wx, [:DATE, :P])
    combine(groupby(wx, :DATE), :P => mean => :P)

    # select!(wx, [:DATE, :MA1])
    # filter!(:MA1 => !ismissing, wx)
    # transform!(wx, :MA1 => ByRow(parse_MA1) => [:SLP, :BP, :q_SLP, :q_BP])
    # filter!(:q_BP => !=("9"), wx)
    # wx.DATE = round.(wx.DATE, Dates.Hour)
    # combine(groupby(wx, :DATE), :BP => mean => :BP)
end

linear_interpolation(Float64.(Dates.value.(wx.DATE)), ustrip.(wx.BP))