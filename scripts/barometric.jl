using DrWatson
@quickactivate "Vortex"

includet(srcdir("weather.jl"))

wx  = barometric_pressure_record(2022)