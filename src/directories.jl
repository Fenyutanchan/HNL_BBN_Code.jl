# Copyright (c) 2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

export data_directory, external_data_directory, output_data_directory
export plot_directory

ensure_directory(directory::AbstractString) = isdir(directory) ? abspath(directory) :
    begin
        ispath(directory) && rm(directory; recursive=true, force=true)
        mkpath(directory; mode=0o755) |> abspath
    end

data_directory = joinpath(@__DIR__, "data") |> ensure_directory

external_data_directory = joinpath(data_directory, "ext") |> ensure_directory
output_data_directory = joinpath(data_directory, "out") |> ensure_directory

plot_directory = joinpath(@__DIR__, "plots") |> ensure_directory
