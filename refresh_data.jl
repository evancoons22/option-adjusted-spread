using HTTP, JSON3, CSV, Dates

const DATA_DIR = joinpath(@__DIR__, "data")

function refresh_series(id; start_date="2020-01-01")
    key = get(ENV, "FRED_API_KEY", "")
    isempty(key) && error("Set FRED_API_KEY first")

    response = HTTP.get(
        "https://api.stlouisfed.org/fred/series/observations";
        query = Dict(
            "api_key"           => key,
            "series_id"         => id,
            "file_type"         => "json",
            "observation_start" => start_date,
            "sort_order"        => "asc",
        ),
    )

    observations = JSON3.read(response.body).observations
    isempty(observations) && error("No observations for $id")

    rows = map(observations) do obs
        value = String(obs.value)
        (
            date = Date(String(obs.date)),
            # FRED reports percent: 4.25 becomes 0.0425.
            rate = value == "." ? missing : parse(Float64, value) / 100,
        )
    end

    mkpath(DATA_DIR)
    path = joinpath(DATA_DIR, "$id.csv")

    # Replace the previous file only after writing succeeds.
    CSV.write(path * ".tmp", rows)
    mv(path * ".tmp", path; force=true)

    println("Saved $id: $(length(rows)) observations")
end

# for id in [
#     "SOFR",
#     "DGS1MO", "DGS3MO", "DGS6MO",
#     "DGS1", "DGS2", "DGS3", "DGS5",
#     "DGS7", "DGS10", "DGS20", "DGS30",
# ]
#     refresh_series(id)
# end

refresh_series("SOFR")
