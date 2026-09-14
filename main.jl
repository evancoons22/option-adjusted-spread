using CSV
using DataFrames
using Distributions
using Dates
using Plots

RFR = 0.05
norm = Normal()
# p_density = pdf(d,0.0)
#

df = CSV.read("./data/UnderlyingOptionsEODCalcs_2023-08-25_cgi_or_historical.csv", DataFrame)
# df = filter(row -> row.underlying_symbol == "SPY", df)
# df = first(df, 1)
# print(describe(df))
# print(df)

df = filter(df) do row
    days = Dates.value(Date(row.expiration) - Date(row.quote_date))
    spot = (row.underlying_bid_1545 + row.underlying_ask_1545) / 2
    row.underlying_symbol == "SPY" &&  # filter for near the money options so vega isn't close to 0
    abs(row.strike / spot - 1) < 0.02 &&
    days >= 30 && 
    row.option_type == "C" &&
    row.bid_1545 > 0 && 
    row.ask_1545 > row.bid_1545
end
# df = first(df, 1)

function Newton(f, df, x_0, tolerance, iters)
    x = x_0
    for _ in 1:iters
        residual = f(x)
        if abs(residual) < tolerance
            return x
        end

        slope = df(x)

        if !isfinite(slope) || slope <= 0
            error("vega is 0 or invalid - cannot continue")
        end

        x = x - residual / slope

        if !isfinite(x) || x <= 0
            error("Newton produced invalid volatility (negative) - try a bracketed solver.")
        end

    end
    if abs(f(x)) < tolerance
        return x
    end
    error("did not converge within $iters iterations.")
end

function C(S_t, t, T, K, sigma, r)
    d_plus = 1 / (sigma * sqrt(T - t)) * (log(S_t / K) + (r + sigma^2 / 2) * (T - t))
    d_minus = d_plus - sigma * sqrt(T - t)
    return cdf(norm, d_plus) * S_t - cdf(norm, d_minus) * K * exp((-1) * r * (T - t))
end

function vega(S_t, t, T, K, sigma, r)
    d_plus = 1 / (sigma * sqrt(T - t)) * (log(S_t / K) + (r + sigma^2 / 2) * (T - t))
    return S_t * pdf(norm, d_plus) * sqrt(T - t)
end

# set up optimization function
# S = (df[1, "underlying_bid_1545"] + df[1, "underlying_ask_1545"]) / 2
# t = 0
# T = Dates.value( Date(df[1, "expiration"]) - Date(df[1, "quote_date"])) / 365
# K = df[1, "strike"] 
# r = RFR
# print("S: ", S, "\nt: ", t, "\nT: ", T, "\nK: ", K, "\nr: ", r, "\n")
# call_price = (df[1, "bid_1545"] + df[1, "ask_1545"]) / 2
# C_test = sigma -> C(S, t, T, K, sigma, r) - call_price
# 
# # set up df
# vega_test = sigma -> vega(S, t, T, K, sigma, r)
# 
# # newton method
# impl_vol = Newton(C_test, vega_test, 0.3, 0.005, 100)
# print(impl_vol)

# Iterate through 

points = Dict()
S = (df[1, "underlying_bid_1545"] + df[1, "underlying_ask_1545"]) / 2
t = 0
r = RFR
rows = nrow(df)
n_invalid = 0
print("\nCalculating volatility surface for strike price $S")
for i = 1:rows
    T = Dates.value( Date(df[i, "expiration"]) - Date(df[i, "quote_date"])) / 365
    K = df[i, "strike"] 
    call_price = (df[i, "bid_1545"] + df[i, "ask_1545"]) / 2
    C_test = sigma -> C(S, t, T, K, sigma, r) - call_price
    vega_test = sigma -> vega(S, t, T, K, sigma, r)

    try
        impl_vol = Newton(C_test, vega_test, 0.3, 0.005, 100)
        points[(T, K)] = impl_vol
    catch err
        # Preserve Ctrl+C so you can still stop the program.
        # err isa InterruptException && rethrow()
        global n_invalid += 1
        # @warn "Newton failed; skipping row" i T K reason=sprint(showerror, err)
        continue
    end

    if i % 20 == 0
        print("\niter $i / $rows")
    end
end
print("\nSkipped $n_invalid invalid results")

Ks = sort(unique([K for ((T, K), vol) in points]))
Ts = sort(unique([T for ((T, K), vol) in points]))

# Rows correspond to expirations; columns correspond to strikes.
# Z = [100 * points[(T, K)] for T in Ts, K in Ks]
# Z = [100 * get(points, (T, K), NaN) for T in Ts, K in Ks]
Z = fill(NaN, length(Ts), length(Ks))

for (i, expiry) in enumerate(Ts)
    # observed = sort([(K, vol) for ((T,K), vol) in points if T == expiry], by = first)
    observed = sort([ (K, vol) for ((T, K), vol) in points if T == expiry && !isnothing(vol) && isfinite(vol) && vol > 0 ]; by=first)
    strikes = first.(observed)
    sigmas = last.(observed)

    for (j, K) in enumerate(Ks)
        idx = searchsortedfirst(strikes, K)

        # if there is a value
        if idx <= length(strikes) && strikes[idx] == K && !isnan(sigmas[idx])
            Z[i,j] = 100 * sigmas[idx]
        # if interp needed
        elseif 1 < idx <= length(strikes)
            left, right = idx - 1, idx
            weight = (K - strikes[left]) / (strikes[right] - strikes[left])
            sigma = (1 - weight) * sigmas[left] + weight * sigmas[right]
            Z[i, j] = 100 * sigma 
        end

    end
end
    


# Z = filter(isfinite, Z)
valid_z = filter(isfinite, vec(Z))
isempty(valid_z) && error("No finite volatility values in Z.")

zmin =  minimum(valid_z) * 0.8
zmax =  maximum(valid_z) * 1.2 
print("\nzmin $zmin zmax $zmax")
figure = surface(
    Ks, Ts, Z;
    xlabel="Strike",
    ylabel="Time to expiry (years)",
    zlabel="Implied volatility (%)",
    color=:viridis,
    zlims=(zmin,zmax),
    label="",
)

savefig(figure, "volatility_surface.png")


# optimization
# want sigma s.t. f(sigma) - C = 0 
# define f as C(S_t, t) = N(d+)S_t - N(d-)Ke^(-r(T-t))
#
# Steps
# 1. write call option formula C(S_t, t)
# 2. write newton optimizer that takes f, df, x_o, tolerance, iters
# 3. Figure out what dC is 
# 4. Input C(S_t, t) - C and figure out df for this formula. 
# 5. get output sigma
# 6. input sigma into formula 
# 7. do this many times to get a surface for a certain expiry date
# 8. with this Surface and interpolation, can price many other options
