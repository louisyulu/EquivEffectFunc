### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 19553487-944a-4e82-80b7-61f977c17f3b
# ╠═╡ show_logs = false
begin
    import Pkg
    Pkg.activate("Project.toml")
	Pkg.status()
end

# ╔═╡ e93b93c1-a4be-444b-ac4a-b6f522aad84d
# ╠═╡ show_logs = false
using Logging, Dierckx, CSV, DataFrames, EquivEffectFunc, WebIO, PlotlyJS, PlutoUI

# ╔═╡ aad7eb01-88ef-4de6-b6b4-09656bf626ca
md"Copyright (c) 2023 Louis Yu Lu, MIT License"

# ╔═╡ cf598f8e-e8dd-4e00-8477-f3e6b88e3f9a
# Make Dispay Area Wider
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	 padding-left: max(50px, 5%);
    	 padding-right: max(50px, 10%);
	}
</style>
"""

# ╔═╡ 8e915fa9-af4b-4b69-988f-8a552578e160
md"""
**Random Samples:**

Sample Count:`  `$(@bind n_samples confirm(Slider(20:100, default=80, show_value=true)))
"""

# ╔═╡ 41452d50-b099-4893-83cd-0b9b7830139f
begin
	x_samples = collect(0.0:1.0:n_samples-1)
	y_samples = 10.0*rand(n_samples) + rand(n_samples) .+ 6.0
	trace1 = scatter(x=x_samples, y=y_samples, name="Samples", mode="markers")
	spl = Spline1D(x_samples, y_samples, k=5, s=0.05)
	xs = collect(0.0:0.1:x_samples[end])
	ys = evaluate(spl, xs)
	trace2 = scatter(x=xs, y=ys, name="Spline", mode="line")
	Plot([trace1, trace2], Layout(height=400))
end

# ╔═╡ 5b56aa71-ad44-4e9d-b128-0cdf4d187436
md"""
**EEF with control points on extrema:**

**Polynomial Order:`  `$(@bind p_order NumberField(0:4, default=3))**
**`  `On Extrema of:`  `$(@bind on_extrema Select([:ys_val => "Function Values", :first_deriv => "Derivatives", :second_deriv => "Second Derivatives"], :first_deriv))**
"""

# ╔═╡ cc0ffc98-1078-4189-8da4-f2e2d0683ee5
let
	trace1 = scatter(x=xs, y=ys, name="Original", mode="line")
	trend, diff = eef_extrema(xs, ys, p_order=p_order, depend_on=on_extrema)
	trace2 = scatter(x=xs, y=trend, name="Trend", mode="line")
	hs = if on_extrema == :ys_val
        ys
    elseif on_extrema == :first_deriv
        deriv(ys, xs)
    elseif on_extrema == :second_deriv
        deriv2(ys, xs)
	end
	trace3 = scatter(x=xs, y=diff, name="Diff", mode="line")
	trace4 = scatter(x=xs, y=hs, name="On extrema", mode="line")
	Plot([trace1, trace2, trace3, trace4], Layout(height=400))
end

# ╔═╡ dd0d9441-49d7-4b52-bed2-334a866836b8
md"""
**EEF with control points on partitions:**

**Polynomial Order:`  `$(@bind p_order2 NumberField(0:4, default=3))**
**`  `On absolute values of:`  `$(@bind depend_on Select([:abs_ys_val => "Function Values", :abs_1st_deriv => "Derivatives", :abs_2nd_deriv => "Second Derivatives"], :abs_2nd_deriv))**
**`  `Split Levels:`  `$(@bind split_levels NumberField(2:7, default=4))**
"""

# ╔═╡ ae2a7045-3ec7-444d-8d3a-af1a435937b1
let
	trace1 = scatter(x=xs, y=ys, name="Original", mode="line")
	trend, diff = eef_partition(xs, ys, split_levels, p_order=p_order2, depend_on=depend_on)
	trace2 = scatter(x=xs, y=trend, name="Trend", mode="line")
	hs = if depend_on == :abs_ys_val
        abs.(ys)
    elseif depend_on == :abs_1st_deriv
        abs.(deriv(ys, xs))
    elseif depend_on == :abs_2nd_deriv
        abs.(deriv2(ys, xs))
	end
	trace3 = scatter(x=xs, y=diff, name="Diff", mode="line")
	trace4 = scatter(x=xs, y=hs, name="On values", mode="line")
	Plot([trace1, trace2, trace3, trace4], Layout(height=400))
end

# ╔═╡ 8232ee28-8ddf-4e17-8009-5bbe08cb90ae
# read stock historical data from csv file
csv_df = reverse(DataFrame(CSV.File("csv/HistoricalData_DOW.csv")))

# ╔═╡ 3bb4faed-fdc9-4214-8b1b-9c1ba6c27f95
md"""
**Stock value EEF with control points on extrema**

**Select column:`  `$(@bind col_name Select(["Close/Last", "Open", :"High", "Low"]))**
**`  `Polynomial Order:`  `$(@bind p_order3 NumberField(0:4, default=3))**
**`  `On Extrema of:`  `$(@bind on_extrema3 Select([:ys_val => "Function Values", :first_deriv => "Derivatives", :second_deriv => "Second Derivatives"], :first_deriv))**
"""

# ╔═╡ 2978e053-c5bd-4e28-8106-336ffb6bc636
let 
	vs = csv_df[!, "$col_name"]
	ds = csv_df[!, "Date"]
	trace1 = scatter(x=ds, y=vs, name="$col_name", mode="line")
	trend, diff = eef_extrema(vs, p_order=p_order3, depend_on=on_extrema3)
	trace2 = scatter(x=ds, y=trend, name="Trend", mode="line")
	hs = if on_extrema3 == :ys_val
        vs
	elseif on_extrema3 == :first_deriv
        deriv(vs)
	elseif on_extrema3 == :second_deriv
        deriv2(vs)
	end
	trace3 = scatter(x=ds, y=diff, name="Diff", mode="line")
	trace4 = scatter(x=ds, y=hs, name="On extrema", mode="line")
	Plot([trace1, trace2, trace3, trace4], Layout(height=400))
end

# ╔═╡ d79c1a8c-6f66-490e-80fe-6a39e05d9bc5
md"""
**Stock Value EEF with control points on partitions:**

**Polynomial Order:`  `$(@bind p_order4 NumberField(0:4, default=3))**
**`  `On absolute values of:`  `$(@bind depend_on4 Select([:abs_ys_val => "Function Values", :abs_1st_deriv => "Derivatives", :abs_2nd_deriv => "Second Derivatives"], :abs_2nd_deriv))**
**`  `Split Levels:`  `$(@bind split_levels4 NumberField(5:9, default=5))**
"""

# ╔═╡ dc874958-41fd-4e47-8fd1-1bdf30174a1f
let 
	vs = csv_df[!, "$col_name"]
	ds = csv_df[!, "Date"]
	trace1 = scatter(x=ds, y=vs, name="$col_name", mode="line")
	trend, diff =  eef_partition(vs, split_levels4, p_order=p_order4, depend_on=depend_on4)
	trace2 = scatter(x=ds, y=trend, name="Trend", mode="line")
	hs = if depend_on4 == :abs_ys_val
        abs.(vs)
    elseif depend_on == :abs_1st_deriv
        abs.(deriv(vs))
    elseif depend_on == :abs_2nd_deriv
        abs.(deriv2(vs))
	end
	trace3 = scatter(x=ds, y=diff, name="Diff", mode="line")
	trace4 = scatter(x=ds, y=hs, name="On values", mode="line")
	Plot([trace1, trace2, trace3, trace4], Layout(height=400))
end

# ╔═╡ Cell order:
# ╟─aad7eb01-88ef-4de6-b6b4-09656bf626ca
# ╠═cf598f8e-e8dd-4e00-8477-f3e6b88e3f9a
# ╠═19553487-944a-4e82-80b7-61f977c17f3b
# ╠═e93b93c1-a4be-444b-ac4a-b6f522aad84d
# ╟─8e915fa9-af4b-4b69-988f-8a552578e160
# ╟─41452d50-b099-4893-83cd-0b9b7830139f
# ╟─5b56aa71-ad44-4e9d-b128-0cdf4d187436
# ╟─cc0ffc98-1078-4189-8da4-f2e2d0683ee5
# ╟─dd0d9441-49d7-4b52-bed2-334a866836b8
# ╟─ae2a7045-3ec7-444d-8d3a-af1a435937b1
# ╠═8232ee28-8ddf-4e17-8009-5bbe08cb90ae
# ╟─3bb4faed-fdc9-4214-8b1b-9c1ba6c27f95
# ╟─2978e053-c5bd-4e28-8106-336ffb6bc636
# ╟─d79c1a8c-6f66-490e-80fe-6a39e05d9bc5
# ╟─dc874958-41fd-4e47-8fd1-1bdf30174a1f
