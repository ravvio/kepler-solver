using GLMakie
using CSV
using DataFrames

# Report

function report_add_line(label, max_eccentricity, axis_avg, axis_max, axis_time)
    df = CSV.read("../build/Release/arm64/output/report_$label.csv", DataFrame)
    dff = filter(cols -> cols.e <= max_eccentricity, df)
    lines!(axis_avg, dff.e, dff.error_avg_rad; label = label)
    lines!(axis_max, dff.e, dff.error_max_rad; label = label)
    lines!(axis_time, dff.e, dff.time_avg_ms; label = label)
end

function report(title, max_eccentricity, labels)
    fig = Figure(size = (800, 1200))
    axis_avg = Axis(fig[1, 1], xlabel = "Eccentricity", ylabel = "Error (rad)", title = "Average Error")
    axis_max = Axis(fig[2, 1], xlabel = "Eccentricity", ylabel = "Error (rad)", title = "Maximum Error")
    axis_time = Axis(fig[3, 1], xlabel = "Eccentricity", ylabel = "Time (ms)", title = "Average Execution Time (ms)")

    for label in labels
        report_add_line(label, max_eccentricity, axis_avg, axis_max, axis_time)
    end

    fig[4, 1] = Legend(fig, [axis_avg, axis_max, axis_time], merge = true,
                       orientation = :vertical, nbanks = 3, tellheight = true)
    save("./report_$title.png", fig)
end

report("hh_e050", 0.5, ["HH_PC_NR_NR", "HH_PC_HA_NR", "HH_PC_OG_NR"])
report("hh_e100", 1, ["HH_PC_NR_NR", "HH_PC_HA_NR", "HH_PC_OG_NR"])

report("ehh_e050", 0.5,["EHH_PC_NR", "EHH_PC_HA", "EHH_PC_OG"])
report("ehh_e100", 1, ["EHH_PC_NR", "EHH_PC_HA", "EHH_PC_OG"])

report("full_e050", 0.5, [
    "HH_PC_NR_NR", "HH_PC_HA_NR", "HH_PC_OG_NR",
    "EHH_PC_NR", "EHH_PC_HA", "EHH_PC_OG"
])
report("full_e100", 1, [
    "HH_PC_NR_NR", "HH_PC_HA_NR", "HH_PC_OG_NR",
    "EHH_PC_NR", "EHH_PC_HA", "EHH_PC_OG"
])

# Orbit

function orbit_add(label, eccentricities, axes)
    df = CSV.read("../build/Release/arm64/output/orbit_$label.csv", DataFrame)
    for (e, axis) in zip(eccentricities, axes)
        dff = filter(cols -> cols.e == e, df)
        scatter!(axis, dff.M, dff.time_ms; label = label, alpha = 0.75, markersize = 5)
    end
end

function orbits(title, labels, eccentricities)
    fig = Figure(size = (800, 400 * length(eccentricities)))

    axes = []
    for (i, e) in enumerate(eccentricities)
        axis = Axis(
            fig[i, 1],
            xlabel = "Mean Anomaly (rad)",
            xticks = (0:π/2:2π, ["0", "π/2", "π", "3π/2", "2π"]),
            ylabel = "Time (ms)",
            yscale = identity,
            title = "e = $e"
        )
        push!(axes, axis)
    end

    for label in labels
        orbit_add(label, eccentricities, axes)
    end

    for axis in axes
        axislegend(axis)
    end
    save("./orbit_$title.png", fig)
end

orbits("hh", ["HH_PC_NR_NR", "HH_PC_HA_NR", "HH_PC_OG_NR"], [0.01, 0.10, 0.50])
orbits("ehh", ["EHH_PC_NR", "EHH_PC_HA", "EHH_PC_OG"], [0.01, 0.10, 0.50])

orbits("halley", ["HH_PC_HA_NR", "EHH_PC_HA"], [0.01, 0.10, 0.50])
