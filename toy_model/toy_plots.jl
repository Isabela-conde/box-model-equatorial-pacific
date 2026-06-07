using Oceananigans
using Oceananigans.Units: hours, days
using CairoMakie
using Printf
using Statistics
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using SeawaterPolynomials: ρ′

function plot_surface_wind_profile(grid, τx_parameters, output_dir)

    Ly = grid.Ly

    ys = range(-Ly/2, Ly/2, length = 300)
    lats = ys ./ 111e3

    u10s = [surface_u10_profile(y, τx_parameters) for y in ys]
    τxs  = [τx_profile(0, y, 0, τx_parameters) for y in ys]

    fig = Figure(size = (700, 400))

    ax1 = Axis(fig[1, 1],
               xlabel = "Surface zonal wind, u₁₀ (m s⁻¹)",
               ylabel = "Latitude",
               title = "Idealised surface zonal wind")

    lines!(ax1, u10s, lats, linewidth = 3)
    vlines!(ax1, 0, color = :gray, linestyle = :dash)

    ax2 = Axis(fig[1, 2],
               xlabel = "Surface stress, τx (m² s⁻²)",
               ylabel = "Latitude",
               title = "Applied Oceananigans flux")

    lines!(ax2, τxs, lats, linewidth = 3)
    vlines!(ax2, 0, color = :gray, linestyle = :dash)

    save(joinpath(output_dir, "surface_wind_profile_check.png"), fig)

    return fig
end

function load_field_time_series(filepath)

    return (u  = FieldTimeSeries(filepath, "u"),
            w  = FieldTimeSeries(filepath, "w"),
            T  = FieldTimeSeries(filepath, "T"),
            S  = FieldTimeSeries(filepath, "S"),
            dTdz = FieldTimeSeries(filepath, "dTdz"))
end

############################################################################################################
##################### video of equatorial cross section of u, w, T, S, and N2 at the equ ###################
############################################################################################################

function plot_equator_cross_section(filepath,
                                    grid,
                                    output_dir;
                                    top_400m = false)

    ts = load_field_time_series(filepath)

    times = ts.w.times
    n = Observable(1)

    yC = Array(ynodes(grid, Center()))
    j_equator = argmin(abs.(yC))

    xC = Array(xnodes(grid, Center()))
    xF = Array(xnodes(grid, Face()))
    zC = Array(znodes(grid, Center()))
    zF = Array(znodes(grid, Face()))

    if top_400m

        kC = findall(zC .>= -400)
        kF = findall(zF .>= -400)

        zC_plot = zC[kC]
        zF_plot = zF[kF]

        zlims = (-400, 0)

        outfile = "equator_cross_section_400m.mp4"

    else

        kC = eachindex(zC)
        kF = eachindex(zF)

        zC_plot = zC
        zF_plot = zF

        zlims = (-grid.Lz, 0)

        outfile = "equator_cross_section_full_depth.mp4"

    end

    uₙ  = @lift Array(interior(ts.u[$n]))[:, j_equator, kC]
    wₙ  = @lift Array(interior(ts.w[$n]))[:, j_equator, kF]
    Tₙ  = @lift Array(interior(ts.T[$n]))[:, j_equator, kC]
    Sₙ  = @lift Array(interior(ts.S[$n]))[:, j_equator, kC]
    dTdzₙ = @lift Array(interior(ts.dTdz[$n]))[:, j_equator, kF]

    fig = Figure(size = (1800, 1000))

    axis_kwargs = (xlabel = "x (m)",
                   ylabel = "z (m)",
                   limits = ((0, grid.Lx), zlims))

    ax_u  = Axis(fig[2, 1]; title = "u at equator",  axis_kwargs...)
    ax_w  = Axis(fig[2, 3]; title = "w at equator",  axis_kwargs...)
    ax_T  = Axis(fig[3, 1]; title = "T at equator",  axis_kwargs...)
    ax_S  = Axis(fig[3, 3]; title = "S at equator",  axis_kwargs...)
    ax_dTdz = Axis(fig[4, 1]; title = "dTdz at equator", axis_kwargs...)

    hm_u  = heatmap!(ax_u,  xF, zC_plot, uₙ;
                     colormap = :balance,
                     colorrange = (-2, 2))

    hm_w  = heatmap!(ax_w,  xC, zF_plot, wₙ;
                     colormap = :balance,
                     colorrange = (-1e-4, 1e-4))

    hm_T  = heatmap!(ax_T,  xC, zC_plot, Tₙ;
                     colormap = :thermal,
                     colorrange = (2, 30))

    hm_S  = heatmap!(ax_S,  xC, zC_plot, Sₙ;
                     colormap = :haline,
                     colorrange = (34.6, 35.4))

    hm_dTdz = heatmap!(ax_dTdz, xC, zF_plot, dTdzₙ;
                     colormap = :thermal,
                     colorrange = (-1, 1))

    Colorbar(fig[2, 2], hm_u)
    Colorbar(fig[2, 4], hm_w)
    Colorbar(fig[3, 2], hm_T)
    Colorbar(fig[3, 4], hm_S)
    Colorbar(fig[4, 2], hm_dTdz)

    title = @lift @sprintf("Equator, y = %.1f m, t = %s",
                           yC[j_equator],
                           prettytime(times[$n]))

    fig[1, 1:4] = Label(fig,
                         title,
                         fontsize = 24,
                         tellwidth = false)

    frames = 1:length(times)

    CairoMakie.record(fig,
                      joinpath(output_dir, outfile),
                      frames,
                      framerate = 8) do i
        n[] = i
    end
end

############################################################################################################
############################### video of y-mean x-z cross section ##########################################
############################################################################################################

function plot_ymean_cross_section(filepath,
                                  grid,
                                  output_dir;
                                  top_400m = false)

    ts = load_field_time_series(filepath)

    times = ts.w.times
    n = Observable(1)

    xC = Array(xnodes(grid, Center()))
    xF = Array(xnodes(grid, Face()))
    zC = Array(znodes(grid, Center()))
    zF = Array(znodes(grid, Face()))

    if top_400m

        kC = findall(zC .>= -400)
        kF = findall(zF .>= -400)

        zC_plot = zC[kC]
        zF_plot = zF[kF]

        zlims = (-400, 0)

        outfile = "ymean_cross_section_400m.mp4"

    else

        kC = eachindex(zC)
        kF = eachindex(zF)

        zC_plot = zC
        zF_plot = zF

        zlims = (-grid.Lz, 0)

        outfile = "ymean_cross_section_full_depth.mp4"

    end

    uₙ  = @lift dropdims(mean(Array(interior(ts.u[$n]))[:, :, kC],  dims = 2), dims = 2)
    wₙ  = @lift dropdims(mean(Array(interior(ts.w[$n]))[:, :, kF],  dims = 2), dims = 2)
    Tₙ  = @lift dropdims(mean(Array(interior(ts.T[$n]))[:, :, kC],  dims = 2), dims = 2)
    Sₙ  = @lift dropdims(mean(Array(interior(ts.S[$n]))[:, :, kC],  dims = 2), dims = 2)
    dTdzₙ = @lift dropdims(mean(Array(interior(ts.dTdz[$n]))[:, :, kF], dims = 2), dims = 2)

    fig = Figure(size = (1800, 1000))

    axis_kwargs = (xlabel = "x (m)",
                   ylabel = "z (m)",
                   limits = ((0, grid.Lx), zlims))

    ax_u  = Axis(fig[2, 1]; title = "y-mean u",  axis_kwargs...)
    ax_w  = Axis(fig[2, 3]; title = "y-mean w",  axis_kwargs...)
    ax_T  = Axis(fig[3, 1]; title = "y-mean T",  axis_kwargs...)
    ax_S  = Axis(fig[3, 3]; title = "y-mean S",  axis_kwargs...)
    ax_dTdz = Axis(fig[4, 1]; title = "y-mean dTdz", axis_kwargs...)

    hm_u  = heatmap!(ax_u,  xF, zC_plot, uₙ;
                     colormap = :balance,
                     colorrange = (-2, 2))

    hm_w  = heatmap!(ax_w,  xC, zF_plot, wₙ;
                     colormap = :balance,
                     colorrange = (-1e-4, 1e-4))

    hm_T  = heatmap!(ax_T,  xC, zC_plot, Tₙ;
                     colormap = :thermal,
                     colorrange = (2, 30))

    hm_S  = heatmap!(ax_S,  xC, zC_plot, Sₙ;
                     colormap = :haline,
                     colorrange = (34.6, 35.4))

    hm_dTdz = heatmap!(ax_dTdz, xC, zF_plot, dTdzₙ;
                     colormap = :thermal,
                     colorrange = (-1, 1))

    Colorbar(fig[2, 2], hm_u)
    Colorbar(fig[2, 4], hm_w)
    Colorbar(fig[3, 2], hm_T)
    Colorbar(fig[3, 4], hm_S)
    Colorbar(fig[4, 2], hm_dTdz)

    title = @lift @sprintf("y-mean fields, t = %s",
                           prettytime(times[$n]))

    fig[1, 1:4] = Label(fig,
                         title,
                         fontsize = 24,
                         tellwidth = false)

    frames = 1:length(times)

    CairoMakie.record(fig,
                      joinpath(output_dir, outfile),
                      frames,
                      framerate = 8) do i
        n[] = i
    end
end

############################################################################################################
######################### video of equatorial band x-z cross section #######################################
############################################################################################################

function plot_equatorial_band_cross_section(filepath,
                                            grid,
                                            output_dir;
                                            y_band = 2.5e5)

    ts = load_field_time_series(filepath)

    times = ts.w.times
    n = Observable(1)

    yC = Array(ynodes(grid, Center()))
    j_band = findall(abs.(yC) .<= y_band)

    xC = Array(xnodes(grid, Center()))
    xF = Array(xnodes(grid, Face()))
    zC = Array(znodes(grid, Center()))
    zF = Array(znodes(grid, Face()))

    uₙ  = @lift dropdims(mean(Array(interior(ts.u[$n]))[:, j_band, :],  dims = 2), dims = 2)
    wₙ  = @lift dropdims(mean(Array(interior(ts.w[$n]))[:, j_band, :],  dims = 2), dims = 2)
    Tₙ  = @lift dropdims(mean(Array(interior(ts.T[$n]))[:, j_band, :],  dims = 2), dims = 2)
    Sₙ  = @lift dropdims(mean(Array(interior(ts.S[$n]))[:, j_band, :],  dims = 2), dims = 2)
    dTdzₙ = @lift dropdims(mean(Array(interior(ts.dTdz[$n]))[:, j_band, :], dims = 2), dims = 2)

    fig = Figure(size = (1800, 1000))

    axis_kwargs = (xlabel = "x (m)",
                   ylabel = "z (m)",
                   limits = ((0, grid.Lx), (-grid.Lz, 0)))

    ax_u  = Axis(fig[2, 1]; title = "band-mean u",  axis_kwargs...)
    ax_w  = Axis(fig[2, 3]; title = "band-mean w",  axis_kwargs...)
    ax_T  = Axis(fig[3, 1]; title = "band-mean T",  axis_kwargs...)
    ax_S  = Axis(fig[3, 3]; title = "band-mean S",  axis_kwargs...)
    ax_dTdz = Axis(fig[4, 1]; title = "band-mean N²", axis_kwargs...)

    hm_u  = heatmap!(ax_u,  xF, zC, uₙ;
                     colormap = :balance,
                     colorrange = (-2, 2))

    hm_w  = heatmap!(ax_w,  xC, zF, wₙ;
                     colormap = :balance,
                     colorrange = (-1e-4, 1e-4))

    hm_T  = heatmap!(ax_T,  xC, zC, Tₙ;
                     colormap = :thermal,
                     colorrange = (2, 30))

    hm_S  = heatmap!(ax_S,  xC, zC, Sₙ;
                     colormap = :haline,
                     colorrange = (34.6, 35.4))

    hm_dTdz = heatmap!(ax_dTdz, xC, zF, dTdzₙ;
                     colormap = :thermal,
                     colorrange = (-1, 1))

    Colorbar(fig[2, 2], hm_u)
    Colorbar(fig[2, 4], hm_w)
    Colorbar(fig[3, 2], hm_T)
    Colorbar(fig[3, 4], hm_S)
    Colorbar(fig[4, 2], hm_dTdz)

    title = @lift @sprintf("Equatorial band mean (|y| ≤ %.0f km), t = %s",
                           y_band / 1000,
                           prettytime(times[$n]))

    fig[1, 1:4] = Label(fig,
                         title,
                         fontsize = 24,
                         tellwidth = false)

    frames = 1:length(times)

    CairoMakie.record(fig,
                      joinpath(output_dir,
                               "equatorial_band_cross_section.mp4"),
                      frames,
                      framerate = 8) do i
        n[] = i
    end
end

############################################################################################################
############################### video of equatorial band top 400m ##########################################
############################################################################################################

function plot_equatorial_band_top_400m(filepath,
                                       grid,
                                       output_dir;
                                       y_band = 2.5e5)

    ts = load_field_time_series(filepath)

    times = ts.w.times
    n = Observable(1)

    yC = Array(ynodes(grid, Center()))
    j_band = findall(abs.(yC) .<= y_band)

    xC = Array(xnodes(grid, Center()))
    xF = Array(xnodes(grid, Face()))
    zC = Array(znodes(grid, Center()))
    zF = Array(znodes(grid, Face()))

    k400C = findall(zC .>= -400)
    k400F = findall(zF .>= -400)

    zC400 = zC[k400C]
    zF400 = zF[k400F]

    uₙ  = @lift dropdims(mean(Array(interior(ts.u[$n]))[:, j_band, k400C],  dims = 2), dims = 2)
    wₙ  = @lift dropdims(mean(Array(interior(ts.w[$n]))[:, j_band, k400F],  dims = 2), dims = 2)
    Tₙ  = @lift dropdims(mean(Array(interior(ts.T[$n]))[:, j_band, k400C],  dims = 2), dims = 2)
    Sₙ  = @lift dropdims(mean(Array(interior(ts.S[$n]))[:, j_band, k400C],  dims = 2), dims = 2)
    dTdzₙ = @lift dropdims(mean(Array(interior(ts.dTdz[$n]))[:, j_band, k400F], dims = 2), dims = 2)

    fig = Figure(size = (1800, 1000))

    axis_kwargs = (xlabel = "x (m)",
                   ylabel = "z (m)",
                   limits = ((0, grid.Lx), (-400, 0)))

    ax_u  = Axis(fig[2, 1]; title = "band-mean u",  axis_kwargs...)
    ax_w  = Axis(fig[2, 3]; title = "band-mean w",  axis_kwargs...)
    ax_T  = Axis(fig[3, 1]; title = "band-mean T",  axis_kwargs...)
    ax_S  = Axis(fig[3, 3]; title = "band-mean S",  axis_kwargs...)
    ax_dTdz = Axis(fig[4, 1]; title = "band-mean dTdz", axis_kwargs...)

    hm_u  = heatmap!(ax_u,  xF, zC400, uₙ;
                     colormap = :balance,
                     colorrange = (-2, 2))

    hm_w  = heatmap!(ax_w,  xC, zF400, wₙ;
                     colormap = :balance,
                     colorrange = (-1e-4, 1e-4))

    hm_T  = heatmap!(ax_T,  xC, zC400, Tₙ;
                     colormap = :thermal,
                     colorrange = (2, 30))

    hm_S  = heatmap!(ax_S,  xC, zC400, Sₙ;
                     colormap = :haline,
                     colorrange = (34.6, 35.4))

    hm_dTdz = heatmap!(ax_dTdz, xC, zF400, dTdzₙ;
                     colormap = :thermal,
                     colorrange = (-1, 1))

    Colorbar(fig[2, 2], hm_u)
    Colorbar(fig[2, 4], hm_w)
    Colorbar(fig[3, 2], hm_T)
    Colorbar(fig[3, 4], hm_S)
    Colorbar(fig[4, 2], hm_dTdz)

    title = @lift @sprintf("Equatorial band mean top 400m (|y| ≤ %.0f km), t = %s",
                           y_band / 1000,
                           prettytime(times[$n]))

    fig[1, 1:4] = Label(fig,
                         title,
                         fontsize = 24,
                         tellwidth = false)

    frames = 1:length(times)

    CairoMakie.record(fig,
                      joinpath(output_dir,
                               "equatorial_band_cross_section_400m.mp4"),
                      frames,
                      framerate = 8) do i
        n[] = i
    end
end

function plot_surface_video(filepath, grid, output_dir)

    ts = (
        u  = FieldTimeSeries(filepath, "u"),
        w  = FieldTimeSeries(filepath, "w"),
        T  = FieldTimeSeries(filepath, "T"),
        S  = FieldTimeSeries(filepath, "S"),
        dTdz = FieldTimeSeries(filepath, "dTdz")
    )

    times = ts.T.times

    xC = Array(xnodes(grid, Center()))
    xF = Array(xnodes(grid, Face()))

    yC = Array(ynodes(grid, Center()))
    yF = Array(ynodes(grid, Face()))

    k_surface_C = size(interior(ts.T[1]), 3)
    k_surface_F = size(interior(ts.w[1]), 3)

    n = Observable(1)

    uₙ = @lift Array(interior(ts.u[$n]))[:, :, k_surface_C]
    wₙ = @lift Array(interior(ts.w[$n]))[:, :, k_surface_F]

    Tₙ = @lift Array(interior(ts.T[$n]))[:, :, k_surface_C]
    Sₙ = @lift Array(interior(ts.S[$n]))[:, :, k_surface_C]
    dTdzₙ = @lift Array(interior(ts.dTdz[$n]))[:, :, k_surface_F]

    fig = Figure(size = (1400, 800))

    axu = Axis(fig[1,1], title = "surface u")
    axw = Axis(fig[1,2], title = "surface w")
    axT = Axis(fig[2,1], title = "surface T")
    axS = Axis(fig[2,2], title = "surface S")
    axN = Axis(fig[3,1], title = "surface N²")

    hmu = heatmap!(axu, xF ./ 1e3, yC ./ 1e3, uₙ,
                   colorrange = (-5, 5),
                   colormap = :balance)

    hmw = heatmap!(axw, xC ./ 1e3, yC ./ 1e3, wₙ,
                   colorrange = (-1e-3, 1e-3),
                   colormap = :balance)

    hmT = heatmap!(axT, xC ./ 1e3, yC ./ 1e3, Tₙ,
                   colorrange = (2, 30),
                   colormap = :thermal)

    hmS = heatmap!(axS, xC ./ 1e3, yC ./ 1e3, Sₙ,
                   colorrange = (34.6, 36.4),
                   colormap = :haline)

    hmdTdz = heatmap!(axN, xC ./ 1e3, yC ./ 1e3,dTdzₙ,
                   colorrange = (-1, 1),
                   colormap = :balance)

    Colorbar(fig[1,3], hmu)
    Colorbar(fig[2,3], hmT)
    Colorbar(fig[3,3], hmdTdz)

    record(fig,
           joinpath(output_dir, "surface_fields.mp4"),
           1:length(times),
           framerate = 8) do i

        n[] = i

    end
end


function plot_surface_wind_boundary_condition(grid, τx_parameters, output_dir;
                                              flip_sign_for_bc = false)

    Ly = grid.Ly

    ys = range(-Ly/2, Ly/2, length = 300)
    lats = ys ./ 111e3

    u10s = [surface_u10_profile(y, τx_parameters) for y in ys]
    τxs  = [τx_profile(0, y, 0, τx_parameters) for y in ys]

    if flip_sign_for_bc
        top_flux = [-τx_profile(0, y, 0, τx_parameters) for y in ys]
    else
        top_flux = [ τx_profile(0, y, 0, τx_parameters) for y in ys]
    end

    implied_u_tendency = -top_flux

    fig = Figure(size = (1200, 700))

    ax1 = Axis(fig[1, 1],
               xlabel = "u₁₀ (m s⁻¹)",
               ylabel = "Latitude",
               title = "Wind speed profile")

    ax2 = Axis(fig[1, 2],
               xlabel = "τx_profile",
               ylabel = "Latitude",
               title = "Raw τx_profile")

    ax3 = Axis(fig[2, 1],
               xlabel = "Top boundary flux",
               ylabel = "Latitude",
               title = "FluxBoundaryCondition value")

    ax4 = Axis(fig[2, 2],
               xlabel = "Implied u tendency sign",
               ylabel = "Latitude",
               title = "- top boundary flux")

    lines!(ax1, u10s, lats, linewidth = 3)
    lines!(ax2, τxs, lats, linewidth = 3)
    lines!(ax3, top_flux, lats, linewidth = 3)
    lines!(ax4, implied_u_tendency, lats, linewidth = 3)

    vlines!(ax1, 0, color = :gray, linestyle = :dash)
    vlines!(ax2, 0, color = :gray, linestyle = :dash)
    vlines!(ax3, 0, color = :gray, linestyle = :dash)
    vlines!(ax4, 0, color = :gray, linestyle = :dash)

    save(joinpath(output_dir, "surface_wind_boundary_condition_check.png"), fig)

    return fig
end

function plot_surface_wind_boundary_condition(grid, τx_parameters, output_dir;
                                              flip_sign_for_bc = false)

    Ly = grid.Ly

    ys = range(-Ly/2, Ly/2, length = 300)
    lats = ys ./ 111e3

    u10s = [surface_u10_profile(y, τx_parameters) for y in ys]
    τxs  = [τx_profile(0, y, 0, τx_parameters) for y in ys]

    if flip_sign_for_bc
        top_flux = [-τx_profile(0, y, 0, τx_parameters) for y in ys]
    else
        top_flux = [ τx_profile(0, y, 0, τx_parameters) for y in ys]
    end

    implied_u_tendency = -top_flux

    fig = Figure(size = (1200, 700))

    ax1 = Axis(fig[1, 1],
               xlabel = "u₁₀ (m s⁻¹)",
               ylabel = "Latitude",
               title = "Wind speed profile")

    ax2 = Axis(fig[1, 2],
               xlabel = "τx_profile",
               ylabel = "Latitude",
               title = "Raw τx_profile")

    ax3 = Axis(fig[2, 1],
               xlabel = "Top boundary flux",
               ylabel = "Latitude",
               title = "FluxBoundaryCondition value")

    ax4 = Axis(fig[2, 2],
               xlabel = "Implied u tendency sign",
               ylabel = "Latitude",
               title = "- top boundary flux")

    lines!(ax1, u10s, lats, linewidth = 3)
    lines!(ax2, τxs, lats, linewidth = 3)
    lines!(ax3, top_flux, lats, linewidth = 3)
    lines!(ax4, implied_u_tendency, lats, linewidth = 3)

    vlines!(ax1, 0, color = :gray, linestyle = :dash)
    vlines!(ax2, 0, color = :gray, linestyle = :dash)
    vlines!(ax3, 0, color = :gray, linestyle = :dash)
    vlines!(ax4, 0, color = :gray, linestyle = :dash)

    save(joinpath(output_dir, "surface_wind_boundary_condition_check.png"), fig)

    return fig
end

function plot_surface_wind_boundary_condition(grid, τx_parameters, output_dir;
    flip_sign_for_bc = false)

    Ly = grid.Ly

    ys = range(-Ly/2, Ly/2, length = 300)
    lats = ys ./ 111e3

    u10s = [surface_u10_profile(y, τx_parameters) for y in ys]
    τxs  = [τx_profile(0, y, 0, τx_parameters) for y in ys]

    if flip_sign_for_bc
        top_flux = [-τx_profile(0, y, 0, τx_parameters) for y in ys]
    else
        top_flux = [ τx_profile(0, y, 0, τx_parameters) for y in ys]
    end

    implied_u_tendency = -top_flux

    fig = Figure(size = (1200, 700))

    ax1 = Axis(fig[1, 1],
               xlabel = "u₁₀ (m s⁻¹)",
               ylabel = "Latitude",
               title = "Wind speed profile")

    ax2 = Axis(fig[1, 2],
               xlabel = "τx_profile",
               ylabel = "Latitude",
               title = "Raw τx_profile")

    ax3 = Axis(fig[2, 1],
               xlabel = "Top boundary flux",
               ylabel = "Latitude",
               title = "FluxBoundaryCondition value")

    ax4 = Axis(fig[2, 2],
               xlabel = "Implied u tendency sign",
               ylabel = "Latitude",
               title = "- top boundary flux")

    lines!(ax1, u10s, lats, linewidth = 3)
    lines!(ax2, τxs, lats, linewidth = 3)
    lines!(ax3, top_flux, lats, linewidth = 3)
    lines!(ax4, implied_u_tendency, lats, linewidth = 3)

    vlines!(ax1, 0, color = :gray, linestyle = :dash)
    vlines!(ax2, 0, color = :gray, linestyle = :dash)
    vlines!(ax3, 0, color = :gray, linestyle = :dash)
    vlines!(ax4, 0, color = :gray, linestyle = :dash)

    save(joinpath(output_dir, "surface_wind_boundary_condition_check.png"), fig)

    return fig
end 

function plot_equator_full_depth(filepath, grid, output_dir)
    plot_equator_cross_section(filepath, grid, output_dir; top_400m = false)
end

function plot_equator_top_400m(filepath, grid, output_dir)
    plot_equator_cross_section(filepath, grid, output_dir; top_400m = true)
end

function plot_ymean_full_depth(filepath, grid, output_dir)
    plot_ymean_cross_section(filepath, grid, output_dir; top_400m = false)
end

function plot_ymean_top_400m(filepath, grid, output_dir)
    plot_ymean_cross_section(filepath, grid, output_dir; top_400m = true)
end

function plot_equatorial_band_full_depth(filepath, grid, output_dir)
    plot_equatorial_band_cross_section(filepath, grid, output_dir)
end

function plot_surface(filepath, grid, output_dir)
    plot_surface_video(filepath, grid, output_dir)
end

function plot_surface_wind_boundary(grid, τx_parameters, output_dir)
    plot_surface_wind_boundary_condition(grid, τx_parameters, output_dir)
end

# Plot the momentum spin-down diagnostics written to *_momentum.jld2.
# U_deep relaxing toward ~0 (and U_deep_abs / KE_deep flattening) confirms the
# deep ocean is spinning down rather than accumulating zonal momentum.
function plot_momentum_diagnostics(ts_file, output_dir)

    U_deep     = FieldTimeSeries(ts_file, "U_deep")
    U_upper    = FieldTimeSeries(ts_file, "U_upper")
    U_deep_abs = FieldTimeSeries(ts_file, "U_deep_abs")
    KE_deep    = FieldTimeSeries(ts_file, "KE_deep")

    times = U_deep.times ./ 86400  # seconds → days

    val(ts) = [interior(ts[n])[1, 1, 1] for n in 1:length(ts.times)]

    fig = Figure(size = (1100, 800))

    ax1 = Axis(fig[1, 1], xlabel = "time (days)",
               ylabel = "∫ u dV  (m⁴ s⁻¹)",
               title = "Volume-integrated zonal momentum")
    lines!(ax1, times, val(U_upper), label = "upper (z ≥ -300 m)")
    lines!(ax1, times, val(U_deep),  label = "deep (z < -300 m)")
    hlines!(ax1, [0.0], color = :black, linestyle = :dash)
    axislegend(ax1, position = :rt)

    ax2 = Axis(fig[2, 1], xlabel = "time (days)",
               ylabel = "deep momentum / energy",
               title = "Deep-ocean spin-down (should flatten)")
    lines!(ax2, times, val(U_deep_abs), label = "∫ |u| dV  (m⁴ s⁻¹)")
    lines!(ax2, times, val(KE_deep),    label = "∫ u² dV  (m⁵ s⁻²)")
    axislegend(ax2, position = :rt)

    save(joinpath(output_dir, "momentum_diagnostics.png"), fig)
    return fig
end

################################################
#### Last-100-day analysis plots (helpers) ####
################################################

const G_EARTH = 9.80665   # m s⁻², matches Oceananigans SeawaterBuoyancy default

# Buoyancy from (Conservative) T and (Absolute) S interior arrays.
# b = -g ρ′ / ρᵣ, with ρ′ the TEOS10 density anomaly and z the cell-centre depth.
function buoyancy_from_TS(T, S, z, eos)
    ρᵣ = eos.reference_density
    b  = similar(T)
    Nx, Ny, Nz = size(T)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        b[i, j, k] = -G_EARTH * ρ′(T[i, j, k], S[i, j, k], z[k], eos) / ρᵣ
    end
    return b
end

# Vertical derivative along a column on a (possibly stretched) z grid.
function ddz(f::AbstractVector, z::AbstractVector)
    n = length(f)
    d = similar(f)
    @inbounds for k in 2:n-1
        d[k] = (f[k+1] - f[k-1]) / (z[k+1] - z[k-1])
    end
    d[1] = (f[2] - f[1]) / (z[2] - z[1])
    d[n] = (f[n] - f[n-1]) / (z[n] - z[n-1])
    return d
end

# Symmetric colour range (robust to outliers) for diverging fields.
function symrange(A; q = 0.995)
    m = quantile(abs.(vec(A)), q)
    m = m == 0 ? 1.0 : m
    return (-m, m)
end

# Mean thermocline depth z(x): depth of maximum stratification db/dz in each
# column, searched within [-depth_max, z_skip] (skip the very surface).
function thermocline_depth(dbdz, z; depth_max = 400.0, z_skip = -20.0)
    Nx = size(dbdz, 1)
    krange = findall(zk -> -depth_max <= zk <= z_skip, z)
    zt = Vector{Float64}(undef, Nx)
    @inbounds for i in 1:Nx
        kloc = krange[argmax(view(dbdz, i, krange))]
        zt[i] = z[kloc]
    end
    return zt
end

################################################
#### y-mean / equatorial-band cross sections, last 100 days ####
################################################
# y_band = nothing → average over all y (full ȳ-mean).
# y_band = 250e3   → average only |y| ≤ 250 km (equatorial-band mean).
function plot_ymean_last100(filepath, output_dir;
                            window_days       = 100,
                            depth_max         = 400.0,
                            reference_density = 1026.0,
                            y_band            = nothing)

    Tts = FieldTimeSeries(filepath, "T")
    Sts = FieldTimeSeries(filepath, "S")
    uts = FieldTimeSeries(filepath, "u")

    grid  = Tts.grid
    times = Tts.times

    xc = xnodes(grid, Center(), Center(), Center())
    yc = ynodes(grid, Center(), Center(), Center())
    zc = znodes(grid, Center(), Center(), Center())
    Nx = length(xc)

    if y_band === nothing
        jsel       = collect(1:length(yc))
        band_label = "ȳ-mean"
        band_tag   = "ymean"
    else
        jsel       = findall(y -> abs(y) <= y_band, yc)
        band_label = @sprintf("eq-band ⟨|y|≤%d km⟩", round(Int, y_band / 1e3))
        band_tag   = @sprintf("eqband%dkm", round(Int, y_band / 1e3))
    end

    t_end  = times[end]
    idxs   = findall(t -> t >= t_end - window_days * days, times)
    i0, i1 = first(idxs), last(idxs)
    nwin   = length(idxs)
    @info "Averaging frames $i0:$i1 ($nwin frames), days " *
          @sprintf("%.0f-%.0f", times[i0] / days, times[i1] / days)

    eos = TEOS10EquationOfState(reference_density = reference_density)

    T0 = interior(Tts[i0])
    Ny = size(T0, 2)
    Nz = size(T0, 3)

    Tmean = zeros(Nx, Ny, Nz)
    bmean = zeros(Nx, Ny, Nz)
    umean = zeros(Nx, Ny, Nz)

    interp_u(u) = size(u, 1) == Nx ? u : 0.5 .* (u[1:end-1, :, :] .+ u[2:end, :, :])

    for n in idxs
        T = interior(Tts[n]); S = interior(Sts[n]); u = interior(uts[n])
        Tmean .+= T
        bmean .+= buoyancy_from_TS(T, S, zc, eos)
        umean .+= interp_u(u)
    end
    Tmean ./= nwin; bmean ./= nwin; umean ./= nwin

    ymean(A) = dropdims(mean(view(A, :, jsel, :), dims = 2), dims = 2)
    Tym = ymean(Tmean)
    bym = ymean(bmean)
    uym = ymean(umean)

    dbdzym = similar(bym)
    for i in 1:Nx
        dbdzym[i, :] = ddz(view(bym, i, :), zc)
    end

    zt = thermocline_depth(dbdzym, zc; depth_max = depth_max)

    xkm = xc ./ 1e3

    panels = [
        (Tym,    "$band_label temperature  ⟨T⟩ (°C)",        :thermal, :seq),
        (bym,    "$band_label buoyancy  ⟨b⟩ (m s⁻²)",        :viridis, :seq),
        (dbdzym, "$band_label stratification  ⟨∂b/∂z⟩ (s⁻²)", :dense,   :seq),
        (uym,    "$band_label zonal velocity  ⟨u⟩ (m s⁻¹)",   :balance, :div),
    ]

    fig = Figure(size = (1200, 1250))

    for (row, (A, ttl, cmap, kind)) in enumerate(panels)
        ax = Axis(fig[row, 1];
                  title  = ttl,
                  xlabel = row == length(panels) ? "x (km)" : "",
                  ylabel = "z (m)")

        if kind == :div
            _, m   = symrange(A)
            levels = range(-m, m, length = 22)   # even edges → odd bands → white band at 0
        else
            lo, hi = quantile(vec(A), 0.005), quantile(vec(A), 0.995)
            levels = range(lo, hi, length = 21)
        end

        cf = contourf!(ax, xkm, zc, A;
                       levels = levels, colormap = cmap,
                       extendlow = :auto, extendhigh = :auto)

        lines!(ax, xkm, zt; color = :red, linewidth = 3.5)

        ylims!(ax, -depth_max, 0)
        Colorbar(fig[row, 2], cf)
    end

    Label(fig[0, 1:2],
          @sprintf("%s, last %d-day mean (days %.0f-%.0f) — red line = thermocline (max ∂b/∂z)",
                   band_label, window_days, times[i0] / days, times[i1] / days);
          fontsize = 18)

    outfile = joinpath(output_dir,
                       @sprintf("%s_last%dd_cross_sections.png", band_tag, window_days))
    save(outfile, fig)
    @info "Saved $outfile"
    return fig
end

################################################
#### Surface fields, last 100 days mean ####
################################################
# Same fields as the surface video (u, w, T, S, dTdz≈N²), as a last-window mean.
function plot_surface_last100(filepath, output_dir; window_days = 100)

    Tts    = FieldTimeSeries(filepath, "T")
    Sts    = FieldTimeSeries(filepath, "S")
    uts    = FieldTimeSeries(filepath, "u")
    wts    = FieldTimeSeries(filepath, "w")
    dTdzts = FieldTimeSeries(filepath, "dTdz")

    grid  = Tts.grid
    times = Tts.times

    xC = xnodes(grid, Center()) ./ 1e3
    xF = xnodes(grid, Face())   ./ 1e3
    yC = ynodes(grid, Center()) ./ 1e3

    t_end = times[end]
    idxs  = findall(t -> t >= t_end - window_days * days, times)
    nwin  = length(idxs)
    @info "Surface mean: frames $(first(idxs)):$(last(idxs)) ($nwin frames), days " *
          @sprintf("%.0f-%.0f", times[first(idxs)] / days, times[last(idxs)] / days)

    # z-face fields (w, dTdz): top face is the surface boundary (≡0), so sample
    # the first face below it (zface=true).
    function mean_surface(tsf; zface = false)
        a1  = interior(tsf[first(idxs)])
        k   = zface ? size(a1, 3) - 1 : size(a1, 3)
        acc = zeros(size(a1, 1), size(a1, 2))
        for n in idxs
            acc .+= interior(tsf[n])[:, :, k]
        end
        return acc ./ nwin
    end

    ū    = mean_surface(uts)
    w̄    = mean_surface(wts;    zface = true)
    T̄    = mean_surface(Tts)
    S̄    = mean_surface(Sts)
    dTdz̄ = mean_surface(dTdzts; zface = true)

    fig = Figure(size = (1200, 1450))

    function panel(row, x, y, A, ttl, cmap, kind)
        ax = Axis(fig[row, 1]; title = ttl,
                  xlabel = row == 5 ? "x (km)" : "", ylabel = "y (km)")
        if kind == :div
            _, m   = symrange(A)
            levels = range(-m, m, length = 22)
        else
            lo, hi = quantile(vec(A), 0.005), quantile(vec(A), 0.995)
            levels = range(lo, hi, length = 21)
        end
        cf = contourf!(ax, x, y, A; levels = levels, colormap = cmap,
                       extendlow = :auto, extendhigh = :auto)
        Colorbar(fig[row, 2], cf)
    end

    panel(1, xF, yC, ū,    "surface ⟨u⟩ (m s⁻¹)",       :balance, :div)
    panel(2, xC, yC, w̄,    "surface ⟨w⟩ (m s⁻¹)",       :balance, :div)
    panel(3, xC, yC, T̄,    "surface ⟨T⟩ (°C)",          :thermal, :seq)
    panel(4, xC, yC, S̄,    "surface ⟨S⟩ (g kg⁻¹)",      :haline,  :seq)
    panel(5, xC, yC, dTdz̄, "surface ⟨∂T/∂z⟩ (°C m⁻¹)",  :balance, :div)

    Label(fig[0, 1:2],
          @sprintf("Surface fields, last %d-day mean (days %.0f-%.0f)",
                   window_days, times[first(idxs)] / days, times[last(idxs)] / days);
          fontsize = 18)

    outfile = joinpath(output_dir, @sprintf("surface_last%dd_mean.png", window_days))
    save(outfile, fig)
    @info "Saved $outfile"
    return fig
end

################################################
#### Surface heat budget (time series) ####
################################################
# Reconstructs the surface heat flux split by mechanism + the actual heat-content
# tendency.  ⚠ PARAMS MUST MATCH THE RUN that produced `filepath`.
function plot_surface_heat_budget(filepath, output_dir;
        ρₒ = 1026.0, cᴾ = 3991.0,
        λ_T = 7e-6, T_warm = 29.5, T_cold = 24.0,
        x_front = 0.6, front_width = 0.18,
        E_base = 1e-8, E_warm = 3e-8, Lᵥ = 2.5e6, ρ_fw = 1000.0,
        Q_solar_Wm2 = 70.0)

    Tts   = FieldTimeSeries(filepath, "T")
    grid  = Tts.grid
    times = Tts.times

    xc = xnodes(grid, Center())
    yc = ynodes(grid, Center())
    Nx = length(xc); Ny = length(yc)
    dx = xc[2] - xc[1]; dy = yc[2] - yc[1]
    dA = dx * dy
    Lx = Nx * dx
    area = Nx * dx * Ny * dy

    latent_coeff = ρ_fw * Lᵥ / (ρₒ * cᴾ)
    warmfrac(x)  = 0.5 * (1 - tanh((x / Lx - x_front) / front_width))
    Ttar = [T_cold + (T_warm - T_cold) * warmfrac(x) for x in xc]
    Evap = [E_base + E_warm * warmfrac(x)            for x in xc]

    Qlat_x    = (ρₒ * cᴾ * latent_coeff) .* Evap
    latent_PW = -sum(Qlat_x) * Ny * dA / 1e15
    solar_PW  =  Q_solar_Wm2 * area / 1e15

    Nt = length(times)
    restore_in  = zeros(Nt)
    restore_out = zeros(Nt)
    for n in 1:Nt
        Tsurf = interior(Tts[n])[:, :, end]
        q = (ρₒ * cᴾ * λ_T) .* (Ttar .- Tsurf)          # +ve = INTO ocean
        restore_in[n]  = sum(max.(q, 0.0)) * dA / 1e15
        restore_out[n] = sum(min.(q, 0.0)) * dA / 1e15
    end
    latent_out  = fill(latent_PW, Nt)
    solar_in    = fill(solar_PW, Nt)
    net_surface = restore_in .+ restore_out .+ latent_out .+ solar_in

    zf = znodes(grid, Face())
    dz = diff(zf)
    HC = zeros(Nt)
    for n in 1:Nt
        T = interior(Tts[n])
        HC[n] = ρₒ * cᴾ * dA * sum(T .* reshape(dz, 1, 1, :))
    end
    dHCdt = zeros(Nt)
    for n in 2:Nt-1
        dHCdt[n] = (HC[n+1] - HC[n-1]) / (times[n+1] - times[n-1])
    end
    dHCdt[1]   = (HC[2]   - HC[1])     / (times[2]   - times[1])
    dHCdt[end] = (HC[end] - HC[end-1]) / (times[end] - times[end-1])
    dHCdt ./= 1e15

    days_ = times ./ days
    fig = Figure(size = (1000, 950))

    ax1 = Axis(fig[1, 1]; xlabel = "time (days)",
               ylabel = "surface heat flux (PW, + into ocean)",
               title  = "Surface heat budget by mechanism — into(+) / out(−)")
    lines!(ax1, days_, solar_in,    color = :darkorange, linewidth = 2, label = "solar IN")
    lines!(ax1, days_, restore_in,  color = :firebrick,  linewidth = 2, label = "restoring IN (cold tongue)")
    lines!(ax1, days_, restore_out, color = :royalblue,  linewidth = 2, label = "restoring OUT")
    lines!(ax1, days_, latent_out,  color = :seagreen,   linewidth = 2, label = "evaporation OUT")
    lines!(ax1, days_, net_surface, color = :black,      linewidth = 3, label = "surface NET")
    hlines!(ax1, [0.0], color = :gray, linestyle = :dash)
    axislegend(ax1, position = :rb)

    ax2 = Axis(fig[2, 1]; xlabel = "time (days)",
               ylabel = "heat tendency (PW)",
               title  = "Ocean heat-content tendency  d/dt∫ρcₚT dV   (→ 0 = equilibrated)")
    lines!(ax2, days_, dHCdt,       color = :purple, linewidth = 3,
           label = "actual dHC/dt (all sources)")
    lines!(ax2, days_, net_surface, color = :black, linewidth = 2, linestyle = :dash,
           label = "surface NET (gap = N/S export)")
    hlines!(ax2, [0.0], color = :gray, linestyle = :dash)
    axislegend(ax2, position = :rt)

    outfile = joinpath(output_dir, "surface_heat_budget.png")
    save(outfile, fig)
    @info "Saved $outfile"
    return fig
end