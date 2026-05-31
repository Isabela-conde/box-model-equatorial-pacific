using Oceananigans
using Oceananigans.Units: hours
using CairoMakie
using Printf
using Statistics

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