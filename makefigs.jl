################################################################################################
#
# Script for producing figure in "Sticky coupling as a control variate for sensitivity analysis"
#
#
#
################################################################################################



using JLD
using Plots
using Statistics
using LaTeXStrings
using GLM
using DataFrames
using Measures

############
# data and save locations
# Define these variables with the directory where the data is located and with 
# the directory to save figures to before running this script
############
lj_cluster_data_dir = 
convex2d_data_dir =
nonconvex2d_data_dir = 
savedir = 

betas = [0.5, 1., 2., 4.]
etas = [0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5]
sticky_mob = zeros(length(betas), length(etas), 2, 200000)
sync_mob = zeros(length(betas), length(etas), 2, 200000)
NEMD_mob = zeros(length(betas), length(etas), 2, 200000) 

sticky_tilt = zeros(length(betas), length(etas), 2, 200000)
sync_tilt = zeros(length(betas), length(etas), 2, 200000)
NEMD_tilt = zeros(length(betas), length(etas), 2, 200000)

for (i, beta) in enumerate(betas), (j, eta) in enumerate(etas)
	sticky_mob[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sticky_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 3]
	sticky_mob[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sticky_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 3]

	sync_mob[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 3]
	sync_mob[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 3]

	NEMD_mob[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 4]
	NEMD_mob[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 4]

	sticky_tilt[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sticky_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 6]
	sticky_tilt[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sticky_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 6]

	sync_tilt[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 6]
	sync_tilt[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 6]

	NEMD_tilt[i, j, 1, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_traj")[:, 7]
	NEMD_tilt[i, j, 2, :] = load(lj_cluster_data_dir*"lj_sync_coupling_distmobilitypressuretilt_sinshear_seed8156733675_L2.0_eta$(eta)_beta$(beta)_T200000.0_dt0.0001_burn100000.0_interval10000.jld", "distresp_avg")[:, 7]

end

begin
	plot()
	let β_ind = 1, last_eta =8
	scatter!(etas, sticky_mob[β_ind, :, 2, end], label = "Sticky Coupling", ms = 12)
	scatter!(etas, sync_mob[β_ind, :, 2, end], label = "Synchronous Coupling", ms = 12)
	scatter!(etas, NEMD_mob[β_ind, :, 2, end], label = "NEMD", ms = 12)

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sticky_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sticky Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[1], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sync_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sync Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[2], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= NEMD_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("NEMD Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[3], lw = 3)
	end

	plot!(title = latexstring("Mobility, Sine Shear, \$ \\beta =\$ $(betas[β_ind])"))
	end

	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	savefig(savedir*"lj_mobility_sinshear_transcoef_fit.png")
end

begin
	plot()
	let β_ind = 4, last_eta =6
	scatter!(etas, sticky_mob[β_ind, :, 2, end], label = "Sticky Coupling", ms = 12)
	scatter!(etas, sync_mob[β_ind, :, 2, end], label = "Synchronous Coupling", ms = 12)
	scatter!(etas, NEMD_mob[β_ind, :, 2, end], label = "NEMD", ms = 12)

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sticky_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sticky Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[1], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sync_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sync Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[2], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= NEMD_mob[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("NEMD Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[3], lw = 3)
	end

	plot!(title = latexstring("Mobility, Sine Shear, \$ \\beta =\$ $(betas[β_ind])"))
	end
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	savefig(savedir*"lj_mobility_sinshear_transcoef_fit4.png")
end


begin
	plot()
	let β_ind = 2, last_eta =6
	scatter!(etas, sticky_tilt[β_ind, :, 2, end], label = "Sticky Coupling", ms = 12)
	scatter!(etas, sync_tilt[β_ind, :, 2, end], label = "Synchronous Coupling", ms = 12)
	scatter!(etas, NEMD_tilt[β_ind, :, 2, end], label = "NEMD", ms = 12)

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sticky_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sticky Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[1], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sync_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sync Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[2], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= NEMD_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("NEMD Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[3], lw = 3)
	end

	plot!(title = latexstring("Tilt, Sine Shear, \$ \\beta =\$ $(betas[β_ind])"))
	end
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	savefig(savedir*"lj_tilt_sinshear_transcoef_fit.png")
end

begin
	plot()
	let β_ind = 3, last_eta =6
	scatter!(etas, sticky_tilt[β_ind, :, 2, end], label = "Sticky Coupling", ms = 12)
	scatter!(etas, sync_tilt[β_ind, :, 2, end], label = "Synchronous Coupling", ms = 12)
	scatter!(etas, NEMD_tilt[β_ind, :, 2, end], label = "NEMD", ms = 12)

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sticky_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sticky Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[1], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= sync_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("Sync Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[2], lw = 3)
	end

	let 
		data = DataFrame(X = etas[1:last_eta], Y= NEMD_tilt[β_ind, 1:last_eta, 2, end])
		beta_1 = lm(@formula(Y~ 0 + X), data).model.pp.beta0[1]
		plot!(etas, beta_1 .* etas, label = latexstring("NEMD Fit line \$\\widehat{\\alpha_R} = \$", round(beta_1; digits = 2)), color = palette(:default)[3], lw = 3)
	end

	plot!(title = latexstring("Tilt, Sine Shear, \$ \\beta =\$ $(betas[β_ind])"))
	end
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	savefig(savedir*"lj_tilt_sinshear_transcoef_fit2.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{mobility}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 0.5\$"))
	plot!(etas, var(NEMD_mob[1, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_mob[1, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_mob[1, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_mobility_varbeta0.5.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{mobility}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 1\$"))
	plot!(etas, var(NEMD_mob[2, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_mob[2, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_mob[2, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_mobility_varbeta1.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{mobility}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 2\$"))
	plot!(etas, var(NEMD_mob[3, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_mob[3, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_mob[3, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_mobility_varbeta2.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{mobility}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 4\$"))
	plot!(etas, var(NEMD_mob[4, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_mob[4, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_mob[4, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_mobility_varbeta4.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{tilt}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 0.5\$"))
	plot!(etas, var(NEMD_tilt[1, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_tilt[1, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_tilt[1, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_tilt_varbeta0.5.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{tilt}}(X_k^{\\eta, \\Delta t}) - R(Y_k^{0, \\Delta t})\$, \$\\beta = 1\$"))
	plot!(etas, var(NEMD_tilt[2, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_tilt[2, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_tilt[2, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_tilt_varbeta1.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{tilt}}(X_k^{\\eta, \\Delta t}) - R_{\\mathrm{tilt}}(Y_k^{0, \\Delta t})\$, \$\\beta = 2\$"))
	plot!(etas, var(NEMD_tilt[3, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_tilt[3, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_tilt[3, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(ylims = (1, maximum(var(sticky_tilt[3, :, 1, :]; dims = 2))*1.05))
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_tilt_varbeta2.png")
end

begin
	plot(title = latexstring("Variance of \$R_{\\mathrm{tilt}}(X_k^{\\eta, \\Delta t}) - R_{\\mathrm{tilt}}(Y_k^{0, \\Delta t})\$, \$\\beta = 4\$"))
	plot!(etas, var(NEMD_tilt[4, :, 1, :]; dims = 2), marker = :circle, label = "NEMD", lw = 3, ms = 12)
	plot!(etas, var(sync_tilt[4, :, 1, :]; dims = 2), marker = :circle, label = "Sync", lw = 3, ms = 12)
	plot!(etas, var(sticky_tilt[4, :, 1, :]; dims = 2), marker = :circle, label = "Sticky", lw = 3, ms = 12)
	plot!(xscale = :log, yscale = :log)
	plot!(ylims = (1, maximum(var(sticky_tilt[4, :, 1, :]; dims = 2))*1.05))
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"lj_tilt_varbeta4.png")
end

begin
	p1 = plot(1.:1.:20000, sticky_mob[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-30, 30))
	plot!(p1, title = latexstring("Mobility, Sticky Coupling \$\\beta = 1\$"))
	plot!(p1, legend = false)
	plot!(p1, legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	p2 = plot(1.:1.:20000, sync_mob[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-30, 30))
	plot!(p2, title = latexstring("Mobility, Synchronous Coupling \$\\beta = 1\$"))
	plot!(p2, legend_title = L"\eta")
	plot!(p2, legendfontsize = 24, legendtitlefontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	plot(p1, p2, layout = (1,2), size = (2400, 800))
	savefig(savedir*"lj_mobility.png")
end

begin
	p1 = plot(1.:1.:20000, sticky_tilt[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-80, 80))
	plot!(p1, title = latexstring("Tilt, Sticky Coupling \$\\beta = 1\$"))
	plot!(p1, legend = false)
	plot!(p1, legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	p2 = plot(1.:1.:20000, sync_tilt[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-80, 80))
	plot!(p2, title = latexstring("Tilt, Synchronous Coupling \$\\beta = 1\$"))
	plot!(p2, legend_title = L"\eta")
	plot!(p2, legendfontsize = 24, legendtitlefontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	plot(p1, p2, layout = (1,2), size = (2400, 800))
	savefig(savedir*"lj_tilt.png")
end

begin
	p1 = plot(1.:1.:20000, sticky_mob[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-30, 30))
	plot!(p1, title = latexstring("Mobility, Sticky Coupling \$\\beta = 1\$"))
	plot!(p1, legend = false)
	plot!(p1, legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	p2 = plot(1.:1.:20000, sync_mob[2, :, 2, 1:20000]' ./ etas', label = etas', ylims = (-30, 30))
	plot!(p2, title = latexstring("Mobility, Synchronous Coupling \$\\beta = 1\$"))
	plot!(p2, legend_title = L"\eta")
	plot!(p2, legendfontsize = 24, legendtitlefontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)

	plot(p1, p2, layout = (2,1), size = (1200, 1600))
	savefig(savedir*"lj_mobility_stacked.png")
end

sticky_convex = zeros(500, 10, 5, 2)
sync_convex = zeros(500, 10, 5, 2)
sticky_nonconvex_resp = zeros(500, 10, 5, 2)
sync_nonconvex_resp = zeros(500, 10, 5, 2)
etas = [0.1, 0.05, 0.025, 0.01, 0.005]
ts = 100:100:1.0e6

for j = 1:5
	for i = 1:5
		seed = 6924031646
		sticky_convex[50*(2i-2)+1:50*(2i-1), :, j, 1] .= load(convex2d_data_dir*string("sticky_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 1]
		sticky_convex[50*(2i-2)+1:50*(2i-1), :,  j, 2] .= load(convex2d_data_dir*string("sticky_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 2]

		sync_convex[50*(2i-2)+1:50*(2i-1), :,  j, 1] .= load(convex2d_data_dir*string("sync_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 1]
		sync_convex[50*(2i-2)+1:50*(2i-1), :,  j, 2] .= load(convex2d_data_dir*string("sync_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 2]
	
		seed = 9909724496
		sticky_convex[50*(2i-1)+1:50*(2i), :,  j, 1] .= load(convex2d_data_dir*string("sticky_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 1]
		sticky_convex[50*(2i-1)+1:50*(2i), :,  j, 2] .= load(convex2d_data_dir*string("sticky_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 2]

		sync_convex[50*(2i-1)+1:50*(2i), :,  j, 1] .= load(convex2d_data_dir*string("sync_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 1]
		sync_convex[50*(2i-1)+1:50*(2i), :,  j, 2] .= load(convex2d_data_dir*string("sync_convex2D_resp_batch", i, "_seed", seed, "_eta", etas[j], "_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld"), "R_avg")[:, 1000:1000:10000, 2]
	end
	sticky_convex[:, :, j, :] /= etas[j]
	sync_convex[:, :, j, :] /= etas[j]
end

sticky_convex_vars = var(sticky_convex, dims=1)[1, :, :, :]
sync_convex_vars = var(sync_convex, dims=1)[1, :, :, :]

for j = 1:5 
	for i = 1:5
		sticky_nonconvex_resp[1+100*(i-1):100*i, :, j, 1] .= load(nonconvex2d_data_dir*"sticky_nonconvex2D_resp_L3.14_batch$(i)_seed9909724496_eta$(etas[j])_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns100.jld", "R_avg")[:, 1000:1000:10000, 1]
		sync_nonconvex_resp[1+100*(i-1):100*i, :, j, 1] .= load(nonconvex2d_data_dir*"sync_nonconvex2D_resp_L3.14_batch$(i)_seed9909724496_eta$(etas[j])_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns100.jld", "R_avg")[:, 1000:1000:10000, 1]
		
		sticky_nonconvex_resp[1+100*(i-1):100*i, :, j, 2] .= load(nonconvex2d_data_dir*"sticky_nonconvex2D_resp_L3.14_batch$(i)_seed9909724496_eta$(etas[j])_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns100.jld", "R_avg")[:, 1000:1000:10000, 2]
		sync_nonconvex_resp[1+100*(i-1):100*i, :, j, 2] .= load(nonconvex2d_data_dir*"sync_nonconvex2D_resp_L3.14_batch$(i)_seed9909724496_eta$(etas[j])_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns100.jld", "R_avg")[:, 1000:1000:10000, 2]

	end
	sticky_nonconvex_resp[:, :, j, :] ./= etas[j]
	sync_nonconvex_resp[:, :, j, :] ./= etas[j]
end

sticky_nonconvex_resp_vars = var(sticky_nonconvex_resp, dims = 1)[1, :, :, :]
sync_nonconvex_resp_vars = var(sync_nonconvex_resp, dims = 1)[1, :, :, :]

begin
	times = 100_000:100_000:1_000_000
	time_strings= [L"T = {%$(i)} \times 10^5" for i in 1:9]
	push!(time_strings, L"T = 10^6")
	p1 = plot()
	for i in [1, 3, 6, 10]
		plot!(p1, 1 ./ etas, sync_convex_vars[i, :, 1, 1]*times[i], label =  time_strings[i], marker = :circle, lw = 3, ms = 12)
	end
	plot!(p1, legend = :bottomright, xticks = (1 ./ etas, etas), xlabel = L"1/\eta", 
	ylabel = latexstring("Variances \$\\times \\, T\$"), title = "Synchronous Coupling in Harmonic Potential", ylims = (0, 1.5))
	plot!(p1, size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"convex_corr_sync.png")

	p2 = plot()
	for i in [1, 3, 6, 10]
		plot!(p2, 1 ./ etas, sticky_convex_vars[i, :, 1, 1]*times[i], label = string(L"T = ", time_strings[i]), marker = :circle, lw = 3, ms = 12)
	end
	plot!(p2, legend = false, xticks = (1 ./ etas, etas), xlabel = L"1/\eta", 
	ylabel = latexstring("Variances \$\\times \\, T\$"), title = "Sticky Coupling in Harmonic Potential")
	plot!(p2, size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 22, margin = 5.0mm)
	savefig(savedir*"convex_corr_sticky.png")
end

begin
	times = 100_000:100_000:1_000_000
	time_strings= [L"T = {%$(i)} \times 10^5" for i in 1:9]
	push!(time_strings, L"T = 10^6")
	p1 = plot()
	for i in [1, 3, 6, 10]
		plot!(p1, 1 ./ etas, sync_nonconvex_resp_vars[i, :, 1, 1]*times[i], label =  time_strings[i], marker = :circle, lw = 3, ms = 12)
	end
	plot!(p1, legend = :bottomright, xticks = (1 ./ etas, etas), xlabel = L"1/\eta", 
	ylabel = latexstring("Variances \$\\times \\, T\$"), title = "Synchronous Coupling in Non-Convex Potential", ylims = (0, 1.5e5))
	plot!(p1, size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 24, margin = 5.0mm)
	savefig(savedir*"nonconvex_dist_sync.png")

	p2 = plot()
	for i in [1, 3, 6, 10]
		plot!(p2, 1 ./ etas, sticky_nonconvex_resp_vars[i, :, 1, 1]*times[i], label =  time_strings[i], marker = :circle, lw = 3, ms = 12)
	end
	plot!(p2, legend = false, xticks = (1 ./ etas, etas), xlabel = L"1/\eta", 
	ylabel = latexstring("Variances \$\\times \\, T\$"), title = "Sticky Coupling in Non-Convex Potential", ylims = (0, 1.5e5))
	plot!(p2, size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 22, margin = 5.0mm)
	savefig(savedir*"nonconvex_dist_sticky.png")
end

begin
	eta = 0.1
	sync_convex_av_traj = load(convex2d_data_dir*"sync_convex2D_resp_batch1_seed6924031646_eta$(eta)_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld", "R_avg")
	sticky_convex_av_traj = load(convex2d_data_dir*"sticky_convex2D_resp_batch1_seed6924031646_eta$(eta)_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld", "R_avg")
	plot()
	plot!(ts, sync_convex_av_traj[1, :, 1]/eta, label = "Synchronous Coupling Estimator", lw = 2)
	plot!(ts, sticky_convex_av_traj[1, :, 1]/eta, label = "Sticky Coupling Estimator", lw = 2)
	plot!(ts, repeat([0.5], length(ts)), label = "True Value", color = :black, lw = 4, ls = :dash)
	plot!(legend = :bottomright, xlabel = "Run Time", title = string("Response Coefficient, η = ", eta))
	plot!(xlims = (0, 5e4))
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 22, margin = 5.0mm)

	savefig(savedir*"convex_cross_convergence_eta$(eta).png")

end

begin
	eta = 0.01
	sync_convex_av_traj = load(convex2d_data_dir*"sync_convex2D_resp_batch1_seed6924031646_eta$(eta)_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld", "R_avg")
	sticky_convex_av_traj = load(convex2d_data_dir*"sticky_convex2D_resp_batch1_seed6924031646_eta$(eta)_beta1.0_burn100000.0_T1.0e6_dt0.005_interval20000_numruns50.jld", "R_avg")
	plot()
	plot!(ts, sync_convex_av_traj[1, :, 1]/eta, label = "Synchronous Coupling Estimator", lw = 2)
	plot!(ts, sticky_convex_av_traj[1, :, 1]/eta, label = "Sticky Coupling Estimator", lw = 2)
	plot!(ts, repeat([0.5], length(ts)), label = "True Value", color = :black, lw = 4, ls = :dash)
	plot!(legend = :bottomright, xlabel = "Run Time", title = string("Response Coefficient, η = ", eta))
	plot!(xlims = (0, 5e4))
	plot!(size = (1200, 800), legendfontsize = 24, tickfontsize = 18, titlefontsize = 26, guidefontsize = 22, margin = 5.0mm)

	savefig(savedir*"convex_cross_convergence_eta$(eta).png")

end