using Distributions
using BAT
using ValueShapes
using IntervalSets


@testset "bat_integrate" begin
    true_log_int(N::Int64) = 0.0
    @testset "normal" begin
        let N=4, s=1, dist=MvNormal(N, s),
            n_samples=10^5, n_chains=5,
            prior_min=-10., prior_max=10.,
            prior_norm = N*log(prior_max-prior_min)

            log_likelihood = params -> LogDVal(log(f(params.a)))
            prior = NamedTupleDist(a = [[prior_min .. prior_max for i in 1:N]...],)
            posterior = PosteriorDensity(log_likelihood, prior)

            f(x::AbstractVector; dist=dist) = pdf(dist, x)

            init = MCMCInitStrategy(
                init_tries_per_chain=80..1280,
                max_nsamples_init=2500,
                max_nsteps_init=2500, max_time_init=1800
            )

            algorithm=MetropolisHastings()

            samples_mcmc, chains_mcmc = bat_sample(posterior, (n_samples, n_chains), algorithm, init=init)

            hmi_data_mcmc = BAT.HMIData(unshaped.(samples_mcmc))
            BAT.hm_integrate!(hmi_data_mcmc)

            int_mcmc = hmi_data_mcmc.integralestimates["cov. weighted result"].final.estimate * exp(prior_norm)
            unc_mcmc = hmi_data_mcmc.integralestimates["cov. weighted result"].final.uncertainty * exp(prior_norm)

            logint_truth = true_log_int(N)

            @test isapprox(int_mcmc, 1, atol=unc_mcmc)
        end
    end
end
