@testset "bayes.jl" begin
    g = BlockMaximaGrid((1,2), 1, [1 0; 1 0])
    c = Chain(g, 4, 2, [0.2, 0.15, 0.04], regional = Bool[1, 1, 0])

    @testset "setinitialvalues!(chain)" begin
        RegionalExtremes.setinitialvalues!(c)

        # Initial values for θ are right
        @test c.θ[1, :, :] ≈ [0.0836 0.0836; -0.3266 -0.3266; 0.0 0.0] atol = 0.001

        # Initial values for κ are right
        @test c.κ[1, :] ≈ [1.0, 1.0] atol = 0.0001

        # Initial values for θlogpdf are right
        @test c.θlogpdf[1, :, :] ≈ [-2.8238 -2.8238; -2.8238 -2.8238; -1.9048 -1.9048] atol = 0.001
    end
    @testset "copylast!(chain, i)" begin
        RegionalExtremes.copylast!(c, 2)

        # last θ iteration is copied
        @test c.θ[1, :, :] == c.θ[2, :, :]

        # last κ iteration is copied
        @test c.κ[1, :] == c.κ[2, :]

        # last θlogpdf iteration is copied
        @test c.θlogpdf[1, :, :] == c.θlogpdf[2, :, :]
    end
    @testset "candidates(chain, i, subset, param)" begin
        candidates = RegionalExtremes.candidates(c, 2, [1], 1)

        # Only param was changed
        @test candidates[1, 1] != c.θ[2, 1, 1]
        @test candidates[2, 1] == c.θ[2, 2, 1]
        @test candidates[3, 1] == c.θ[2, 3, 1]
    end
    @testset "MHupdatechain!(chain, i, subset, param, c, logpdfP)" begin
        # With logpdf 0, all the candidates are accepted
        c.θlogpdf[2, 1, :] = [0, 0]
        candidates = [Inf Inf; Inf Inf; Inf Inf]
        lodpdf = [0.0, 0.0]
        RegionalExtremes.MHupdatechain!(c, 2, [1,2], 1, candidates, lodpdf)

        @test c.θacc[2, 1, :] == Bool[1, 1]
        @test c.θ[2, 1, :] == candidates[1, :]
        @test c.θlogpdf[2, 1, :] == lodpdf

        # With logpdf -Inf, all the candidates are rejected
        c.θlogpdf[2, 1, :] = [0, 0]
        candidates = [-Inf -Inf; -Inf -Inf; -Inf -Inf]
        lodpdf = [-Inf, -Inf]
        RegionalExtremes.MHupdatechain!(c, 2, [1,2], 1, candidates, lodpdf)

        @test c.θacc[2, 1, :] == Bool[0, 0]
        @test c.θ[2, 1, :] != candidates[1, :]
        @test c.θlogpdf[2, 1, :] != lodpdf
    end
    @testset "MHlocal!(chain, i, p)" begin
        currenti = 2
        RegionalExtremes.copylast!(c, currenti)

        allbutcurrenti = [1, 3, 4]
        allbutcurrentp = [2 , 3]

        initialθ = copy(c.θ)
        initialθacc = copy(c.θacc)
        initialθlogpdf = copy(c.θlogpdf)

        RegionalExtremes.MHlocal!(c, 2, 1)

        # All iterations but current did not change
        @test initialθ[allbutcurrenti, :, :] == c.θ[allbutcurrenti, :, :]
        @test initialθacc[allbutcurrenti, :, :] == c.θacc[allbutcurrenti, :, :]
        @test initialθlogpdf[allbutcurrenti, :, :] == c.θlogpdf[allbutcurrenti, :, :]

        # All parameters but current did not change
        @test initialθ[currenti, allbutcurrentp, :] == c.θ[currenti, allbutcurrentp, :]
        @test initialθacc[currenti, allbutcurrentp, :] == c.θacc[currenti, allbutcurrentp, :]
        @test initialθlogpdf[currenti, allbutcurrentp, :] == c.θlogpdf[currenti, allbutcurrentp, :]
    end
    @testset "MHregional!(chain, i, p)" begin
        currenti = 2
        RegionalExtremes.copylast!(c, currenti)

        allbutcurrenti = [1, 3, 4]
        allbutcurrentp = [2 , 3]

        initialθ = copy(c.θ)
        initialθacc = copy(c.θacc)
        initialθlogpdf = copy(c.θlogpdf)

        RegionalExtremes.MHregional!(c, 2, 1)

        # All iterations but current did not change
        @test initialθ[allbutcurrenti, :, :] == c.θ[allbutcurrenti, :, :]
        @test initialθacc[allbutcurrenti, :, :] == c.θacc[allbutcurrenti, :, :]
        @test initialθlogpdf[allbutcurrenti, :, :] == c.θlogpdf[allbutcurrenti, :, :]

        # All parameters but current did not change
        @test initialθ[currenti, allbutcurrentp, :] == c.θ[currenti, allbutcurrentp, :]
        @test initialθacc[currenti, allbutcurrentp, :] == c.θacc[currenti, allbutcurrentp, :]
        @test initialθlogpdf[currenti, allbutcurrentp, :] == c.θlogpdf[currenti, allbutcurrentp, :]
    end
    @testset "MHrealisation!(chain, i)" begin
        currenti = 2
        RegionalExtremes.copylast!(c, currenti)

        allbutcurrenti = [1, 3, 4]

        initialθ = copy(c.θ)
        initialθacc = copy(c.θacc)
        initialθlogpdf = copy(c.θlogpdf)

        RegionalExtremes.MHrealisation!(c, 2)

        # All iterations but current did not change
        @test initialθ[allbutcurrenti, :, :] == c.θ[allbutcurrenti, :, :]
        @test initialθacc[allbutcurrenti, :, :] == c.θacc[allbutcurrenti, :, :]
        @test initialθlogpdf[allbutcurrenti, :, :] == c.θlogpdf[allbutcurrenti, :, :]
    end
    @testset "κrealisation!(chain, i)" begin
        currenti = 2
        RegionalExtremes.copylast!(c, currenti)

        allbutcurrenti = [1, 3, 4]

        initialκ = copy(c.κ)

        RegionalExtremes.κrealisation!(c, 2)

        # All iterations but current did not change
        @test initialκ[allbutcurrenti, :] == c.κ[allbutcurrenti, :]

        # Current iteration changed
        @test initialκ[currenti, :] != c.κ[currenti, :]
    end
    @testset "run!(chain)" begin
        # Stationary
        κs = [50.0, 50.0, 3000.0]
        g, θ, x = genGrid((10, 12), 1, κs..., 300)
        c = Chain(g, 1000, 500, [0.2, 0.15, 0.15], regional = [true, true, false])

        RegionalExtremes.run!(c)

        # GEV params land in quantile
        np = sum([length(v) for v in c.index])
        n = size(c.g.y, 1)
        for p in 1:np
            inrange = zeros(Bool, n)
            for i in 1:n
                q1 = quantile!(c.θ[c.burnin:end, p, i], 0.025)
                q2 = quantile!(c.θ[c.burnin:end, p, i], 0.975)

                inrange[i] = q1 <= θ[i, p] && θ[i, p] <= q2
            end
            @test mean(inrange) >= 0.9
        end

        # κ land in quantile
        for p in 1:size(c.κ, 2)
            q1 = quantile!(c.κ[c.burnin:end, p], 0.025)
            q2 = quantile!(c.κ[c.burnin:end, p], 0.975)
            @test q1 <= κs[p]
            @test κs[p] <= q2
        end

        # Non-stationary
        κs = [50.0, 10.0, 50.0, 3000.0]
        g, θ, x = genGrid((10, 12), 1, κs[1], κs[3], κs[4], 300, evμ = [(RegionalExtremes.rdm, κs[2])])
        c = Chain(g, 1000, 500, [0.2, 0.3, 0.15, 0.15], regional = Bool[1, 1, 1, 0], μcov = x[1])

        RegionalExtremes.run!(c)

        # GEV Values land in quantile
        np = sum([length(v) for v in c.index])
        n = size(c.g.y, 1)
        for p in 1:np
            inrange = zeros(Bool, n)
            for i in 1:n
                q1 = quantile!(c.θ[c.burnin:end, p, i], 0.025)
                q2 = quantile!(c.θ[c.burnin:end, p, i], 0.975)

                inrange[i] = q1 <= θ[i, p] && θ[i, p] <= q2
            end
            @test mean(inrange) >= 0.9
        end

        # κ values land in quantile
        for p in 1:size(c.κ, 2)
            q1 = quantile!(c.κ[c.burnin:end, p], 0.025)
            q2 = quantile!(c.κ[c.burnin:end, p], 0.975)
            @test q1 <= κs[p]
            @test κs[p] <= q2
        end
    end
end
