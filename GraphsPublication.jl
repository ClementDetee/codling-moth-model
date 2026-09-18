###############################################################################################################
#### GraphsPublication — Figures de publication : compétition larvaire et bifurcations                  ####
####                                                                                                     ####
#### Ce script génère les figures destinées à la publication, pour les 4 modèles de parasitoïdie.       ####
####                                                                                                     ####
#### Graphiques disponibles (activer/désactiver avec les flags true/false en tête de fichier) :         ####
####   • Comp_FluxLarve       : proportion de larves atteignant le stade nymphal selon leur densité     ####
####                            (courbe analytique — exp(-alphaL·L·tauL), sans simulation)              ####
####   • Allee_LogX_SansMorta : densité larvaire selon l'intensité Allee [deltaL = 0]   — axe X log    ####
####   • Allee_LogX_AvecMorta : densité larvaire selon l'intensité Allee [deltaL = 0.1] — axe X log    ####
####   • MortaL_LinX_Allee0   : densité larvaire selon la mortalité larvaire [Allee = 0]               ####
####   • MortaL_LinX_Allee10  : densité larvaire selon la mortalité larvaire [Allee = 10]              ####
####   • MortaL_LinX_Allee100 : densité larvaire selon la mortalité larvaire [Allee = 100]             ####
###############################################################################################################

using Revise
includet(joinpath(@__DIR__, "../fonction/MD_dynamique.jl"))

using .module_analyse
using DifferentialEquations
using Plots
using ColorSchemes
using Statistics

Plots.default(size=(750, 450))

# ────────────────────────────────────────────────
# Graphiques à afficher (true / false)
# ────────────────────────────────────────────────
Comp_FluxLarve        = true    # Flux larves→nymphes selon densité larvaire (analytique)
Allee_LogX_SansMorta  = false    # Densité larvaire selon Intensité Allee [deltaL = 0.0]  — X log
Allee_LogX_AvecMorta  = false    # Densité larvaire selon Intensité Allee [deltaL = 0.1]  — X log
MortaL_LinX_Allee0    = false    # Densité larvaire selon Mortalité larvaire [Allee = 0]   — X linéaire
MortaL_LinX_Allee10   = false    # Densité larvaire selon Mortalité larvaire [Allee = 10]  — X linéaire
MortaL_LinX_Allee100  = false    # Densité larvaire selon Mortalité larvaire [Allee = 100] — X linéaire

# ────────────────────────────────────────────────
# Modèles
# ────────────────────────────────────────────────
model0 = module_analyse.without_parasitoid_model
model1 = module_analyse.eggs_parasitoid_model
model2 = module_analyse.larval_parasitoid_model
model3 = module_analyse.nymphal_parasitoid_model

module_analyse.ALLEE_TYPE[] = "Test"

# ────────────────────────────────────────────────
# Paramètres de base
# ────────────────────────────────────────────────
DicoBase = Dict(
    "deltaE"  => 0.0,    "tauU"   => 4.0,
    "tauL"    => 16.0,   "tauM"   => 8.0,
    "tauJ"    => 16.0,   "alphaL" => 0.0001,
    "deltaL"  => 0.0,    "deltaN" => 0.0,
    "deltaA"  => 0.1,    "deltaP" => 0.1,
    "a"       => 0.01,   "q"      => 70,
    "H"       => 0.25,   "mi"     => 1,
    "allee"   => 0.0,
)

ordre_para = ("q","a","alphaL","deltaE","deltaL","deltaN","deltaA","deltaP",
              "tauU","tauL","tauM","tauJ","H","mi","allee")

function fct_lags(d)
    [d["tauU"], d["tauU"]+d["tauL"], d["tauL"],
     d["tauU"]+d["tauL"]+d["tauM"],
     d["tauL"]+d["tauM"], d["tauM"],
     d["tauJ"], d["tauJ"], d["tauU"], d["tauL"]]
end

# ────────────────────────────────────────────────
# Conditions initiales : 300 adultes, 0.1 parasitoïdes
# ────────────────────────────────────────────────
function fct_hini_0(d, type)
    type == "h" && return [0.0, 0.0, 0.0, 0.0, -d["tauL"]*d["alphaL"]]
    [0.0, 0.0, 0.0, 300.0, -d["tauL"]*(d["alphaL"]*0 + d["deltaL"])]
end

function fct_hini_1(d, type)
    type == "h" && return [0.0,0.0,0.0,0.0,0.0,0.0,
        -d["tauU"]*(d["a"]*0   + d["deltaE"]),
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"])]
    [0.0,0.0,0.0, 20.0, 0.0, 0.1,
        -d["tauU"]*(d["a"]*0.1 + d["deltaE"]),
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"])]
end

function fct_hini_2(d, type)
    type == "h" && return [0.0,0.0,0.0,0.0,0.0,0.0,
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"] + d["a"]*0),
        -d["tauJ"]*(d["alphaL"]*0 + d["deltaL"])]
    [0.0,0.0,0.0, 300.0, 1.0, 0.1,
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"] + d["a"]*0.1),
        -d["tauJ"]*(d["alphaL"]*0 + d["deltaL"])]
end

function fct_hini_3(d, type)
    type == "h" && return [0.0,0.0,0.0,0.0,0.0,0.0,
        -d["tauM"]*(d["a"]*0   + d["deltaN"]),
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"])]
    [0.0,0.0,0.0, 300.0, 0.0, 0.1,
        -d["tauM"]*(d["a"]*0.1 + d["deltaN"]),
        -d["tauL"]*(d["alphaL"]*0 + d["deltaL"])]
end

para_models = [model0, model1, model2, model3]
para_hinis  = [fct_hini_0, fct_hini_1, fct_hini_2, fct_hini_3]
para_nb     = [0, 1, 2, 3]
para_labels = ["Sans para.", "Para. œufs", "Para. larves", "Para. nymphes"]

# ────────────────────────────────────────────────
# Couleurs (cohérentes sur tous les graphiques)
# ────────────────────────────────────────────────
let pal = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions = collect(range(0, 1; length=12))
    indices   = [2, 4, 6, 10, 1, 3, 8, 11]
    cs        = [RGB{Float64}(get(pal, positions[i])) for i in indices]
    cs[3]     = RGB{Float64}(0.0, 0.55, 0.15)
    global colors_pub = cs[1:4]
end

# ────────────────────────────────────────────────
# Paramètres temporels
# ────────────────────────────────────────────────
t_fin    = 20000.0
t_mesure = 10000.0
tspan    = (0.0, t_fin)
pdt      = 0.01

seuil_para_extinction = 0.01

# ────────────────────────────────────────────────
# Simulation
# ────────────────────────────────────────────────
function run_sim(model, fct_hini, dico, tspan, t_mes; var_idx=2)
    p    = tuple([dico[k] for k in ordre_para]...)
    lags = fct_lags(dico)
    u0   = fct_hini(dico, "u")
    h0   = fct_hini(dico, "h")
    h_f(p, t) = h0
    prob = DDEProblem(model, u0, h_f, tspan, p; constant_lags=lags)
    sol  = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=50.0,
                 abstol=1e-6, reltol=1e-6, maxiters=Int(1e9))
    idx  = findall(x -> x >= tspan[2] - t_mes, sol.t)
    isempty(idx) && (idx = eachindex(sol.t))
    vals  = [max(0.0, sol[var_idx, i]) for i in idx]
    moy   = mean(vals); mn = minimum(vals); mx = maximum(vals)
    para_pres = false
    if model !== model0
        pfin      = [sol[6, i] for i in idx]
        para_pres = mean(pfin) > seuil_para_extinction
    end
    return moy, mn, mx, para_pres
end

# ────────────────────────────────────────────────
# Tracé avec segmentation pointillé/plein + ribbon
# – légende : trait plein coloré par modèle (toujours)
# – une entrée noire pointillée commune = "Absence de parasitoïde"
# ────────────────────────────────────────────────
const EPS_LOG = 0.1

function plot_segments!(plt, x, moy, mn, mx, nb_m, para_arr, c, lbl; lw=2)
    moy_c = max.(moy, EPS_LOG)
    mn_c  = max.(mn,  EPS_LOG)
    mx_c  = max.(mx,  EPS_LOG)

    if nb_m == 0
        Plots.plot!(plt, x, mx_c, fillrange=mn_c,
                    label="", color=c, fillalpha=0.13, linewidth=0)
        Plots.plot!(plt, x, moy_c,
                    label=lbl, color=c, linestyle=:solid, linewidth=lw)
        return
    end

    # Entrée de légende toujours en trait plein
    Plots.plot!(plt, [], [], label=lbl, color=c, linewidth=lw, linestyle=:solid)

    i = 1; n = length(x)
    while i <= n
        cur = para_arr[i]
        j   = i + 1
        while j <= n && para_arr[j] == cur; j += 1; end
        seg   = i:min(j, n)
        style = cur ? :solid : :dot
        Plots.plot!(plt, x[seg], mx_c[seg], fillrange=mn_c[seg],
                    label="", color=c, fillalpha=0.13, linewidth=0)
        Plots.plot!(plt, x[seg], moy_c[seg],
                    label="", color=c, linestyle=style, linewidth=lw)
        i = j
    end
end

# ────────────────────────────────────────────────
# Valeurs X
# ────────────────────────────────────────────────
valeurs_allee_pub  = 10 .^ range(log10(0.1), log10(5000.0), length=40)
valeurs_deltaL_pub = collect(range(0.0, 1.0, length=60))


########################################
# Flux larves→nymphes selon densité larvaire — analytique
########################################

if Comp_FluxLarve

    tauL_c   = DicoBase["tauL"]
    alphaL_c = DicoBase["alphaL"]
    L_range  = collect(range(0, 5000.0, length=500))

    plt_comp = Plots.plot(
        title          = "Proportion de larves atteignant le stade nymphal",
        xlabel         = "Densité larvaire (L)",
        ylabel         = "Proportion survivant à la compétition",
        grid           = false,
        titlefontsize  = 14,
        guidefontsize  = 12,
        tickfontsize   = 10,
        ylims          = (0.0, 1.05),
        legend         = false,
    )

    survie = exp.(-(alphaL_c .* L_range) .* tauL_c)
    Plots.plot!(plt_comp, L_range, survie, color=colors_pub[1], linewidth=2)

    display(plt_comp)
    println("✓ Comp_FluxLarve terminé")
end


########################################
# Densité larvaire selon Intensité Allee [deltaL = 0.0] — X log
########################################

if Allee_LogX_SansMorta

    plt_1 = Plots.plot(
        title          = "Densité larvaire selon Intensité Allee",
        xlabel         = "Intensité de l'effet Allee (log)",
        ylabel         = "Densité larvaire (log)",
        grid           = false,
        xscale         = :log10,
        yscale         = :log10,
        titlefontsize  = 9,
        guidefontsize  = 8,
        legendfontsize = 7,
        legend         = :topright,
    )
    Plots.plot!(plt_1, [], [], label="Absence de parasitoïde",
                color=:black, linestyle=:dot, linewidth=2)

    for (ip, (model, hini, nb_m)) in enumerate(zip(para_models, para_hinis, para_nb))
        moyennes = Float64[]; mins = Float64[]; maxs = Float64[]; para_arr = Bool[]
        for allee_val in valeurs_allee_pub
            dico           = copy(DicoBase)
            dico["allee"]  = allee_val
            dico["deltaL"] = 0.0
            moy, mn, mx, pp = run_sim(model, hini, dico, tspan, t_mesure)
            push!(moyennes, moy); push!(mins, mn); push!(maxs, mx); push!(para_arr, pp)
        end
        plot_segments!(plt_1, valeurs_allee_pub, moyennes, mins, maxs,
                       nb_m, para_arr, colors_pub[ip], para_labels[ip])
    end

    Plots.xlims!(plt_1, 0.1, 5000.0)
    display(plt_1)
    println("✓ Allee_LogX_SansMorta terminé")
end


########################################
# Densité larvaire selon Intensité Allee [deltaL = 0.1] — X log
########################################

if Allee_LogX_AvecMorta

    plt_2 = Plots.plot(
        title          = "Densité larvaire selon Intensité Allee  [deltaL = 0.1]",
        xlabel         = "Intensité de l'effet Allee (log)",
        ylabel         = "Densité larvaire (log)",
        grid           = false,
        xscale         = :log10,
        yscale         = :log10,
        titlefontsize  = 9,
        guidefontsize  = 8,
        legendfontsize = 7,
        legend         = :topright,
    )
    Plots.plot!(plt_2, [], [], label="Absence de parasitoïde",
                color=:black, linestyle=:dot, linewidth=2)

    for (ip, (model, hini, nb_m)) in enumerate(zip(para_models, para_hinis, para_nb))
        moyennes = Float64[]; mins = Float64[]; maxs = Float64[]; para_arr = Bool[]
        for allee_val in valeurs_allee_pub
            dico           = copy(DicoBase)
            dico["allee"]  = allee_val
            dico["deltaL"] = 0.1
            moy, mn, mx, pp = run_sim(model, hini, dico, tspan, t_mesure)
            push!(moyennes, moy); push!(mins, mn); push!(maxs, mx); push!(para_arr, pp)
        end
        plot_segments!(plt_2, valeurs_allee_pub, moyennes, mins, maxs,
                       nb_m, para_arr, colors_pub[ip], para_labels[ip])
    end

    Plots.xlims!(plt_2, 0.1, 5000.0)
    display(plt_2)
    println("✓ Allee_LogX_AvecMorta terminé")
end


########################################
# Densité larvaire selon Mortalité larvaire [Allee = 0] — X linéaire
########################################

if MortaL_LinX_Allee0

    plt_3 = Plots.plot(
        title          = "Densité larvaire selon Mortalité larvaire",
        xlabel         = "Mortalité larvaire (deltaL)",
        ylabel         = "Densité larvaire (log)",
        grid           = false,
        yscale         = :log10,
        titlefontsize  = 9,
        guidefontsize  = 8,
        legendfontsize = 7,
        legend         = :topright,
    )
    Plots.plot!(plt_3, [], [], label="Absence de parasitoïde",
                color=:black, linestyle=:dot, linewidth=2)

    for (ip, (model, hini, nb_m)) in enumerate(zip(para_models, para_hinis, para_nb))
        moyennes = Float64[]; mins = Float64[]; maxs = Float64[]; para_arr = Bool[]
        for deltaL_val in valeurs_deltaL_pub
            dico           = copy(DicoBase)
            dico["allee"]  = 0.0
            dico["deltaL"] = deltaL_val
            moy, mn, mx, pp = run_sim(model, hini, dico, tspan, t_mesure)
            push!(moyennes, moy); push!(mins, mn); push!(maxs, mx); push!(para_arr, pp)
        end
        plot_segments!(plt_3, valeurs_deltaL_pub, moyennes, mins, maxs,
                       nb_m, para_arr, colors_pub[ip], para_labels[ip])
    end

    Plots.xlims!(plt_3, 0.0, 1.0)
    display(plt_3)
    println("✓ MortaL_LinX_Allee0 terminé")
end


########################################
# Densité larvaire selon Mortalité larvaire [Allee = 10] — X linéaire
########################################

if MortaL_LinX_Allee10

    plt_4 = Plots.plot(
        title          = "Densité larvaire selon Mortalité larvaire  [Allee = 10]",
        xlabel         = "Mortalité larvaire (deltaL)",
        ylabel         = "Densité larvaire (log)",
        grid           = false,
        yscale         = :log10,
        titlefontsize  = 9,
        guidefontsize  = 8,
        legendfontsize = 7,
        legend         = :topright,
    )
    Plots.plot!(plt_4, [], [], label="Absence de parasitoïde",
                color=:black, linestyle=:dot, linewidth=2)

    for (ip, (model, hini, nb_m)) in enumerate(zip(para_models, para_hinis, para_nb))
        moyennes = Float64[]; mins = Float64[]; maxs = Float64[]; para_arr = Bool[]
        for deltaL_val in valeurs_deltaL_pub
            dico           = copy(DicoBase)
            dico["allee"]  = 10.0
            dico["deltaL"] = deltaL_val
            moy, mn, mx, pp = run_sim(model, hini, dico, tspan, t_mesure)
            push!(moyennes, moy); push!(mins, mn); push!(maxs, mx); push!(para_arr, pp)
        end
        plot_segments!(plt_4, valeurs_deltaL_pub, moyennes, mins, maxs,
                       nb_m, para_arr, colors_pub[ip], para_labels[ip])
    end

    Plots.xlims!(plt_4, 0.0, 1.0)
    display(plt_4)
    println("✓ MortaL_LinX_Allee10 terminé")
end


########################################
# Densité larvaire selon Mortalité larvaire [Allee = 100] — X linéaire
########################################

if MortaL_LinX_Allee100

    plt_5 = Plots.plot(
        title          = "Densité larvaire selon Mortalité larvaire  [Allee = 100]",
        xlabel         = "Mortalité larvaire (deltaL)",
        ylabel         = "Densité larvaire (log)",
        grid           = false,
        yscale         = :log10,
        titlefontsize  = 9,
        guidefontsize  = 8,
        legendfontsize = 7,
        legend         = :topright,
    )
    Plots.plot!(plt_5, [], [], label="Absence de parasitoïde",
                color=:black, linestyle=:dot, linewidth=2)

    for (ip, (model, hini, nb_m)) in enumerate(zip(para_models, para_hinis, para_nb))
        moyennes = Float64[]; mins = Float64[]; maxs = Float64[]; para_arr = Bool[]
        for deltaL_val in valeurs_deltaL_pub
            dico           = copy(DicoBase)
            dico["allee"]  = 100.0
            dico["deltaL"] = deltaL_val
            moy, mn, mx, pp = run_sim(model, hini, dico, tspan, t_mesure)
            push!(moyennes, moy); push!(mins, mn); push!(maxs, mx); push!(para_arr, pp)
        end
        plot_segments!(plt_5, valeurs_deltaL_pub, moyennes, mins, maxs,
                       nb_m, para_arr, colors_pub[ip], para_labels[ip])
    end

    Plots.xlims!(plt_5, 0.0, 1.0)
    display(plt_5)
    println("✓ MortaL_LinX_Allee100 terminé")
end
