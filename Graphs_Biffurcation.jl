###############################################################################################################
#### Graphs_Biffurcation — Diagrammes de bifurcation : effet Allee, mortalités et fécondité               ####
####                                                                                                     ####
#### Ce script génère des diagrammes de bifurcation explorant la densité de population à l'équilibre   ####
#### en fonction de différents paramètres, pour les 4 modèles de parasitoïdie.                          ####
####                                                                                                     ####
#### Graphiques disponibles (contrôlés par les flags true/false en tête de fichier) :                  ####
####   • Ggen       : bifurcation générique (X = Allee, mortalité ou fécondité q)                      ####
####   • Ggen_fusion: plusieurs mortalités côte à côte                                                  ####
####   • Gcomp      : compétition larvaire selon l'intensité Allee ou la mortalité                     ####
####   • Gq0adu     : fécondité effective q₀ selon l'intensité Allee (analytique)                      ####
####   • Gadu_q0    : densité adulte selon la fécondité effective (analytique)                         ####
####   • Gcrit      : frontière d'extinction (densité initiale vs intensité Allee critique)             ####
###############################################################################################################


## Graphique générique ##
afficher_Ggen       = true                   # true = activer le graphique générique
X_variable          = "morta"                # "morta" = mortalité en X  |  "allee" = allee en X  |  "q" = fécondité q en X (de 100 à 0)
morta_affiche       = "deltaL"                 # paramètre de mortalité : "deltaE", "deltaL", "deltaN", "deltaA"
para_type           = [3]                    # modèles : [0]=sans para  [1]=oeufs  [2]=larves  [3]=nymphes
variable_y_Ggen     = "lar"                  # variable en Y : "euf"=oeufs | "lar"=larves | "nym"=nymphes | "adu"=adultes | "q0"=fécondité effective q₀=q·N/(N+allee) | "q0xu4"=production totale d'oeufs q₀×N
ymax_Ggen           = Inf    #  [0,300]            # limite Y (Inf = auto)
xmax_Ggen           = 0.5 #  zoom X : Inf = auto | pour ÉTENDRE au-delà de 150 → modifier valeurs_allee_Ggen | xmax_Ggen sert uniquement à ZOOMER (réduire)
x_log               = false              # true = axe X en échelle logarithmique
afficher_para_Ggen  = false                  # true = afficher aussi la courbe parasitoïde           # bornes Y du graphique de pente, comme ymax_Ggen (ex: [0,100])

afficher_Ggen_fusion = false                  # true = générer un Ggen par mortalité et les fusionner côte à côte
mortas_fusion        = ["deltaE", "deltaL", "deltaN"]   # mortalités à fusionner
layout_fusion        = (1, 3)               # disposition : (lignes, colonnes)
afficher_Gcomp      = false                   # allee en X, Y = competition larvaire (alphaL*L), courbes = valeurs_morta_Ggen
afficher_comp_lar   = false                  # true = graphique compagnon : densité larvaire pour la 1ère valeur de la liste
afficher_comp_para  = false                   # true = courbe supplémentaire : mortalité compétition des larves parasitées (modèle larvaire uniquement)
valeurs_seuils_Ggen = [10]                    # valeurs d'allee utilisées comme courbes (si X_variable="morta")
valeurs_allee_Ggen  = x_log ? 10 .^ range(0.0, log10(Float64(xmax_Ggen)), length=60) :
                              range(0.0, Float64(xmax_Ggen), length=60)   # plage continue d'allee pour l'axe X    (si X_variable="allee")
valeurs_morta_Ggen  = [0]                # 0.25  valeurs de mortalité utilisées comme courbes (si X_variable="allee")
valeurs_q_Ggen      = range(70, 0.0, length=40)            # plage de fécondité q pour l'axe X, de 70 à 0 (si X_variable="q")
precision = 1
coef_directeur      = false                  # true = afficher aussi la pente des courbes (%) : -100=verticale ↓, 0=horizontale, 100=verticale ↑
ymax_pente          = [-100, 100] 
afficher_Gq0adu     = false                   # true = q0 vs intensité Allee, courbes = valeurs d'adultes N (analytique)
valeurs_adu_Gq0     = [10, 100, 500, 1000]     # densités adultes N utilisées comme courbes dans afficher_Gq0adu
afficher_Gadu_q0     = false                   # true = densité adulte en Y, q0 en X, courbes = intensités Allee (analytique)
valeurs_allee_Gadu   = [0, 10, 250, 5000]     # intensités Allee utilisées comme courbes dans afficher_Gadu_q0
afficher_Gcrit       = false                   # true = frontière d'extinction : uStart en X, intensité Allee critique en Y
valeurs_ustart_Gcrit = range(0, 500.0, length=20)    # plage de densités initiales (X)
allee_sweep_Gcrit    = range(0.0, 5000.0, length=60)  # plage Allee balayée pour trouver le seuil de crash (Y)
morta_Gcrit          = "deltaE"              # paramètre de mortalité varié entre les courbes de Gcrit ("deltaE","deltaL","deltaN","deltaA") ou "" pour une seule courbe
valeurs_morta_Gcrit  = [0.0, 0.1, 0.2]      # valeurs du paramètre morta_Gcrit utilisées comme courbes
# Sécurité : les paramètres de mortalité sont bornés entre 0 et 1, donc un zoom > 10 n'a pas de sens
xmax_Ggen = (X_variable == "morta" && (xmax_Ggen isa AbstractVector ? xmax_Ggen[end] : xmax_Ggen) > 10) ? 2 : xmax_Ggen

# ymax_Ggen / xmax_Ggen = valeur (borne haute, borne basse = 0) ou [bas, haut] (les deux bornes)
ylims_Ggen(y) = y isa AbstractVector ? (y[1], y[2]) : (0.0, y)
xlims_Ggen(x) = x isa AbstractVector ? (x[1], x[2]) : (0.0, x)


include(joinpath(@__DIR__, "../fonction/MD_dynamique.jl"))

using .module_analyse
using DifferentialEquations
import Plots
using Plots: plot, plot!, ylims!, hline!, xlims!, RGB, font
using Statistics: mean
using ColorSchemes
Plots.default(size=(560, 400))   # largeur réduite de 25% (600×400 → 450×400)

DicoDesParas = Dict(
    "deltaE" => 0.0 ,  "tauU" => 4.0,  "tauL" => 16.0, "tauM" => 8.0,
    "tauJ"   => 16.0, "alphaL" => 0.0001, "deltaL" => 0.0, "deltaN" => 0.0,
    "deltaA" => 0.1,  "deltaP" => 0.1, "a" => 0.01, "q" => 70,
    "H" => 0.25, "mi" => 1,
    "allee" => 0,
)   

t_fin    = 10000.0   # temps final de simulation
t_mesure = 5000.0   # fenêtre de mesure : bifurcation calculée sur [t_fin - t_mesure, t_fin]
tspan    = (0.0, t_fin)
pdt      = 0.001

module_analyse.ALLEE_TYPE[] = "Test"

ordre_para = ("q", "a", "alphaL", "deltaE", "deltaL", "deltaN", "deltaA", "deltaP",
              "tauU", "tauL", "tauM", "tauJ", "H", "mi", "allee")

function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
            (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
            (dico_para["tauL"] + dico_para["tauM"]), dico_para["tauM"], dico_para["tauJ"],
            dico_para["tauJ"], dico_para["tauU"], dico_para["tauL"]]
end

function fct_hini_0(dico_para,type) 
    if type == "h"
        hStart = 0.0
        hPara = 0.0
        return [hStart, hStart, hStart, hStart,- dico_para["tauL"] * dico_para["alphaL"]]
    else 
        uStart = 20
        uPara = 0
        return [0, 0, 0, uStart, - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    end 
end

function fct_hini_1(dico_para,type) 
    if type == "h"
        hStart = 0.0
        hPara = 0.0
        return [hStart, hStart, hStart, hStart, hPara, hPara, - dico_para["tauU"] * (dico_para["a"] * hPara + dico_para["deltaE"]), 
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    else 
        uStart = 300
        uPara = 0.1
        return [0, 0, 0, uStart, 0, uPara, - dico_para["tauU"] * (dico_para["a"] * uPara + dico_para["deltaE"]), 
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    end 
end 

function fct_hini_2(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0
        return [hStart, hStart, hStart, hStart, hPara, hPara, 
        - dico_para["tauL"] * ((hPara + hStart) * dico_para["alphaL"] + dico_para["deltaL"] + 
        dico_para["a"] * hPara),
        - dico_para["tauJ"] * ((hPara + hStart) * dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        uPara = 0.1
        uStart = 300
        return [0, 0, 0, uStart, 0, uPara, 
        - dico_para["tauL"] * ((0 + 0) * dico_para["alphaL"] + dico_para["deltaL"] + 
        dico_para["a"]* uPara),
        - dico_para["tauJ"] * ((0 + 0) * dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

function fct_hini_3(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0

        return [hStart, hStart, hStart, hStart, hPara, hPara, 
        - dico_para["tauM"] * (dico_para["a"] * hPara + dico_para["deltaN"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * hStart + dico_para["deltaL"])]
    else 
        uStart = 300
        uPara = 0.1
        return [0, 0, 0, uStart, 0, uPara, 
        - dico_para["tauM"] * (dico_para["a"] * uPara + dico_para["deltaN"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    end 
end

models         = [module_analyse.without_parasitoid_model,
                  module_analyse.eggs_parasitoid_model,
                  module_analyse.larval_parasitoid_model,
                  module_analyse.nymphal_parasitoid_model]
fct_hinis      = [fct_hini_0, fct_hini_1, fct_hini_2, fct_hini_3]
labels_modeles = ["Sans parasitoïde", "Para. oeufs", "Para. larves", "Para. nymphes"]
nb_modeles     = [0, 1, 2, 3]

palette   = ColorSchemes.colorschemes[:twelvebitrainbow]
positions = collect(range(0, 1; length=12))
indices   = [2, 4, 6, 10, 1, 3, 8, 11]
colors    = [get(palette, positions[i]) for i in indices]
colors[3] = RGB{Float64}(0.0, 0.55, 0.15)

seuil_para_extinction = 0.001

function run_simulation(model, fct_hini, dico, ordre_para, fct_lags,
                        tspan, t_fin, pdt, nb_modele, var_idx=2)
    p       = tuple([dico[k] for k in ordre_para]...)
    lags    = fct_lags(dico)
    u0      = fct_hini(dico, "u")
    h0      = fct_hini(dico, "h")
    h_fct(p, t) = h0

    prob = DDEProblem(model, u0, h_fct, tspan, p; constant_lags=lags)
    sol  = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=1.0,
                 abstol=1e-6, reltol=1e-6)

    idx_fin = findall(x -> x >= tspan[2] - t_fin, sol.t)
    isempty(idx_fin) && (idx_fin = eachindex(sol.t))
    if isempty(idx_fin)
        return 0.0, 0.0, 0.0, false, 0.0, 0.0, 0.0
    end
    larves_fin = [max(0.0, sol[var_idx, i]) for i in idx_fin]
    moy_larves = mean(larves_fin)
    min_larves = minimum(larves_fin)
    max_larves = maximum(larves_fin)

    moy_para = 0.0
    min_para = 0.0
    max_para = 0.0
    para_present = false
    if nb_modele != 0
        para_fin     = [sol[6, i] for i in idx_fin]
        moy_para     = mean(para_fin)
        min_para     = minimum(para_fin)
        max_para     = maximum(para_fin)
        para_present = moy_para > seuil_para_extinction
    end

    return moy_larves, min_larves, max_larves, para_present, moy_para, min_para, max_para
end

########################################
# GRAPHIQUE GÉNÉRIQUE
########################################

if afficher_Ggen

    morta_ranges = Dict("deltaE" => range(0.0, 1.0, length=60),
                        "deltaL" => range(0.0, 1.0, length=60),
                        "deltaN" => range(0.0, 1.0, length=60),
                        "deltaA" => range(0.0, 1.0, length=60))
    morta_labels = Dict("deltaE" => "Mortalité des œufs (deltaE)",
                        "deltaL" => "Mortalité larvaire (deltaL)",
                        "deltaN" => "Mortalité nymphale (deltaN)",
                        "deltaA" => "Mortalité adulte (deltaA)")
    para_models  = [module_analyse.without_parasitoid_model,
                    module_analyse.eggs_parasitoid_model,
                    module_analyse.larval_parasitoid_model,
                    module_analyse.nymphal_parasitoid_model]
    para_hinis   = [fct_hini_0, fct_hini_1, fct_hini_2, fct_hini_3]
    para_nb      = [0, 1, 2, 3]
    para_labels  = ["Sans para.", "Para. œufs", "Para. larves", "Para. nymphes"]
    para_styles  = [:solid, :solid, :solid, :solid]

    var_map_g    = Dict("euf" => 1, "lar" => 2, "nym" => 3, "adu" => 4)
    var_labels_g = Dict("euf" => "Densité des œufs", "lar" => "Densité larvaire",
                        "nym" => "Densité nymphale", "adu" => "Densité des adultes",
                        "q0"    => "Fécondité effective q₀ = q·N*/(N*+intensité Allee)",
                        "q0xu4" => "Production totale d'œufs q₀×N*")
    var_idx_g = variable_y_Ggen ∈ ("q0", "q0xu4") ? 4 : var_map_g[variable_y_Ggen]

    # Selon X_variable, on détermine ce qui va en X et ce qui varie entre courbes
    if X_variable == "morta"
        x_max_eff = xlims_Ggen(xmax_Ggen)[2] < Inf ? xlims_Ggen(xmax_Ggen)[2] : maximum(morta_ranges[morta_affiche])
        valeurs_x       = collect(range(0.0, x_max_eff, length=length(morta_ranges[morta_affiche])))
        valeurs_courbes = valeurs_seuils_Ggen
        xlabel_str      = morta_labels[morta_affiche]
        curve_label(v)  = "intensité Allee = $v"
        set_dico!(d, xval, cval) = (d[morta_affiche] = xval; d["allee"] = cval)
    elseif X_variable == "q"
        valeurs_x       = collect(valeurs_q_Ggen)
        x_max_eff       = abs(maximum(valeurs_q_Ggen) - minimum(valeurs_q_Ggen))
        xlabel_str      = "Fécondité q"
        if isempty(morta_affiche)
            valeurs_courbes = [0.0]
            curve_label(v)  = get(var_labels_g, variable_y_Ggen, variable_y_Ggen)
            set_dico!(d, xval, cval) = (d["q"] = xval)
        else
            valeurs_courbes = valeurs_morta_Ggen
            curve_label(v)  = "$(morta_affiche) = $v"
            set_dico!(d, xval, cval) = (d["q"] = xval; d[morta_affiche] = cval)
        end
    else  # "allee"
        x_max_eff = xlims_Ggen(xmax_Ggen)[2] < Inf ? xlims_Ggen(xmax_Ggen)[2] : maximum(valeurs_allee_Ggen)
        valeurs_x = collect(range(0.0, x_max_eff, length=length(valeurs_allee_Ggen)))
        xlabel_str      = "Intensité Allee"
        if isempty(morta_affiche)
            valeurs_courbes = [0.0]
            curve_label(v)  = get(var_labels_g, variable_y_Ggen, variable_y_Ggen)
            set_dico!(d, xval, cval) = (d["allee"] = xval)
        else
            valeurs_courbes = valeurs_morta_Ggen
            curve_label(v)  = "$(morta_affiche) = $v"
            set_dico!(d, xval, cval) = (d["allee"] = xval; d[morta_affiche] = cval)
        end
    end

    println("=== Graphique générique : $X_variable en X, para_type=$para_type ===")

    para_type_str = join([para_labels[pt + 1] for pt in para_type], ", ")
    plt_gen = plot(title  = "$(var_labels_g[variable_y_Ggen]) selon $xlabel_str\n$para_type_str",
                   xlabel = xlabel_str,
                   ylabel = var_labels_g[variable_y_Ggen],
                   grid   = false)

    courbes_data  = []

    for (ip, pt) in enumerate(para_type)
        model    = para_models[pt + 1]
        fct_hini = para_hinis[pt + 1]
        nb_m     = para_nb[pt + 1]
        pstyle   = para_styles[ip]

        for (idx_s, cval) in enumerate(valeurs_courbes)
            moyennes     = Float64[]
            mins         = Float64[]
            maxs         = Float64[]
            moyennes_p   = Float64[]
            mins_p       = Float64[]
            maxs_p       = Float64[]
            para_present = Bool[]

            for xval in valeurs_x
                dico = copy(DicoDesParas)
                set_dico!(dico, xval, cval)

                moy, mn, mx, para, moy_p, mn_p, mx_p = run_simulation(
                    model, fct_hini, dico, ordre_para,
                    fct_lags, tspan, t_mesure, pdt, nb_m, var_idx_g)
                if variable_y_Ggen == "q0"
                    moy = module_analyse._allee_q(dico["q"], moy, dico["allee"])
                    mn  = module_analyse._allee_q(dico["q"], mn,  dico["allee"])
                    mx  = module_analyse._allee_q(dico["q"], mx,  dico["allee"])
                elseif variable_y_Ggen == "q0xu4"
                    moy = module_analyse._allee_q(dico["q"], moy, dico["allee"]) * moy
                    mn  = module_analyse._allee_q(dico["q"], mn,  dico["allee"]) * mn
                    mx  = module_analyse._allee_q(dico["q"], mx,  dico["allee"]) * mx
                end
                push!(moyennes, moy); push!(mins, mn); push!(maxs, mx)
                push!(para_present, para)
                if afficher_para_Ggen
                    push!(moyennes_p, moy_p); push!(mins_p, mn_p); push!(maxs_p, mx_p)
                end
            end

            multi_para   = length(para_type) > 1
            multi_courbe = length(valeurs_courbes) > 1
            c   = (multi_para && !multi_courbe) ? colors[ip] : colors[idx_s]
            lbl = if multi_para && multi_courbe
                "$(para_labels[pt+1]) — $(curve_label(cval))"
            elseif multi_para
                para_labels[pt+1]
            elseif multi_courbe
                curve_label(cval)
            else
                curve_label(cval)
            end

            push!(courbes_data, (valeurs_x, moyennes, c, lbl, pstyle))

            plot!(plt_gen, [], [], color=c, linewidth=2, linestyle=pstyle, label=lbl)

            i = 1
            while i <= length(valeurs_x)
                present = para_present[i]
                j = i
                while j <= length(valeurs_x) && para_present[j] == present
                    j += 1
                end
                seg = i:min(j, length(valeurs_x))  # +1 point pour connecter au segment suivant
                plot!(plt_gen, valeurs_x[seg], moyennes[seg],
                      ribbon=(moyennes[seg] .- mins[seg], maxs[seg] .- moyennes[seg]),
                      label="", color=c, fillalpha=0.2, linewidth=2,
                      linestyle=(nb_m == 0 || present) ? pstyle : :dot)
                i = j
            end

            if afficher_para_Ggen
                moyennes_p_plot = copy(moyennes_p)
                last_para = findlast(v -> v > seuil_para_extinction, moyennes_p_plot)
                if !isnothing(last_para) && last_para < length(valeurs_x)
                    moyennes_p_plot[last_para + 1] = moyennes[last_para + 1]
                end
                plot!(plt_gen, valeurs_x, moyennes_p_plot,
                      label="", color=c, fillalpha=0.0, linewidth=2, linestyle=:dot)
            end
        end
    end

    multi_courbe_leg = length(valeurs_courbes) > 1
    if length(para_type) > 1 && multi_courbe_leg
        # plusieurs para_type ET plusieurs courbes → légende styles (couleur=morta, style=para_type)
        if any(pt -> pt != 0, para_type)
            plot!(plt_gen, [], [], linestyle=:solid, color=:black, linewidth=2, label="Hôte")
        end
        if afficher_para_Ggen
            plot!(plt_gen, [], [], linestyle=:dot, color=:grey, linewidth=2, label="Parasitoïde")
        end
        for (ip, pt) in enumerate(para_type)
            plot!(plt_gen, [], [], linestyle=para_styles[ip], color=:black, linewidth=2, label=para_labels[pt + 1])
        end
    elseif afficher_para_Ggen
        plot!(plt_gen, [], [], linestyle=:dot, color=:grey, linewidth=2, label="Parasitoïde (pointillés)")
    end

    ylims!(plt_gen, ylims_Ggen(ymax_Ggen)...)
    if X_variable == "q"
        Plots.xlims!(plt_gen, minimum(valeurs_q_Ggen), maximum(valeurs_q_Ggen))
        Plots.xflip!(plt_gen)
    else
        xlo, xhi = xlims_Ggen(xmax_Ggen)
        Plots.xlims!(plt_gen, x_log ? max(xlo, 1e0) : xlo, xhi)
    end
    x_log && Plots.plot!(plt_gen, xscale=:log10)
    display(plt_gen)
    println("Graphique générique terminé")

    # Graphique de la pente (%) des courbes : -100 = verticale ↓, 0 = horizontale, 100 = verticale ↑
    if coef_directeur
        all_y = vcat([my for (_, my, _, _, _) in courbes_data]...)
        y_global_min, y_global_max = extrema(all_y)
        y_range = y_global_max - y_global_min
        y_range = y_range > 0 ? y_range : 1.0
        x_range = x_max_eff > 0 ? x_max_eff : 1.0

        plt_pente = plot(title  = "Pente — $(var_labels_g[variable_y_Ggen]) selon $xlabel_str\n$para_type_str",
                         xlabel = xlabel_str,
                         ylabel = "Pente (%)",
                         grid   = false)

        for (vx, my, c, lbl, pstyle) in courbes_data
            x_mid      = (vx[1:end-1] .+ vx[2:end]) ./ 2
            slope_norm = (diff(my) ./ y_range) ./ (diff(vx) ./ x_range)
            pente_pct  = (2 / π) .* atan.(slope_norm) .* 100

            plot!(plt_pente, x_mid, pente_pct, color=c, linewidth=2, linestyle=pstyle, label=lbl)
        end

        Plots.hline!(plt_pente, [0.0], color=:black, linestyle=:dash, label="")
        ylims!(plt_pente, ylims_Ggen(ymax_pente)...)
        xlo_p, xhi_p = xlims_Ggen(xmax_Ggen)
        Plots.xlims!(plt_pente, x_log ? max(xlo_p, 1e-3) : xlo_p, xhi_p)
        x_log && Plots.plot!(plt_pente, xscale=:log10)
        display(plt_pente)
        println("Graphique pente terminé")
    end

end # afficher_Ggen

########################################
# FUSION Ggen — plusieurs mortalités côte à côte
########################################

if afficher_Ggen_fusion

    morta_labels_f = Dict("deltaE" => "Mortalité des œufs (deltaE)",
                          "deltaL" => "Mortalité larvaire (deltaL)",
                          "deltaN" => "Mortalité nymphale (deltaN)",
                          "deltaA" => "Mortalité adulte (deltaA)")
    var_map_f    = Dict("euf" => 1, "lar" => 2, "nym" => 3, "adu" => 4)
    var_labels_f = Dict("euf" => "Densité des œufs", "lar" => "Densité larvaire",
                        "nym" => "Densité nymphale", "adu" => "Densité des adultes")
    var_idx_f    = var_map_f[variable_y_Ggen]
    ylabel_f     = var_labels_f[variable_y_Ggen]

    # Modèles à afficher dans la fusion (reprend para_type du Ggen principal)
    para_models_f = [module_analyse.without_parasitoid_model,
                     module_analyse.eggs_parasitoid_model,
                     module_analyse.larval_parasitoid_model,
                     module_analyse.nymphal_parasitoid_model]
    para_hinis_f  = [fct_hini_0, fct_hini_1, fct_hini_2, fct_hini_3]
    para_labels_f = ["Sans para.", "Para. œufs", "Para. larves", "Para. nymphes"]

    # X = intensité Allee, courbes = valeurs_morta_Ggen
    x_max_eff_f = xlims_Ggen(xmax_Ggen)[2] < Inf ? xlims_Ggen(xmax_Ggen)[2] : maximum(valeurs_allee_Ggen)
    valeurs_x_f = collect(range(0.0, x_max_eff_f, length=length(valeurs_allee_Ggen)))

    plts_fusion  = []
    n_total      = length(mortas_fusion)
    idx_milieu   = ceil(Int, n_total / 2)   # subplot du milieu → xlabel
    morta_titres = Dict("deltaE" => "Mortalité des œufs",
                        "deltaL" => "Mortalité larvaire",
                        "deltaN" => "Mortalité nymphale",
                        "deltaA" => "Mortalité adulte")

    for (idx_f, morta) in enumerate(mortas_fusion)

        show_ylabel = idx_f == 1
        show_xlabel = idx_f == idx_milieu

        plt_f = plot(title          = morta_titres[morta],
                     xlabel         = show_xlabel ? "Intensité de l'effet Allee" : "",
                     ylabel         = show_ylabel ? ylabel_f : "",
                     grid           = false,
                     titlefontsize  = 10,
                     guidefontsize  = 9,
                     legendfontsize = 6,
                     legend         = :topright)

        for (idx_m, pt) in enumerate(para_type)
            model = para_models_f[pt + 1]
            hini  = para_hinis_f[pt + 1]
            lbl_m = para_labels_f[pt + 1]

            for (idx_c, morta_val) in enumerate(valeurs_morta_Ggen)
                moyennes = Float64[]; mins = Float64[]; maxs = Float64[]

                for allee_val in valeurs_x_f
                    dico           = copy(DicoDesParas)
                    dico["allee"]  = allee_val
                    dico[morta]    = morta_val
                    lags = fct_lags(dico)
                    p    = tuple([dico[k] for k in ordre_para]...)
                    u0   = hini(dico, "u"); h0 = hini(dico, "h")
                    prob = DDEProblem(model, u0, (p2,t)->h0, (0.0,t_fin), p; constant_lags=lags)
                    sol  = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=50.0,
                                 abstol=1e-6, reltol=1e-6, maxiters=Int(1e9))
                    idx  = findall(x -> x >= t_fin - t_mesure, sol.t)
                    isempty(idx) && (idx = eachindex(sol.t))
                    vals = [max(0.0, sol[var_idx_f, i]) for i in idx]
                    push!(moyennes, mean(vals)); push!(mins, minimum(vals)); push!(maxs, maximum(vals))
                end

                lbl = length(para_type) > 1 ? "$lbl_m — $morta = $morta_val" : "$morta = $morta_val"
                plot!(plt_f, valeurs_x_f, moyennes,
                      ribbon=(moyennes .- mins, maxs .- moyennes),
                      label=lbl, color=colors[idx_c],
                      linewidth=2, fillalpha=0.15)
            end
        end

        push!(plts_fusion, plt_f)
        println("✓ Fusion — $morta terminé")
    end

    n_cols  = layout_fusion[2]
    n_rows  = layout_fusion[1]
    titre_global = "$ylabel_f selon l'intensité de l'effet Allee\npour différentes valeurs de mortalité"
    display(plot(plts_fusion..., layout=layout_fusion, size=(550*n_cols, 450*n_rows),
                 plot_title=titre_global, plot_titlefontsize=11))

end # afficher_Ggen_fusion

########################################
# GRAPHIQUE COMPÉTITION : allee en X
# Y = compétition larvaire (alphaL * L_moy)
# courbes = valeurs_morta_Ggen
########################################

if afficher_Gcomp

    para_models_c = [module_analyse.without_parasitoid_model,
                     module_analyse.eggs_parasitoid_model,
                     module_analyse.larval_parasitoid_model,
                     module_analyse.nymphal_parasitoid_model]
    para_hinis_c  = [fct_hini_0, fct_hini_1, fct_hini_2, fct_hini_3]
    para_nb_c     = [0, 1, 2, 3]
    para_labels_c = ["Sans para.", "Para. œufs", "Para. larves", "Para. nymphes"]

    # X_variable contrôle l'axe X du Gcomp comme pour le Ggen
    if X_variable == "morta"
        x_max_eff_c  = xlims_Ggen(xmax_Ggen)[2] < Inf ? xlims_Ggen(xmax_Ggen)[2] : maximum(morta_ranges[morta_affiche])
        valeurs_x_c  = collect(range(0.0, x_max_eff_c, length=length(morta_ranges[morta_affiche])))
        valeurs_cbc  = valeurs_seuils_Ggen
        xlabel_c     = morta_labels[morta_affiche]
        clabel_c(v)  = "intensité Allee = $v"
        setdico_c!(d, xv, cv) = (d[morta_affiche] = xv; d["allee"] = cv)
    else  # "allee"
        x_max_eff_c  = xlims_Ggen(xmax_Ggen)[2] < Inf ? xlims_Ggen(xmax_Ggen)[2] : maximum(valeurs_allee_Ggen)
        valeurs_x_c  = collect(range(0.0, x_max_eff_c, length=length(valeurs_allee_Ggen)))
        xlabel_c     = "Intensité Allee"
        if isempty(morta_affiche)
            valeurs_cbc  = [0.0]
            clabel_c(v)  = "Nb de larves mourant de la compétition"
            setdico_c!(d, xv, cv) = (d["allee"] = xv)
        else
            valeurs_cbc  = valeurs_morta_Ggen
            clabel_c(v)  = "$(morta_affiche) = $v"
            setdico_c!(d, xv, cv) = (d["allee"] = xv; d[morta_affiche] = cv)
        end
    end

    println("=== Graphique competition : $X_variable en X, courbes = $(X_variable == "allee" ? morta_affiche : "allee") ===")

    para_type_str_c = join([para_labels_c[pt + 1] for pt in para_type], ", ")
    plt_comp = plot(title  = "Compétition larvaire selon $xlabel_c\n$para_type_str_c",
                    xlabel = xlabel_c,
                    ylabel = "Compétition larvaire (alphaL * L²)",
                    grid   = false)

    alphaL_val  = DicoDesParas["alphaL"]
    lar_moy_ref = Float64[]
    lar_min_ref = Float64[]
    lar_max_ref = Float64[]

    for (ip, pt) in enumerate(para_type)
        model    = para_models_c[pt + 1]
        fct_hini = para_hinis_c[pt + 1]
        nb_m     = para_nb_c[pt + 1]

        pstyle_c = para_styles[ip]

        for (idx_s, cval) in enumerate(valeurs_cbc)
            moyennes     = Float64[]
            mins         = Float64[]
            maxs         = Float64[]
            moy_lp_arr   = Float64[]
            para_present = Bool[]

            for xval in valeurs_x_c
                dico = copy(DicoDesParas)
                setdico_c!(dico, xval, cval)

                moy, mn, mx, para, _, _, _ = run_simulation(
                    model, fct_hini, dico, ordre_para,
                    fct_lags, tspan, t_mesure, pdt, nb_m, 2)
                push!(moyennes, alphaL_val * moy^2)
                push!(mins,     alphaL_val * mn^2)
                push!(maxs,     alphaL_val * mx^2)
                push!(para_present, para)
                if idx_s == 1
                    push!(lar_moy_ref, moy)
                    push!(lar_min_ref, mn)
                    push!(lar_max_ref, mx)
                end
                if afficher_comp_para && nb_m == 2
                    moy_lp, _, _, _, _, _, _ = run_simulation(
                        model, fct_hini, dico, ordre_para,
                        fct_lags, tspan, t_mesure, pdt, nb_m, 5)
                    push!(moy_lp_arr, alphaL_val * (moy + moy_lp) * moy_lp)
                end
            end

            c   = colors[idx_s]
            lbl = length(para_type) > 1 ?
                  "$(para_labels_c[pt+1]) — $(clabel_c(cval))" :
                  clabel_c(cval)

            plot!(plt_comp, valeurs_x_c, moyennes,
                  ribbon=(moyennes .- mins, maxs .- moyennes),
                  label=lbl, color=c, fillalpha=0.2, linewidth=2, linestyle=pstyle_c)

            if afficher_comp_para && nb_m == 2 && !isempty(moy_lp_arr)
                lbl_lp = (length(para_type) > 1 ? "$(para_labels_c[pt+1]) — " : "") *
                         "Lp parasitées" *
                         (length(valeurs_cbc) > 1 ? " — $(clabel_c(cval))" : "")
                plot!(plt_comp, valeurs_x_c, moy_lp_arr,
                      label=lbl_lp, color=c, linewidth=2, linestyle=:dash)
            end
        end
    end

    ylims!(plt_comp, ylims_Ggen(ymax_Ggen)...)
    Plots.xlims!(plt_comp, xlims_Ggen(xmax_Ggen)...)
    display(plt_comp)
    println("Graphique competition termine")

    if afficher_comp_lar
        cval_ref = valeurs_cbc[1]
        lbl_ref  = clabel_c(cval_ref)
        plt_lar = plot(title  = "Densité larvaire selon $xlabel_c ($lbl_ref)",
                       xlabel = xlabel_c,
                       ylabel = "Densité larvaire",
                       grid   = false)
        plot!(plt_lar, valeurs_x_c, lar_moy_ref,
              ribbon=(lar_moy_ref .- lar_min_ref, lar_max_ref .- lar_moy_ref),
              label=lbl_ref, color=colors[1], fillalpha=0.2, linewidth=2)
        ylims!(plt_lar, ylims_Ggen(ymax_Ggen)...)
        Plots.xlims!(plt_lar, xlims_Ggen(xmax_Ggen)...)
        display(plt_lar)
        println("Graphique comp_lar terminé")
    end

end # afficher_Gcomp

########################################
# GRAPHIQUE q0 vs densité adulte (courbe analytique)
# X = N (densité adulte, de 0 à xmax_Ggen)
# Y = q0 = q·N/(N+seuil)
# courbes = valeurs_seuils_Ggen (différentes intensités Allee)
# aucune simulation — formule directe
########################################

if afficher_Gq0adu

    q_val    = DicoDesParas["q"]
    A_vals   = collect(range(0.0, xlims_Ggen(xmax_Ggen)[2], length=300))
    n_adu    = length(valeurs_adu_Gq0)
    palette_adu = ColorSchemes.colorschemes[:twelvebitrainbow]
    colors_adu  = [get(palette_adu, p) for p in range(0.1, 0.9; length=max(n_adu, 2))]

    plt_q0adu = plot(title  = "Fécondité effective selon l'intensité Allee",
                     xlabel = "Intensité Allee",
                     ylabel = "Fécondité des adultes au temps t",
                     grid   = false)

    for (idx_s, N) in enumerate(valeurs_adu_Gq0)
        q0_vals = [iszero(A) ? q_val : q_val * N / (N + A) for A in A_vals]
        c   = colors_adu[idx_s]
        lbl = "Densité d'adultes = $N"
        plot!(plt_q0adu, A_vals, q0_vals, color=c, linewidth=2, label=lbl)
    end

    ylims!(plt_q0adu, 0.0, q_val * 1.05)
    Plots.xlims!(plt_q0adu, xlims_Ggen(xmax_Ggen)...)
    display(plt_q0adu)
    println("Graphique q0 vs intensité Allee terminé")

end # afficher_Gq0adu

########################################
# GRAPHIQUE q0 en X, densité adulte en Y (courbe analytique)
# X = fécondité effective q0 (0 → q)
# Y = densité adulte N
# courbes = intensités Allee
########################################

if afficher_Gadu_q0

    q_val2   = DicoDesParas["q"]
    q0_vals  = collect(range(0.0, q_val2 * 0.999, length=300))
    n_allee  = length(valeurs_allee_Gadu)
    palette_a2  = ColorSchemes.colorschemes[:twelvebitrainbow]
    colors_a2   = [get(palette_a2, p) for p in range(0.1, 0.9; length=max(n_allee, 2))]

    plt_adu_q0 = plot(title  = "Fécondité effective selon la densité adulte",
                      xlabel = "Densité d'adultes N",
                      ylabel = "Fécondité des adultes au temps t",
                      grid   = false)

    N_range = collect(range(0.0, 1000.0, length=300))

    for (idx_a, allee_val) in enumerate(valeurs_allee_Gadu)
        if iszero(allee_val)
            q0_line = fill(q_val2, length(N_range))
        else
            q0_line = [q_val2 * N / (N + allee_val) for N in N_range]
        end
        plot!(plt_adu_q0, N_range, q0_line,
              color=colors_a2[idx_a], linewidth=2, label="Allee = $allee_val")
    end

    Plots.xlims!(plt_adu_q0, 0.0, 1000.0)
    Plots.ylims!(plt_adu_q0, 0.0, q_val2 * 1.05)
    display(plt_adu_q0)
    println("Graphique densité adulte vs q0 terminé")

end # afficher_Gadu_q0

########################################
# GRAPHIQUE CRITIQUE : frontière d'extinction
# Y = uStart (densité adulte initiale)
# X = intensité Allee critique (premier seuil où la pop s'effondre)
# courbes = valeurs_morta_Ggen
########################################

if afficher_Gcrit

    println("=== Graphique critique : frontière d'extinction uStart vs Allee ===")

    plt_crit = plot(title  = "Frontière d'extinction",
                    xlabel = "Densité adulte initiale (uStart)",
                    ylabel = "Intensité Allee critique",
                    xlims  = (0, 500),
                    ylims  = (0, 5000),
                    grid   = false)

    valeurs_crit = isempty(morta_Gcrit) ? [NaN] : valeurs_morta_Gcrit

    for (idx_s, cval) in enumerate(valeurs_crit)
        A_crit_vals  = Float64[]
        ustart_valid = Float64[]

        for uStart_val in valeurs_ustart_Gcrit
            fct_hini_ust = (dico_para, type) -> begin
                if type == "h"
                    return fct_hini_0(dico_para, "h")
                else
                    u0 = fct_hini_0(dico_para, "u")
                    u0[4] = uStart_val
                    return u0
                end
            end

            survived_prev = true
            A_crit = NaN
            for allee_val in allee_sweep_Gcrit
                dico = copy(DicoDesParas)
                dico["allee"] = allee_val
                if !isempty(morta_Gcrit)
                    dico[morta_Gcrit] = cval
                end

                moy, _, _, _, _, _, _ = run_simulation(
                    module_analyse.without_parasitoid_model,
                    fct_hini_ust, dico, ordre_para,
                    fct_lags, tspan, t_mesure, pdt, 0, 4)

                survives = moy > seuil_para_extinction
                if !survives && survived_prev
                    A_crit = allee_val
                    break
                end
                survived_prev = survives
            end

            println("  uStart=$uStart_val → A_crit=$A_crit")
            if !isnan(A_crit)
                push!(A_crit_vals, A_crit)
                push!(ustart_valid, uStart_val)
            end
        end

        c   = colors[idx_s]
        lbl = isempty(morta_Gcrit) ? "Frontière d'extinction" : "$(morta_Gcrit) = $cval"
        plot!(plt_crit, ustart_valid, A_crit_vals,
              color=c, linewidth=2, label=lbl)
    end

    display(plt_crit)
    println("Graphique critique terminé")

end # afficher_Gcrit
