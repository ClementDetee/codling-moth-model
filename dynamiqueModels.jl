###############################################################################################################
#### dynamiqueModels — Visualisation des dynamiques temporelles                                          ####
####                                                                                                     ####
#### Ce script simule et trace les trajectoires de densité au cours du temps pour les 4 modèles :       ####
####   sans parasitoïde, parasitoïde des œufs, des larves ou des nymphes.                               ####
####                                                                                                     ####
#### Modes disponibles :                                                                                 ####
####   • Normal     : une simulation par modèle sélectionné                                             ####
####   • Stade      : décomposition par stade de développement (oeufs / larves / nymphes / adultes)     ####
####   • Multi-courbes : un modèle, N valeurs d'un même paramètre (Allee ou mortalité)                  ####
####                                                                                                     ####
#### Graphiques auxiliaires optionnels :                                                                 ####
####   • Compétition larvaire, fécondité effective q0, taux de croissance per capita                    ####
###############################################################################################################

using Revise
includet(joinpath(@__DIR__, "../fonction/MD_dynamique.jl"))

using .module_analyse
using DifferentialEquations
using Plots
Plots.default(size=(450, 400))   
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using ColorSchemes



#########################
## Créaction du modèle ##
#########################

model0 = module_analyse.without_parasitoid_model
model1 = module_analyse.eggs_parasitoid_model
model2 = module_analyse.larval_parasitoid_model
model3 = module_analyse.nymphal_parasitoid_model

####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
DicoDesParasGlobalNul = Dict(
    "deltaE" => 0.0,
    "tauU" => 4.0,
    "tauL" => 16.0,
    "tauM" => 8.0,
    "tauJ" => 16.0,
    "alphaL" => 0.0001,
    "deltaL" => 0.25,
    "deltaN" => 0.0,
    "deltaA" => 0.1,
    "deltaP" =>  0.1,
    "a" => 0.01,
    "q" => 70,
    "H" => 0.25,
    "mi" => 1,
    "allee" => 10,
)
##Palette de comande##

modeles_affiches      = [4]    # ex. [1] sans parasitoïde, [1,2,3,4] tous
parasitoides_affiches = [4]    # indices des modèles dont afficher le parasitoïde ([] pour aucun)
variable_y            = "nym"  # axe Y hôte : "euf", "lar", "nym", "adu"
afficher_competition  = false   # graphique séparé : taux de mortalité larvaire (deltaL + alphaL*u[2]) * u[2]
afficher_comp_para    = false    # true = ajouter la courbe de mortalité des larves parasitées (alphaL*(L+Lp)*Lp), modèle larvaire uniquement
log_scale_y           = false   # true = ordonnées en échelle logarithmique (log10)
stade_para_morta      = ""      # stade parasité dont afficher les morts : "euf", "lar", "nym", ou "" pour désactiver
afficher_derivee      = false   # true = graphique de la dérivée + impression du premier t où elle devient négative
afficher_q0           = false   # true = graphique de q0 (fécondité effective Allee) en fonction du temps
afficher_q0xu4        = false  # true = graphique de q0×u[4] (production totale d'œufs) en fonction du temps
courbes_affich        = ""      # "" = mode normal | "allee" | "morta" | "stade"
morta_courbe          = "deltaL"  # paramètre de mortalité si courbes_affich = "morta" : "deltaE","deltaL","deltaN","deltaA"
valeurs_allee_dyn     = [5600,5700]   # valeurs d'allee si courbes_affich = "allee"
valeurs_morta_dyn     = [0]     # valeurs de mortalité si courbes_affich = "morta"
stade_affich          = [2]  # stades à afficher si courbes_affich = "stade" : 1=oeufs 2=larves 3=nymphes 4=adultes
ylim = Inf # scalaire (borne haute, basse = 0) ou [bas, haut]
xlimgr = 2000

display_value         = "max"   # "max" ou "mean"
t_display_min         = 3000.0  # début de la fenêtre de mesure pour display_value
t_display_max         = 5000.0  # fin de la fenêtre de mesure pour display_value

t_total = Inf  # durée totale de simulation (Inf = utiliser xlimgr[2] automatiquement)

ylims_dyn(y) = y isa AbstractVector ? (y[1], y[2]) : (0.0, y)
xlims_dyn(x) = x isa AbstractVector ? (x[1], x[2]) : (20, x)

module_analyse.ALLEE_TYPE[] = "Test"
tspan = (0.0, isinf(t_total) ? xlims_dyn(xlimgr)[2] : t_total)
alg = "fixe"
pdt = 0.0001


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "a", "alphaL", "deltaE", "deltaL", "deltaN", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ", "H", "mi", "allee")

## Fonction pertant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
   (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
   (dico_para["tauL"] + dico_para["tauM"]),dico_para["tauM"],dico_para["tauJ"],
   dico_para["tauJ"],dico_para["tauU"],dico_para["tauL"]]
end

# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (DicoDesParasGlobalNul,fct_lags, ordre_para)


## Fonction permettant de choisir les valeurs initiales 

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
        return [0, 0, 0, uStart, 1, uPara, - dico_para["tauL"] * ((0 + 0) * dico_para["alphaL"] + dico_para["deltaL"] + 
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


# Structurer les paramètres pour les inserer dans le modèle


# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
    # int des derniers pas de temps observé
para_run_0 = (fct_hini_0,tspan,"fixe",xlims_dyn(xlimgr)[2])
para_run_1 = (fct_hini_1,tspan,"fixe",xlims_dyn(xlimgr)[2])
para_run_2 = (fct_hini_2,tspan,"fixe",xlims_dyn(xlimgr)[2])
para_run_3 = (fct_hini_3,tspan,"fixe",xlims_dyn(xlimgr)[2])

# para save 
    # numéro du modèle 
    # Bool save ou non 
    # Save Y or N   
para_save_0 = ("Without parasitoid",false,false,0)   
para_save_1 = ("Eggs parasitoid",false,false,1) 
para_save_2 = ("Larval parasitoid",false,false,2) 
para_save_3 = ("Nymphal parasitoid",false,false,3) 


models_all     = [model0, model1, model2, model3]
para_runs_all  = [para_run_0, para_run_1, para_run_2, para_run_3]
para_saves_all = [para_save_0, para_save_1, para_save_2, para_save_3]

rez = [[],[],[],[]]
if isempty(courbes_affich)
    # Mode normal : une simulation par modèle
    @time rez[1] = module_analyse.runDynamics(model0, para, para_run_0, para_save_0);
    @time rez[2] = module_analyse.runDynamics(model1, para, para_run_1, para_save_1);
    @time rez[3] = module_analyse.runDynamics(model2, para, para_run_2, para_save_2);
    @time rez[4] = module_analyse.runDynamics(model3, para, para_run_3, para_save_3);
elseif courbes_affich == "stade" || courbes_affich == "stade+"
    # Mode stade : une simulation, stades tracés séparément ou additionnés
    idx_m = modeles_affiches[1]
    rez[idx_m] = module_analyse.runDynamics(models_all[idx_m], para, para_runs_all[idx_m], para_saves_all[idx_m])
else
    # Mode multi-courbes : un modèle, N valeurs d'un paramètre
    idx_m = modeles_affiches[1]
    valeurs   = courbes_affich == "allee" ? valeurs_allee_dyn : valeurs_morta_dyn
    rez_multi = []
    for val in valeurs
        dico = copy(DicoDesParasGlobalNul)
        courbes_affich == "allee" ? (dico["allee"] = val) : (dico[morta_courbe] = val)
        para_v = (dico, fct_lags, ordre_para)
        push!(rez_multi, module_analyse.runDynamics(
            models_all[idx_m], para_v, para_runs_all[idx_m], para_saves_all[idx_m]))
    end
    # Stocker dans rez pour les graphes secondaires (dernière valeur)
    rez[idx_m] = rez_multi[end]
end



# para_plot
    # string du titre du plot
    # strig du titre de l'axe des X
    # strig du titre de l'axe des Y
    # int de la dynamique à suivre 
    # vecteur de string du label de chaque courbes


var_map_dyn    = Dict("euf" => 1, "lar" => 2, "nym" => 3, "adu" => 4)
var_labels_dyn = Dict("euf" => "Densité des œufs", "lar" => "Densité larvaire",
                      "nym" => "Densité nymphale", "adu" => "Densité des adultes")

pop_labels_dyn  = Dict("euf" => "Population œufs", "lar" => "Population larvaire", "nym" => "Population nymphale", "adu" => "Population adultes")
model_names_all = ["Sans parasitoïde", "Parasitoïde des œufs", "Parasitoïde larvaire", "Parasitoïde nymphal"]
if isempty(courbes_affich)
    has_para_curves = !isempty(filter(x -> x != 1, parasitoides_affiches))
    multi_m = length(modeles_affiches) > 1 || has_para_curves
    labels_affiches = length(modeles_affiches) > 1 ? [model_names_all[i] for i in modeles_affiches] :
                                                     [var_labels_dyn[variable_y] for i in modeles_affiches]
    para_plot = [" ", "Temps", var_labels_dyn[variable_y], var_map_dyn[variable_y], labels_affiches]
    plt = @time module_analyse.traceDynamics(rez[modeles_affiches], para_plot)
elseif courbes_affich == "stade" || courbes_affich == "stade+"
    stage_names = Dict(1 => "Œufs", 2 => "Larves", 3 => "Nymphes", 4 => "Adultes")
    palette_s   = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions_s = collect(range(0, 1; length=12))
    colors_s    = [get(palette_s, positions_s[i]) for i in [2,4,6,10]]
    idx_m = modeles_affiches[1]
    plt = Plots.plot(title=" ", xlabel="Temps", ylabel="Densité de population",
                     grid=false, legend=:topright, legendfont=font(8))
    if courbes_affich == "stade"
        for (i, stade) in enumerate(stade_affich)
            Plots.plot!(plt, rez[idx_m][2], [u[stade] for u in rez[idx_m][1]],
                        label=stage_names[stade], color=colors_s[i], linewidth=2)
        end
    else  # "stade+"
        summed = [sum(u[s] for s in stade_affich) for u in rez[idx_m][1]]
        lbl    = join([stage_names[s] for s in stade_affich], " + ")
        Plots.plot!(plt, rez[idx_m][2], summed, label=lbl, color=colors_s[1], linewidth=2)
    end
else
    valeurs = courbes_affich == "allee" ? valeurs_allee_dyn : valeurs_morta_dyn
    labels_c = courbes_affich == "allee" ? ["intensité Allee = $v" for v in valeurs] :
                                           ["$(morta_courbe) = $v" for v in valeurs]
    para_plot = [" ", "Temps", var_labels_dyn[variable_y], var_map_dyn[variable_y], labels_c]
    plt = @time module_analyse.traceDynamics(rez_multi, para_plot)
end

# Ajout des courbes de parasitoïdes (u[6])
if !isempty(parasitoides_affiches) && isempty(courbes_affich)
    palette   = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions = collect(range(0, 1; length=12))
    indices_p = [1, 3, 5, 8]
    colors_p  = [get(palette, positions[i]) for i in indices_p]
    labels_para = ["Parasitoïde des œufs", "Parasitoïde larvaire", "Parasitoïde nymphal"]
    for (i, idx_m) in enumerate(filter(x -> x != 1, parasitoides_affiches))
        Plots.plot!(plt, rez[idx_m][2], [u[6] for u in rez[idx_m][1]],
                    label=labels_para[idx_m - 1], color=colors_p[i],
                    linewidth=2, linestyle=:dash)
    end
end

Plots.xlims!(plt, xlims_dyn(xlimgr)...)
if log_scale_y
    Plots.plot!(plt, yscale=:log10)
    Plots.ylabel!(plt, var_labels_dyn[variable_y] * " (log)")
    Plots.ylims!(plt, 0.01, Inf)
    log_ticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5]
    log_labels = ["0.01", "0.1", "1", "10", "100", "1 000", "10 000", "100 000"]
    Plots.yticks!(plt, log_ticks, log_labels)
else
    Plots.ylims!(plt, ylims_dyn(ylim)...)
end
let t = module_analyse.ALLEE_TYPE[]
    models_txt = join([model_names_all[i] for i in modeles_affiches], " & ")
    allee_txt = if t == "none"
        "Temps\n$models_txt"
    elseif t == "Wittmann"
        "Temps\n$models_txt — intensité Allee = $(DicoDesParasGlobalNul["allee"])"
    elseif t == "test" || t == "Test"
        "Temps\n$models_txt — intensité Allee = $(DicoDesParasGlobalNul["allee"])"
    elseif t == "Wu"
        "Temps\n$models_txt — intensité Allee = $(DicoDesParasGlobalNul["allee"])"
    elseif t == "yue"
        "Temps\n$models_txt — intensité Allee = $(DicoDesParasGlobalNul["allee"])"
    else
        "Temps\n$models_txt — intensité Allee = $(DicoDesParasGlobalNul["allee"])"
    end
    Plots.xlabel!(plt, allee_txt)
end

display(plt)

var_idx_dyn = var_map_dyn[variable_y]
summarize(vals) = display_value == "max" ? maximum(vals) : mean(vals)
label_dv = display_value == "max" ? "max" : "mean"

function window_vals(rez_m, idx)
    t_vals = rez_m[2]
    mask = findall(t -> t_display_min <= t <= t_display_max, t_vals)
    isempty(mask) ? [rez_m[1][end][idx]] : [rez_m[1][i][idx] for i in mask]
end

if isempty(courbes_affich)
    println("=== $(uppercase(display_value)) values (t=$(t_display_min)–$(t_display_max)) ===")
    for idx_m in modeles_affiches
        v = window_vals(rez[idx_m], var_idx_dyn)
        println("  $(model_names_all[idx_m]) — $(pop_labels_dyn[variable_y])  min: $(minimum(v))  ($label_dv): $(summarize(v))")
    end
    if !isempty(parasitoides_affiches)
        for idx_m in filter(x -> x != 1, parasitoides_affiches)
            v_p = window_vals(rez[idx_m], 6)
            println("  $(model_names_all[idx_m]) — Parasitoid  min: $(minimum(v_p))  ($label_dv): $(summarize(v_p))")
        end
    end
end

if stade_para_morta != "" && isempty(courbes_affich)
    # idx_m du modèle correspondant au stade parasité demandé
    # 2 = eggs parasitoid (u[5] = oeufs parasités)
    # 3 = larval parasitoid (u[5] = larves parasitées)
    # 4 = nymphal parasitoid (u[5] = nymphes parasitées)
    stage_model_map = Dict("euf" => 2, "lar" => 3, "nym" => 4)
    stage_label_map = Dict("euf" => "œufs parasités", "lar" => "larves parasitées", "nym" => "nymphes parasitées")
    target_idx = stage_model_map[stade_para_morta]

    if target_idx in modeles_affiches
        t_vals    = Float64.(rez[target_idx][2])
        mort_rate = if stade_para_morta == "euf"
            [DicoDesParasGlobalNul["deltaE"] * u[5] for u in rez[target_idx][1]]
        elseif stade_para_morta == "lar"
            [(DicoDesParasGlobalNul["alphaL"] * (u[2] + u[5]) + DicoDesParasGlobalNul["deltaL"]) * u[5]
             for u in rez[target_idx][1]]
        else  # nym
            [DicoDesParasGlobalNul["deltaN"] * u[5] for u in rez[target_idx][1]]
        end
        val_morta = summarize(mort_rate)
        println("=== Parasitized mortality ($(stage_label_map[stade_para_morta])) ===")
        println("  $(model_names_all[target_idx])  min: $(minimum(mort_rate))  ($label_dv): $val_morta")
    else
        println("⚠ stade_para_morta = \"$stade_para_morta\" mais le modèle correspondant (idx $target_idx) n'est pas dans modeles_affiches.")
    end
end

# Graphique de la mortalité larvaire par concurrence
if afficher_competition
    alphaL = DicoDesParasGlobalNul["alphaL"]
    deltaL = DicoDesParasGlobalNul["deltaL"]

    palette_c  = ColorSchemes.colorschemes[:twelvebitrainbow]
    n_curves_c = !isempty(courbes_affich) && courbes_affich ∈ ("allee","morta") ? length(rez_multi) : length(modeles_affiches)
    colors_c   = [get(palette_c, p) for p in range(0.1, 0.9; length=max(n_curves_c, 2))]

    labels_all = ["Sans parasitoïde", "Parasitoïde des œufs", "Parasitoïde larvaire", "Parasitoïde nymphal"]

    plt_comp = Plots.plot(title  = "Larval mortality rate — competition",
                          xlabel = "Time",
                          ylabel = "Concurence larvaire : alphaL·(L+Lp) · L",
                          grid   = false)

    if !isempty(courbes_affich) && courbes_affich ∈ ("allee", "morta")
        valeurs_c  = courbes_affich == "allee" ? valeurs_allee_dyn : valeurs_morta_dyn
        labels_c   = courbes_affich == "allee" ? ["intensité Allee = $v" for v in valeurs_c] :
                                                  ["$(morta_courbe) = $v" for v in valeurs_c]
        idx_m = modeles_affiches[1]
        for (i, rz) in enumerate(rez_multi)
            t_vals    = rz[2]
            mort_vals = [let L = u[2], Lp = (idx_m == 3 ? u[5] : 0.0)
                             DicoDesParasGlobalNul["alphaL"] * (L + Lp) * L
                         end for u in rz[1]]
            Plots.plot!(plt_comp, t_vals, mort_vals,
                        label=labels_c[i], color=colors_c[i], linewidth=2)
        end
        if afficher_comp_para && modeles_affiches[1] == 3
            for (i, rz) in enumerate(rez_multi)
                t_vals     = rz[2]
                mort_para  = [DicoDesParasGlobalNul["alphaL"] * (u[2] + u[5]) * u[5] for u in rz[1]]
                Plots.plot!(plt_comp, t_vals, mort_para,
                            label="Lp parasitées — $(labels_c[i])", color=colors_c[i],
                            linewidth=2, linestyle=:dot)
            end
        end
    else
        for (i, idx_m) in enumerate(modeles_affiches)
            t_vals    = rez[idx_m][2]
            mort_vals = [let L = u[2], Lp = (idx_m == 3 ? u[5] : 0.0)
                             DicoDesParasGlobalNul["alphaL"] * (L + Lp) * L
                         end for u in rez[idx_m][1]]
            Plots.plot!(plt_comp, t_vals, mort_vals,
                        label=labels_all[idx_m], color=colors_c[i], linewidth=2)
        end
        if afficher_comp_para && 3 ∈ modeles_affiches
            i_para = findfirst(==(3), modeles_affiches)
            t_vals    = rez[3][2]
            mort_para = [DicoDesParasGlobalNul["alphaL"] * (u[2] + u[5]) * u[5] for u in rez[3][1]]
            Plots.plot!(plt_comp, t_vals, mort_para,
                        label="Larves parasitées (compétition)", color=colors_c[i_para],
                        linewidth=2, linestyle=:dot)
        end
    end

    Plots.xlims!(plt_comp, xlims_dyn(xlimgr)...)
    if log_scale_y
        Plots.plot!(plt_comp, yscale=:log10)
        Plots.ylims!(plt_comp, 1e-3, ylims_dyn(ylim)[2])
        log_ticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5]
        log_labels = ["0.001", "0.01", "0.1", "1", "10", "100", "1 000", "10 000", "100 000"]
        Plots.yticks!(plt_comp, log_ticks, log_labels)
    else
        Plots.ylims!(plt_comp, ylims_dyn(ylim)...)
    end
    display(plt_comp)

    println("=== Mean competition values ===")
    for idx_m in modeles_affiches
        val_mean_comp = mean(let L = u[2], Lp = (idx_m == 3 ? u[5] : 0.0)
                                DicoDesParasGlobalNul["alphaL"] * (L + Lp) * L
                            end for u in rez[idx_m][1])
        println("  $(model_names_all[idx_m]) : $val_mean_comp")
    end
end

# Taux de croissance per capita (dN/dt)/N vs N (adultes)
if afficher_derivee && isempty(courbes_affich)
    palette_d   = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions_d = collect(range(0, 1; length=12))
    indices_d   = [4, 6, 10, 2]
    colors_d    = [get(palette_d, positions_d[i]) for i in indices_d]

    plt_deriv = Plots.plot(title  = "Taux de croissance per capita des adultes",
                           xlabel = "Adult pop (N)",
                           ylabel = "(dN/dt) / N",
                           grid   = false)
    Plots.hline!(plt_deriv, [0.0], color=:black, linewidth=1, linestyle=:dash, label="")

    println("=== Densité adulte au croisement zéro (seuil de décroissance) ===")
    for (i, idx_m) in enumerate(modeles_affiches)
        adu_vals = [u[4] for u in rez[idx_m][1]]
        t_vals   = Float64.(rez[idx_m][2])
        dN       = diff(adu_vals) ./ diff(t_vals)
        N_mid    = (adu_vals[1:end-1] .+ adu_vals[2:end]) ./ 2
        r_vals   = dN ./ max.(N_mid, 1e-10)   # évite division par zéro

        Plots.plot!(plt_deriv, N_mid, r_vals,
                    label=model_names_all[idx_m], color=colors_d[i], linewidth=2)

        # Premier croisement de zéro (passage positif → négatif)
        idx_cross = findfirst(i -> r_vals[i] >= 0 && r_vals[i+1] < 0, 1:length(r_vals)-1)
        if idx_cross !== nothing
            N_cross = N_mid[idx_cross]
            println("  $(model_names_all[idx_m]) — N* ≈ $N_cross")
        else
            println("  $(model_names_all[idx_m]) — pas de croisement détecté")
        end
    end

    Plots.xlims!(plt_deriv, 0.0, 10.0)
    display(plt_deriv)
end

# Fécondité effective q0 = q * DA / (DA + seuil_allee) en fonction du temps
if afficher_q0
    q_val = DicoDesParasGlobalNul["q"]

    palette_q  = ColorSchemes.colorschemes[:twelvebitrainbow]
    n_curves_q = !isempty(courbes_affich) && courbes_affich ∈ ("allee","morta") ? length(rez_multi) : length(modeles_affiches)
    colors_q   = [get(palette_q, p) for p in range(0.1, 0.9; length=max(n_curves_q, 2))]

    plt_q0 = Plots.plot(title  = "Fécondité effective q0 en fonction du temps",
                        xlabel = "Time",
                        ylabel = "q0 = q · DA(t) / (DA(t) + seuil)",
                        grid   = false)

    if !isempty(courbes_affich) && courbes_affich ∈ ("allee", "morta")
        valeurs_c = courbes_affich == "allee" ? valeurs_allee_dyn : valeurs_morta_dyn
        labels_c  = courbes_affich == "allee" ? ["intensité Allee = $v" for v in valeurs_c] :
                                                 ["$(morta_courbe) = $v" for v in valeurs_c]
        for (i, rz) in enumerate(rez_multi)
            seuil  = courbes_affich == "allee" ? valeurs_c[i] : DicoDesParasGlobalNul["allee"]
            t_vals  = rz[2]
            q0_vals = [let a = u[4]
                           iszero(seuil) ? q_val : q_val * a / (a + seuil)
                       end for u in rz[1]]
            Plots.plot!(plt_q0, t_vals, q0_vals,
                        label=labels_c[i], color=colors_q[i], linewidth=2)
        end
    else
        seuil = DicoDesParasGlobalNul["allee"]
        for (i, idx_m) in enumerate(modeles_affiches)
            t_vals  = rez[idx_m][2]
            q0_vals = [let a = u[4]
                           iszero(seuil) ? q_val : q_val * a / (a + seuil)
                       end for u in rez[idx_m][1]]
            Plots.plot!(plt_q0, t_vals, q0_vals,
                        label=model_names_all[idx_m], color=colors_q[i], linewidth=2)
        end
    end

    Plots.xlims!(plt_q0, xlims_dyn(xlimgr)...)
    Plots.ylims!(plt_q0, ylims_dyn(ylim)...)
    display(plt_q0)
end

# Production totale d'œufs q0×u[4] en fonction du temps
if afficher_q0xu4
    q_val = DicoDesParasGlobalNul["q"]

    palette_q4  = ColorSchemes.colorschemes[:twelvebitrainbow]
    n_curves_q4 = !isempty(courbes_affich) && courbes_affich ∈ ("allee","morta") ? length(rez_multi) : length(modeles_affiches)
    colors_q4   = [get(palette_q4, p) for p in range(0.1, 0.9; length=max(n_curves_q4, 2))]

    plt_q0xu4 = Plots.plot(title  = "Production totale d'œufs q0×u[4]",
                            xlabel = "Time",
                            ylabel = "q0 × u[4] (œufs / temps)",
                            grid   = false)

    if !isempty(courbes_affich) && courbes_affich ∈ ("allee", "morta")
        valeurs_c = courbes_affich == "allee" ? valeurs_allee_dyn : valeurs_morta_dyn
        labels_c  = courbes_affich == "allee" ? ["intensité Allee = $v" for v in valeurs_c] :
                                                 ["$(morta_courbe) = $v" for v in valeurs_c]
        for (i, rz) in enumerate(rez_multi)
            seuil  = courbes_affich == "allee" ? valeurs_c[i] : DicoDesParasGlobalNul["allee"]
            t_vals = rz[2]
            vals   = [let a = u[4]
                          q0 = iszero(seuil) ? q_val : q_val * a / (a + seuil)
                          q0 * a
                      end for u in rz[1]]
            Plots.plot!(plt_q0xu4, t_vals, vals,
                        label=labels_c[i], color=colors_q4[i], linewidth=2)
        end
    else
        seuil = DicoDesParasGlobalNul["allee"]
        for (i, idx_m) in enumerate(modeles_affiches)
            t_vals = rez[idx_m][2]
            vals   = [let a = u[4]
                          q0 = iszero(seuil) ? q_val : q_val * a / (a + seuil)
                          q0 * a
                      end for u in rez[idx_m][1]]
            Plots.plot!(plt_q0xu4, t_vals, vals,
                        label=model_names_all[idx_m], color=colors_q4[i], linewidth=2)
        end
    end

    Plots.xlims!(plt_q0xu4, xlims_dyn(xlimgr)...)
    Plots.ylims!(plt_q0xu4, ylims_dyn(ylim)...)
    display(plt_q0xu4)
end

