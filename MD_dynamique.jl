###############################################################################################################
#### MD_dynamique — Module central du modèle ravageur/parasitoïdes (Carpocapse)                         ####
####                                                                                                     ####
#### Ce module définit :                                                                                 ####
####   • Les 4 modèles DDE (sans parasitoïde, parasitoïde des œufs / des larves / des nymphes)          ####
####   • L'effet Allee sur la fécondité (types : "none", "Test", "Wittmann", "Wu", "yue")               ####
####   • runDynamics      : intégration numérique d'une dynamique temporelle (RK4 à pas fixe)           ####
####   • traceDynamics    : tracé des trajectoires de densité en fonction du temps                      ####
####   • simulate_bifurcation_Final : balayage paramétrique pour diagrammes de bifurcation              ####
####   • plot_bifurcation_Final     : affichage des diagrammes de bifurcation                           ####
###############################################################################################################

### Code générale ----
module module_analyse

# Fonction exporté 
export runDynamics, traceDynamics, without_parasitoid_model, eggs_parasitoid_model, larval_parasitoid_model, nymphal_parasitoid_model, simulate_bifurcation_Final, plot_bifurcation_Final

# Import packages
using DifferentialEquations, Plots, NumericalIntegration, CSV, Trapz, FFTW
using Trapz, Statistics, MAT, DelimitedFiles, Dates, MAT, DataFrames, ColorSchemes

# ─────────────────────────────────────────────────────────────────────────────
# Type d'effet Allee appliqué à la fécondité q
#   Valeurs possibles : "none", "Wittmann", "test", "wu", "yue"
#   Pour changer depuis un script : module_analyse.ALLEE_TYPE[] = "Wittmann"
# ─────────────────────────────────────────────────────────────────────────────
const ALLEE_TYPE = Ref{String}("Test")

# Paramètres de l'effet Allee (à renseigner dans le dictionnaire de paramètres) :
#   seuil_allee : seuil de l'effet Allee (Wittmann, Test, Wu, yue)

function _allee_q(q, u4, seuil_allee)
    type = ALLEE_TYPE[]
    if type == "none"
        return q
    elseif type == "Wittmann"
        return q * exp(- seuil_allee / (u4))
    elseif type == "Test"
        iszero(seuil_allee) && return q
        return ((u4 / (u4 + seuil_allee)) * q)
    elseif type == "Wu"
        iszero(seuil_allee) && return q
        return ((u4 / (u4 + seuil_allee)) * q)
    elseif type == "yue"
        return q * (1.0 - seuil_allee / (u4 + 1.0))
    else
        error("ALLEE_TYPE inconnu : \"$type\". Valeurs valides : \"none\", \"Wittmann\", \"Test\", \"Wu\", \"yue\"")
    end
end


### Création de la fonction de dymanique ----

# Paramètres 
#   modele : fonction du modèle à retard 
#           
#   para   : tuple d'un tuple et d'un vecteur : 
#               tuple des paramètres du modèle 
#               vecteur des lags utilisé dans le modèle 
#
#   para_run : tuple contenant : 
#               vecteur des valeurs à 0
#               vecteur des valeurs historiques
#               tuple du début et fin de la simul
#               la méthode de résolution 
#               taile du pas de temps (facultatif)
#               
#   para_save: tuple des élément de localisation 
#               numéro du modèle 
#               Bool présence de densité dépendence sur les larves 
#               Bool save ou non
        
#
function runDynamics(modele, para, para_run, para_save)
    
    # Extraire les paramètres de sauvgarde
    presence_dd = para_save[2]
    YN_save = para_save[3]
    nb_modele = para_save[4]

    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    p = fixed_params[ordre_para[1]]
    for (h) in 2:length(ordre_para)
        p = (p...,fixed_params[ordre_para[h]])
    end
 

    # Extraire les paramètres de simulation

    u0 = para_run[1](fixed_params,"u")
    h0 = para_run[1](fixed_params,"h")
    h(p,t) = h0
    tspan = para_run[2]
    alg = para_run[3]
    t_plot = para_run[4]
    if length(para_run) >= 5
        pdt = para_run[5]
    else 
        pdt = 0.0001
    end


    # Fixer l'horraire
    heure = string(now())

    # Lag total 
    lags = lags_fct(fixed_params)

    # Résoudre le problème selon le solveur choisis 
    if alg == "maison"
        # Définir le problème et le résoudre
        t_values, u_values, eq_type = solveur_maison.solve_dde_rk2_first(modele, lags, tspan, u0, pdt, p, h0)
    
    elseif alg == "fixe"
        # Définir le problème et le résoudre
        prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
        freq_save = 0.1

        sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save, abstol=1e-9, reltol=1e-9, maxiters=Int(1e9))
     
        # Extraire les solutions
        t_values = sol.t # Les temps où la solution est évaluée
        u_values = sol.u # Les valeurs des solutions correspondantes
    end

    t_values = collect(t_values)

    # Séléctionner les valeurs tracé
    if t_values[end] > t_plot
        t_last_values = [x for x in t_values if x > (t_values[end] - t_plot)]
        u_last_values = u_values[end-length(t_last_values)+1:end] 
    else 
        t_last_values = t_values
        u_last_values = u_values
    end 
    
    t_lite = t_last_values
    u_lite = u_last_values

    return [u_lite, t_lite]
end


### Création de la fonction de traceDynamics ----

# Paramètres 
#   rez   : Vector{Vector{Any}} contenant les résultats des dynamiques 
#               
#   para_plot : tuple des paramètres pour l'affichage des courbes
#               tuple du numéro des variables tracé
#               tuple du nom des varaibles
#               taille du segment de la dynamique tracé

function traceDynamics(rez, para_plot)

    plot_title = para_plot[1]
    x_name = para_plot[2]
    y_name = para_plot[3]
    curve_plot = para_plot[4]
    labels = para_plot[5]

    # Créer les couleurs dynamiquement selon le nombre de courbes
    palette = ColorSchemes.colorschemes[:twelvebitrainbow]
    n = length(rez)
    colors = [RGB{Float64}(get(palette, p)) for p in range(0.25, 0.9; length=max(n, 2))]

    # Créer le plot vide
    plt = plot(title=plot_title, xlabel = x_name, ylabel = y_name, lw=4, grid=false, titlefont = font(15), legend=:topright, legendfont = font(8))

    # Ajouter les courbes
    for i in 1:length(rez)
            plot!(plt, rez[i][2],  [u[curve_plot] for u in rez[i][1]], 
            label=labels[i], color=colors[i], linewidth = 2)
    end

    display(plt)  
    
    return(plt)
end 






### Création de la fonction diagramme de bifurcation Final  ----

# Paramètres 
#
#   modele : fonction : fonction permettant la simulation numérique des différents modèles
#
#   para   : tuple d'un tuple et d'un vecteur : 
#               dictionnaire des paramètres du modèle et de leur valeur
#               vecteur des string des lags utilisé dans le modèle 
#               ordre des para dans p
#
#   para_run : tuple contenant : 
#               fonction créant les valeurs des varialbes à 0 et en historique
#               tuple du début et fin de la simul
#      
#   para_bifucation : tuple contenant : 
#               noms des paramètres qui serviront d'axe à la bifurcation
#               min et max du paramètre              
#               nb d'itération de calcule 
#               int des X dernier pas de temps étudié
#          
#         
#   para_save : tuple des élément de localisation 
#               numéro du modèle 
#               Bool save ou non 
#               Save Y or N
#   max : valeur maximum prose par la simulation
#                  
#   type : string du type de bifurcation étudié
#
#   solveur_type : String type de solveur ("default", "maison")
function simulate_bifurcation_Final(modele, para, para_run, para_bifucation, para_save, max, type, popobs=2)

    if max isa Float64 
    else
        max=1
    end

    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    # Extraire les paramètres de simulation
    tspan = para_run[2]
    if length(para_run) >= 3
        pdt = para_run[3]
    else 
        pdt = 0.0001
    end

    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[4]
    YN_save = para_save[3]
    
    # Extraire les paramètres pour la création de la bifurcation
     axe_range = range(para_bifucation[2][1],para_bifucation[2][2],para_bifucation[3])


    t_stud = para_bifucation[4]


    # Préparation du vecteur contenant les résultats totaux
    bifurcation_rez_tot = Vector{Vector{Vector{Float64}}}()
    bifurcation_rez = Vector{Vector{Float64}}()

    #for (o,morta) in enumerate(para_bifucation[1])
    for (o,morta) in enumerate(para_bifucation[5])

        if type[1] == "1D"
            if morta == "Without Mortality"
                fixed_params["deltaE"]=0
                fixed_params["deltaL"]=0
                fixed_params["deltaN"]=0 
            elseif morta == "Egg Mortality"
                fixed_params["deltaE"]=0.1
                fixed_params["deltaL"]=0
                fixed_params["deltaN"]=0
            elseif morta == "Larvae Mortality"
                fixed_params["deltaE"]=0
                fixed_params["deltaL"]=0.1
                fixed_params["deltaN"]=0
            elseif morta == "Nymphae Mortality"
                fixed_params["deltaE"]=0
                fixed_params["deltaL"]=0
                fixed_params["deltaN"]=0.1
            end
        end 
        
        # Fixer le nom des axes 
        #axe_name = string(morta)
        if type[1] == "bifurcation_mortalite"
            axe_name = string(para_bifucation[1][o])
        else
            axe_name = string(para_bifucation[1][1])
        end

        # Préparation du vecteur contenant les résultats locaux
        bifurcation_rez = Vector{Vector{Float64}}()


        # Boucle créatrice de la bifurcation
        for (i, axe) in enumerate(axe_range)
            
            # Créer les paramètres dynamiques pour ce run
            p = Dict()
            
            for (key, val) in fixed_params
                p[key] = val
            end

            p[axe_name] = axe

            # Transformer les dico en tuple
            lags = lags_fct(p)

            p_tuple = p[ordre_para[1]]
            for (h) in 2:length(ordre_para)
                p_tuple = (p_tuple...,p[ordre_para[h]])
            end
            p_dico = p
            p = p_tuple 

            # Définir les valeur des varaibles en 0 et leur historique
            u0 = para_run[1](p_dico,"u")
            h0 = para_run[1](p_dico,"h")
            h(p,t) = h0


            # Définir le problème et le résoudre
            prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
            freq_save = 0.1

            sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save,abstol=1e-9, reltol=1e-9)
        
            # Extraire les solutions
            t_values = sol.t # Les temps où la solution est évaluée
            u_values = sol.u # Les valeurs des solutions correspondantes
            

            # Sélectionner les dernières unités de temps 
            t_values = collect(t_values)
                if t_values[end] > t_stud
                t_last_values = [x for x in t_values if x > (t_values[end] - t_stud)]
                u_last_values = u_values[end-length(t_last_values)+1:end] 
            else 
                t_last_values = t_values
                u_last_values = u_values
            end
            
            # Garder min/max/mean
            MaxMinMean = [0.0,0.0,0.0,0.0]
            tot_ops = [v[popobs] for v in u_last_values]
            MaxMinMean[1] = maximum(tot_ops)/max
            MaxMinMean[2] = minimum(tot_ops)/max
            MaxMinMean[3] = mean(tot_ops)/max

            # Ajouter la valeur moyenne du para 
            if nb_modele != 0 && mean([v[6] for v in u_last_values]) > 0.001
                MaxMinMean[4] = 1.0 
            end 

            # Inserer les points pour cette mortalité dans la sortie total
            push!(bifurcation_rez, MaxMinMean)
                 
        end  
        push!(bifurcation_rez_tot, bifurcation_rez)
    end 

    return bifurcation_rez_tot
end


### Création de la fonction diagramme de bifurcation Final  ----

# Paramètres 
#   bifurcation_rez_tot   : données utilisé pour le graph
#      
#   para_bifucation : tuple contenant : 
#               noms des paramètres qui serviront d'axe à la bifurcation
#               min et max du paramètre              
#               nb d'itération de calcule 
#               int des X dernier pas de temps étudié
#          
#         
#   para_save : tuple des élément de localisation 
#               numéro du modèle 
#               Bool save ou non 
#               Save Y or N

function plot_bifurcation_Final(bifurcation_rez_tot, para_bifucation, para_save, para_plot_panel, para_plot)

    # Extraire les paramètres de sauvgarde
    nom_modele = para_save[1]
    YN_save = para_save[3]
    nb_modele = para_save[4]

    # Extraire les paramètres pour la création de la bifurcation
    axe_range = range(para_bifucation[2][1],para_bifucation[2][2],para_bifucation[3])



    # Tracer le diagramme de bifurcation
    plt = plot(title="",
    xlabel=para_plot[2],
    ylabel=para_plot[3],
    grid=false)

 
    if para_bifucation[6] == "log"
        Plots.plot!(yticks = (log10.([1,10,100,1000,10000]),["1","10","100","1000","10000"]))
    elseif para_bifucation[6] == "ratio"
        Plots.plot!(yticks = ([0,20,40,60,80,100],["0%","20%","40%","60%","80%","100%"])) 
        end 

    if para_plot[1] isa Vector{Float64}
        Plots.plot!(ylims=(para_plot[1][1],para_plot[1][2]))
    end 

    labels = para_bifucation[5]
    palette = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions = collect(range(0, 1; length=12))
    indices = [2, 4, 6, 10]
    colors = [RGB{Float64}(get(palette, positions[i])) for i in indices]


    #for (o,morta) in enumerate(para_bifucation[1])
    for (o,morta) in enumerate(para_bifucation[5])

        # Afficher la zone min/max,
        min_vals = [bifurcation_rez_tot[o][i][1] for (i, _) in enumerate(axe_range)]
        max_vals = [bifurcation_rez_tot[o][i][2] for (i, _) in enumerate(axe_range)]

        # Ajuster les valeur en log pour éviter les -inf
        if para_bifucation[6] == "log"
            min_vals = log10.(min_vals .+ 1)
            max_vals = log10.(max_vals .+ 1)
        elseif para_bifucation[6] == "ratio"
            min_vals = exp.(min_vals)*100
            max_vals = exp.(max_vals)*100
        end 


        Plots.plot!(
            plt,
            axe_range,
            max_vals,
            fillrange = min_vals,
            label = para_bifucation[5][o],
            color = colors[o],
            alpha = 0.6,
            linewidth = 0,
            fillalpha = 0.15,
        )

        if para_bifucation[6] == "ratio"
            Plots.plot!(
                legend = false
            )
        end 


        # Tracer les lignes
        vect_line = [v[3] for v in bifurcation_rez_tot[o]] 
        
        # Ajuster les valeur en log pour éviter les -inf
        if para_bifucation[6] == "log"
            vect_line = log10.(vect_line .+ 1) 
        elseif para_bifucation[6] == "ratio"
            vect_line = exp.(vect_line)*100 
        end 

        presence_para = [v[4] for v in bifurcation_rez_tot[o]] 
        if para_bifucation[6] == "ratio"
            presence_para = ones(length(presence_para))
        end 
        x_line = axe_range

        # Trouver l'indice de changement
        change_idx = findfirst(presence_para .!= presence_para[1])

        if isnothing(change_idx)
            # Pas de changement
            style = presence_para[1] == 0.0 ? :dashdotdot : :solid
            plot!(x_line, vect_line, linestyle=style, linewidth=2, color=colors[o], label = false)
        else
            # Première partie, jusqu'au point juste avant le changement
            plot!(x_line[1:change_idx-1], vect_line[1:change_idx-1],
                linestyle=presence_para[1] == 0.0 ? :dashdotdot : :solid, linewidth=2, color=colors[o], label = false)
            
            # Deuxième partie, à partir du point de changement
            plot!(x_line[change_idx:end], vect_line[change_idx:end],
                linestyle=presence_para[change_idx] == 0.0 ? :dashdotdot : :solid, linewidth=2, color=colors[o], label = false)
        end
    end

    if para_plot_panel[1] == "panel"
        plot!(legend=false)
        plot!(plot_titlefontsize=1)
        plot!(xguidefontsize=10)
        plot!(yguidefontsize=10)

        if nb_modele == 0
            plot!(xlabel = "")
        elseif nb_modele == 1
            plot!(xlabel = "")
            plot!(ylabel = "")
        elseif nb_modele == 2

        elseif nb_modele == 3
            plot!(ylabel = "")
        end
    end 


    display(plt)
    
end

### Création du modèle de dynamique ravageur/parasitoïdes des oeufs ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoïdes, alphaL = intensité de la compétition larvaire,
#           delta = taux de mortalité (E = des œufs, L = des larves, N = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des œufs, L = des larves, M = des nymphes, J = des parasitoïdes juvéniles)
#           H et mi non implémentés
#
#   t   : float : temps actuel

function without_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaE, deltaL, deltaN, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi, seuil_allee = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[5]
    hist5 = h(p, t - tauL)[2]

    q0 = _allee_q(q, u[4],  seuil_allee)
    q1 = _allee_q(q, hist1, seuil_allee)
    q2 = _allee_q(q, hist2, seuil_allee)
    q3 = _allee_q(q, hist3, seuil_allee)
    

    du[1] = q0 * u[4] - deltaE * u[1] - q1 * hist1 * exp(- deltaE * tauU)
    du[2] = q1 * hist1 * exp(- deltaE * tauU)   -
        q2 * hist2 * exp(-deltaE * tauU) * exp(u[5])   -   (deltaL + alphaL * u[2]) * u[2]
    du[3] = q2 * hist2 * exp(-deltaE * tauU) * exp(u[5]) - deltaN * u[3] -
        q3 * hist3 * exp(- deltaE * tauU) * exp(-deltaN * tauM) * exp(hist4)
    du[4] = q3 * hist3 * exp(- deltaE * tauU) * exp(-deltaN * tauM) * exp(hist4)   -   deltaA * u[4]
    du[5] = - alphaL * (u[2] - hist5)
end




### Création du modèle de dynamique ravageur/parasitoïdes des oeufs ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoïdes, alphaL = intensité de la compétition larvaire,
#           delta = taux de mortalité (E = des œufs, L = des larves, N = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des œufs, L = des larves, M = des nymphes, J = des parasitoïdes juvéniles)
#           H et mi non implémentés
#
#   t   : float : temps actuel

function eggs_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaE, deltaL, deltaN, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi, seuil_allee = p

    hist1  = h(p, t - tauU)[4]
    hist2  = h(p, t - (tauU + tauL))[4]
    hist3  = h(p, t - tauL)[7]
    hist4  = h(p, t - (tauU + tauL + tauM))[4]
    hist5  = h(p, t - (tauL + tauM))[7]
    hist6  = h(p, t - tauM)[8]
    hist7  = h(p, t - tauJ)[1]
    hist8  = h(p, t - tauJ)[6]
    hist9  = h(p, t - tauU)[6]
    hist10 = h(p, t - tauL)[2]

    q0 = _allee_q(q, u[4],  seuil_allee)
    q1 = _allee_q(q, hist1, seuil_allee)
    q2 = _allee_q(q, hist2, seuil_allee)
    q4 = _allee_q(q, hist4, seuil_allee)

    du[1] = q0 * u[4] - a * u[1] * u[6] - deltaE * u[1] - q1 * hist1 * exp(u[7])
    du[2] = q1 * hist1 * exp(u[7]) - (alphaL * u[2] + deltaL) * u[2] - q2 * hist2 * exp(hist3) * exp(u[8])
    du[3] = q2 * hist2 * exp(hist3) * exp(u[8]) - deltaN * u[3] - q4 * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaN)
    du[4] = q4 * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaN)   -   deltaA * u[4]
    du[5] = a * u[1] * u[6] - deltaE * u[5] - a * hist7 * hist8 * exp(- deltaE * tauJ)
    du[6] = a * hist7 * hist8 * exp(- deltaE * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist9)
    du[8] = - alphaL * (u[2] - hist10)
end





### Création du modèle de dynamique ravageur/parasitoïdes des larves ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoïdes, alphaL = intensité de la compétition larvaire,
#           delta = taux de mortalité (E = des œufs, L = des larves, N = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des œufs, L = des larves, M = des nymphes, J = des parasitoïdes juvéniles)
#           H et mi non implémentés
#
#   t   : float : temps actuel

function larval_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaE, deltaL, deltaN, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi, seuil_allee = p

    hist1  = h(p, t - tauU)[4]
    hist2  = h(p, t - (tauU + tauL))[4]
    hist3  = h(p, t - (tauU + tauL + tauM))[4]
    hist4  = h(p, t - tauM)[7]
    hist5  = h(p, t - tauJ)[2]
    hist6  = h(p, t - tauJ)[6]
    hist7  = h(p, t - tauL)[6]
    hist8  = h(p, t - tauL)[2]
    hist9  = h(p, t - tauL)[5]
    hist10 = h(p, t - tauJ)[2]
    hist11 = h(p, t - tauJ)[5]

    q0 = _allee_q(q, u[4],  seuil_allee)
    q1 = _allee_q(q, hist1, seuil_allee)
    q2 = _allee_q(q, hist2, seuil_allee)
    q3 = _allee_q(q, hist3, seuil_allee)

    du[1] = q0 * u[4] - deltaE * u[1] - q1 * hist1 * exp(- deltaE * tauU)
    du[2] = q1 * hist1 * exp(- deltaE * tauU) - a * u[2] * u[6] -
       (alphaL * (u[2] + u[5]) + deltaL) * u[2] - q2 * hist2 * exp(- deltaE * tauU) * exp(u[7])
    du[3] = q2 * hist2 * exp(- deltaE * tauU) * exp(u[7]) - deltaN * u[3] -
        q3 * hist3 * exp(- deltaE * tauU) * exp(hist4) * exp(- tauM * deltaN)
    du[4] = q3 * hist3 * exp(- deltaE * tauU) * exp(hist4) * exp(- tauM * deltaN)   -   deltaA * u[4]
    du[5] = a * u[2] * u[6] - (alphaL * (u[2] + u[5]) + deltaL) * u[5] - a * hist5 * hist6 * exp(u[8])
    du[6] = a * hist5 * hist6 * exp(u[8]) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) - alphaL * ((u[2] + u[5]) - (hist8 + hist9))
    du[8] = - alphaL * ((u[2] + u[5])-(hist10 + hist11))
end



### Création du modèle de dynamique ravageur/parasitoïdes des nymphes ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoïdes, alphaL = intensité de la compétition larvaire,
#           delta = taux de mortalité (E = des œufs, L = des larves, N = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des œufs, L = des larves, M = des nymphes, J = des parasitoïdes juvéniles)
#           H et mi non implémentés
#
#   t   : float : temps actuel

function nymphal_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaE, deltaL, deltaN, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi, seuil_allee = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[8]
    hist5 = h(p, t - tauJ)[3]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauM)[6]
    hist8 = h(p, t - tauL)[2]

    q0 = _allee_q(q, u[4],  seuil_allee)
    q1 = _allee_q(q, hist1, seuil_allee)
    q2 = _allee_q(q, hist2, seuil_allee)
    q3 = _allee_q(q, hist3, seuil_allee)

    du[1] = q0 * u[4] - deltaE * u[1] - q1 * hist1 * exp(- deltaE * tauU)
    du[2] = q1 * hist1 * exp(- deltaE * tauU) - (alphaL * u[2] + deltaL) * u[2] -
        q2 * hist2 * exp(- deltaE * tauU) * exp(u[8])
    du[3] = q2 * hist2 * exp(- deltaE * tauU) * exp(u[8]) - a * u[3] * u[6] - deltaN * u[3] -
        q3 * hist3 * exp(- deltaE * tauU) * exp(hist4) * exp(u[7])
    du[4] = q3 * hist3 * exp(- deltaE * tauU) * exp(hist4) * exp(u[7])   -   deltaA * u[4]
    du[5] = a * u[3] * u[6] - deltaN * u[5] - a * hist5 * hist6 * exp(- deltaN * tauJ)
    du[6] =  a * hist5 * hist6 * exp(- deltaN * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7)
    du[8] = - alphaL * (u[2] - hist8)
end


end