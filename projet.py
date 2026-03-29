# Bibliothèques utilisées dans tout le projet
import math
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.optimize import linear_sum_assignment


# ==============================================================================
# PARTIE 1 — Génération des profils de vote
# ==============================================================================

# Génère un profil dans A^n avec niveau de polarisation contrôlé.
# polarisation ∈ [0,1] : 0 = consensus, 1 = bipolaire maximal
# p : probabilité de bruit par coordonnée (loi binomiale)
# Retourne une liste de n bulletins (listes de 0/1 de longueur m)
def generation_profile_approbation(n, m, polarisation, p=0.1):

    profile = []
    vote = [np.random.randint(0, 2) for _ in range(m)]
    vote_opp = [1 - x for x in vote]

    if polarisation == 0:
        return [vote.copy() for _ in range(n)]

    if polarisation == 1:
        return [vote.copy() for _ in range(n // 2)] + [vote_opp.copy() for _ in range(n // 2)]

    # fonction qui prends en argument une liste de vote et modifie aléatoirement un nombre de vote selon une loi binomial avec comme paramètre p
    def ajouter_bruit(base):
        copie = base.copy()
        nb_modifs = np.random.binomial(m, p)
        indices = np.random.choice(m, nb_modifs, replace=False)
        for i in indices:
            copie[i] = 1 - copie[i]
        return copie

    # pour les n électeur on décide de qu'elle "côté de la polarisation" il va être puis on appliquqe ajouter_bruit
    for _ in range(n):
        r = np.random.rand()

        if r < polarisation / 2:
            profile.append(ajouter_bruit(vote_opp))

        elif r < polarisation:
            profile.append(ajouter_bruit(vote))

        else:
            profile.append(vote.copy())

    return profile


# Génère un profil dans L^n avec niveau de polarisation contrôlé.
# scale, spread : paramètres du bruit local sur le classement
# Retourne une liste de n classements (permutations de {0,...,m-1})
def generation_profile_ordretotal(n, m, polarisation, p=0.1, scale=1, spread=1):

    # classement de base
    vote = np.random.permutation(m).tolist()
    vote_opp = vote[::-1]

    # cas extrêmes
    if polarisation == 0:
        return [vote.copy() for _ in range(n)]

    if polarisation == 1:
        return [vote.copy() for _ in range(n // 2)] + [vote_opp.copy() for _ in range(n // 2)]

    def ajouter_bruit_ordre(base):
        copie = base.copy()
        nb_modifs = np.random.binomial(m, p)
        for _ in range(nb_modifs):
            # tirage autour du milieu
            i = int(np.random.normal(loc=m // 2, scale=scale))

            # premier indice
            i = max(0, min(m - 2, i))  # m-2 pour éviter le dépassement avec j

            # déplacement local
            diff = np.random.randint(1, spread + 1)
            j = min(m - 1, i + diff)

            copie[i], copie[j] = copie[j], copie[i]

        return copie

    profile = []

    for _ in range(n):
        r = np.random.rand()

        if r < polarisation / 2:
            profile.append(ajouter_bruit_ordre(vote_opp))

        elif r < polarisation:
            profile.append(ajouter_bruit_ordre(vote))

        else:
            profile.append(vote.copy())

    return profile


# ==============================================================================
# PARTIE 2 — Mesure φ²
# ==============================================================================

# Calcule d_{c_k, c_l}(p) pour toutes les paires de candidates.
# Retourne un dict {(k, j): valeur} pour chaque paire k < j
def calcul_d_approbation(p):
    d = {}

    m = len(p[0])
    for k in range(m):
        for j in range(k+1, m):
            n_kj = 0
            n_jk = 0
            for a in p:
                if a[k] == 1 and a[j] == 0: n_kj += 1
                elif a[k] == 0 and a[j] == 1: n_jk += 1
            d[(k, j)] = abs(n_kj - n_jk)
    return d


# Même calcul pour un profil par ordre total.
# p[i] : classement (liste), p[i][r] = candidate au rang r
def calcul_d_ordre(p):
    d = {}
    n = len(p)
    m = len(p[0])
    rang = []
    for a in p:
        rank = {candidate: i for i, candidate in enumerate(a)}
        rang.append(rank)
    for k in range(m):
        for j in range(k+1, m):
            n_kj = sum([1 for rank in rang if rank[k] < rank[j]])
            n_jk = n - n_kj
            d[(k, j)] = abs(n_kj - n_jk)
    return d


# Calcule φ²(p) pour un profil par approbation.
# Retourne la moyenne normalisée de (n - d_{c_k c_l}) sur toutes les paires.
def phi2_approbation(p):
    n = len(p)
    m = len(p[0])
    d = calcul_d_approbation(p)
    tot = 0
    for i in d.values():
        tot += (n - i) / n
    return tot / (math.factorial(m) / (math.factorial(m - 2) * 2))


# Même mesure pour un profil par ordre total.
def phi2_ordre(p):
    n = len(p)
    m = len(p[0])
    d = calcul_d_ordre(p)
    tot = 0
    for i in d.values():
        tot += (n - i) / n
    return tot / (math.factorial(m) / (math.factorial(m - 2) * 2))


# ==============================================================================
# PARTIE 3 — Distances et mesures φ_dH / φ_dS
# ==============================================================================

# Distance de Hamming : nombre de positions différentes entre deux bulletins a1, a2 ∈ A.
def distance_hamming(a1, a2):
    return np.sum(np.array(a1) != np.array(a2))


# Distance de Spearman : somme des écarts de rang entre deux ordres totaux.
# s1, s2 : classements sous forme de liste (indice = rang, valeur = candidate)
def distance_spearman(s1, s2):
    rank = {c: i for i, c in enumerate(s1)}
    rank2 = {c: i for i, c in enumerate(s2)}
    return sum(abs(rank[c] - rank2[c]) for c in rank)


# Calcule u1*(p) pour des votes par approbation : min sum d_H(a, a') sur a ∈ A.
# Exploite la décomposition coordonnée par coordonnée (valeur majoritaire par position).
# Retourne la somme des erreurs minimales par coordonnée.
def u1_approbation(p):
    n = len(p)
    m = len(p[0])
    u1 = 0
    for j in range(m):
        n_j1 = sum(a[j] for a in p)
        n_j0 = n - n_j1
        u1 += min(n_j0, n_j1)
    return u1


# Calcule u1*(p) pour des votes par ordre total via couplage parfait de poids min.
# Retourne la valeur optimale (entier).
def u1_ordre(p):
    m = len(p[0])

    # matrice des coûts
    W = np.zeros((m, m), dtype=int)

    ranks = []
    for a in p:
        rank = {candidate: i for i, candidate in enumerate(a)}
        ranks.append(rank)

    # matrice de coûts W[i,k] = coût de placer la candidate i au rang k
    for i in range(m):
        for k in range(m):
            W[i, k] = sum(abs(k - rank[i]) for rank in ranks)

    # problème de couplage parfait
    row_ind, col_ind = linear_sum_assignment(W)

    return int(W[row_ind, col_ind].sum())


# Calcule le bulletin consensus d'un cluster (votes par approbation) : valeur majoritaire par position.
# Retourne une liste de 0/1 de longueur m.
def centroide_approbation(cluster):
    n = len(cluster)
    m = len(cluster[0])
    centroide = []
    for j in range(m):
        n_j1 = sum(a[j] for a in cluster)
        if n_j1 >= n - n_j1:
            centroide.append(1)
        else:
            centroide.append(0)
    return centroide


# Calcule le classement consensus d'un cluster (votes par ordre total) via couplage parfait.
# Retourne un classement sous forme de liste (indice = rang, valeur = candidate).
def centroide_ordre(cluster):
    m = len(cluster[0])

    # on construit les rangs de chaque votante du cluster
    ranks = []
    for a in cluster:
        rank = {candidate: i for i, candidate in enumerate(a)}
        ranks.append(rank)

    # matrice de coûts W[i,k] = coût de placer la candidate i au rang k
    W = np.zeros((m, m), dtype=int)
    for i in range(m):
        for k in range(m):
            W[i, k] = sum(abs(k - rank[i]) for rank in ranks)

    row_ind, col_ind = linear_sum_assignment(W)

    # col_ind[i] = le rang attribué à la candidate i ; donc ordre[rang] = candidate, ce qui correspond à notre format de profil
    ordre = [0] * m
    for i in range(m):
        ordre[col_ind[i]] = i
    return ordre


# Estime ũ2*(p) par k-means (k=2) pour des votes par approbation.
# nb_restarts : relances pour limiter l'impact des optima locaux
# Retourne le meilleur coût obtenu (min sur toutes les relances).
def u2_approbation(p, nb_restarts=20):

    if len(set(tuple(a) for a in p)) < 2:
        return 0

    best = float('inf')  # on initialise best à l'infini pour pouvoir faire min dessus

    for _ in range(nb_restarts):
        a1 = random.choice(p).copy()
        a2 = random.choice([a for a in p if a != a1]).copy()

        while True:
            cluster_a1 = []
            cluster_a2 = []
            for a in p:
                if distance_hamming(a, a1) <= distance_hamming(a, a2):
                    cluster_a1.append(a)
                else:
                    cluster_a2.append(a)

            if not cluster_a1 or not cluster_a2:
                break

            a1_new = centroide_approbation(cluster_a1)
            a2_new = centroide_approbation(cluster_a2)

            if a1_new == a1 and a2_new == a2:
                break

            a1, a2 = a1_new, a2_new

        if cluster_a1 and cluster_a2:
            val = sum(distance_hamming(a, a1) for a in cluster_a1) + sum(distance_hamming(a, a2) for a in cluster_a2)
            best = min(best, val)  # on garde le min sur toutes les relances

    return best


# Même estimation pour des votes par ordre total.
def u2_ordre(p, nb_restarts=20):
    # même structure que u2_approbation
    if len(set(tuple(a) for a in p)) < 2:
        return 0

    best = float('inf')

    for _ in range(nb_restarts):
        a1 = random.choice(p).copy()
        a2 = random.choice([a for a in p if a != a1]).copy()

        cluster_a1 = []
        cluster_a2 = []

        while True:
            cluster_a1 = []
            cluster_a2 = []
            for a in p:
                if distance_spearman(a, a1) <= distance_spearman(a, a2):
                    cluster_a1.append(a)
                else:
                    cluster_a2.append(a)

            if not cluster_a1 or not cluster_a2:
                break

            a1_new = centroide_ordre(cluster_a1)
            a2_new = centroide_ordre(cluster_a2)

            if a1_new == a1 and a2_new == a2:
                break

            a1, a2 = a1_new, a2_new

        if cluster_a1 and cluster_a2:
            val = sum(distance_spearman(a, a1) for a in cluster_a1) + sum(distance_spearman(a, a2) for a in cluster_a2)
            best = min(best, val)

    return best


# Calcule φ_dH(p) = (2 / n*m) * (u1*(p) - ũ2*(p)) pour un profil par approbation.
def phi_dH(p):
    n = len(p)
    m = len(p[0])
    return (2 / (n * m)) * (u1_approbation(p) - u2_approbation(p))


# Calcule φ_dS(p) = (4 / n*m²) * (u1*(p) - ũ2*(p)) pour un profil par ordre total.
def phi_dS(p):
    n = len(p)
    m = len(p[0])
    return (4 / (n * m**2)) * (u1_ordre(p) - u2_ordre(p))


# ==============================================================================
# GRAPHIQUES (questions 6 et 15)
# ==============================================================================

# Trace l'évolution de φ² en fonction du paramètre de polarisation
# pour les deux types de votes (approbation et ordre total)
def plot_phi2(n=1000, m=10, nb_points=20):
    polarisations = np.linspace(0, 1, nb_points)
    phi2_approbation_val = []

    for pol in polarisations:
        profile = generation_profile_approbation(n=n, m=m, polarisation=pol, p=0.1)
        phi2_approbation_val.append(phi2_approbation(profile))

    plt.plot(polarisations, phi2_approbation_val, marker='o')
    plt.title('Évolution de φ2 en fonction de la polarisation (approbation)')
    plt.xlabel('Polarisation')
    plt.ylabel('φ2')
    plt.grid()
    plt.show()

    phi2_ordre_val = []
    for pol in polarisations:
        profile = generation_profile_ordretotal(n=n, m=m, polarisation=pol, p=0.1, scale=1, spread=1)
        phi2_ordre_val.append(phi2_ordre(profile))

    plt.plot(polarisations, phi2_ordre_val, marker='o')
    plt.title('Évolution de φ2 en fonction de la polarisation (ordre total)')
    plt.xlabel('Polarisation')
    plt.ylabel('φ2')
    plt.grid()
    plt.show()


# Trace l'évolution de φ_dH et φ_dS en fonction du paramètre de polarisation.
# n, m : taille des profils ; nb_tests : répétitions par valeur (résultat moyenné)
def experiment(n=50, m=10, nb_tests=10):
    polarisations = np.linspace(0, 1, 20)

    phi_H_vals = []
    phi_S_vals = []

    for pol in polarisations:
        phi_H_mean = 0
        phi_S_mean = 0

        for _ in range(nb_tests):
            # génération des profils
            p_H = generation_profile_approbation(n, m, pol)
            p_S = generation_profile_ordretotal(n, m, pol)

            # calculs
            phi_H_mean += phi_dH(p_H)
            phi_S_mean += phi_dS(p_S)

        phi_H_vals.append(phi_H_mean / nb_tests)
        phi_S_vals.append(phi_S_mean / nb_tests)

    # tracé
    plt.figure()
    plt.plot(polarisations, phi_H_vals, label="φ_dH (approbation)")
    plt.plot(polarisations, phi_S_vals, label="φ_dS (ordre)")
    plt.xlabel("Polarisation")
    plt.ylabel("Mesure φ")
    plt.title("Évolution de la polarisation")
    plt.legend()
    plt.grid()

    plt.show()


# ==============================================================================
# MAIN — Démonstration minimale du script
# ==============================================================================

def main():
    print("=== Génération de profils ===")
    p_app = generation_profile_approbation(n=10, m=4, polarisation=0.5)
    print("Profil approbation (n=10, m=4, pol=0.5) :")
    for bulletin in p_app:
        print(" ", bulletin)

    p_ord = generation_profile_ordretotal(n=10, m=4, polarisation=0.5)
    print("\nProfil ordre total (n=10, m=4, pol=0.5) :")
    for classement in p_ord:
        print(" ", classement)

    print("\n=== Mesure φ² ===")
    print("φ²(approbation) :", phi2_approbation(p_app))
    print("φ²(ordre total) :", phi2_ordre(p_ord))

    print("\n=== Distances ===")
    a1 = p_app[0]
    a2 = p_app[1]
    print(f"d_H({a1}, {a2}) =", distance_hamming(a1, a2))

    s1 = p_ord[0]
    s2 = p_ord[1]
    print(f"d_S({s1}, {s2}) =", distance_spearman(s1, s2))

    print("\n=== u1* ===")
    print("u1*(approbation) :", u1_approbation(p_app))
    print("u1*(ordre total) :", u1_ordre(p_ord))

    print("\n=== φ_dH et φ_dS ===")
    print("φ_dH(approbation) :", phi_dH(p_app))
    print("φ_dS(ordre total) :", phi_dS(p_ord))

    print("\n=== Graphiques ===")
    print("Tracé φ² (approbation et ordre total)...")
    plot_phi2(n=500, m=6)

    print("Tracé φ_dH et φ_dS (expérience complète)...")
    experiment(n=30, m=5, nb_tests=5)


if __name__ == "__main__":
    main()
