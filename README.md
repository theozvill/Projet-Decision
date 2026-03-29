# Projet — Analyse de la Polarisation

Projet L3 Informatique (2025–2026), cours *Fondements Mathématiques pour l'Aide à la Décision*.

Ce script Python implémente les fonctions développées dans le notebook du projet sur la **mesure de la polarisation des préférences électorales**.

---

## Contenu du fichier `projet.py`

Le fichier est organisé en quatre parties :

### Partie 1 — Génération des profils de vote

- `generation_profile_approbation(n, m, polarisation, p=0.1)` — Génère un profil de `n` bulletins d'approbation sur `m` candidates. Le paramètre `polarisation` (entre 0 et 1) contrôle le niveau de polarisation : 0 = consensus pur, 1 = deux groupes totalement opposés. Un bruit binomial de paramètre `p` est ajouté pour simuler la variabilité réelle.
- `generation_profile_ordretotal(n, m, polarisation, p=0.1, scale=1, spread=1)` — Même logique pour des votes par ordre total (classements). Le bruit est introduit par des échanges locaux tirés autour de la position médiane du classement.

### Partie 2 — Mesure φ²

- `calcul_d_approbation(p)` — Calcule les valeurs `d_{c_k, c_l}(p)` pour toutes les paires de candidates, pour un profil par approbation. Retourne un dictionnaire `{(k, j): valeur}`.
- `calcul_d_ordre(p)` — Même calcul pour un profil par ordre total.
- `phi2_approbation(p)` — Calcule la mesure φ² pour un profil par approbation.
- `phi2_ordre(p)` — Calcule la mesure φ² pour un profil par ordre total.

### Partie 3 — Distances et mesures φ_dH / φ_dS

- `distance_hamming(a1, a2)` — Distance de Hamming entre deux bulletins d'approbation.
- `distance_spearman(s1, s2)` — Distance de Spearman entre deux classements (ordres totaux).
- `u1_approbation(p)` — Calcule `u1*(p)` pour des votes par approbation en exploitant la valeur majoritaire par coordonnée (complexité O(nm)).
- `u1_ordre(p)` — Calcule `u1*(p)` pour des votes par ordre total via un couplage parfait de poids minimum (`scipy.optimize.linear_sum_assignment`).
- `centroide_approbation(cluster)` — Calcule le bulletin consensus d'un cluster (valeur majoritaire par position).
- `centroide_ordre(cluster)` — Calcule le classement consensus d'un cluster via couplage parfait.
- `u2_approbation(p, nb_restarts=20)` — Estime `ũ2*(p)` par l'algorithme k-means (k=2) pour des votes par approbation. Plusieurs relances aléatoires sont effectuées pour limiter l'impact des optima locaux.
- `u2_ordre(p, nb_restarts=20)` — Même estimation pour des votes par ordre total.
- `phi_dH(p)` — Calcule φ_dH(p) = (2 / n·m) · (u1*(p) − ũ2*(p)) pour un profil par approbation.
- `phi_dS(p)` — Calcule φ_dS(p) = (4 / n·m²) · (u1*(p) − ũ2*(p)) pour un profil par ordre total.

### Graphiques

- `plot_phi2(n, m, nb_points)` — Trace l'évolution de φ² en fonction du paramètre de polarisation pour les deux types de votes (questions 6).
- `experiment(n, m, nb_tests)` — Trace l'évolution de φ_dH et φ_dS en fonction du paramètre de polarisation, en moyennant sur plusieurs profils générés (question 15).

### Fonction principale

- `main()` — Démontre une utilisation minimale du script : génération de profils, calcul des mesures et affichage des graphiques.

---

## Prérequis

- Python 3.8 ou supérieur

---

## Dépendances

| Bibliothèque | Usage |
|---|---|
| `numpy` | Génération aléatoire, calculs vectoriels |
| `scipy` | Résolution du couplage parfait (`linear_sum_assignment`) |
| `matplotlib` | Tracé des graphiques |
| `math` | Calcul de factorielles pour la normalisation de φ² |
| `random` | Sélection aléatoire des centroïdes initiaux dans k-means |

---

## Installation des dépendances

Si les bibliothèques ne sont pas déjà installées :

```bash
pip install numpy scipy matplotlib
```

---

## Lancer le script

```bash
python projet.py
```

L'exécution du `main()` enchaîne les étapes suivantes :

1. Génération d'un profil d'approbation et d'un profil par ordre total avec polarisation intermédiaire.
2. Affichage des bulletins générés.
3. Calcul et affichage des mesures φ², d_H, d_S, u1*, φ_dH, φ_dS sur ces profils.
4. Tracé des courbes d'évolution de φ² en fonction de la polarisation (deux graphiques séparés).
5. Tracé de la courbe d'évolution de φ_dH et φ_dS (graphique unique avec légende).

Les graphiques s'affichent à l'écran via `matplotlib`. Fermer chaque fenêtre pour passer au suivant.

---

## Remarques importantes

### Temps de calcul

- `u1_ordre` et `u2_ordre` sont coûteux sur des grands profils car ils construisent une matrice m×m pour chaque appel à `linear_sum_assignment`. Avec `nb_restarts=20` dans `u2_ordre`, le calcul peut être lent si `n` et `m` sont grands.
- Pour les graphiques (`experiment`), réduire `n`, `m` ou `nb_tests` accélère significativement l'exécution.

### scipy — `linear_sum_assignment`

La fonction `linear_sum_assignment` de `scipy.optimize` résout le **problème d'affectation linéaire** (couplage parfait de poids minimum dans un graphe biparti). Elle est utilisée pour calculer `u1_ordre` et `centroide_ordre`. Sa complexité est O(m³).

### numpy — génération aléatoire

Les fonctions de génération utilisent `np.random` (non seedé par défaut). Les résultats varient à chaque exécution. Pour des résultats reproductibles, ajouter `np.random.seed(42)` en début de script ou au début de `main()`.

### matplotlib — affichage interactif

Les graphiques utilisent `plt.show()`, qui bloque l'exécution jusqu'à la fermeture de la fenêtre. En environnement non interactif (serveur, notebook), remplacer par `plt.savefig("nom.png")`.
