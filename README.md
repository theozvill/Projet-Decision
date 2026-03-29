Projet : Mesure de la polarisation dans des profils électoraux

------------------------------------------------------------
Description
------------------------------------------------------------

Ce projet a pour objectif d’étudier la polarisation dans des profils électoraux.
On considère deux types de votes :

- Votes par approbation (vecteurs binaires)
- Votes par ordres totaux (permutations)

On implémente différentes méthodes permettant :

- de générer des profils avec un niveau de polarisation contrôlé
- de définir des distances entre votes (Hamming et Spearman)
- de calculer des mesures de dispersion (u1* et u2*)
- de définir et analyser une mesure de polarisation φ
- d’étudier expérimentalement le comportement de cette mesure

------------------------------------------------------------
Contenu du projet
------------------------------------------------------------

Le projet contient les éléments suivants :

1. Génération de profils
   - generation_profile_approbation
   - generation_profile_ordretotal

2. Distances
   - distance_hamming
   - distance_spearman

3. Calcul de u1*
   - u1_approbation (optimisé en O(nm))
   - u1_ordre (résolution via problème d’affectation)

4. Calcul de u2*
   - u2_approbation (heuristique de type k-médian)
   - u2_ordre (heuristique similaire)

5. Calcul de la mesure de polarisation
   - phi_dH
   - phi_dS

6. Expérimentation
   - tracé de l’évolution de φ en fonction du paramètre de polarisation

------------------------------------------------------------
Choix algorithmiques
------------------------------------------------------------

- Pour u1* en approbation : utilisation d’un vote majoritaire coordonnée par coordonnée.
- Pour u1* en ordre total : réduction à un problème d’affectation résolu avec scipy.optimize.linear_sum_assignment.
- Pour u2* : utilisation d’une heuristique de type k-médian (initialisation aléatoire + itérations).
- Pour φ : application directe de la formule théorique.

------------------------------------------------------------
Complexité
------------------------------------------------------------

- u1_approbation : O(nm)
- u1_ordre : O(nm^2 + m^3)
- u2 (heuristique) : dépend du nombre d’itérations
- φ : coût dominé par u1 et u2

------------------------------------------------------------
Utilisation
------------------------------------------------------------

1. Générer un profil :
   - approbation : generation_profile_approbation(n, m, polarisation)
   - ordre total : generation_profile_ordretotal(n, m, polarisation)

2. Calculer φ :
   - phi_dH(p) pour approbation
   - phi_dS(p) pour ordres

3. Lancer l’expérimentation :
   - utiliser la fonction experiment pour tracer les courbes

------------------------------------------------------------
Dépendances
------------------------------------------------------------

Le projet nécessite les bibliothèques suivantes :

- numpy
- scipy
- matplotlib
- random

Installation possible avec :
pip install numpy scipy matplotlib

------------------------------------------------------------
Remarques
------------------------------------------------------------

- Le calcul de u2* est approché par une heuristique, ce qui peut introduire de légères variations.
- Les résultats expérimentaux sont moyennés pour limiter l’effet du hasard.

------------------------------------------------------------
Auteur
------------------------------------------------------------

Théo EL ZOGHBI-VILLETTE
Iliam SOUAMI
Alessandro UNSWORTH
