Documentation du fichier : profiles.dat
---------------------------------------

Jean-François Mahfouf (23/11/2023)


Ce fichier contient 1535 profils verticaux issus du modèle du CEPMMT dans une configuration
à 60 niveaux verticaux et une troncature TL511 (probablement CY23R4) que j'avais dû créer
en 2001 juste avant de partir au Canada (voir publication Mahfouf et al. (2005) au QJRMS). 
Comme cela fait plus de 20 ans, diverses incertitudes demeurent sur le contenu exact de ce fichier.

Il s'agit probablement de prévisions à courte échéance (6 h). Les profils ont été échantillonnés
sur 8 bandes de latitudes autour de l'équateur entre 1.23°N et -1.23°S (tous les 40 km correspondant
à l'espacement de la grille de Gauss). Sur chaque cercle de latitude sont échantillonnées environ
200 longitudes (sur les 1023 de la grille de Gauss). Cet échantillonnage n'est pas régulier et
doit donc correspondre à des points où le schéma de convection profonde a été activé. 
Le pas de temps de ce modèle était de 900 s (pour estimer les variables évoluées à partir des tendances)

Latitude 1 : nombre de points = 207
Latitude 2 : nombre de points = 214
Latitude 3 : nombre de points = 195
Latitude 4 : nombre de points = 188
Latitude 5 : nombre de points = 191
Latitude 6 : nombre de points = 186
Latitude 7 : nombre de points = 187
Latitude 8 : nombre de points = 167

Le fichier est structuré ainsi pour chaque profil :

1) Header : chaine de 20 caractères indiquant le numéro du profil (1 à 1535)

2) 5 réels : Latitude (degrés), Longitude (degrés), Indice terre/mer (1 terre et 0 pour mer),
Température de surface (K), Pression de surface (Pa)

3) 4 réels dont la signification m'échappe, tous les profils ont les mêmes valeurs :
(2.71) - (-39.72) - (-0.039) - (0.007)

4) 7 tableaux de 60 valeurs (profils verticaux) ordonnés depuis le sommet du modèle vers la surface
- Colonne 1 : Température (K)
- Colonne 2 : Humidité spécifique (kg/kg)
- Colonne 3 : Pression aux niveaux entiers (Pa)
- Colonne 4 : Pression aux demi-niveaux (Pa) - la pression de surface est à un demi-niveau
- Colonne 5 : Vitesse verticale (m/s)
- Colonne 6 : Vitesse du vent horizontal - composante zonale (m/s)
- Colonne 7 : Vitesse du vent horizontal - composante méridienne (m/s)

