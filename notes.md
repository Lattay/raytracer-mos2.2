# Intro

Raytracing/Lancer de rayon:
Envoyer des rayons dans une scène pour faire un rendu.
Par oposition à la "rasterization".

## Rasterization
Etant donné une géométrie de tris, on projette les tris sur le plan de l'écran
pour les tracer.

Limitations de la rasterization:
Difficile d'ajouter les phénomènes physiques (réflexions, ombres, distance
focale de la caméra...)

Ex: Pour construire les ombres on créer une caméra à la position de chaque
source de lumière pour construire des cartes de profondeures pour chaque
source?
On utilse alors ces cartes de profondeures pour calculer les ombres portées.

## Raytracer

On lance un rayon qui passe par le centre de la caméra et le pixel observé.
On observe le point d'intersection avec la géométrie et on en déduit la
couleur.
Les différents effets sont obtenu en modifiant les trajectoires des rayons.

# Algo

Rayon: origine C et direction u

Intersection avec une sphère de rayon R et de centre O
```
{P / P = C + t.u & ||P - O||^2 = R^2}

{t / t^2.||u||^2 + 2t<u | C - O> + ||C - O||^2 - R^2 = 0
```



