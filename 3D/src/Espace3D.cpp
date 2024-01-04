#include "../include/Espace3D.hpp"

/*
* Ici on va représenté tout ce qu'il faut pour notre espace. Autrement dit, notre espace est constitué d'un ensemble de particule
* Tous ces particules sont contenue dans un grand cube
* Le but est d'avoir une structure de sorte que chaque particule soit contenu dans un sous cube
*/


/*
* Une boite dans l'espace tridimensionnel est definie par deux ensembles de coordonnées par convention
* On a d'une part les coordonnées du coin inférieur gauche de la boite (x, y, z)
* Et les point2 de la boite d'autre part (largeur, hauteur, profondeur)*/
Boite::Boite(): point1(0,0,0),point2(1,1,1) {}

Boite::Boite(PosParticule3D point1, PosParticule3D point2): Boite(point1,point2){}