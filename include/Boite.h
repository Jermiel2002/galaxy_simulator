#ifndef BOITE_H
#define BOITE_H

#pragma once

#include "Particule3D.h"

/*
 * Ici, on va représenté tout ce qu'il faut pour notre espace. Autrement dit, notre espace est constitué d'un ensemble de particule.
 * Tous ces particules sont contenues dans un grand cube.
 * Le but est d'avoir une structure de sorte que chaque particule soit contenu dans un sous cube.
 */

/**Avant l'OctreeNode, il faut une boite 3D dont les point2 représentent l'espace qu'occupe une particule dans l'espace 3D*/

/*
 * Une boite dans l'espace tridimensionnel est definie par deux ensembles de coordonnées par convention
 * On a d'une part les coordonnées du coin inférieur gauche de la boite (x, y, z)
 * Et les point2 de la boite d'autre part (largeur, hauteur, profondeur)*/
class Boite
{
public:
    // constructeurs
    Boite();
    Boite(PosParticule3D point1, PosParticule3D point2);

    // fonctions
    // bool inBoite(PosParticule3D p, Boite b);
    // bool intersect(Boite b1, Boite b2);

    // variables
    PosParticule3D point1;
    PosParticule3D point2;

    PosParticule3D GetMin() const;
    PosParticule3D GetMax() const;
    // ParticuleData p; // une boite contient normalement des infos sur sa particule
};

#endif // BOITE_H