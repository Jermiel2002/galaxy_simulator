#ifndef _ESPACE3D_H 
#define _ESPACE3D_H 

#pragma once

#include <iostream>
#include "math.h"
#include <stdbool.h>
#include "Particule3D.hpp"


/*
* Ici on va représenté tout ce qu'il faut pour notre espace. Autrement dit, notre espace est constitué d'un ensemble de particule
* Tous ces particules sont contenue dans un grand cube
* Le but est d'avoir une structure de sorte que chaque particule soit contenu dans un sous cube
*/


/*
* Une boite dans l'espace tridimensionnel est definie par deux ensembles de coordonnées par convention
* On a d'une part les coordonnées du coin inférieur gauche de la boite (x, y, z)
* Et les dimensions de la boite d'autre part (largeur, hauteur, profondeur)*/
struct Boite
{
    //constructeurs
    Boite(double x, double y, double z, double w, double h, double d);
    Boite();

    //fonctions
    bool inBoite(PosParticule3D p, Boite b);
    bool intersect(Boite b1, Boite b2);

    //variables
    double x, y, z;
    double w,h,d;
};

struct Octree
{
    //variables
    Boite cube;
    int size;
    int id;
    /*
        La ligne suivante déclare un tableau tridimensionnel de pointeurs vers des structures Octree. Cela représente les huit sous-arbres (octants) de l'arbre actuel. Chaque dimension du tableau peut avoir deux valeurs, indiquant ainsi les huit subdivisions possibles de l'espace tridimensionnel. Un pointeur à nullptr indique qu'une subdivision particulière n'est pas encore créée.
    */
    Octree *substrees[2][2][2];

    //fonction
    bool is_divided;
    /*newOctree qui crée et initialise un nouvel Octree en allouant de la mémoire dynamique pour la structure Octree et en initialisant ses propriétés. */
    Octree* newOctree(Boite b);
};

#endif
