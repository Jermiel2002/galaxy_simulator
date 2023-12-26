#include "../include/Octree.hpp"

//--- Standard includes --------------------------------------------------------
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>

//------------------------------------------------------------------------------
// static variables

/*s_theta est une variable statique qui représente le paramètre 0 de l'algo de barnes hut. Il est utilisé pour
* décider si un noeud doit être approximé comme une particule unique lors du calcul des forces.*/
double Octree::s_theta = 0.9;

/**s_renegades est une tableau dynamique contenant des particules d'un noeud*/
std::vector<ParticuleData> Octree::s_renegades;

/*s_stat est une structure de type Octree::DebugStat, qui est une structure utilisé pour stocker des statistiques de débogage.
il contient seulement un membre _nNumCalc qui représente le nombre total de calculs effectués lors de la simulation*/
Octree::DebugStat Octree::s_stat = {0};

double Octree::s_gamma = 0; //constante gravitationnelle

/*s_soft représente le paramètre de softening. Cela ajoute une petite constante pour éviter les singularités
lorsque deux particules sont très proches. La valeur 0.1*0.1 est utilisée ici, représentant environ 3 années lumière.*/
double Octree::s_soft = 0.1 * 0.1;

//codage du constructeur

/*créer un noeud (une feuille) de l'octree*/
Octree::Octree(const PosParticule3D &min, const PosParticule3D &max, Octree *parent): 
_particle(), 
_mass(0), 
_cm(), 
_min(min), 
_max(max),
_center(/*a calculer*/),
_parent(parent),
_num(0),
_bSubdivided(false)
{
    _octreeNode[0] = _octreeNode[1] = _octreeNode[2] = _octreeNode[3] = _octreeNode[4] = _octreeNode[5] = _octreeNode[6] = _octreeNode[7] = nullptr;
}

bool Octree::IsRoot() const{
    return _parent == nullptr;//renvoie true si le noeud n'a pas de parent
}

bool Octree::IsExternal() const //vérifie si un noeud est une feuille de l'octree
{
    return _octreeNode[0] == nullptr &&
           _octreeNode[1] == nullptr &&
           _octreeNode[2] == nullptr &&
           _octreeNode[4] == nullptr &&
           _octreeNode[5] == nullptr &&
           _octreeNode[6] == nullptr &&
           _octreeNode[7] == nullptr  ;
}

/*-------------Accès aux propriété spécifiques d'un noeud-----------------*/
bool Octree::WasTooClose() const //utilisée pour déterminer si le nœud actuel a été subdivisé en raison d'une proximité trop grande lors du calcul des forces gravitationnelles.
{
    return _bSubdivided;
}

const PosParticule3D &Octree::GetMin() const//représente les coins de la boite englobant du nœud (la boîte englobante).
{
    return _min;
}


const PosParticule3D &Octree::GetMax() const
{
    return _max;
}


const PosParticule3D &Octree::GetCenterOfMass() const //Le centre de masse en 3D
{
    return _cm;
}


double Octree::GetTheta() const
{
    return s_theta;
}

void Octree::SetTheta(double theta)
{
    s_theta = theta;
}


int Octree::StatGetNumCalc() const
{
    return s_stat._nNumCalc;
}

/** \brief Returns the number of particles not assigned to any node. */
int Octree::GetNumRenegades() const
{
    return s_renegades.size();
}


/** \brief Returns the number of particles inside this node. */
int Octree::GetNum() const
{
    return _num;
}
/*-------------------------------------------------------------------*/

