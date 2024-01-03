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
double OctreeNode::s_theta = 0.9;

/**s_renegades est une tableau dynamique contenant des particules d'un noeud*/
std::vector<ParticuleData> OctreeNode::s_renegades;

/*s_stat est une structure de type OctreeNode::DebugStat, qui est une structure utilisé pour stocker des statistiques de débogage.
il contient seulement un membre _nNumCalc qui représente le nombre total de calculs effectués lors de la simulation*/
OctreeNode::DebugStat OctreeNode::s_stat = {0};

double OctreeNode::s_gamma = 0; //constante gravitationnelle

/*s_soft représente le paramètre de softening. Cela ajoute une petite constante pour éviter les singularités
lorsque deux particules sont très proches. La valeur 0.1*0.1 est utilisée ici, représentant environ 3 années lumière.*/
double OctreeNode::s_soft = 0.1 * 0.1;
//--------------------------------------------------------------------------------------------------------------------------------

//codage du constructeur

/*Representons un octant : un octant ou cube regroupe une ou plusieurs particules

*-------------calcul du centre de la boite (center)---------------------
* Il faut savoir que Le centre du cube est situé à égale distance le long de chaque axe à partir du sommet inférieur gauche.
* ajouter la moitié de la largeur(dimension x) à la coordonnée (point1 x) du coin inférieur gauche pour obtenir la coordonnée x du centre, 
* faire de même pour les coordonnées y et z en ajoutant la moitié de la hauteur et de la profondeur respectivement.

*--------------Comment repérer l'espace qu'occupe un octant ?-----------------
* Pour un quadtree on utilise les coordonnées min et max respectivement des coins inférieur gauche et supérieur droit
* Pour un OctreeNode, on va utiliser une boite car pour retrouver l'espace qu'il occupe, il suffit de connaitre les coordonnées d'un sommet et les point2
*/
OctreeNode::OctreeNode(const Boite &reper_boite, OctreeNode *parent = nullptr):
_particle(), 
_mass(0), 
_cm(), 
_boite(reper_boite),
_center((reper_boite.point1.x + (reper_boite.point2.x/2)), (reper_boite.point1.y + (reper_boite.point2.y/2)), (reper_boite.point1.z + (reper_boite.point2.z/2))),
_parent(parent),
_num(0),
_bSubdivided(false)
{
    //garantie d'initialisation correcte des fils.
    _noeudFils[0] = _noeudFils[1] = _noeudFils[2] = _noeudFils[3] = _noeudFils[4] = _noeudFils[5] = _noeudFils[6] = _noeudFils[7] = nullptr;
}

//definition des méthodes d'un noeud

bool OctreeNode::IsRoot() const{
    return _parent == nullptr;//renvoie true si le noeud n'a pas de parent
}

bool OctreeNode::IsExternal() const //vérifie si un noeud est une feuille de l'OctreeNode
{
    return _noeudFils[0] == nullptr &&
           _noeudFils[1] == nullptr &&
           _noeudFils[2] == nullptr &&
           _noeudFils[4] == nullptr &&
           _noeudFils[5] == nullptr &&
           _noeudFils[6] == nullptr &&
           _noeudFils[7] == nullptr  ;
}

bool OctreeNode::WasTooClose() const //utilisée pour déterminer si le nœud actuel a été subdivisé en raison d'une proximité trop grande lors du calcul des forces gravitationnelles.
{
    return _bSubdivided;
}

const Boite &OctreeNode::GetBoite() const//représente les point2 de la boite englobant du nœud (la boîte englobante).
{
    return _boite;
}

const PosParticule3D &OctreeNode::GetCenterOfMass() const //renvoie une référence cste vers le centre de masse d'un OctreeNode 3D
{
    return _cm;
}


double OctreeNode::GetTheta() const
{
    return s_theta;
}

void OctreeNode::SetTheta(double theta)
{
    s_theta = theta;
}


int OctreeNode::StatGetNumCalc() const
{
    return s_stat._nNumCalc;
}

/** \brief Returns the number of particles not assigned to any node. */
int OctreeNode::GetNumRenegades() const
{
    return s_renegades.size();
}


/** \brief Returns the number of particles inside this node. */
int OctreeNode::GetNum() const
{
    return _num;
}

/*-------------------------------------------------------------------*/

/*
la fonction OctreeNode::StatReset() 
cette fonction est destinée à réinitialiser les statistiques et les drapeaux de subdivision associés à un arbre Barnes-Hut, dans le contexte d'une simulation physique. 
La partie qui réinitialise les drapeaux de subdivision est réalisée de manière récursive pour tous les nœuds de l'arbre.
Cette opération est utilisée pour marquer tous les nœuds comme non subdivisés, indiquant potentiellement que l'arbre va être reconstruit ou mis à jour pour la prochaine phase de calculs.
*/
void OctreeNode::StatReset()
{
    if(!IsRoot())
        throw std::runtime_error("Only the root node may reset statistics data.");
    
    s_stat._nNumCalc = 0;

    struct ResetSubdivideFlags
    {
        ResetSubdivideFlags(OctreeNode *pRoot)
        {
            ResetFlag(pRoot);
        }

        void ResetFlag(OctreeNode *pNode)
        {
            pNode->_bSubdivided = false;
            for (int i = 0; i < 8; ++i)
            {
                if(pNode->_noeudFils[i])
                    ResetFlag(pNode->_noeudFils[i]);
            }
        }
    } ResetFlagNow(this);
}

/**
 * Voyons la racine de notre OctreeNode comme un grang cube contenant 8 autres sous cubes
 * en fonction du nombre de particules aléatoire dans l'espace considéré.
 * La fonction-ci redessine un grand cube vide qui occupe un espace dont les coordonnées sont fourni
 * Autrement dit, elle met à jour ses coordonnées et ses propriétés, et prépare le nœud pour une nouvelle utilisation dans la simulation physique.
*/
void OctreeNode::Reset(const Boite &reperBoite)
{
    if(!IsRoot())
        throw std::runtime_error("Only the root node may reset the tree.");

    //Supprime les sous-noeuds existants et les initialise à nullptr
    for (int i = 0; i < 8; ++i)
    {
        delete _noeudFils[i];
        _noeudFils[i] =nullptr;
    }

    //Met à jour les coordonnées de la boite englobante
    _boite = reperBoite;

    //Calcule et met à jour le centre de la boite englobante
    _center = PosParticule3D((reperBoite.point1.x + (reperBoite.point2.x/2)), (reperBoite.point1.y + (reperBoite.point2.y/2)), (reperBoite.point1.z + (reperBoite.point2.z/2)));

    //Réinitialise le nombre de particules dans le noeud à 0
    _num = 0;

    //Réinitialise la masse totale du noeud à 0
    _mass = 0;

    //Réinitialise le centre de masse à PosParticule3D(0,0,0)
    _cm = PosParticule3D(0,0,0);

    //Efface la liste de "renegades" (des particules qui ne sont pas correctement associé à un noeud)
    s_renegades.clear();
}

/**
 * 3- Définir une fonction pour insérer une particule dans l'arbre
 * 4- Créer une fonction d'insertion principale
*/

/*Fonction pour déterminer la boite (le quadrant) d'une particule*/
OctreeNode::BoiteAParticule OctreeNode::GetTypeBoite(PosParticule3D &p) const
{
    if(p.x <= _center.x) //SWD,NWD,SWU,NWU
    {
        if(p.y <= _center.y)//SWD,SWU
        {
            if(p.z <= _center.z)
            {
                return BoiteAParticule::SWD;
            }
            else
            {
                return BoiteAParticule::SWU;
            }
        }
        else//NWD,NWU
        {
            if(p.z <= _center.z)
            {
                return BoiteAParticule::NWD;
            }
            else
            {
                return BoiteAParticule::NWU;
            }
        }
    }
    else //SED,NED,SEU,NEU
    {
        if(p.y <= _center.y)//SED,SEU
        {
            if(p.z <= _center.z)
            {
                return BoiteAParticule::SED;
            }
            else
            {
                return BoiteAParticule::SEU;
            }
        }
        else//NED,NEU
        {
            if(p.z <= _center.z)
            {
                return BoiteAParticule::NED;
            }
            else
            {
                return BoiteAParticule::NEU;
            }
        }
    }

     if (p.x < _boite.point1.x || p.y < _boite.point1.y || p.z < _boite.point1.z || p.x > _boite.point2.x || p.y > _boite.point2.y || p.z > _boite.point2.z) //les coordonnée du particule depasse celle de la boite
    {
        std::stringstream ss;
        ss << "Can't determine Octant!\n"
           << "particle  : "
           << "(" << p.x << ", " << p.y << ")\n"
           << "OctantMin   : "
           << "(" << _boite.point1.x << ", " << _boite.point1.y << ", " << _boite.point1.z << ")\n"
           << "OctantMax   : "
           << "(" << _boite.point2.x << ", " << _boite.point2.y << ", " << _boite.point2.z << ")\n"
           << "Center: "
           << "(" << _center.x << ", " << _center.y << _center.z << ")\n";
        throw std::runtime_error(ss.str().c_str());
    }
    else
    {
        throw std::runtime_error("Can't determine Octant!");
    }
}

/**Fonction pour subdiviser un noeud
 * En gros, cette fonction crée les sous noeuds d'un noeud, divisant ainsi l'espace en huit parties égales
 * Cette fonction prend un paramètre boitP qui indique l'octant (la boite) dans lequel un nouveau noeud doit être créé
 * 
 * --> Si c'est une boite sud ouest bas (SWD), les coordonnées du coin inférieur gauche de la boite sont celle du coin inférieur gauche
 * de la boite englobante et les coordonnées représentant les point2, sont les coordonnées du centre de la boite
 * 
 * --> Si c'est une boite nord ouest bas (NWD) les coordonnées du coin inférieur gauche de la boite sont celle du coin inférieur gauche
 * de la boite englobante et les coordonnées représentant les point2, sont les coordonnées du centre de la boite
*/
OctreeNode *OctreeNode::CreateOctreeNodeNode(BoiteAParticule boiteP)
{
    switch (boiteP)
    {
        case SWD:  
            return &OctreeNode(Boite(_boite.point1,_center),this);
        case SWU:
            return &OctreeNode(Boite(PosParticule3D(_boite.point1.x,_boite.point1.y,_center.z),PosParticule3D(_center.x,_center.y,_boite.point2.z)),this);
        case NWD:
            return &OctreeNode(Boite(PosParticule3D(_boite.point1.x,_center.y,_boite.point1.z),PosParticule3D(_center.x,_boite.point2.y,_center.z)),this);
        case NWU:
            return &OctreeNode(Boite(PosParticule3D(_boite.point1.x,_center.y,_center.z),PosParticule3D(_center.x,_boite.point2.y,_boite.point2.z)),this);
        case SED:
            return &OctreeNode(Boite(PosParticule3D(_center.x,_boite.point1.y,_boite.point1.z),PosParticule3D(_boite.point2.x,_center.y,_center.z)), this);
        case SEU:
            return &OctreeNode(Boite(PosParticule3D(_center.x,_boite.point1.y,_center.z),PosParticule3D(_boite.point2.x,_center.y,_boite.point2.z)),this);
        case NED:
            return &OctreeNode(Boite(PosParticule3D(_center.x,_center.y,_boite.point1.z),PosParticule3D(_boite.point2.x,_boite.point2.y,_center.z)),this);
        case NEU:
            return &OctreeNode(Boite(_center,_boite.point2), this);
        default:
        {
            std::stringstream ss;
            ss << "Can't determine octant!\n";
            throw std::runtime_error(ss.str().c_str());
        }
    }
};

/*La fonction DumpNode semble être une fonction utilisée à des fins de débogage (debugging). 
Elle affiche des informations sur le contenu du nœud de l'arbre Barnes-Hut, y compris son quadrant, 
son niveau et les détails sur les particules qu'il contient.*/
void OctreeNode::DumpNode(int oct, int level)
{
    std::string space;
    for (int i = 0; i < level; ++i)
        space += "  ";

    std::cout << space << "Octant" << oct << ": ";
    std::cout << space << "(num=" << _num << "; ";
    std::cout << space << "mass=" << _mass << ";";
    std::cout << space << "cx=" << _cm.x << ";";
    std::cout << space << "cy=" << _cm.y << ")\n";

    for (int i = 0; i < 8; ++i)
    {
        if (_noeudFils[i])
        {
            _noeudFils[i]->DumpNode(i, level + 1);
        }
    }
}


/**
 * La fonction OctreeNode::ComputeMassDistribution() est responsable du calcul de la masse totale de toutes les particules contenues
 * dans les noeuds enfants de chaque noeud de l'arbre ainsi que leur centre de masse. La distribution de masse
 * doit être calculé pour chaque noeud de l'arbre
*/
void OctreeNode::ComputeMassDistribution()
{
    if (_num == 1)
    {
        //on récupère les états de la particule actuel
        EtatParticule *ps = _particle._pState;
        EtatAuxiliaire *pa = _particle._pAuxState;
        assert(ps);
        assert(pa);

        _mass = pa->masse;
        _cm = *(ps->pos);
    }
    else
    {
        _mass = 0;
        _cm = PosParticule3D();

        for (int i=0; i<8; ++i)
        {
            if(_noeudFils[i])
            {
                _noeudFils[i]->ComputeMassDistribution();
                _mass += _noeudFils[i]->_mass;
                _cm.x += (_noeudFils[i]->_particle._pState->pos->x) * _noeudFils[i]->_mass;
                _cm.y += (_noeudFils[i]->_particle._pState->pos->y) * _noeudFils[i]->_mass;
                _cm.z += (_noeudFils[i]->_particle._pState->pos->z) * _noeudFils[i]->_mass;
            }
        }
        _cm.x /= _mass;
        _cm.y /= _mass;
        _cm.z /= _mass;
    }
}

/**
 * La fonction CalcAcc calcule l'accélération gravitationnelle entre deux particules en utilisant une approximation de la loi 
 * de gravité en 3D. Elle prend en compte la distance entre les particules, évite les singularités et retourne le vecteur d'accélération
 * résultant.
 * 
 * Autrement dit, elle calcul l'accélération due à la gravité exercée par la particule p2, sur la
 * particule p1.
*/
PosParticule3D OctreeNode::CalcAcc(const ParticuleData &p1, const ParticuleData &p2) const
{
    PosParticule3D acc;

    //on vérifie si les deux particules p1 et p2 référencent le mpeme objet en mémoire
    if(&p1 == &p2)
        return acc;
    //on crée des références aux coordonnées de la particule p1 d'une part et de la particule p2 d'autre part pour plus de lisibilité
    const double &x1(p1._pState->pos->x), &y1(p1._pState->pos->y), &z1(p1._pState->pos->z);
    const double &x2(p2._pState->pos->x), &y2(p2._pState->pos->y), &z2(p2._pState->pos->z);

    //on crée une référence à la masse de p2
    const double &m2(p2._pAuxState->masse);

    //calcul de la distance entre les particules
    double dist = sqrt(
        (x1 - x2)*(x1 - x2) +
        (y1 - y2)*(y1 - y2) +
        (z1 - z2)*(z1 - z2) + s_soft);

    //calcul de l'accélération gravitationnelle
    if(dist > 0)
    {
        double k = s_gamma * m2 / (dist * dist * dist);

        acc.x += k * (x2 - x1);
        acc.y += k * (y2 - y1);
        acc.z += k * (z2 - z1);
    }
    else
    {
        //cas où les particules sont au même endroit
        acc.x = acc.y = acc.z = 0;
    }

    return acc;
}


/**
 * La fonction CalcTreeForce calcule la force exercée par un noeud de l'arbre(l'octree actuel est ses enfants) sur une particule donnée p1.
*/
PosParticule3D OctreeNode::CalcTreeForce(const ParticuleData &p) const
{
    PosParticule3D acc;

    //variables pour la distance entre la particule et le centre de masse du noeud (dist), une cste de proportionnalité (k) et la largeur de la boite englobante du noeud (cote)
    double dist(0), k(0), cote(0);

    /*cas où le noeud contient une seule particule, la fonction calcule directement la force gravitationnelle
    entre cette particule et la particule p et on incrémente le nombre de calcul*/
    if(_num == 1)
    {
        acc = CalcAcc(p,_particle);
        s_stat._nNumCalc++;
    }
    else //cas où le noeud est subdivisé
    {
        dist = sqrt(
            (p._pState->pos->x - _cm.x)*(p._pState->pos->x - _cm.x) +
            (p._pState->pos->y - _cm.y)*(p._pState->pos->y - _cm.y) +
            (p._pState->pos->z - _cm.z)*(p._pState->pos->z - _cm.z)
        );
        cote = _boite.point2.x - _boite.point1.x;

        /*On vérifie si le rapport entre la largeur de la boite englobante et la distance
        ens inférieur ou égal à un seuil s_theta. Si c'est la cas, le noeud est traité comme une particule unique, et la force gravitationnelle est calculées directement
        comme dans le cas précédent. Sinon, le noeud est considérée comme une entité composite et la force est calculée récursivement en appelant la fonction CalcTreeForce sur 
        ses huits noeuds enfants*/
        if(cote/dist <= s_theta)
        {
            _bSubdivided = false;
            k = s_gamma * _mass / (dist * dist * dist);
            acc.x = k * (_cm.x - p._pState->pos->x);
            acc.y = k * (_cm.y - p._pState->pos->y);
            acc.z = k * (_cm.z - p._pState->pos->z);

            s_stat._nNumCalc++;
        }
        else
        {
            _bSubdivided = true;
            PosParticule3D buf;
            for (int q = 0; q < 8; ++q)
            {
                if(_noeudFils[q])
                {
                    buf = _noeudFils[q]->CalcTreeForce(p);
                    acc.x += buf.x;
                    acc.y += buf.y;
                    acc.z += buf.z;
                }
            }
        }
    }
    return acc;
}


/** \brief Cette fonction calcule la force totale agissant sur une particule donnée p1.
 * Cette force est calcule en combinant la force provenant de l'arbre Barne Hut (calctreeforce)
 * et la force provenant des particules renégates
 */
PosParticule3D OctreeNode::CalcForce(const ParticuleData &p) const
{
    //calculate the force from the barnes hut tree to the particle p
    PosParticule3D acc = CalcTreeForce(p);

    //calculate the force from particles not in the barnes hut tree on particle p
    if (s_renegades.size())
    {
        for(std::size_t i = 0; i < s_renegades.size(); ++i)
        {
            PosParticule3D buf = CalcAcc(p,s_renegades[i]);
            acc.x += buf.x;
            acc.y += buf.y;
            acc.z += buf.z;
        }
    }
    return acc;
}


/**Fonction pour insérer une particule dans l'arbre
 * Cette fonction insère une particule dans l'arbre en parcourant récursivement les sous noeuds
*/
void OctreeNode::InsertParticule(const ParticuleData &newParticule, int level)
{
    //Récuperons les coordonnées de la particule à insérer
    /*const*/ PosParticule3D &p1 = *(newParticule._pState->pos);

    //vérifions que la particule se trouve dans à l'intérieur des limites de la boite
    //englobante du noeud de l'arbre
    if((p1.x < _boite.point1.x || p1.x > _boite.point2.x) || (p1.y < _boite.point1.y || p1.y > _boite.point2.y) || (p1.z < _boite.point1.z || p1.z > _boite.point2.z))
    {
        std::stringstream ss;
        ss << "Particle position ("<<p1.x<<","<<p1.y<<","<<p1.z<<")"
        << "is outside tree node ("
        <<"min.x="<<_boite.point1.x<<","
        <<"max.x="<<_boite.point2.x<<","
        <<"min.y="<<_boite.point1.y<<","
        <<"max.y="<<_boite.point2.y<<","
        <<"min.z="<<_boite.point1.z<<","
        <<"max.z="<<_boite.point2.z<<",";
        throw std::runtime_error(ss.str());
    }

    if(_num > 1)
    {
        //determinons le quadrant approprié pour la new particule
        BoiteAParticule cubeDeLaParticule = GetTypeBoite(p1);
        if(!_noeudFils[cubeDeLaParticule])
            _noeudFils[cubeDeLaParticule] = CreateOctreeNodeNode(cubeDeLaParticule);
        _noeudFils[cubeDeLaParticule]->InsertParticule(newParticule, level + 1);

    }
    else if (_num == 1)
    {
        assert(IsExternal() || IsRoot());

        PosParticule3D &p2 = *(_particle._pState->pos);

        //const EtatParticule &p2 = *(_particle._pState);

        /**Si la nouvelle particule a les mêmes coordonnées que la particule déjà présente dans le noeud, cela est considéré comme impossible
         * et la nouvelle particule est ajoutée au vecteur de renegates. Sinon le noeud doit etre subdivisé
        */
        if((p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z))
        {
            s_renegades.push_back(newParticule);
        }
        else
        {
            BoiteAParticule cubeDeLaParticule = GetTypeBoite(p2);
            if(!_noeudFils[cubeDeLaParticule])
                _noeudFils[cubeDeLaParticule] == CreateOctreeNodeNode(cubeDeLaParticule);
            _noeudFils[cubeDeLaParticule]->InsertParticule(_particle,level+1);
            _particle.Reset();

            cubeDeLaParticule = GetTypeBoite(p1);
            if(!_noeudFils[cubeDeLaParticule])
                _noeudFils[cubeDeLaParticule] = CreateOctreeNodeNode(cubeDeLaParticule);
            _noeudFils[cubeDeLaParticule]->InsertParticule(newParticule,level+1);
        }
    }
    else if (_num == 0)
    {
        _particle = newParticule;
    }
    _num++;
};

OctreeNode::~OctreeNode()
{
    for (int i = 0; i < 8; i++)
        delete _noeudFils[i];
}
