#include "BHTree.h"

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
double BHTreeNode::s_theta = 0.9;

/**s_renegades est une tableau dynamique contenant des particules d'un noeud
 * Dans l'algorithme de Barnes-Hut, qui est utilisé pour accélérer le calcul des forces gravitationnelles entre particules dans des simulations physiques,
 * chaque particule est associée à un nœud de l'arbre de subdivision hiérarchique. Cependant, il peut arriver que certaines particules ne soient pas correctement
 * associées à un nœud en raison de diverses conditions ou erreurs. La liste s_renegades est probablement utilisée pour suivre ces particules qui ne sont pas
 * correctement associées à un nœud dans le cadre de l'algorithme.
 */
std::Espace3D<ParticleData> BHTreeNode::s_renegades;

/*s_stat est une structure de type BHTreeNode::DebugStat, qui est une structure utilisé pour stocker des statistiques de débogage.
il contient seulement un membre _nNumCalc qui représente le nombre total de calculs effectués lors de la simulation*/
BHTreeNode::DebugStat BHTreeNode::s_stat = {0};

double BHTreeNode::s_gamma = 0; // gravitational constant is set from the outside

/*s_soft représente le paramètre de softening. Cela ajoute une petite constante pour éviter les singularités
lorsque deux particules sont très proches. La valeur 0.1*0.1 est utilisée ici, représentant environ 3 années lumière.*/
double BHTreeNode::s_soft = 0.1 * 0.1; // approx. 3 light year

/*
dans le contexte de l'algorithme Barnes-Hut, ces coordonnées (min et max) peuvent être utilisées pour délimiter l'espace associé à un nœud de l'arbre.
on crée ainsi une sorte de boite rectangulaire qui contient un noeud

Comment les coordonnées du centre sont calculé?
pour le x : MINx représente la coordonnées minimale sur l'axe x: autrement dis, la distance du Minx par rapport à l'origine.
            MAXx - MINx représente la distance entre les deux coordonnées x des coin inférieur et supérieur du rectangle
            On divise par 2 pour avoir le centre
*/
BHTreeNode::BHTreeNode(const Vec2D &min,
                       const Vec2D &max,
                       BHTreeNode *parent)
    : _particle() // il n'y a aucune particule à l'intérieur
      ,
      _mass(0) // évident que la masse soit null
      ,
      _cm(), _min(min)
      // min et max définissent les coordonnées de la boite qui englobe la particule
      ,
      _max(max), _center(min.x + (max.x - min.x) / 2.0, min.y + (max.y - min.y) / 2.0), _parent(parent) // en créant un noeud s'il a un parent c'est ici qu'il faut le notifié ou jamais
      ,
      _num(0) // on garde une trace du nombre de particule dans ce noeud
      ,
      _bSubdivided(false) // permet de savoir si le noeud a été subdivisé
{
    _quadNode[0] = _quadNode[1] = _quadNode[2] = _quadNode[3] = nullptr; // les quatres fils du noeud qu'on vient de construire
}

bool BHTreeNode::IsRoot() const
{
    return _parent == nullptr;
}

bool BHTreeNode::IsExternal() const
{
    return _quadNode[0] == nullptr &&
           _quadNode[1] == nullptr &&
           _quadNode[2] == nullptr &&
           _quadNode[3] == nullptr;
}

bool BHTreeNode::WasTooClose() const
{
    return _bSubdivided;
}

const Vec2D &BHTreeNode::GetMin() const
{
    return _min;
}

const Vec2D &BHTreeNode::GetMax() const
{
    return _max;
}

const Vec2D &BHTreeNode::GetCenterOfMass() const
{
    return _cm;
}

double BHTreeNode::GetTheta() const
{
    return s_theta;
}

void BHTreeNode::SetTheta(double theta)
{
    s_theta = theta;
}

int BHTreeNode::StatGetNumCalc() const
{
    return s_stat._nNumCalc;
}

/** \brief Returns the number of particles not assigned to any node.
 * La méthode GetNumRenegades retourne le nombre de particules qui n'ont pas été assignées à un nœud spécifique de l'arbre. Dans le contexte de cet arbre Barnes-Hut,
 * lorsqu'une particule ne peut pas être insérée dans un nœud particulier (par exemple, parce que le nœud est déjà occupé par une autre particule ou que la particule
 * se trouve en dehors des limites du nœud), cette particule est ajoutée à un vecteur appelé s_renegades.
 * Ainsi, cette méthode renvoie la taille de ce vecteur, ce qui représente le nombre de particules non assignées.
 */
int BHTreeNode::GetNumRenegades() const
{
    return s_renegades.size();
}

/** \brief Returns the number of particles inside this node. */
int BHTreeNode::GetNum() const
{
    return _num;
}

/*
la fonction BHTreeNode::StatReset() réinitialise certaines statistiques associées à un nœud particulier dans un arbre de subdivision hiérarchique.
Elle vérifie d'abord que la fonction est appelée sur le nœud racine, puis réinitialise les statistiques et les drapeaux de subdivision pour tous les
nœuds de l'arbre.
*/
void BHTreeNode::StatReset()
{
    if (!IsRoot())
        throw std::runtime_error("Only the root node may reset statistics data.");

    s_stat._nNumCalc = 0;

    struct ResetSubdivideFlags
    {
        ResetSubdivideFlags(BHTreeNode *pRoot)
        {
            ResetFlag(pRoot);
        }

        void ResetFlag(BHTreeNode *pNode)
        {
            pNode->_bSubdivided = false;
            for (int i = 0; i < 4; ++i)
            {
                if (pNode->_quadNode[i])
                    ResetFlag(pNode->_quadNode[i]);
            }
        }
    } ResetFlagNow(this);
}

// Réinitialise un noeud
void BHTreeNode::Reset(const Vec2D &min, const Vec2D &max)
{
    if (!IsRoot())
        throw std::runtime_error("Only the root node may reset the tree.");

    // Supprime les sous-nœuds existants et les initialise à nullptr
    for (int i = 0; i < 4; ++i)
    {
        delete _quadNode[i];
        _quadNode[i] = nullptr;
    }

    // Met à jour les coordonnées minimales et maximales de la boîte englobante
    _min = min;
    _max = max;

    // Calcule et met à jour le centre de la boîte englobante
    _center = Vec2D(min.x + (max.x - min.x) / 2.0,
                    min.y + (max.y - min.y) / 2.0);

    // Réinitialise le nombre de particules dans le nœud à 0
    _num = 0;

    // Réinitialise la masse totale du nœud à 0
    _mass = 0;

    // Réinitialise le centre de masse à Vec2D(0, 0)
    _cm = Vec2D(0, 0);

    // Efface la liste de "renegades" (potentiellement des particules qui ne sont pas correctement associées à un nœud)
    s_renegades.clear();
}

BHTreeNode::EQuadrant BHTreeNode::GetQuadrant(double x, double y) const
{
    if (x <= _center.x && y <= _center.y)
    {
        return SW;
    }
    else if (x <= _center.x && y >= _center.y)
    {
        return NW;
    }
    else if (x >= _center.x && y >= _center.y)
    {
        return NE;
    }
    else if (x >= _center.x && y <= _center.y)
    {
        return SE;
    }
    else if (x > _max.x || y > _max.y || x < _min.x || y < _min.y)
    {
        std::stringstream ss;
        ss << "Can't determine quadrant!\n"
           << "particle  : "
           << "(" << x << ", " << y << ")\n"
           << "quadMin   : "
           << "(" << _min.x << ", " << _min.y << ")\n"
           << "quadMax   : "
           << "(" << _max.x << ", " << _max.y << ")\n"
           << "quadCenter: "
           << "(" << _center.x << ", " << _center.y << ")\n";
        throw std::runtime_error(ss.str().c_str());
    }
    else
    {
        throw std::runtime_error("Can't determine quadrant!");
    }
}

/*
 * En résumé, cette fonction est responsable de créer un nouveau nœud pour un quadrant spécifié lors de la subdivision d'un nœud existant dans l'arbre
 * de subdivision hiérarchique de Barnes-Hut.
 */
BHTreeNode *BHTreeNode::CreateQuadNode(EQuadrant eQuad)
{
    switch (eQuad)
    {
    case SW:
        return new BHTreeNode(_min, _center, this);
    case NW:
        return new BHTreeNode(Vec2D(_min.x, _center.y),
                              Vec2D(_center.x, _max.y),
                              this);
    case NE:
        return new BHTreeNode(_center, _max, this);
    case SE:
        return new BHTreeNode(Vec2D(_center.x, _min.y),
                              Vec2D(_max.x, _center.y),
                              this);
    default:
    {
        std::stringstream ss;
        ss << "Can't determine quadrant!\n";
        /*
                   << "particle  : " << "(" << x          << ", " << y          << ")\n"
                     << "quadMin   : " << "(" << _min.x    << ", " << _min.y    << ")\n"
                     << "quadMax   : " << "(" << _max.x    << ", " << _max.y    << ")\n"
                     << "quadCenter: " << "(" << _center.x << ", " << _center.y << ")\n";
        */
        throw std::runtime_error(ss.str().c_str());
    }
    }
}

/**
 * La fonction BHTreeNode::ComputeMassDistribution() est responsable du calcul de la distribution de masse dans le nœud et de la mise à jour du centre d
 * e masse (_cm) en conséquence. Elle est conçue pour être appelée récursivement sur les nœuds de l'arbre de subdivision hiérarchique de Barnes-Hut afin
 * de propager les calculs de masse et de centre de masse à tous les niveaux de l'arbre.
 */
void BHTreeNode::ComputeMassDistribution()
{

    if (_num == 1)
    {
        PODState *ps = _particle._pState;
        PODAuxState *pa = _particle._pAuxState;
        assert(ps);
        assert(pa);

        _mass = pa->mass;
        _cm = Vec2D(ps->x, ps->y);
    }
    else
    {
        _mass = 0;
        _cm = Vec2D(0, 0);

        for (int i = 0; i < 4; ++i)
        {
            if (_quadNode[i])
            {
                _quadNode[i]->ComputeMassDistribution();
                _mass += _quadNode[i]->_mass;
                _cm.x += _quadNode[i]->_cm.x * _quadNode[i]->_mass;
                _cm.y += _quadNode[i]->_cm.y * _quadNode[i]->_mass;
            }
        }

        _cm.x /= _mass;
        _cm.y /= _mass;
    }
}

/** \brief Calculate the accelleration caused by gravitaion of p2 on p1. */
Vec2D BHTreeNode::CalcAcc(const ParticleData &p1, const ParticleData &p2) const
{
    Vec2D acc;

    if (&p1 == &p2)
        return acc;

    // assign references to the variables in a readable form
    const double &x1(p1._pState->x),
        &y1(p1._pState->y);
    const double &x2(p2._pState->x),
        &y2(p2._pState->y),
        &m2(p2._pAuxState->mass);

    double r = sqrt((x1 - x2) * (x1 - x2) +
                    (y1 - y2) * (y1 - y2) + s_soft);
    if (r > 0)
    {
        double k = s_gamma * m2 / (r * r * r);

        acc.x += k * (x2 - x1);
        acc.y += k * (y2 - y1);
    } // if distance is greater zero
    else
    {
        // two particles on the same spot is physical nonsense!
        // nevertheless it may happen. I just need to make sure
        // there is no singularity...
        acc.x = acc.y = 0;
    }

    return acc;
}

/**
 * La fonction BHTreeNode::CalcForce(const ParticleData &p1) const semble calculer la force exercée sur une particule (p1) en prenant en compte à la
 * fois les forces calculées à partir de l'arbre de Barnes-Hut (CalcTreeForce(p1)) et les forces provenant de particules qui ne sont pas correctement
 * associées à un nœud de l'arbre (s_renegades).
 */
Vec2D BHTreeNode::CalcForce(const ParticleData &p1) const
{
    // calculate the force from the barnes hut tree to the particle p1
    Vec2D acc = CalcTreeForce(p1);

    // calculate the force from particles not in the barnes hut tree on particle p
    if (s_renegades.size())
    {
        for (std::size_t i = 0; i < s_renegades.size(); ++i)
        {
            Vec2D buf = CalcAcc(p1, s_renegades[i]);
            acc.x += buf.x;
            acc.y += buf.y;
        }
    }

    return acc;
}

/**  \brief Compute the force acting from this node and it's child
            to a particle p.
    La fonction BHTreeNode::CalcTreeForce(const ParticleData &p1) const semble calculer la force gravitationnelle
    exercée sur la particule p1 en utilisant l'approximation de l'algorithme de Barnes-Hut. Cette approximation permet de réduire le coût
    computationnel en considérant les forces exercées par des groupes de particules plutôt que par chaque particule individuellement.
*/
Vec2D BHTreeNode::CalcTreeForce(const ParticleData &p1) const
{
    Vec2D acc;

    double r(0), k(0), d(0);
    if (_num == 1)
    {
        acc = CalcAcc(p1, _particle);
        s_stat._nNumCalc++;
    }
    else
    {
        r = sqrt((p1._pState->x - _cm.x) * (p1._pState->x - _cm.x) +
                 (p1._pState->y - _cm.y) * (p1._pState->y - _cm.y));
        d = _max.x - _min.x;
        if (d / r <= s_theta)
        {
            _bSubdivided = false;
            k = s_gamma * _mass / (r * r * r);
            acc.x = k * (_cm.x - p1._pState->x);
            acc.y = k * (_cm.y - p1._pState->y);

            // keep track of the number of calculations
            s_stat._nNumCalc++;
        }
        else
        {

            _bSubdivided = true;
            Vec2D buf;
            for (int q = 0; q < 4; ++q)
            {
                if (_quadNode[q])
                {
                    //          const PODState &state = *(p1._pState);
                    buf = _quadNode[q]->CalcTreeForce(p1);
                    acc.x += buf.x;
                    acc.y += buf.y;
                } // if node exists
            }     // for all child nodes
        }
    }

    return acc;
}

/**
 * La fonction BHTreeNode::DumpNode(int quad, int level) semble être une fonction de débogage utilisée pour afficher des informations sur le contenu
 * d'un nœud de l'arbre de Barnes-Hut. Elle parcourt récursivement les sous-nœuds de l'arbre à partir du nœud actuel et affiche différentes informations
 * telles que le quadrant, le nombre de particules dans le nœud, la masse totale du nœud, et les coordonnées du centre de masse.
 */
void BHTreeNode::DumpNode(int quad, int level)
{
    std::string space;
    for (int i = 0; i < level; ++i)
        space += "  ";

    std::cout << space << "Quadrant " << quad << ": ";
    std::cout << space << "(num=" << _num << "; ";
    std::cout << space << "mass=" << _mass << ";";
    std::cout << space << "cx=" << _cm.x << ";";
    std::cout << space << "cy=" << _cm.y << ")\n";

    for (int i = 0; i < 4; ++i)
    {
        if (_quadNode[i])
        {
            _quadNode[i]->DumpNode(i, level + 1);
        }
    }
}

/**
 * La fonction BHTreeNode::Insert(const ParticleData &newParticle, int level) semble être responsable de l'insertion d'une nouvelle particule dans
 * l'arbre de Barnes-Hut. Elle suit le principe de subdivision de l'arbre en quadrants pour organiser les particules de manière hiérarchique. Cette
 * fonction gère les cas où le nœud est vide (_num == 0), contient une seule particule (_num == 1), ou contient plusieurs particules (_num > 1).
 */
void BHTreeNode::Insert(const ParticleData &newParticle, int level)
{
    const PODState &p1 = *(newParticle._pState);
    if ((p1.x < _min.x || p1.x > _max.x) || (p1.y < _min.y || p1.y > _max.y))
    {
        std::stringstream ss;
        ss << "Particle position (" << p1.x << ", " << p1.y << ") "
           << "is outside tree node ("
           << "min.x=" << _min.x << ", "
           << "max.x=" << _max.x << ", "
           << "min.y=" << _min.y << ", "
           << "max.y=" << _max.y << ")";
        throw std::runtime_error(ss.str());
    }

    if (_num > 1)
    {
        EQuadrant eQuad = GetQuadrant(p1.x, p1.y);
        if (!_quadNode[eQuad])
            _quadNode[eQuad] = CreateQuadNode(eQuad);

        _quadNode[eQuad]->Insert(newParticle, level + 1);
    }
    else if (_num == 1)
    {
        assert(IsExternal() || IsRoot());

        const PODState &p2 = *(_particle._pState);

        // This is physically impossible: There are
        // two bodies at the exact same coordinates. In these
        // cases do not add the second body and place
        // it in the renegade Espace3D.
        if ((p1.x == p2.x) && (p1.y == p2.y))
        {
            s_renegades.push_back(newParticle);
        }
        else
        {
            // There is already a particle
            // subdivide the node and relocate that particle
            EQuadrant eQuad = GetQuadrant(p2.x, p2.y);
            if (_quadNode[eQuad] == nullptr)
                _quadNode[eQuad] = CreateQuadNode(eQuad);
            _quadNode[eQuad]->Insert(_particle, level + 1);
            _particle.Reset();

            eQuad = GetQuadrant(p1.x, p1.y);
            if (!_quadNode[eQuad])
                _quadNode[eQuad] = CreateQuadNode(eQuad);
            _quadNode[eQuad]->Insert(newParticle, level + 1);
        }
    }
    else if (_num == 0)
    {
        _particle = newParticle;
    }

    _num++;
}

BHTreeNode::~BHTreeNode()
{
    for (int i = 0; i < 4; ++i)
        delete _quadNode[i];
}
