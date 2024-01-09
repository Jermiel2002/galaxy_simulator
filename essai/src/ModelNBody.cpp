#include "ModelNBody.hpp"

//--- Standard includes --------------------------------------------------
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <omp.h>

using namespace std;

ModelNBody::ModelNBody():
                            _pInitial(nullptr),
                            _pAux(nullptr),
                            _reperBoite(),
                            _root(OctreeNode(Boite(PosParticule3D(),PosParticule3D()))),
                            _pointMin(_reperBoite.point1),
                            _pointMax(_reperBoite.point2),
                            _center(),
                            _roi(1),
                            _num(0),
                            _bVerbose(true)
{}

ModelNBody::~ModelNBody()
{
    delete _pInitial;
    delete _pAux;
}

void ModelNBody::SetROI(double roi)
{
    _roi = roi;
}

/*double ModelNBody::GetSuggestedTimeStep() const{
    return _timeStep;
}

double ModelNBody::GetROI() const
{
    return _roi;
}
*/
PosParticule3D ModelNBody::GetCenterOfMass() const
{
    const PosParticule3D &cm3d = _root.GetCenterOfMass();
    return cm3d;
}
/*
const PosParticule3D &ModelNBody::GetCamDir() const
{
    return _camDir;
}

const PosParticule3D &ModelNBody::GetCamPos() const
{
    return _camPos;
}

double *ModelNBody::GetInitialState()
{
    return reinterpret_cast<double *>(_pInitial);
}
*/
/**
 *  GetOrbitalVelocity calcule la vitesse orbitale nécessaire pour maintenir une orbite circulaire entre deux particules dans un système n-body, 
 *  en utilisant la formule de la troisième loi de Kepler.
*/
void ModelNBody::GetOrbitalVelocity(const ParticuleData &p1, const ParticuleData &p2)
{
    double x1 = p1._pState->pos.x, y1 = p1._pState->pos.y, z1 = p1._pState->pos.z, m1 = p1._pAuxState->masse;
    double x2 = p2._pState->pos.x, y2 = p2._pState->pos.y, z2 = p2._pState->pos.z;

    //Calcul de la distance entre les deux particules
    double r[3], dist;
    r[0] = x1 - x2;
    r[1] = y1 - y2;
    r[2] = z1 - z2;

    //disantance in parsec
    dist = sqrt((r[0] * r[0]) + (r[1] * r[1]) + (r[2] * r[2]));

    //Calcul de la vitesse orbitale nécessaire
    double v = sqrt(gamma_1 * m1/dist);

    //Calcul du vecteur de vitesse perpendiculaire au vecteur de distance 
    double &vx = p2._pState->vitesse.x, &vy = p2._pState->vitesse.y, &vz = p2._pState->vitesse.z;
    vx = (r[1]/dist) * v;
    vy = (-r[0]/dist) * v;
    vz = (-r[2]/dist) * v;
}

void ModelNBody::ResetDim(int num, double stepsize)
{
    _num = num;
    //SetDim(_num * 8);

    delete _pInitial;
    _pInitial = new EtatParticule[num];

    delete _pAux;
    _pAux = new EtatAuxiliaire[num];

    _timeStep = stepsize;

    //reset bounding box and center
    _pointMax.x = _pointMax.y = _pointMax.z = std::numeric_limits<double>::min();
    _pointMin.x = _pointMin.y = _pointMin.z = std::numeric_limits<double>::max();
    _center = PosParticule3D(0,0,0);
}

/*Init initialise un modèle de simulation N-Body avec des particules en utilisant un algorithme spécifique*/
/*void ModelNBody::Init()
{
    //Reset model size
    ResetDim(5000, 100000);

    double mass = 0; // utilisée pour stocker la masse totale du système.

    //Une double boucle (for) est utilisée pour initialiser les particules.
    int ct = 0;
    ParticuleData blackHole/*particule spéciale trou noir*///, macho[10]/*tableau de particules*/;
/*
    /**
     * Les trois boucle imbriquées sont utilisées pour initialiser les particules en fonction des
     * indices k, l. Ces indices controlent la configuration spatiale des particules dans le modèle
     * 
     * k controle le nombre d'anneuax concentriques. Prenons au total 40.
     * l controle la position des particules à l'intérieur de chaque anneau
     * 
     * Chaque groupe correspond à un anneau concentrique de particules autour du centre du modèle
     * Autrement dit, les particules sont positionnées de manière à former un cercle ou un anneau avec le centre du 
     * modèle comme point central
     * Plus précisement, chaque anneau est constitué de particules qui se trouvent à une distance fixe du centre du modèle
     * Les particules à l'intérieur d'un même anneau sont positionnées de manière équidistante les unes des autres sur le cercle ou l'anneau
    */
    /*for(int k = 0; k < 40; ++k)
    {
        for (int l = 0; l < 100; ++l)
        {
            if (ct >= _num)//si le nombre particules inséré plus grand que le nombre de particules défini(_num) dans notre modèle
                goto hell;

            /*À chaque itération des boucles, les références st et st_aux sont utilisées pour accéder aux états initial et auxiliaire d'une particule.*/
            /*EtatParticule &st = _pInitial[ct];
            EtatAuxiliaire &st_aux = _pAux[ct];

            if (ct == 0)//La première particule est un trou noir spécial avec une masse importante, positionné au centre.
            {
                blackHole._pState = &st;
                blackHole._pAuxState = &st_aux;

                /*La particule 0 est un trou noir (objet blackHole), positionné au centre (x = 0, y = 0, z = 0) avec une masse élevée (st_aux.mass = 1000000).*/
                /*st.pos.x = st.pos.y = st.pos.z = 0;
                st.vitesse.x = st.vitesse.y = st.vitesse.z = 0;
                st_aux.masse = 1000000; // 431000;   // 4.31 Millionen Sonnenmassen
            }
            else if (ct == 1)
            {
                macho[0]._pState = &st;
                macho[0]._pAuxState = &st_aux;

                st_aux.masse = blackHole._pAuxState->masse / 10.0;
                st.pos.x = 5000;
                st.pos.y = 5000;
                st.pos.z = 5000;

                GetOrbitalVelocity(blackHole, ParticuleData(&st,&st_aux));
            }
            else if (ct == 2)
            {
                macho[1]._pState = &st;
                macho[1]._pAuxState = &st_aux;

                st_aux.masse = blackHole._pAuxState->masse / 10.0;
                st.pos.x = -5000;
                st.pos.y = -5000;
                st.pos.z = -5000;

                GetOrbitalVelocity(blackHole, ParticuleData(&st,&st_aux));
            }
            else
            {
                /**
                 * Ici on va utilisé les coordonnées sphériques pour positionner les particules. Les coordonnées spériques nécessite trois paramètres:
                 * la distance radiale(r), l'angle azimutal(theta) et l'angle polaire(phi)
                */
                /*st_aux.masse = 0.76 + 100 * ((double)rand() / RAND_MAX);
                double rad = 1200 + k*100;//La distance radiale de la particule courante par rapport au centre est calculée en fonction de l'indice kk
                double theta = 2 * M_PI * l /100;
                double phi = M_PI /2.0 ;//angle polaire (90 degrés pour rester dans le plan)
                st.pos.x = rad * sin(phi) * cos(theta);
                st.pos.y = rad * sin(phi) * sin(theta);
                st.pos.z = rad * cos(phi);
                GetOrbitalVelocity(blackHole, ParticuleData(&st,&st_aux));
            }

            // determine the size of the area including all particles
            _pointMax.x = std::max(_pointMax.x,st.pos.x);
            _pointMax.y = std::max(_pointMax.y, st.pos.y);
            _pointMax.z = std::max(_pointMax.z,st.pos.z);
            _pointMin.x = std::max(_pointMin.x,st.pos.x);
            _pointMin.y = std::max(_pointMin.y, st.pos.y);
            _pointMin.z = std::max(_pointMin.z,st.pos.z);

            _center.x += st.pos.x * st_aux.masse;
            _center.y += st.pos.y * st_aux.masse;
            _center.z += st.pos.z * st_aux.masse;
            mass += st_aux.masse;
            ++ct;
        }
    }*/
    /**
     * Si le nombre particules inséré plus grand que le nombre de particules défini(_num) dans notre modèle alors on redéfinit
     * la région d'intérêt (_ROI) et ajuster la position des particules dans cette région*/
 /*   hell:
        //compute the center of mass
        _center.x /= mass;
        _center.y /= mass;
        _center.z /= mass;

        // The Barnes Hut algorithm needs square shaped quadrants.
        // calculate the height of the square including all particles (and a bit more space)
        _roi = 1.5 * std::max(std::max(_pointMax.x - _pointMin.x, _pointMax.y - _pointMin.y),_pointMax.z - _pointMin.z);

        // compute the center of the region including all particles
        _pointMin.x = _center.x - _roi;
        _pointMax.x = _center.x + _roi;
        _pointMin.y = _center.y - _roi;
        _pointMax.y = _center.y + _roi;
        _pointMin.z = _center.z - _roi;
        _pointMax.z = _center.z + _roi;
        
    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin   = " << _pointMin.x << ", ymin=" << _pointMin.y << ", zmin=" << _pointMin.z << "\n";
    std::cout << "  xmax   = " << _pointMax.y << ", ymax=" << _pointMax.y << ", zmax=" << _pointMax.z << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  center = " << _center.x << ", cy  =" << _center.y << _center.z << "\n";
    std::cout << "  roi    = " << _roi << "\n";        
}*/

/**
 * Dans une sphère en trois dimensions (3D), les coordonnées sont généralement spécifiées en utilisant les coordonnées sphériques. Les coordonnées sphériques comprennent trois composantes : la distance radiale rr, l'angle polaire θθ (mesuré à partir de l'axe zz), et l'angle azimutal ϕϕ (mesuré dans le plan xyxy).

Les conversions entre les coordonnées cartésiennes (x,y,z)(x,y,z) et les coordonnées sphériques (r,θ,ϕ)(r,θ,ϕ) sont données par les équations suivantes :

x=rsin⁡(θ)cos⁡(ϕ)x=rsin(θ)cos(ϕ)
y=rsin⁡(θ)sin⁡(ϕ)y=rsin(θ)sin(ϕ)
z=rcos⁡(θ)z=rcos(θ)

Dans ce système de coordonnées, rr est la distance radiale de l'origine au point, θθ est l'angle polaire (mesuré à partir de l'axe zz), et ϕϕ est l'angle azimutal (mesuré dans le plan xyxy).

L'angle polaire (θθ) est généralement mesuré par rapport à l'axe zz dans un système de coordonnées sphériques. Si l'angle azimutal est mesuré autour de l'axe zz dans le plan xyxy, l'angle polaire peut être généré de manière aléatoire entre 00 et ππ (ou 00 à 180180 degrés). Voici comment vous pourriez le calculer :

double theta=π⋅((double)rand()/RAND_MAX)double theta=π⋅((double)rand()/RAND_MAX)

Cela vous donnera un angle polaire aléatoire compris entre 00 et ππ.
*/
/*
void ModelNBody::InitCollision()//initialise les particules du modèle pour simuler une collision entre deux trous noirs (black holes)
{
    ResetDim(5000, 100);

    ParticuleData blackHole;
    ParticuleData blackHole2;

    for (int i = 0; i < _num; ++i)
    {
        EtatParticule &st = _pInitial[i];
        EtatAuxiliaire &st_aux = _pAux[i];

        if(i == 0)//premiere particule, elle est spéciale et représente une trace qui n'est pas partie de la simulation principale.
        {
            blackHole._pState = &st;
            blackHole._pAuxState = &st_aux;

            st.pos.x = st.pos.y = st.pos.z = 0;
            st.vitesse.x = st.vitesse.y = st.vitesse.z = 0;
            st_aux.masse = 1000000; // 431000;   // 4.31 Millionen Sonnenmassen
        }
        else if(i < 4000)//initalisation des particules dans la région externe, en leur attribuant des orbites aléatoires
        {
            const double rad = 10;
            double r = 0.1 + .8 * (rad * ((double)rand() / RAND_MAX));
            double theta = 2.0 * M_PI * ((double)rand() / RAND_MAX);
            double phi = M_PI * ((double)rand()/RAND_MAX);//angle polaire 3D (angle vertical)

            st_aux.masse = 0.03 + 20 * ((double)rand() / RAND_MAX);//entre [0.03,20]

            st.pos.x = r * sin(phi)*cos(theta);
            st.pos.y = r * sin(phi)*sin(theta);
            st.pos.z = r * cos(phi);


            GetOrbitalVelocity(blackHole, ParticuleData(&st, &st_aux));
        }
        else if (i == 4000)//une deuxième particule est créée avec une masse importante et une orbite autour du centre
        {
            blackHole2._pState = &st;
            blackHole2._pAuxState = &st_aux;

            st.pos.x = st.pos.y = st.pos.z = 10;
            st_aux.masse = 100000;
            GetOrbitalVelocity(blackHole, blackHole2);
            blackHole2._pState->vitesse.x *= 0.9;
            blackHole2._pState->vitesse.y *= 0.9;
            blackHole2._pState->vitesse.z *= 0.9;
        }
        else//Ces particules sont générées dans une région interne et sont influencées par le deuxième trou noir. Leur vitesse est ajustée en fonction de la vitesse du trou noir.
        {
            const double rad = 3;
            double r = 0.1 + .8 * (rad * ((double)rand() / RAND_MAX));
            double theta = 2.0 * M_PI * ((double)rand() / RAND_MAX);
            double phi = M_PI * ((double)rand()/RAND_MAX);//angle polaire 3D

            st_aux.masse = 0.03 + 20 * ((double)rand() / RAND_MAX);

            st.pos.x = blackHole2._pState->pos.x + r * sin(phi)*cos(theta);
            st.pos.y = blackHole2._pState->pos.y + r * sin(phi)*sin(theta);
            st.pos.z = blackHole2._pState->pos.x + r * cos(phi);

            GetOrbitalVelocity(blackHole2, ParticuleData(&st, &st_aux));

            st.vitesse.x += blackHole2._pState->vitesse.x;
            st.vitesse.y += blackHole2._pState->vitesse.y;
            st.vitesse.z += blackHole2._pState->vitesse.z;
        }

            // determine the size of the area including all particle
            _pointMax.x = std::max(_pointMax.x,st.pos.x);
            _pointMax.y = std::max(_pointMax.y, st.pos.y);
            _pointMax.z = std::max(_pointMax.z,st.pos.z);
            _pointMin.x = std::max(_pointMin.x,st.pos.x);
            _pointMin.y = std::max(_pointMin.y, st.pos.y);
            _pointMin.z = std::max(_pointMin.z,st.pos.z);
    }

    double l = 1.05 * std::max(std::max(_pointMax.x - _pointMin.x, _pointMax.y - _pointMin.y),_pointMax.z - _pointMin.z);
    _roi = l * 1.5;

    //centre de mass de la region incluant tous les particules
    PosParticule3D c(_pointMin.x + (_pointMax.x - _pointMin.x) / 2.0, _pointMin.y + (_pointMax.y - _pointMin.y) / 2.0, _pointMin.z + (_pointMax.z - _pointMin.z) / 2.0);

    _pointMin.x = c.x - l/2.0;
    _pointMax.x = c.x + l/2.0;
    _pointMin.y = c.y - l/2.0;
    _pointMax.y = c.y + l/2.0;
    _pointMin.z = c.z - l/2.0;
    _pointMax.z = c.z + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << _pointMin.x << ", ymin=" << _pointMin.y << ", zmin=" << _pointMin.z << "\n";
    std::cout << "  xmax=" << _pointMax.x << ", ymax=" << _pointMax.y << ", zmax=" << _pointMax.z << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x << ", cy  =" << c.y << ", cz  =" << c.z << "\n";
    std::cout << "  l  =" << l << "\n";
}*/

/**
 * Init3Body initialise le modèle pour simuler un système à trois corps. Chacun des trois corps est positionné avec des masses spécifiques
 * et des positions initiales. Elle ajuste egalement la taille du modèle en fonction de ces corps
*/
/*
void ModelNBody::Init3Body()
{
    ResetDim(3, .5);//réinitialisation de la taille du modèle
    _root.SetTheta(0.9);
    EtatParticule *st(nullptr);
    EtatAuxiliaire *st_aux(nullptr);

    st = &_pInitial[0];//st pointe vers la première particule du modèle
    st_aux = &_pAux[0];
    st->pos.x = 1;
    st->pos.y = 3;
    st->pos.z = 5;
    st->vitesse.x = st->vitesse.y = st->vitesse.z = 0;
    st_aux->masse = 3;

    st = &_pInitial[1];
    st_aux = &_pAux[1];
    st->pos.x = -2;
    st->pos.y = 1;
    st->pos.z = 3;
    st->vitesse.x = st->vitesse.y = st->vitesse.z = 0;
    st_aux->masse = 4;

    st = &_pInitial[2];
    st_aux = &_pAux[2];
    st->pos.x = 1;
    st->pos.y = -1;
    st->pos.z = 2;
    st->vitesse.x = st->vitesse.y = st->vitesse.z = 0;
    st_aux->masse = 5;

    // determine the size of the area including all particles
    for (int i = 0; i < _num; ++i)
    {
        EtatParticule &st = _pInitial[i];
        _pointMax.x = std::max(_pointMax.x,st.pos.x);
        _pointMax.y = std::max(_pointMax.y,st.pos.y);
        _pointMax.z = std::max(_pointMax.z,st.pos.z);

        _pointMin.y = std::max(_pointMax.y,st.pos.x);
        _pointMin.y = std::max(_pointMin.y,st.pos.y);
        _pointMin.y = std::max(_pointMin.y,st.pos.z);
    }

    double l = 1.05 * std::max(std::max(_pointMax.x - _pointMin.x, _pointMax.y - _pointMin.y),_pointMax.z - _pointMin.z);
    _roi = l*1.5;

    PosParticule3D c(_pointMin.x + (_pointMax.x - _pointMin.x) / 2.0, _pointMin.y + (_pointMax.y - _pointMin.y) / 2.0, _pointMin.z + (_pointMax.z - _pointMin.z) / 2.0);

    _pointMin.x = c.x - l/2.0;
    _pointMax.x = c.x + l/2.0;
    _pointMin.y = c.y - l/2.0;
    _pointMax.y = c.y + l/2.0;
    _pointMin.z = c.z - l/2.0;
    _pointMax.z = c.z + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << _pointMin.x << ", ymin=" << _pointMin.y << ", zmin=" << _pointMin.z << "\n";
    std::cout << "  xmax=" << _pointMax.x << ", ymax=" << _pointMax.y << ", zmax=" << _pointMax.z << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x << ", cy  =" << c.y << ", cz  =" << c.z << "\n";
    std::cout << "  l  =" << l << "\n";
}*/

//void ModelNBody::CalcBHArea(const ParticuleData &data)
//{
    /*
    _pointMax.x = _pointMax.y = std::numeric_limits<double>::min();
    _pointMin.x = _pointMin.y = std::numeric_limits<double>::max();

    for(int i=0; i<_num; ++i)
    {
        EtatParticule &s = data._pState[i];
        _pointMax.x = std::max(_pointMax.x,s.pos.x);
        _pointMax.y = std::max(_pointMax.y,s.pos.y);
        _pointMax.z = std::max(_pointMax.z,s.pos.z);

        _pointMin.y = std::max(_pointMax.y,s.pos.x);
        _pointMin.y = std::max(_pointMin.y,s.pos.y);
        _pointMin.y = std::max(_pointMin.y,s.pos.z);   
    }

    double l = 1.05 * std::max(std::max(_pointMax.x - _pointMin.x, _pointMax.y - _pointMin.y),_pointMax.z - _pointMin.z);
    _roi = l*1.5;

    PosParticule3D c(_pointMin.x + (_pointMax.x - _pointMin.x) / 2.0, _pointMin.y + (_pointMax.y - _pointMin.y) / 2.0, _pointMin.z + (_pointMax.z - _pointMin.z) / 2.0);

    _pointMin.x = c.x - l/2.0;
    _pointMax.x = c.x + l/2.0;
    _pointMin.y = c.y - l/2.0;
    _pointMax.y = c.y + l/2.0;
    _pointMin.z = c.z - l/2.0;
    _pointMax.z = c.z + l/2.0;*/
//}

/*
* \brief Build the barnes hut tree by adding all particles that are inside
           the region of interest.
* responsable de la construction de l'arbre Barnes-Hut en ajoutant toutes les particules qui se trouvent à l'intérieur de la région d'intérêt (roi).
*/
void ModelNBody::BuiltTree(const ParticuleData &all)
{
    // Reset the octree, make sure only particles inside the roi
    // are handled. The renegade ones may live long and prosper
    // outside my simulation
    _root.Reset(Boite(PosParticule3D(_center.x -_roi, _center.y -_roi, _center.z -_roi),PosParticule3D(_center.x + _roi, _center.y + _roi, _center.z + _roi)));

    //build the octree
    int ct = 0;
    for(int i = 0; i < _num; ++i)
    {
        try
        {
            {
                /**Pour chaque particule, une instance de ParticleData est créée en utilisant les données de position et auxiliaires de la particule actuelle.
                 * "all" est un objet de type ParticleData qui est passé à la fonction BuiltTree en tant que paramètre. Il contient toutes les données d'état (_pState) et auxiliaires (_pAuxState) de toutes les particules de la simulation.
                */
                ParticuleData p(&(all._pState[i]), &(all._pAuxState[i]));

                //insert the particle, but only if its inside the roi
                _root.InsertParticule(p,0);
                ++ct;
            }

        }
        catch(const std::exception& exc)
        {
            std::cout << exc.what() << "\n";
            std::cout << "Particle " << i << " (" << all._pState->pos.x << ", " << all._pState->pos.y << ", " << all._pState->pos.z << ") is outside the roi (skipped).\n";
            std::cout << "  roi size   =   " << _roi << "\n";
            std::cout << "  roi center = (" << _center.x << ", " << _center.y << ", " << _center.z << ")\n";
        }     
    }
    std::cout << ct << " particles added succesfully\n";

    // compute masses and center of mass on all scales of the tree
    _root.ComputeMassDistribution();

    if (_bVerbose)//si _bVerbose est true, un affichage détaillé de l'arbre est réalisé en utilisant la méthode _root.DumpNode. Cela peut être utile pour le débogage et la visualisation de la structure de l'arbre.
    {
        std::cout << "Tree Dump\n";
        std::cout << "---------\n";
        _root.DumpNode(-1,0);
        std::cout << "\n\n";
    }

    //Update the center of mass
    _center = _root.GetCenterOfMass();
}

const EtatAuxiliaire *ModelNBody::GetAuxState() const
{
    return _pAux;
}

OctreeNode *ModelNBody::GetRootNode()
{
    return &_root;
}

int ModelNBody::GetN() const
{
    return _num;
}

double ModelNBody::GetTheta() const
{
    return _root.GetTheta();
}

void ModelNBody::SetVerbose(bool bVerbose)
{
    _bVerbose = bVerbose;
}


void ModelNBody::SetTheta(double theta)
{
    _root.SetTheta(theta);
}

/**
 * ModelNBody::Eval, evalue l'accélération des particules à un instant donné dans le temps.
*/
/*void ModelNBody::Eval(double *a_state, double a_time, double *a_deriv)
{
    /**Conversion des pointeurs d'état et de dérivées
     * Les lignes ci-dessous convertissent les pointeurs bruts a_state et a_deriv en pointeurs typés (PODState et PODDeriv). 
     * Ensuite, elle crée une instance de ParticleData appelée all en utilisant les données d'état de toutes les particules et les données auxiliaires de la classe ModelNBody.
    */
   /* EtatParticule *pState = reinterpret_cast<EtatParticule *>(a_state);
    MajEtat *pDeriv = reinterpret_cast<MajEtat *> (a_deriv);
    ParticuleData all(pState,_pAux);

    CalcBHArea(all);//calcule l'aire nécessaire pour la méthode de Barnes-Hut en utilisant les données de toutes les particules.
    BuiltTree(all);

    #pragma omp parallel for
        for (int i = 1; i < _num; ++i)
        {
            ParticuleData p(&pState[i],&_pAux[i]);
            PosParticule3D acc = _root.CalcForce(p);
            pDeriv[i].acceleration.x = acc.x;
            pDeriv[i].acceleration.y = acc.y;
            pDeriv[i].acceleration.z = acc.z;
            pDeriv[i].vitesse.x = pState[i].vitesse.x;
            pDeriv[i].vitesse.y = pState[i].vitesse.y;
            pDeriv[i].vitesse.z = pState[i].vitesse.z;
        }

    // Particle 0 is calculated last, because the statistics
    // data relate to this particle. They would be overwritten
    // otherwise 
    _root.StatReset();
    ParticuleData p(&pState[0], &_pAux[0]);
    PosParticule3D acc = _root.CalcForce(p);
    pDeriv[0].acceleration.x = acc.x;
    pDeriv[0].acceleration.y = acc.y;
    pDeriv[0].acceleration.z = acc.z;
    pDeriv[0].vitesse.x = pState[0].vitesse.x;
    pDeriv[0].vitesse.y = pState[0].vitesse.y;
    pDeriv[0].vitesse.z = pState[0].vitesse.z;  

    // Save Boites for camera orientations
    //  m_camDir.x = pState[0].x - pState[4000].x;
    //  m_camDir.y = pState[0].y - pState[4000].y;
    //m_camDIr.z = pState[0].z - pState[4000].y;
    _camPos.x = _root.GetCenterOfMass().x;
    _camPos.y = _root.GetCenterOfMass().y;
    _camPos.z = _root.GetCenterOfMass().z;
}

bool ModelNBody::IsFinished(double *state)
{
    return false;
}*/