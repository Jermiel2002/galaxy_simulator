#include "../include/Particule3D.hpp"
#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>

using namespace std;

//codage des constructeurs pour identifier les positions d'une particule dans l'espace
PosParticule3D::PosParticule3D(): x(nullptr), y(nullptr), z(nullptr){};
PosParticule3D::PosParticule3D(double *xCoord, double *yCoord, double *zCoord): x(xCoord), y(yCoord), z(zCoord)
{
    assert(x);
    assert(y);
    assert(z);
}

//codage des constructeurs pour représenter l'état d'une particule dans l'espace
EtatParticule::EtatParticule(): pos(nullptr), vitesse(nullptr), acceleration(nullptr){};
EtatParticule::EtatParticule(PosParticule3D *ppos,PosParticule3D *vvitesse, PosParticule3D *aacceleration): pos(ppos),vitesse(vvitesse),acceleration(aacceleration)
{
    assert(pos);
    assert(vitesse);
    assert(acceleration);
};

//codage des constructeurs pour représenter l'état auxiliaire d'une particule dans l'espace
EtatAuxiliaire::EtatAuxiliaire():masse(nullptr){};
EtatAuxiliaire::EtatAuxiliaire(double *weight): masse(weight)
{
    assert(masse);
};

//codage des constructeurs pour représenter les données liées à une particule
ParticuleData::ParticuleData(): _pState(nullptr), _pAuxState(nullptr){};
ParticuleData::ParticuleData(EtatParticule *pS, EtatAuxiliaire *pA): _pState(pS), _pAuxState(pA)
{
    assert(_pState);
    assert(_pAuxState);
};
ParticuleData::ParticuleData(ParticuleData const& autre): _pState(autre._pState), _pAuxState(autre._pAuxState)
{};
ParticuleData &ParticuleData::operator=(const ParticuleData &ref)
{
    if (this != &ref)
    {
        _pState = ref._pState;
        _pAuxState = ref._pAuxState;
    }

    return *this;
};

void ParticuleData::Reset()
{
    _pState = nullptr;
    _pAuxState = nullptr;
}


bool ParticuleData::IsNull() const
{
    return _pState && _pAuxState;
}

bool PosParticule3D::ParticuleEgaux(PosParticule3D p1, PosParticule3D p2) 
    {
        if(p1.x == p2.x && p1.y == p2.y && p1.z == p2.z)
            return true;
        else
            return false;
    }

double PosParticule3D::distEntrParticule(PosParticule3D p1, PosParticule3D p2)
    {
        return sqrt((p1.x - p2.x)*(p1.x-p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
    }

/*PosParticule3D PosParticule3D::fusionParticule(PosParticule3D p1, PosParticule3D p2)
    {
        return PosParticule3D((p1.x+p2.x) , p1.y+p2.y , p1.z+p2.z);
    }*/

/*PosParticule3D PosParticule3D::MultiplierParticule(PosParticule3D p, double s)
    {
        return PosParticule3D(p.x * s, p.y * s, p.z * s);
    }*/