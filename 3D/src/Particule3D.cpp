#include "../include/Particule3D.hpp"
#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>

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
EtatParticule::EtatParticule(): pos(nullptr), vitesse(nullptr), accelearation(nullptr){};
EtatParticule::EtatParticule(PosParticule3D *ppos,PosParticule3D *vvitesse, PosParticule3D *aaccelearation): pos(ppos),vitesse(vvitesse),accelearation(aaccelearation)
{
    assert(pos);
    assert(vitesse);
    assert(accelearation);
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