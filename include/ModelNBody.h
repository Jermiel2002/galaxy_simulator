#ifndef _MODEL_N_BODY_H
#define _MODEL_N_BODY_H

#include <cmath>

#include "Constants.h"
#include "Octree.h"
#include "Boite.h"
#include "IModel.h"

/** \brief
 * Le problème des N-corps est un problème classique en physique où l'on cherche à simuler le mouvement de N particules sous l'influence mutuelle de leurs forces gravitationnelles.
 * La classe ModelNBody est conçue pour représenter un modèle traitant le problème des N-corps (n-body problem)
 * ModelNBody est une classe héritant de la classe IModel
 */
class ModelNBody final : public IModel
{
public:
    ModelNBody();
    virtual ~ModelNBody();

    void Init();
    void InitCollision();
    void Init3Body();

    virtual void Eval(double *state, double time, double *deriv) override;
    virtual bool IsFinished(double *state) override;
    virtual double *GetInitialState() override;

    double GetSuggestedTimeStep() const;
    double GetTheta() const;
    OctreeNode *GetRootNode();
    const EtatAuxiliaire *GetAuxState() const;
    int GetN() const;

    const PosParticule3D &GetCamDir() const;
    const PosParticule3D &GetCamPos() const;
    PosParticule3D GetCenterOfMass() const;

    void SetTheta(double theta);
    void SetVerbose(bool bVerbose);
    void SetROI(double roi);
    double GetROI() const;

private:
    void CalcBHArea(const ParticuleData &p);
    void BuiltTree(const ParticuleData &p);
    void GetOrbitalVelocity(const ParticuleData &p1, const ParticuleData &p2);
    void ResetDim(int num, double stepsize);

    EtatParticule *_pInitial; ///< The initial state
    EtatAuxiliaire *_pAux;    ///< Auxilliary state information

    OctreeNode _root; ///< The root node of the barnes hut tree
    Boite *_reperBoite;
    PosParticule3D _pointMin; /// point inférieur gauche coté avant du cube, of the bounding box containing all particles
    PosParticule3D _pointMax; /// point supérieur droit coté arrière du cube, of the bounding box containing all particles
    PosParticule3D _center;   ///< The center of the simulation, the barnes hut tree is centered at this point
    PosParticule3D _camDir;   ///< Direction of the camera
    PosParticule3D _camPos;   ///< Position of the camera
    double _roi;
    double _timeStep;

    static constexpr double gamma_1 = Constants::Gamma / (Constants::ParsecInMeter * Constants::ParsecInMeter * Constants::ParsecInMeter) * Constants::MassOfSun * (365.25 * 86400) * (365.25 * 86400);

    int _num;
    bool _bVerbose;
};
#endif//_MODEL_N_BODY_H