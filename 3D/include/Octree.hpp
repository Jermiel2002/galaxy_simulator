#ifndef OCTREE_HPP
#define OCTREE_HPP

//--- Standard includes --------------------------------------------------------
#include <vector>

//--- Implementation -----------------------------------------------------------
#include "Particule3D.hpp"
#include "Espace3D.hpp"

/**
 * Les calculs effectués dans d'autres parties de la classe, tels que le calcul des forces gravitationnelles, 
 * devraient être adaptés pour travailler avec des coordonnées en trois dimensions. Les méthodes CalcAcc et CalcTreeForce ainsi que d'autres parties du code nécessiteraient des ajustements pour prendre en compte les changements dimensionnels.
*/
//------------------------------------------------------------------------------
class Octree
{
    public:
        struct DebugStat
        {
            int _nNumCalc;///< Total number of calculations for estimating the force
        };

        //constructeurs
        Octree(const PosParticule3D &min, const PosParticule3D &max, Octree *parent = nullptr);
        ~Octree();//on libère la mémoire du noeud

        //fonctions
        void Reset(const PosParticule3D &min, const PosParticule3D &max);
        bool IsRoot() const;
        bool IsExternal() const;
        bool WasTooClose() const;

        void StatReset();
        int StatGetNumCalc() const;

        int GetNumRenegades() const;
        int GetNum() const;
        const PosParticule3D &GetCenterOfMass() const;
        const PosParticule3D &GetMin() const;
        const PosParticule3D &GetMax() const;

        double GetTheta() const;
        void SetTheta(double Theta);

        void Insert(const ParticuleData &newParticule, int level);

        /*EQuadrant GetQuadrant(double x, double y) const;
        BHTreeNode *CreateQuadNode(EQuadrant eQuad);*/

        void ComputeMassDistribution();

        PosParticule3D CalcForce(const ParticuleData &p) const;
        void DumpNode(int quad, int level);

    public:
        Octree *_octreeNode[8];
    
    private:
        PosParticule3D CalcAcc(const ParticuleData &p1, const ParticuleData &p2) const;
        PosParticule3D CalcTreeForce(const ParticuleData &p) const;

            /** \brief Data for the particle.

            Only valid if this is a leaf node.
            */
        ParticuleData _particle;
        double _mass;              ///< Mass of all particles inside the node
        PosParticule3D _cm;                 ///< Center of Mass
        PosParticule3D _min;                ///< Upper left edge of the node
        PosParticule3D _max;                ///< Lower right edge of the node
        PosParticule3D _center;             ///< Center of the node
        Octree *_parent;       ///< The parent node
        int _num;                  ///< The number of particles in this node
        mutable bool _bSubdivided; ///< True if this node is too close to use the approximation for the force calculation

        static double s_theta;
        static std::vector<ParticuleData> s_renegades;


    public:
        static double s_gamma;

    
    private:
        static double s_soft;
        static DebugStat s_stat;
};

#endif