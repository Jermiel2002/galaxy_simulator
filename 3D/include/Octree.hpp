#pragma once
//--- Implementation -----------------------------------------------------------
#include "../include/Boite.hpp"
#include <vector>
#include "../include/Constants.hpp"

//------------------------------------------------------------------------------
/**A quoi ressemble un noeud de notre OctreeNode ? c'est ce à quoi répond la class OctreeNode*/
class OctreeNode
{
    public:

        /** \brief Enumeration for the Octrants. */
        enum BoiteAParticule
        {
            SWD = 0,//Octant Sud-Ouest Bas
            SWU,//Octant Sud-Ouest Haut
            NWD,//Octant Nord-Ouest Bas
            NWU,//Octant Nord-Ouest Haut
            SED,//Octant Sud-Est Bas
            SEU,//Octant Sud-Est Haut
            NED,//Octant Nord-Est Bas
            NEU,//Octant Nord-est haut
            NONE
        };

        //constructeurs
        OctreeNode(const Boite &reper_boite, OctreeNode *parent = nullptr);//on crée un noeud
        ~OctreeNode();//on libère la mémoire du noeud en le supprimant

        //fonctions
        void Reset(const Boite &reperBoite);
        bool IsRoot() const;
        bool IsExternal() const;
        bool WasTooClose() const;

        void StatReset();
        int StatGetNumCalc() const;

        int GetNumRenegades() const;
        int GetNum() const;
        const PosParticule3D &GetCenterOfMass() const;
        Boite GetBoite() const;
        const ParticuleData &GetParticule() const;
        void SetBoite(Boite boite);

        double GetTheta() const;
        void SetTheta(double Theta);

        BoiteAParticule GetTypeBoite(PosParticule3D const& p) const;
        OctreeNode *CreateOctreeNode(BoiteAParticule boiteP);
        void InsertParticule(const ParticuleData &newParticule, int level);

        void ComputeMassDistribution();
        void GetOrbitalVelocity(const ParticuleData &p1, const ParticuleData &p2);


        PosParticule3D CalcForce(const ParticuleData &p) const;
        void DumpNode(int oct, int level);

        struct DebugStat
        {
            int _nNumCalc;///< Total number of calculations for estimating the force
        };
    
    public:
        OctreeNode *_noeudFils[8];//un noeud est constitué de huit fils

    private://l'utilisateur ne doit pas modifié ces variables; pour les consulté ou modifié si nécessaire, on créera une fonction plus tard
        PosParticule3D CalcAcc(const ParticuleData &p1, const ParticuleData &p2) const;
        PosParticule3D CalcTreeForce(const ParticuleData &p) const;
        //ParticuleData _particle;
        double _mass;              ///< Mass of all particles inside the node
        PosParticule3D _cm;                 ///< Center of Mass
        Boite _boite;             //La boite contient les coordonnées de l'espace cubique qu'occupe une particule, la boite englobante  
        PosParticule3D _center;             ///< Center of the node
        OctreeNode *_parent;       ///< The parent node
        int _num;                  ///< The number of particles in this node
        mutable bool _bSubdivided; ///< True if this node is too close to use the approximation for the force calculation; "mutable" = peut être modifié même si l'objet est déclaré comme const

        static double s_theta;
        static std::vector<ParticuleData> s_renegades;


    //public:
        //static constexpr double gamma_1 = Constants::Gamma / (Constants::ParsecInMeter * Constants::ParsecInMeter * Constants::ParsecInMeter) * Constants::MassOfSun * (365.25 * 86400) * (365.25 * 86400);

  
    private:
        static double s_soft;
        static DebugStat s_stat;
};
