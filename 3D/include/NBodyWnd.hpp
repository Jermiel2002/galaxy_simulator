#ifndef _NBODYWND_H
#define _NBODYWND_H

#include <stdint.h>
#include <fstream>

#include "SDLWnd.hpp"
#include "Octree.hpp"
#include "ModelNBody.hpp"
#include "IIntegrator.hpp"

/** \brief Main window of th n-body simulation. */
class NBodyWnd final : public SDLWindow
{
public:

    NBodyWnd(int sz, std::string caption);
    virtual ~NBodyWnd();

    virtual void Render() override;
    virtual void Update() override;
    virtual void OnProcessEvents(uint8_t type) override;

    void Init(int num);

private:

    /**
     * L'énumération DisplayState est utilisée dans la classe NBodyWnd pour représenter différents états d'affichage pour la simulation n-body. 
     * Chaque option peut être activée ou désactivée individuellement en les combinant à l'aide d'opérations binaires (par exemple, en utilisant l'opérateur | pour activer plusieurs options simultanément).
    */
    enum DisplayState : unsigned int
    {
        dspNONE = 0,//Aucun état d'affichage en particulier
        dspAXIS = 1 << 0, //Activer cette option pour afficher les axes
        dspBODIES = 1 << 1, //Activer cette option pour afficher les particules du système.
        dspSTAT = 1 << 2, //Activer cette option pour afficher des statistiques sur la simulation.
        dspTREE = 1 << 3, //Activer cette option pour afficher l'arbre de Barnes-Hut.
        dspTREE_COMPLETE = 1 << 4, //Activer cette option pour afficher l'arbre de Barnes-Hut complet.
        dspCENTER_OF_MASS = 1 << 5,//Activer cette option pour afficher le centre de masse du système.
        dspPAUSE = 1 << 6,//Activer cette option pour mettre la simulation en pause.
        dspVERBOSE = 1 << 7,//Activer cette option pour afficher les détails de manière verbeuse.
        dspHELP = 1 << 8,//Activer cette option pour afficher une aide.
        dspARROWS = 1 << 9,//Activer cette option pour afficher des flèches.
        dspROI = 1 << 10//Activer cette option pour afficher la région d'intérêt (ROI).
    };

    NBodyWnd(const NBodyWnd &orig);

    void DrawBodies();
    void DrawStat();
    void DrawTree();
    void DrawHelp();
    void DrawROI();
    void DrawCenterOfMass();
    void DrawNode(OctreeNode *pNode, int level);

    ModelNBody *_pModel;
    IIntegrator *_pSolver;

    int _camOrient;  ///< Index of the camera orientation to use
    uint32_t _flags; ///< The display flags
    bool _bDumpState = false;
    bool _bDumpImage = false;
    std::ofstream _outfile;
};

#endif
