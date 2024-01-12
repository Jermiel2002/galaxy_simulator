#include "../include/NBodyWnd.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>
#include <omp.h>

#include "../include/IntegratorRK4.hpp"
#include "../include/IntegratorRK5.hpp"
#include "../include/IntegratorADB6.hpp"

NBodyWnd::NBodyWnd(int sz, std::string caption)
    :SDLWindow(sz, sz, 30.0, caption)
    ,_pModel(nullptr)
    ,_pSolver(nullptr)
    ,_camOrient(0)
    ,_flags(dspBODIES | /*dspTREE |*/ dspAXIS | dspSTAT | dspVERBOSE)
{}

NBodyWnd::~NBodyWnd()
{
    delete _pModel;
    delete _pSolver;
}

void NBodyWnd::Init(int num)
{
    // Create the n-body model class
    delete _pModel;
    _pModel = new ModelNBody();

    // assign model to the solver and set the integration step width
    delete _pSolver;

    //  _pSolver = new IntegratorADB4(_pModel, 5);
//    _pSolver = new IntegratorRK5(_pModel, _pModel->GetSuggestedTimeStep());
    _pSolver = new IntegratorADB6(_pModel, _pModel->GetSuggestedTimeStep());
    _pSolver->SetInitialState(_pModel->GetInitialState());

    if (_bDumpState)
    {
        std::string sName = std::string("traces_barnes_hut_") + _pSolver->GetID() + ".dat";
        _outfile.open(sName.c_str());
    }

    // OpenGL initialization
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    SetCameraOrientation(PosParticule3D(0, 1, 0));

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void NBodyWnd::Update()
{
    static int ct = 0;

    if (_flags & dspPAUSE)
        return;

    _pSolver->SingleStep();
    ++ct;

    if (_bDumpState && ct % 1000000 == 0)
    {
        int num = _pModel->GetDim() / 8;
        EtatParticule *state = reinterpret_cast<EtatParticule *>(_pSolver->GetState());
        _outfile << _pSolver->GetTime() << ", ";
        for (int i = 0; i < num; ++i)
        {
            _outfile << state[i].pos.x << ", "
                    << state[i].pos.y << ", "
                    << state[i].pos.z << ", ";
        }
        _outfile << std::endl;
    }

    if (_bDumpImage && ct % 2000 == 0)
    {
        SaveToTGA();
    }
}

void NBodyWnd::Render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    PosParticule3D orient;
    switch (_camOrient)
    {
    case 0:
        orient.x = 0;
        orient.y = 1;
        orient.z = 0;
        break;

    case 1:
        orient = _pModel->GetCamDir();
        break;
    }
    SetCameraOrientation(orient);

    if (_flags & dspAXIS)
    {
        // Draw axis at position of the first particle
        const PosParticule3D &cm = _pModel->GetCenterOfMass();
        DrawAxis(PosParticule2D(cm.x, cm.y));
    }

    if (_flags & dspTREE)
        DrawTree();

    if (_flags & dspBODIES)
        DrawBodies();

    if (_flags & dspROI)
        DrawROI();

    if (_flags & dspSTAT)
        DrawStat();

    if (_flags & dspHELP)
        DrawHelp();

    SDL_GL_SwapBuffers();
}

void NBodyWnd::DrawBodies()
{
    assert(_pSolver);

    EtatParticule *state = reinterpret_cast<EtatParticule *>(_pSolver->GetState());
    //  const PODAuxState *state_aux  = _pModel->GetAuxState();

    glColor3f(1, 1, 1);
    glPointSize(2); // state_aux[i].mass/10);
    glBegin(GL_POINTS);

    for (int i = 0; i < _pModel->GetN(); ++i)
    {
        glVertex3f(state[i].pos.x, state[i].pos.y, state[i].pos.z);
    }

    glEnd();
}

void NBodyWnd::DrawStat()
{
    double x0 = 10, y0 = 20, dy = 20;
    int line = 0;

    // acquire pointer to the root node of the barnes hut tree from the model
    OctreeNode *pRoot = _pModel->GetRootNode();

    const PosParticule3D &camPos = GetCamPos(),
                &camLookAt = GetCamLookAt();
    glColor3f(1, 1, 1);
    TextOut(x0, y0 + dy * line++, "Number of bodies (outside tree): %d (%d)", pRoot->GetNum(), pRoot->GetNumRenegades());
    TextOut(x0, y0 + dy * line++,0, "Theta: %2.1f", pRoot->GetTheta());
    TextOut(x0, y0 + dy * line++,0, "FPS: %d", GetFPS());
    TextOut(x0, y0 + dy * line++,0, "Time: %2.1f y", _pSolver->GetTime());
    TextOut(x0, y0 + dy * line++,0, "Camera: %2.2f, %2.2f, %2.2f", camPos.x, camPos.y, camPos.z);
    TextOut(x0, y0 + dy * line++,0, "LookAt: %2.2f, %2.2f, %2.2f", camLookAt.x, camLookAt.y, camLookAt.z);
    TextOut(x0, y0 + dy * line++,0, "Field of view: %2.2f pc", GetFOV());
    TextOut(x0, y0 + dy * line++,0, "Calculations: %d", pRoot->StatGetNumCalc());
    TextOut(x0, y0 + dy * line++,0, "Solver: %s", _pSolver->GetID().c_str());
}

void NBodyWnd::DrawHelp()
{
    double x0 = 10, y0 = 20, dy = 20;
    int line = 0;
    PosParticule3D p;

    glColor3f(1, 1, 1);
    TextOut(x0, y0 + dy * line++,0, "Keyboard commands");
    TextOut(x0, y0 + dy * line++,0, "a     - display axis");
    TextOut(x0, y0 + dy * line++,0, "t     - display Barnes Hut tree");
    TextOut(x0, y0 + dy * line++,0, "s     - display statistic data");
    TextOut(x0, y0 + dy * line++,0, "c     - display center of mass");
    TextOut(x0, y0 + dy * line++,0, "b     - display particles");
    TextOut(x0, y0 + dy * line++,0, "h     - display help text");
    TextOut(x0, y0 + dy * line++,0, "0..1  - Set camera orientation");
    TextOut(x0, y0 + dy * line++,0, "pause - halt simulation");
}

void NBodyWnd::DrawROI()
{
    const PosParticule3D &cm = _pModel->GetCenterOfMass();

    double l = GetFOV() / 20;
    glColor3f(1, 0, 0);

    glPushMatrix();
    glTranslatef(cm.x, cm.y, cm.z);
    // cross at the center of mass
    glBegin(GL_LINES);
    glVertex3f(-l, 0, 0);
    glVertex3f(l, 0, 0);
    glVertex3f(0, -l, 0);
    glVertex3f(0, l, 0);
    glEnd();

    // region of interest
    l = _pModel->GetROI() / 2.0;
    glBegin(GL_LINE_STRIP);
    glVertex3f(-l, l, 0);
    glVertex3f(l, l, 0);
    glVertex3f(l, -l, 0);
    glVertex3f(-l, -l, 0);
    glVertex3f(-l, l, 0);
    glEnd();

    glPopMatrix();
}

void NBodyWnd::DrawTree()
{
    struct DrawTree
    {
        enum EWhat
        {
            COMPLETE, ///< Display all nodes
            APPROX,   ///< Display only the nodes as they where used for the force calculation
        };

        DrawTree(OctreeNode *pRoot, EWhat what, int flags, int sz)
            : displayFlags(flags)
        {
            DrawNode(pRoot, 0, what, sz);
        }

        void DrawNode(OctreeNode *pNode, int level, EWhat what, int sz)
        {
            assert(pNode);

            /**
             * déterminer la couleur à utiliser pour dessiner le nœud en fonction de son niveau (level) et du type de dessin (what).
            */
            double col = 1 - level * 0.2;
            switch (what)
            {
            case COMPLETE:
                glColor3f(col, 1, col);
                break;
            case APPROX:
                glColor3f(0, 1, 0);
                break;
            }

            // Draw node cube
            if (what == COMPLETE ||
                (what == APPROX && !pNode->WasTooClose()))
            {
                const Boite cube = pNode->GetBoite();
                const PosParticule3D &min = cube.point1;
                const PosParticule3D &max = cube.point2;

                //glBegin(GL_LINE_STRIP);
                glBegin(GL_QUADS);

                //front face
                glVertex3f(min.x, min.y, min.z);
                glVertex3f(max.x, min.y, min.z);
                glVertex3f(max.x, max.y, min.z);
                glVertex3f(min.x, max.y, min.z);

                //Back face
                glVertex3f(min.x,min.y,max.z);
                glVertex3f(max.x,min.y,max.z);
                glVertex3f(max.x,max.y,max.z);
                glVertex3f(min.x,max.y,max.z);

                // Left face
                glVertex3f(min.x, min.y, min.z);
                glVertex3f(min.x, max.y, min.z);
                glVertex3f(min.x, max.y, max.z);
                glVertex3f(min.x, min.y, max.z);

                // Right face
                glVertex3f(max.x, min.y, min.z);
                glVertex3f(max.x, max.y, min.z);
                glVertex3f(max.x, max.y, max.z);
                glVertex3f(max.x, min.y, max.z);

                // Top face
                glVertex3f(min.x, max.y, min.z);
                glVertex3f(max.x, max.y, min.z);
                glVertex3f(max.x, max.y, max.z);
                glVertex3f(min.x, max.y, max.z);

                // Bottom face
                glVertex3f(min.x, min.y, min.z);
                glVertex3f(max.x, min.y, min.z);
                glVertex3f(max.x, min.y, max.z);
                glVertex3f(min.x, min.y, max.z);

                glEnd();
            
                // Draw a cross at the center of mass is the corresponding flag is set
                if (displayFlags & dspCENTER_OF_MASS && !pNode->IsExternal())
                {
                    double len = (double)sz / 50 * std::max(1 - level * 0.2, 0.1);
                    glPointSize(4);
                    glColor3f(col, 1, col);

                    const PosParticule3D cm = pNode->GetCenterOfMass();
                    glBegin(GL_LINES);
                    glVertex3f(cm.x - len, cm.y, cm.z);
                    glVertex3f(cm.x + len, cm.y, cm.z);
                    glEnd();
                    glBegin(GL_LINES);
                    glVertex3f(cm.x, cm.y - len, cm.z);
                    glVertex3f(cm.x, cm.y + len, cm.z);
                    glBegin(GL_LINES);
                    glVertex3f(cm.x, cm.y, cm.z - len);
                    glVertex3f(cm.x, cm.y, cm.z + len);
                    glEnd();

                }
            }

            // If the node was not subdivided in the force calculation
            // dont draw its subnodes either
            if (what != COMPLETE && !pNode->WasTooClose())
                return;

            for (int i = 0; i < 8; ++i)
            {
                if (pNode->_noeudFils[i])
                    DrawNode(pNode->_noeudFils[i], level + 1, what, sz);
            }
        } // DrawTree::DrawNode

        int displayFlags;
    }; // struct DrawTree

    OctreeNode *pRoot = _pModel->GetRootNode();
    if ((_flags & dspTREE) && (_flags & dspTREE_COMPLETE))
        DrawTree DrawComplete(pRoot, DrawTree::COMPLETE, _flags, GetFOV());
    else if ((_flags & dspTREE) && !(_flags & dspTREE_COMPLETE))
        DrawTree DrawFar(pRoot, DrawTree::APPROX, _flags, GetFOV());
}

void NBodyWnd::DrawNode(OctreeNode *pNode, int level)
{
    assert(pNode);
    double len = 0.01 * std::max(1 - level * 0.2, 0.1);

    double col = 1 - level * 0.2;

    if (pNode->WasTooClose())
        glColor3f(1, col, col);
    else
        glColor3f(col, 1, col);

    const Boite cube = pNode->GetBoite();
    const PosParticule3D &min = cube.point1;
    const PosParticule3D &max = cube.point2;

    glBegin(GL_QUADS);

    // Front face
    glVertex3f(min.x, min.y, min.z);
    glVertex3f(max.x, min.y, min.z);
    glVertex3f(max.x, max.y, min.z);
    glVertex3f(min.x, max.y, min.z);

    // Back face
    glVertex3f(min.x, min.y, max.z);
    glVertex3f(max.x, min.y, max.z);
    glVertex3f(max.x, max.y, max.z);
    glVertex3f(min.x, max.y, max.z);

    // Left face
    glVertex3f(min.x, min.y, min.z);
    glVertex3f(min.x, max.y, min.z);
    glVertex3f(min.x, max.y, max.z);
    glVertex3f(min.x, min.y, max.z);

    // Right face
    glVertex3f(max.x, min.y, min.z);
    glVertex3f(max.x, max.y, min.z);
    glVertex3f(max.x, max.y, max.z);
    glVertex3f(max.x, min.y, max.z);

    // Top face
    glVertex3f(min.x, max.y, min.z);
    glVertex3f(max.x, max.y, min.z);
    glVertex3f(max.x, max.y, max.z);
    glVertex3f(min.x, max.y, max.z);

    // Bottom face
    glVertex3f(min.x, min.y, min.z);
    glVertex3f(max.x, min.y, min.z);
    glVertex3f(max.x, min.y, max.z);
    glVertex3f(min.x, min.y, max.z);

    glEnd();

    if (_flags & dspCENTER_OF_MASS && !pNode->IsExternal())
    {
        PosParticule3D cm = pNode->GetCenterOfMass();

        glPointSize(4);
        glColor3f(col, 1, col);
        glBegin(GL_LINES);
        glVertex3f(cm.x - len, cm.y, cm.z);
        glVertex3f(cm.x + len, cm.y, cm.z);
        glVertex3f(cm.x, cm.y - len, cm.z);
        glVertex3f(cm.x, cm.y + len, cm.z);
        glVertex3f(cm.x, cm.y, cm.z - len);
        glVertex3f(cm.x, cm.y, cm.z + len);
        glEnd();
    }

    if (!pNode->WasTooClose())
        return;

    for (int i = 0; i < 8; ++i)
    {
        if (pNode->_noeudFils[i])
        {
            DrawNode(pNode->_noeudFils[i], level + 1);
        }
    }
}

void NBodyWnd::OnProcessEvents(uint8_t type)
{
    switch (type)
    {
        /* commented: does not work
          case SDL_MOUSEBUTTONDOWN:
          {
            // Place a tracer particle at the mouse coordinates
            if (!_pSolver)
              break;

            PODState *state = reinterpret_cast<PODState *>(_pSolver->GetState());
            Vec3D p = GetOGLPos(_event.button.x,
                                _event.button.y);
            state[0].x = p.x;
            state[0].y = p.y;
            SetCamera(p, p, Vec3D(0, 1, 0));

            // the solver may need to rest its temporary arrays. I can't just
            // overwrite part of its data because ADB schemes will go mad
            // if i change a particles position. without restarting the engine
            _pSolver->SetInitialState(reinterpret_cast<double *>(state));
            _pSolver->SingleStep();
          }
          break;
        */
    case SDL_KEYDOWN:
        switch (_event.key.keysym.sym)
        {
        case SDLK_1:
            _camOrient = 0;
            break;

        case SDLK_2:
            _camOrient = 1;
            break;

        case SDLK_a:
            std::cout << "Display:  Toggling axis " << ((_flags & dspAXIS) ? "off" : "on") << "\n";
            _flags ^= dspAXIS;
            break;

        case SDLK_b:
            std::cout << "Display:  Toggling bodies " << ((_flags & dspBODIES) ? "off" : "on") << "\n";
            _flags ^= dspBODIES;
            break;

        case SDLK_t:
        {
            if (!(_flags & dspTREE))
            {
                _flags ^= dspTREE;
                std::cout << "Display:  Tree cells used in force calculation\n";
            }
            else if ((_flags & dspTREE) && !(_flags & dspTREE_COMPLETE))
            {
                _flags ^= dspTREE_COMPLETE;
                std::cout << "Display:  Complete tree\n";
            }
            else if ((_flags & dspTREE) && (_flags & dspTREE_COMPLETE))
            {
                _flags &= ~(dspTREE | dspTREE_COMPLETE);
                std::cout << "Display:  No tree\n";
            }
        }
        break;

        case SDLK_c:
            std::cout << "Display:  Center of mass " << ((_flags & dspCENTER_OF_MASS) ? "off" : "on") << "\n";
            _flags ^= dspCENTER_OF_MASS;
            break;

        case SDLK_h:
            std::cout << "Display:  Help text" << ((_flags & dspHELP) ? "off" : "on") << "\n";
            _flags ^= dspHELP;
            _flags &= ~dspSTAT;
            break;

        case SDLK_PAUSE:
            std::cout << "Simulation:  pause " << ((_flags & dspPAUSE) ? "off" : "on") << "\n";
            _flags ^= dspPAUSE;
            break;

        case SDLK_v:
            std::cout << "Simulation:  verbose mode " << ((_flags & dspVERBOSE) ? "off" : "on") << "\n";
            _flags ^= dspVERBOSE;
            if (_pModel)
                _pModel->SetVerbose(_flags & dspVERBOSE);
            break;

        case SDLK_s:
            std::cout << "Display:  statistics " << ((_flags & dspSTAT) ? "off" : "on") << "\n";
            _flags ^= dspSTAT;
            break;

        case SDLK_f:
            std::cout << "Display:  force indicator " << ((_flags & dspARROWS) ? "off" : "on") << "\n";
            _flags ^= dspARROWS;
            break;

        case SDLK_r:
            std::cout << "Display:  region of interest " << ((_flags & dspROI) ? "off" : "on") << "\n";
            _flags ^= dspROI;
            break;

        case SDLK_p:
            SaveToTGA();
            break;

        case SDLK_y:
            _pModel->SetTheta(_pModel->GetTheta() + 0.1);
            break;

        case SDLK_x:
            _pModel->SetTheta(std::max(_pModel->GetTheta() - 0.1, 0.1));
            break;

        case SDLK_KP_PLUS:
            ScaleAxis(0.9);
            SetCameraOrientation(PosParticule3D(0, 1, 0));
            break;

        case SDLK_KP_MINUS:
            ScaleAxis(1.1);
            SetCameraOrientation(PosParticule3D(0, 1, 0));
            break;

        default:
            break;
        }

        break;
    }
}