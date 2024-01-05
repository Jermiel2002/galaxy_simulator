#ifndef _SDL_WINDOW_H
#define _SDL_WINDOW_H

#include <string>

#include <SDL/SDL.h>
#include <SDL/SDL_gfxPrimitives.h>
#include <SDL/SDL_opengl.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "../include/Particule3D.hpp"
/** \brief Basic infrastructure for grafical output using SDL/OpenGL */
class SDLWindow
{
public:
    SDLWindow(int width, int height, double axisLen, const std::string &caption);
    virtual ~SDLWindow();
    void MainLoop();
    void ExitMainLoop();
    void SetCaption(const std::string &caprion);
    int GetWidth() const;
    int GetHeight() const;
    virtual void Render() = 0;
    virtual void Update() = 0;
    /**
     * En programmation orientée objet, une méthode protégée (protected method) est une méthode d'une classe qui est accessible à la fois par la classe elle-même et par ses classes dérivées (sous-classes).
     * Cela signifie que la méthode ne peut pas être appelée depuis l'extérieur de la classe, mais elle peut être utilisée dans les classes qui héritent de cette classe.
     */
protected:
    virtual void PollEvents();
    virtual void OnProcessEvents(uint8_t type);

    //----------------------------------
    // Camera
    //----------------------------------
    const PosParticule3D &GetCamPos() const;
    const PosParticule3D &GetCamOrient() const;
    const PosParticule3D &GetCamLookAt() const;
    void SetCameraOrientation(const PosParticule3D &orientation);
    void SetCamera(const PosParticule3D &pos, const PosParticule3D &lookAt, const PosParticule3D &orient);
    void AdjustCamera();

    //-------------------------------------
    // Basic graphics functionality
    //-------------------------------------
    void DrawAxis(const PosParticule3D &origin);
    int GetFPS() const;
    void SaveToTGA(const std::string &sName);
    void SaveToTGA(int idx = -1);

    //----------------------------------------
    // misc
    //---------------------------------------------
    /**
     * Ces méthodes fournissent des fonctionnalités graphiques de base telles que dessiner un axe, obtenir le nombre d'images par seconde, et sauvegarder l'écran au format TGA.
     */
    void ScaleAxis(double scale);
    double GetFOV() const;
    SDL_Surface *Surface();
    SDL_Event _event;

    static void InitFont();
    static void KillFont();
    static void TextOut(const char *fmt, ...);
    static void TextOut(const PosParticule3D position, const char *fmt, ...);
    static PosParticule3D GetOGLPos(int x, int y);

    static GLuint s_fontBase;

private:
    void InitGL();
    double _fov; // Length of an axis
    int _width;  // width of the window in pixel
    int _height; // Height of the window in pixel
    int _fps;
    int _idxSnapshot;

    PosParticule3D _camPos;    // Position de la camera
    PosParticule3D _camLookAt; // Point at which the camera is aimed
    PosParticule3D _camOrient; // orientation of the camera

    SDL_Surface *_pScreen;

    volatile bool _bRunning;
};

#endif