#ifndef _IMODEL_H
#define _IMODEL_H

#include <string>

/**
 * Cette interface définit un ensemble de méthodes que tout modèle devant être simulé doit implémenter. L'interface fournit une structure de base pour représenter un modèle,
 * permettant ainsi une certaine abstraction et une utilisation générique dans un simulateur ou un cadre de simulation.
 *
 * Virtual signifie que la classe ou la methode doit être redéfinie ou implémenter dans
 * les classes dérivées
 */
class IModel
{
public:
    IModel(const std::string &name, unsigned dim = 1) : _dim(dim), _name(name) {}

    virtual ~IModel() {}

    unsigned GetDim() const noexcept { return _dim; }
    void SetDim(unsigned dim) noexcept { _dim = dim; }

    std::string GetName() const noexcept { return _name; }

    /*pour évaluer l'état du modèle à un certain moment time en fonction de l'état actuel state et peut mettre à jour les dérivées deriv_in.*/
    virtual void Eval(double *state, double time, double *deriv_in) = 0;
    /*Cette méthode détermine si le modèle a atteint un état final ou s'il doit continuer à être simulé. Elle renvoie true si le modèle est considéré comme terminé, sinon false.*/
    virtual bool IsFinished(double *state) = 0;
    /*Cette méthode renvoie un pointeur vers l'état initial du modèle. Elle est utilisée pour initialiser la simulation.*/
    virtual double *GetInitialState() = 0;

private:
    unsigned _dim;
    std::string _name;
};

#endif //_IMODEL_H