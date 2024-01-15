#ifndef _PARTICULE3D_HPP
#define _PARTICULE3D_HPP

#pragma pack(push, 1)

class PosParticule2D
{
public:
    PosParticule2D(double _x = 0, double _y = 0);
    double x;
    double y;
};

//Particule3D représente les coordonnées d'une particule ou un corps dans l'espace 3D, peut également etre utilisé pour d'autre objet necessitant une position 3D
class PosParticule3D{
    public:
        //constructeurs
        PosParticule3D();
        PosParticule3D(double x, double y, double z);

     /*   //fonctions
        bool ParticuleEgaux(PosParticule3D p1, PosParticule3D p2); //verifie que deux particules ont la même position

        double distEntrParticule(PosParticule3D p1, PosParticule3D p2); //calcul de la distance euclidienne entre deux point1s tridimensionnels dans l'espace

        //PosParticule3D fusionParticule(PosParticule3D p1, PosParticule3D p2);

        //PosParticule3D MultiplierParticule(PosParticule3D p, double s);*/

        //variables
        double x;
        double y;
        double z;

};

/*
* EtatParticule : représentation de l'état d'une particule (position et vitesses)
* Pour une particule à un instant t donnée, l'état peut être décrit par ces valeurs. L'état d'une particule peut changer au cours du temps
* en raison des forces appliquées, des collisions...
*/
class EtatParticule
{
    public:
    EtatParticule();
    EtatParticule(PosParticule3D pos, PosParticule3D vitesse);
    PosParticule3D pos;
    PosParticule3D vitesse;
};

class EtatAuxiliaire
{
    public:
    EtatAuxiliaire();
    EtatAuxiliaire(double masse);
    double masse;
};

//class de mis à jour de l'état d'une particule au fil du temps contenant les dérivées de l'état d'une particule
struct MajEtat
{
    PosParticule3D vitesse;
    PosParticule3D acceleration;
};
#pragma pack(pop)

/* ParticleData définit une classe qui encapsule des données liées à une particule*/
class ParticuleData
{
    public:
        //constructeurs
        ParticuleData(); //constructeur par defaut à initialisé dans le .cpp
        ParticuleData(EtatParticule *pState, EtatAuxiliaire *pAuxState);//Initialise une particule connaissant son etat et son etatAux
        ParticuleData(const ParticuleData &autre);//constructeur de copie simple qui initialise une nouvelle particule en copiant les valeurs des attributs d'une autre paticule existante
        ParticuleData &operator = (const ParticuleData &ref); //Opérateur d'affectation utiliser qd on souhaite copier les valeurs d'une particules vers une autre existante

        //fonctions
        void Reset(); //réinitialise l'état de l'objet ParticleData
        bool IsNull() const; //Vérifie si les variables de la particule sont nulles
        //variables
        EtatParticule *_pState; //Une particule est caractérisé par son état
        EtatAuxiliaire *_pAuxState; //Et son état auxiliaire
};

#endif