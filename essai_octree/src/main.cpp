#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <omp.h>

#include "../include/Octree.hpp"
#include "../include/Boite.hpp"
//#include "../include/ModelNBody.hpp"

int main()
{
    try
    {
        //On crée la boite englobante
        Boite boite_englobante = Boite(PosParticule3D(0,0,0),PosParticule3D(4,4,4));
        OctreeNode oct = OctreeNode(boite_englobante);//création de l'octree symbolisant un noeud (un cube)

        std::cout << "Coordonnées de la boite englobante\n" << 
        "point inf avant gauche : (" << oct.GetBoite().point1.x << ", " << oct.GetBoite().point1.y << ", " << oct.GetBoite().point1.z << ")" 
        " \npoint sup arrière droit : (" << oct.GetBoite().point2.x << ", " << oct.GetBoite().point2.y << ", " << oct.GetBoite().point2.z << ")" << '\n';

        /*TEST DES POSITIONS DANS LE CUBE. On suppose un cube dont le coin inférieur gauche de la face avant : (0,0,0)
        et le coin supérieur droit de la face arrière : (4,4,4)*/
        PosParticule3D test1(1,1,1);//sous cube inférieur gauche
        PosParticule3D test2(3,1,1);//sous cube inférieur droit
        PosParticule3D test3(1,3,1);//sous cube supérieur gauche
        PosParticule3D test4(3,3,1);//sous cube supérieur droit
        PosParticule3D test5(1,1,3);//sous cube inférieur arrière gauche
        PosParticule3D test6(3,1,3);//sous cube inférieur arrière droit
        PosParticule3D test7(1,3,3);//sosu cube supérieur arrière gauche
        PosParticule3D test8(3,3,3);//sous cube supérieur arrière droit
        /*Test que chaque position appartient bien à un sous cube
        std::cout << oct.GetTypeBoite(test1) << "\n";
        std::cout << oct.GetTypeBoite(test2) << "\n";
        std::cout << oct.GetTypeBoite(test3) << "\n";
        std::cout << oct.GetTypeBoite(test4) << "\n";
        std::cout << oct.GetTypeBoite(test5) << "\n";
        std::cout << oct.GetTypeBoite(test6) << "\n";
        std::cout << oct.GetTypeBoite(test7) << "\n";
        std::cout << oct.GetTypeBoite(test8) << "\n";*/

        /*CREATION DE PARTICULES AYANT CES POSITIONS*/
        EtatParticule etat1;
        etat1.pos = test1;
        etat1.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat1;
        auxEtat1.masse = 40000;//le premier particule est au milieu elle a une grosse masse
        ParticuleData particule1(&etat1,&auxEtat1);

        EtatParticule etat2;
        etat2.pos = test2;
        etat2.vitesse = PosParticule3D(23,22,2);
        EtatAuxiliaire auxEtat2;
        auxEtat2.masse = 500;
        ParticuleData particule2(&etat2,&auxEtat2);

        EtatParticule etat3;
        etat3.pos = test3;
        etat3.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat3;
        auxEtat3.masse = 2632;
        ParticuleData particule3(&etat3,&auxEtat3);

        EtatParticule etat4;
        etat4.pos = test4;
        etat4.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat4;
        auxEtat4.masse = 4220;
        ParticuleData particule4(&etat4,&auxEtat4);

        EtatParticule etat5;
        etat5.pos = test5;
        etat5.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat5;
        auxEtat5.masse = 260;
        ParticuleData particule5(&etat5,&auxEtat5);

        EtatParticule etat6;
        etat6.pos = test6;
        etat6.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat6;
        auxEtat6.masse = 262;
        ParticuleData particule6(&etat6,&auxEtat6);

        EtatParticule etat7;
        etat7.pos = test7;
        etat7.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat7;
        auxEtat7.masse = 827;
        ParticuleData particule7(&etat7,&auxEtat7);

        EtatParticule etat8;
        etat8.pos = test8;
        etat8.vitesse = PosParticule3D(4,7,3);
        EtatAuxiliaire auxEtat8;
        auxEtat8.masse = 272;
        ParticuleData particule8(&etat8,&auxEtat8);

        /*INSERTION DES PARTICULES CREER*/
        oct.InsertParticule(particule1,0);
        oct.InsertParticule(particule2,0);
        oct.InsertParticule(particule3,0);
        oct.InsertParticule(particule4,0);
        oct.InsertParticule(particule5,0);
        oct.InsertParticule(particule6,0);
        oct.InsertParticule(particule7,0);
        oct.InsertParticule(particule8,0);

        /*AFFICHAGE DE L'ARBRE*/
        oct.DumpNode(0,0);
        oct.ComputeMassDistribution();
        std::cout << "le centre de masse a pour coordonnées : (" << oct.GetCenterOfMass().x << ", " << oct.GetCenterOfMass().y << ", " << oct.GetCenterOfMass().z << ")\n";
        //std::cout << "la masse totale est :" << oct.GetMass();
        /*
        Pour diviser le cube en 8 sous-cubes équivalents, nous pouvons utiliser les coordonnées médianes sur chaque axe. Les coordonnées du centre du cube sont simplement la moyenne des coordonnées des coins opposés. Voici les coordonnées des coins des sous-cubes et 8 points à l'intérieur de chaque sous-cube :

    Sous-cube inférieur gauche :
        Coins : (0,0,0) à (2,2,2)
        Centre : (1,1,1)
        Points à l'intérieur :
            (0.5, 0.5, 0.5)
            (1, 0.5, 0.5)
            (0.5, 1, 0.5)
            (1, 1, 0.5)
            (0.5, 0.5, 1)
            (1, 0.5, 1)
            (0.5, 1, 1)
            (1, 1, 1)

    Sous-cube inférieur droit :
        Coins : (2,0,0) à (4,2,2)
        Centre : (3,1,1)
        Points à l'intérieur :
            (2.5, 0.5, 0.5)
            (3, 0.5, 0.5)
            (2.5, 1, 0.5)
            (3, 1, 0.5)
            (2.5, 0.5, 1)
            (3, 0.5, 1)
            (2.5, 1, 1)
            (3, 1, 1)

    Sous-cube supérieur gauche :
        Coins : (0,2,0) à (2,4,2)
        Centre : (1,3,1)
        Points à l'intérieur :
            (0.5, 2.5, 0.5)
            (1, 2.5, 0.5)
            (0.5, 3, 0.5)
            (1, 3, 0.5)
            (0.5, 2.5, 1)
            (1, 2.5, 1)
            (0.5, 3, 1)
            (1, 3, 1)

    Sous-cube supérieur droit :
        Coins : (2,2,0) à (4,4,2)
        Centre : (3,3,1)
        Points à l'intérieur :
            (2.5, 2.5, 0.5)
            (3, 2.5, 0.5)
            (2.5, 3, 0.5)
            (3, 3, 0.5)
            (2.5, 2.5, 1)
            (3, 2.5, 1)
            (2.5, 3, 1)
            (3, 3, 1)

    Sous-cube inférieur arrière gauche :
        Coins : (0,0,2) à (2,2,4)
        Centre : (1,1,3)
        Points à l'intérieur :
            (0.5, 0.5, 2.5)
            (1, 0.5, 2.5)
            (0.5, 1, 2.5)
            (1, 1, 2.5)
            (0.5, 0.5, 3)
            (1, 0.5, 3)
            (0.5, 1, 3)
            (1, 1, 3)

    Sous-cube inférieur arrière droit :
        Coins : (2,0,2) à (4,2,4)
        Centre : (3,1,3)
        Points à l'intérieur :
            (2.5, 0.5, 2.5)
            (3, 0.5, 2.5)
            (2.5, 1, 2.5)
            (3, 1, 2.5)
            (2.5, 0.5, 3)
            (3, 0.5, 3)
            (2.5, 1, 3)
            (3, 1, 3)

    Sous-cube supérieur arrière gauche :
        Coins : (0,2,2) à (2,4,4)
        Centre : (1,3,3)
        Points à l'intérieur :
            (0.5, 2.5, 2.5)
            (1, 2.5, 2.5)
            (0.5, 3, 2.5)
            (1, 3, 2.5)
            (0.5, 2.5, 3)
            (1, 2.5, 3)
            (0.5, 3, 3)
            (1, 3, 3)

    Sous-cube supérieur arrière droit :
        Coins : (2,2,2) à (4,4,4)
        Centre : (3,3,3)
        Points à l'intérieur :
            (2.5, 2.5, 2.5)
            (3, 2.5, 2.5)
            (2.5, 3, 2.5)
            (3, 3, 2.5)
            (2.5, 2.5, 3)
            (3, 2.5, 3)
            (2.5, 3, 3)
            (3, 3, 3)
        */

        return 0;
    } catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}