#include <cstdlib>
#include <iostream>

#include "Octree.hpp"
#include "Boite.hpp"

int main()
{
    try
    {
        //Create an instance of OctreeNode
        OctreeNode rootNode(Boite(PosParticule3D(0,0,0), PosParticule3D(100,100,100)),nullptr);

        /*------------------------------------------------------------------------------*/
        //Créer particule 1
        PosParticule3D position1(1,2,3);

        EtatParticule etat1;
        etat1.pos = position1;
        etat1.vitesse = PosParticule3D(4,7,3);

        EtatAuxiliaire auxEtat1;
        auxEtat1.masse = 40;

        ParticuleData particule1(&etat1,&auxEtat1);

        //Créer particule 2
        PosParticule3D position2(3,4,5);
        //std::cout << "position créer\n";

        EtatParticule etat2;
        etat2.pos = position2;
        etat2.vitesse = PosParticule3D(4,6,8);

        EtatAuxiliaire auxEtat2;
        auxEtat2.masse = 70;
        //std::cout << "etats créer\n";

        ParticuleData particule2(&etat2,&auxEtat2);
        //std::cout << "Particule créée \n";

        //Créer particule 3
        PosParticule3D position3(5,8,3);
        //std::cout << "position créer\n";

        EtatParticule etat3;
        etat3.pos = position3;
        etat3.vitesse = PosParticule3D(23,6,12);

        EtatAuxiliaire auxEtat3;
        auxEtat2.masse = 80;
        //std::cout << "etats créer\n";

        ParticuleData particule3(&etat3,&auxEtat3);
        //std::cout << "Particule créée \n";
        
        //Insert some particles
        rootNode.InsertParticule(particule1,0);
        rootNode.InsertParticule(particule2,0);
        rootNode.InsertParticule(particule3,0);

        //Particules parfaitement insérée

        //Display information about the OctreeNode
        rootNode.DumpNode(0,0);

        return 0;
    } catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}