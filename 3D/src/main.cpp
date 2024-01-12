#include <cstdlib>
#include <iostream>

#include "../include/NBodyWnd.hpp"


int main(int argc, char **argv)
{
    try
    {
        NBodyWnd wndMain(700, "NBody Simulation (Barnes Hut algorithm)");

        // Define simulation size
        wndMain.Init(/*50*/5);
        wndMain.MainLoop();
    }
    catch (std::exception &exc)
    {
        std::cout << "Fatal error: " << exc.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "Fatal error: unexpected exception" << std::endl;
    }

    return (EXIT_SUCCESS);
}
