#include <iostream>
#include <SFML/Graphics.hpp>
#include "Fluid.h"

int main() {
    std::cout << "2D PIC/FLIP Simulation.\n" << std::endl;
    Fluid fluid;
    // define viewing window
    sf::RenderWindow window(sf::VideoMode(LENGTH, HEIGHT), "2D PIC/FLIP Simulator");
    // define circles that act as particles
    std::vector<sf::CircleShape> circles(fluid.numParticles);
    // set the initial radius and color of the particles
    for (auto& circle : circles)    {
        circle.setRadius(fluid.particleRadius);
        circle.setFillColor(sf::Color::Blue);
    }
    std:: cout << fluid.numParticles;
    for (int i = 0; i < fluid.numParticles; i++)    {
        circles[i].setPosition(fluid.particleXPositions[i], fluid.particleYPositions[i]);
    }

    // define main game loop
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)    {
                window.close();
            }
        }
        window.clear(sf::Color::White);
        // draw the particles on screen
        for (auto& circle: circles) {
            window.draw(circle);
        }
        // update the display
        window.display();
        // update fluid state
        fluid.simulateFluid();
        for (int i = 0; i < fluid.numParticles; i++)    {
            circles[i].setPosition(fluid.particleXPositions[i], fluid.particleYPositions[i]);
        }
    }

    return 0;
}
