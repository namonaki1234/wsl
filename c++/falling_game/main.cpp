#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>

const int BOARD_WIDTH = 10;
const int BOARD_HEIGHT = 20;
const int TILE_SIZE = 30;

int main() {
    // 1. 盤面データ (0:空, 1:ブロック)
    std::vector<std::vector<int>> board(BOARD_HEIGHT, std::vector<int>(BOARD_WIDTH, 0));

    // 操作中のブロック位置
    int blockX = 5;
    int blockY = 0;

    // 2. ウィンドウ作成
    sf::RenderWindow window(sf::VideoMode({BOARD_WIDTH * TILE_SIZE, BOARD_HEIGHT * TILE_SIZE}), "Falling Game B+A");
    window.setFramerateLimit(60);

    sf::Clock clock;

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        // --- A: ロジック（時間経過で落下） ---
        if (clock.getElapsedTime().asSeconds() >= 1.0f) {
            blockY++;
            if (blockY >= BOARD_HEIGHT) blockY = 0; // ループ
            clock.restart();

            // ターミナルに表示（Aの確認）
            system("clear");
            std::cout << "Terminal Debug View (A):" << std::endl;
            for (int y = 0; y < BOARD_HEIGHT; ++y) {
                std::cout << "|";
                for (int x = 0; x < BOARD_WIDTH; ++x) {
                    if (y == blockY && x == blockX) std::cout << "[]";
                    else std::cout << "  ";
                }
                std::cout << "|" << std::endl;
            }
        }

        // --- B: グラフィック（ウィンドウ描画） ---
        window.clear(sf::Color::Black);

        sf::RectangleShape shape(sf::Vector2f(TILE_SIZE - 2.f, TILE_SIZE - 2.f));
        shape.setFillColor(sf::Color::Cyan);
        shape.setPosition({(float)blockX * TILE_SIZE, (float)blockY * TILE_SIZE});

        window.draw(shape);
        window.display();
    }

    return 0;
}