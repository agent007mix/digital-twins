#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Структура для хранения информации о зоне проникновения
struct PenetrationZone {
    float radius;
    float resistance;
};

// Структура для хранения информации о пласте
struct Layer {
    float thickness;
    float resistance;
    PenetrationZone zone;
};

// Структура для хранения информации о скважине и пластах
struct Medium {
    float boreholeRadius;
    float boreholeResistance;
    float upperResistance;
    float lowerResistance;
    vector<Layer> layers;
};

// Структура для хранения информации о зонде
struct Probe {
    float AM;
    float MN;
    float current;
    int measurementsCount;
    vector<float> depths;
    vector<float> potentialDifferences;
};

// Структура для хранения данных сетки
struct Grid {
    float hz;
    float hr;
    float max_r;
    int n_vm;
};

// Функция чтения данных из grid.txt
Grid ReadGrid(const string& filename) {
    ifstream file(filename);
    Grid grid;

    file >> grid.hz;
    file >> grid.hr;
    file >> grid.max_r;
    file >> grid.n_vm;

    file.close();
    return grid;
}

// Чтение medium.txt
Medium ReadMedium(const string& filename) {
    ifstream file(filename);
    Medium medium;

    file >> medium.boreholeRadius;
    file >> medium.boreholeResistance;
    file >> medium.upperResistance;
    file >> medium.lowerResistance;

    int N;
    file >> N; // Число пластов

    for (int i = 0; i < N; ++i) {
        Layer layer;
        file >> layer.thickness;
        file >> layer.resistance;
        file >> layer.zone.radius;
        file >> layer.zone.resistance;

        if (layer.zone.radius == 0 && layer.zone.resistance == 0) {
            layer.zone.radius = -1; // Если зона проникновения отсутствует
        }

        medium.layers.push_back(layer);
    }

    file.close();
    return medium;
}

// Функция для создания и заполнения сетки сопротивлениями
vector<vector<float>> CreateResistivityGrid(const Medium& medium, const Grid& grid) {
    int r_nodes = static_cast<int>(grid.max_r / grid.hr) + 1; // количество узлов по радиусу
    int z_nodes = static_cast<int>((medium.upperResistance + medium.lowerResistance) / grid.hz) + grid.n_vm * 2;

    
    vector<vector<float>> resistivityGrid(z_nodes, vector<float>(r_nodes));

    // Заполняем сетку соответствующими значениями сопротивлений
    for (int z = 0; z < z_nodes; ++z) {
        for (int r = 0; r < r_nodes; ++r) {
            float radius = r * grid.hr;

            // Внутри скважины
            if (radius < medium.boreholeRadius) {
                resistivityGrid[z][r] = medium.boreholeResistance;
            }
            else {
                
                if (z < medium.upperResistance) {
                    resistivityGrid[z][r] = medium.upperResistance;
                }
                else if (z >= z_nodes - grid.n_vm) {
                    resistivityGrid[z][r] = medium.lowerResistance;
                }
                else {
                    
                    int layer_idx = z % medium.layers.size();
                    const Layer& layer = medium.layers[layer_idx];

                    
                    if (radius <= layer.zone.radius && layer.zone.radius > 0) {
                        resistivityGrid[z][r] = layer.zone.resistance;
                    }
                    else {
                        resistivityGrid[z][r] = layer.resistance;
                    }
                }
            }
        }
    }
    return resistivityGrid;
}
// Функция для записи сетки в бинарный файл
void WriteResistivityGridBinary(const vector<vector<float>>& grid, const string& filename) {
    ofstream outfile(filename, ios::binary);
    for (const auto& row : grid) {
        outfile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(float));
    }

    outfile.close();
}

int main() {
    // Чтение данных
    Medium medium = ReadMedium("medium.txt");
    Grid grid = ReadGrid("grid.txt");

    // Создание сетки сопротивления
    vector<vector<float>> resistivityGrid = CreateResistivityGrid(medium, grid);

    // Запись сетки в файл
    WriteResistivityGridBinary(resistivityGrid, "medium.raw");

    return 0;
}
