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

    // Параметры сетки
    float hz;   // Шаг по вертикали
    float hr;   // Шаг по радиусу
    float maxR; // Максимальный радиус
    int enclosingNodesCount; // Число узлов для вмещающих пластов
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

// Чтение данных из medium.txt
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

// Чтение данных о зондажных измерениях
vector<Probe> ReadProbes(const string& filename) {
    ifstream file(filename);
    vector<Probe> probes;

    int m; // Число зондов
    file >> m;

    for (int i = 0; i < m; ++i) {
        Probe probe;
        file >> probe.AM;
        file >> probe.MN;
        file >> probe.current;
        file >> probe.measurementsCount;

        // Считывание глубин для зонда
        ifstream depthFile("probe_" + to_string(i + 1) + ".txt");
        for (int j = 0; j < probe.measurementsCount; ++j) {
            float depth;
            depthFile >> depth;
            probe.depths.push_back(depth);
        }
        depthFile.close();

        // Считывание разностей потенциалов
        ifstream dataFile("probe_data_" + to_string(i + 1) + ".txt");
        for (int j = 0; j < probe.measurementsCount; ++j) {
            float potentialDifference;
            dataFile >> potentialDifference;
            probe.potentialDifferences.push_back(potentialDifference);
        }
        dataFile.close();

        probes.push_back(probe);
    }

    file.close();
    return probes;
}

// Чтение параметров сетки из файла grid.txt
void ReadGrid(const string& filename, Medium& medium) {
    ifstream file(filename);
    file >> medium.hz;
    file >> medium.hr;
    file >> medium.maxR;
    file >> medium.enclosingNodesCount;
    file.close();
}

// Дискретизация среды и заполнение значениями сопротивлений
vector<float> DiscretizeMedium(const Medium& medium) {
    // Количество узлов по R и Z
    size_t nodesR = static_cast<size_t>(medium.maxR / medium.hr);
    size_t nodesZ = 0;

    // Подсчет общего числа узлов по Z
    for (const Layer& layer : medium.layers) {
        nodesZ += static_cast<size_t>(layer.thickness / medium.hz);
    }

    // Общий размер сетки, включая вмещающие пласты
    nodesZ += 2 * medium.enclosingNodesCount;

    vector<float> grid(nodesZ * nodesR, 0);

    // Внесение сопротивлений вмещающих пластов
    for (size_t z = 0; z < medium.enclosingNodesCount; ++z) {
        for (size_t r = 0; r < nodesR; ++r) {
            grid[z * nodesR + r] = medium.upperResistance;
            grid[(nodesZ - z - 1) * nodesR + r] = medium.lowerResistance;
        }
    }
    // Внесение сопротивлений для каждого пласта
    size_t zStart = medium.enclosingNodesCount;
    for (const Layer& layer : medium.layers) {
        size_t layerNodesZ = static_cast<size_t>(layer.thickness / medium.hz);
        for (size_t z = zStart; z < zStart + layerNodesZ; ++z) {
            for (size_t r = 0; r < nodesR; ++r) {
                float radius = r * medium.hr;
                if (radius < layer.zone.radius) {
                    grid[z * nodesR + r] = layer.zone.resistance;
                }
                else {
                    grid[z * nodesR + r] = layer.resistance;
                }
            }
        }
        zStart += layerNodesZ;
    }

    // Внесение сопротивления бурового раствора в пределах радиуса скважины
    for (size_t z = 0; z < nodesZ; ++z) {
        for (size_t r = 0; r < nodesR; ++r) {
            if (r * medium.hr < medium.boreholeRadius) {
                grid[z * nodesR + r] = medium.boreholeResistance;
            }
        }
    }

    return grid;
}

// Запись дискретизированной среды в бинарный файл medium.raw
void WriteMediumBinary(const vector<float>& grid, const string& filename) {
    ofstream outfile(filename, ios::binary);
    outfile.write(reinterpret_cast<const char*>(grid.data()), grid.size() * sizeof(float));
    outfile.close();
}

// Запись данных зонда в бинарный файл
void WriteProbeDataBinary(const Probe& probe, const string& filename) {
    ofstream outfile(filename, ios::binary);

    for (size_t i = 0; i < probe.depths.size(); ++i) {
        float depth = probe.depths[i];
        float potentialDifference = probe.potentialDifferences[i];
        outfile.write(reinterpret_cast<char*>(&depth), sizeof(depth));
        outfile.write(reinterpret_cast<char*>(&potentialDifference), sizeof(potentialDifference));
    }

    outfile.close();
}

int main() {
    // Чтение данных из файлов
    Medium medium = ReadMedium("medium.txt");
    ReadGrid("grid.txt", medium);
    vector<Probe> probes = ReadProbes("probes.txt");

    // Дискретизация среды и запись в файл
    vector<float> grid = DiscretizeMedium(medium);
    WriteMediumBinary(grid, "medium.raw");

    // Запись данных зондов в бинарные файлы
    for (size_t i = 0; i < probes.size(); ++i) {
        WriteProbeDataBinary(probes[i], "probe_" + to_string(i + 1) + ".raw");
    }
    return 0;
}