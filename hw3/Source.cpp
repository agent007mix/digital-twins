#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// ��������� ��� �������� ���������� � ���� �������������
struct PenetrationZone {
    float radius;
    float resistance;
};

// ��������� ��� �������� ���������� � ������
struct Layer {
    float thickness;
    float resistance;
    PenetrationZone zone;
};

// ��������� ��� �������� ���������� � �������� � �������
struct Medium {
    float boreholeRadius;
    float boreholeResistance;
    float upperResistance;
    float lowerResistance;
    vector<Layer> layers;
};

// ��������� ��� �������� ���������� � �����
struct Probe {
    float AM;
    float MN;
    float current;
    int measurementsCount;
    vector<float> depths;
    vector<float> potentialDifferences;
};

// ��������� ��� �������� ������ �����
struct Grid {
    float hz;
    float hr;
    float max_r;
    int n_vm;
};

// ������� ������ ������ �� grid.txt
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

// ������ medium.txt
Medium ReadMedium(const string& filename) {
    ifstream file(filename);
    Medium medium;

    file >> medium.boreholeRadius;
    file >> medium.boreholeResistance;
    file >> medium.upperResistance;
    file >> medium.lowerResistance;

    int N;
    file >> N; // ����� �������

    for (int i = 0; i < N; ++i) {
        Layer layer;
        file >> layer.thickness;
        file >> layer.resistance;
        file >> layer.zone.radius;
        file >> layer.zone.resistance;

        if (layer.zone.radius == 0 && layer.zone.resistance == 0) {
            layer.zone.radius = -1; // ���� ���� ������������� �����������
        }

        medium.layers.push_back(layer);
    }

    file.close();
    return medium;
}

// ������� ��� �������� � ���������� ����� ���������������
vector<vector<float>> CreateResistivityGrid(const Medium& medium, const Grid& grid) {
    int r_nodes = static_cast<int>(grid.max_r / grid.hr) + 1; // ���������� ����� �� �������
    int z_nodes = static_cast<int>((medium.upperResistance + medium.lowerResistance) / grid.hz) + grid.n_vm * 2;

    
    vector<vector<float>> resistivityGrid(z_nodes, vector<float>(r_nodes));

    // ��������� ����� ���������������� ���������� �������������
    for (int z = 0; z < z_nodes; ++z) {
        for (int r = 0; r < r_nodes; ++r) {
            float radius = r * grid.hr;

            // ������ ��������
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
// ������� ��� ������ ����� � �������� ����
void WriteResistivityGridBinary(const vector<vector<float>>& grid, const string& filename) {
    ofstream outfile(filename, ios::binary);
    for (const auto& row : grid) {
        outfile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(float));
    }

    outfile.close();
}

int main() {
    // ������ ������
    Medium medium = ReadMedium("medium.txt");
    Grid grid = ReadGrid("grid.txt");

    // �������� ����� �������������
    vector<vector<float>> resistivityGrid = CreateResistivityGrid(medium, grid);

    // ������ ����� � ����
    WriteResistivityGridBinary(resistivityGrid, "medium.raw");

    return 0;
}
