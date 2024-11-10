#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <fstream>
#include <stdexcept>

using namespace std;
class IGridInterpolation {
public:
    virtual double interpolate(const vector<double>& grid, double x, double step) const = 0;
    virtual ~IGridInterpolation() = default;
};

class NearestNeighborInterpolation : public IGridInterpolation {
public:
    double interpolate(const vector<double>& grid, double x, double step) const override {
        int idx = static_cast<int>(round(x / step));
        if (idx < 0) idx = 0;
        if (idx >= grid.size()) idx = grid.size() - 1;
        return grid[idx];
    }
};

class LinearInterpolation : public IGridInterpolation {
public:
    double interpolate(const vector<double>& grid, double x, double step) const override {
        int idx = static_cast<int>(floor(x / step));
        if (idx < 0) return grid[0];
        if (idx >= grid.size() - 1) return grid.back();

        double t = (x - idx * step) / step;
        return (1 - t) * grid[idx] + t * grid[idx + 1];
    }
};

class Grid {
public:
    Grid(double step, int nodes, const vector<double>& values)
        : step_(step), nodes_(nodes), values_(values) {}

    static Grid loadFromFile(const string& filename, double step) {
        ifstream file(filename, ios::binary);
        if (!file) {
            throw runtime_error("Cannot open file for reading");
        }

        vector<double> values;
        double value;
        while (file.read(reinterpret_cast<char*>(&value), sizeof(double))) {
            values.push_back(value);
        }
        file.close();

        int nodes = values.size();
        return Grid(step, nodes, values);
    }

    void saveToRawFile(const string& filename) const {
        ofstream file(filename, ios::binary);
        if (!file) {
            throw runtime_error("Cannot open file for writing");
        }
        file.write(reinterpret_cast<const char*>(values_.data()), values_.size() * sizeof(double));
        if (!file) {
            throw runtime_error("Error writing to file");
        }
    }

    Grid Upscale(int k, const IGridInterpolation& interpolator) const {
        double newStep = step_ / k;
        int newNodes = nodes_ * k;
        vector<double> newValues(newNodes);

        for (int i = 0; i < newNodes; ++i) {
            double x = i * newStep;
            newValues[i] = interpolator.interpolate(values_, x, step_);
        }
        return Grid(newStep, newNodes, newValues);
    }

    Grid Downscale(int k, const IGridInterpolation& interpolator) const {
        double newStep = step_ * k;
        int newNodes = nodes_ / k;
        vector<double> newValues(newNodes);

        for (int i = 0; i < newNodes; ++i) {
            double x = i * newStep;
            newValues[i] = interpolator.interpolate(values_, x, step_);
        }
        return Grid(newStep, newNodes, newValues);
    }

    void print() const {
        for (double value : values_) {
            cout << value << " ";
        }
        cout << endl;
    }

private:
    double step_;
    int nodes_;
    vector<double> values_;
};

// Генерация функции f(x) = 10 * sin(x) + sin(10 * x) на сетке с заданным шагом
vector<double> generateFunction(double step, int nodes) {
    vector<double> values(nodes);
    for (int i = 0; i < nodes; ++i) {
        double x = i * step;
        values[i] = 10 * sin(x) + sin(10 * x);
    }
    return values;
}

int main() {
    // Параметры сетки
    double step = 0.1;
    int nodes = 300;

    // Генерация значений функции
    vector<double> values = generateFunction(step, nodes);

    // Создание исходной сетки
    Grid originalGrid(step, nodes, values);

    originalGrid.saveToRawFile("original_grid.raw");

    // Загрузка сетки из файла
   
    //Grid loadedGrid = Grid::loadFromFile("grid1.raw", step);
    Grid loadedGrid = Grid::loadFromFile("grid2.raw", step);
    //Grid loadedGrid = Grid::loadFromFile("original_grid.raw", step);
    cout << "Loaded grid: ";
    loadedGrid.print();

    // Интерполяторы
    NearestNeighborInterpolation nearestNeighbor;
    LinearInterpolation linear;

    
    Grid upscaled_Linear = loadedGrid.Upscale(4, linear);
    Grid upscaled_nearestNeighbor = loadedGrid.Upscale(4, nearestNeighbor);
    Grid downscaled_Linear = loadedGrid.Downscale(4, linear);
    Grid downscaled_nearestNeighbor = loadedGrid.Downscale(4, nearestNeighbor);

    
    upscaled_Linear.saveToRawFile("upscaled_Linear.raw");
    upscaled_nearestNeighbor.saveToRawFile("upscaled_nearestNeighbor.raw");
    downscaled_Linear.saveToRawFile("downscaled_Linear.raw");
    downscaled_nearestNeighbor.saveToRawFile("downscaled_nearestNeighbor.raw");

    cout << "Results saved to .raw files successfully." << endl;

    return 0;
}
