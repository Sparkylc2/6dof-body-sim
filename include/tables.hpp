#pragma once
#include <span>
#include <string>
#include <vector>

// --- lookup tables ---
// for anything given as data instead of a formula (thrust curves, CD vs mach, etc)
// everything is linear interpolation and held flat past the ends (no extrapolation), and the
// breakpoints have to be strictly increasing
namespace BSIM {

// --- views ---
// these don't own anything, so the same code can run on arrays from anywhere (the matlab bridge
// points them straight at matlab's arrays)
// they don't check anything either, the tables below check once when they're made

// y(xq), x and y the same length
double lookup1(std::span<const double> x, std::span<const double> y, double xq);

// z(xq, yq), z is nx by ny and column major, so z[i + nx * j] is the value at (x[i], y[j])
// (the same layout as a matlab matrix with x down the rows and y along the columns)
double lookup2(std::span<const double> x, std::span<const double> y, std::span<const double> z,
               double xq, double yq);

// throws std::invalid_argument if the breakpoints are empty or not strictly increasing, or the
// values are the wrong size
void check_table(std::span<const double> x, std::span<const double> y);
void check_table(std::span<const double> x, std::span<const double> y, std::span<const double> z);

// --- tables ---
// versions with param owning
struct Table1D {
    std::vector<double> x, y;

    Table1D() = default;
    Table1D(std::vector<double> x, std::vector<double> y);

    double operator()(double xq) const { return lookup1(x, y, xq); }
};

struct Table2D {
    std::vector<double> x, y, z; // z is column major, see lookup2

    Table2D() = default;
    Table2D(std::vector<double> x, std::vector<double> y, std::vector<double> z);

    double operator()(double xq, double yq) const { return lookup2(x, y, z, xq, yq); }
};

// --- csv files ---
// blank lines and lines starting with # are skipped, and so is a first line that isn't numbers
// (so a header is fine). throws std::runtime_error if the file can't be read or parsed
//
// 1d: two columns, x then y
//      mach, cd
//      0.0, 0.45
//      0.8, 0.47
//
// 2d: the first row is the y breakpoints (the first cell is just a label), every row after that
// is an x breakpoint then the values along y
//      mach\alpha, 0.0, 0.1, 0.2
//      0.0,        0.45, 0.48, 0.55
//      0.8,        0.47, 0.50, 0.58
Table1D load_table1d(const std::string &path);
Table2D load_table2d(const std::string &path);

} // namespace BSIM
