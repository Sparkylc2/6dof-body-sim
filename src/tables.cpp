#include "tables.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace BSIM {

namespace {

// the segment [x[k], x[k+1]] that xq is in and how far along it is (t in [0, 1])
// past the ends it keeps the end segment with t clamped, so the value is flat
struct Segment {
    std::size_t k = 0;
    double t = 0.0;
};

Segment find_segment(std::span<const double> x, double xq) {
    const std::size_t n = x.size();

    if (n < 2 || xq <= x.front()) {
        return {0, 0.0};
    }

    if (xq >= x.back()) {
        return {n - 2, 1.0};
    }

    // first breakpoint strictly above xq, so x[k] <= xq < x[k+1]
    const std::size_t k =
        static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), xq) - x.begin()) - 1;

    return {k, (xq - x[k]) / (x[k + 1] - x[k])};
}

void check_breakpoints(std::span<const double> x, const char *name) {
    if (x.empty()) {
        throw std::invalid_argument(std::string("table: no ") + name + " breakpoints");
    }

    for (std::size_t i = 1; i < x.size(); i++) {

        if (!(x[i] > x[i - 1])) {
            throw std::invalid_argument(std::string("table: ") + name +
                                        " breakpoints aren't strictly increasing");
        }
    }
}

// --- csv ---
// splits a line on commas (and whitespace round them) into numbers,
// false if any cell isn't a number
bool parse_row(const std::string &line, std::vector<double> &out) {
    out.clear();
    std::stringstream ss(line);
    std::string cell;

    while (std::getline(ss, cell, ',')) {
        std::size_t used = 0;
        try {
            out.push_back(std::stod(cell, &used));
        } catch (const std::exception &) {
            return false;
        }

        // anything other than whitespace after the number means it isn't one
        for (std::size_t i = used; i < cell.size(); i++) {
            if (!std::isspace(static_cast<unsigned char>(cell[i]))) {
                return false;
            }
        }
    }

    return !out.empty();
}

bool skip_line(const std::string &line) {
    const auto first = line.find_first_not_of(" \t\r");
    return first == std::string::npos || line[first] == '#';
}

// every row of the file as text, without the skipped lines
std::vector<std::string> read_lines(const std::string &path) {
    std::ifstream file(path);

    if (!file) {
        throw std::runtime_error("table: can't open " + path);
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(file, line)) {
        if (!skip_line(line))
            lines.push_back(line);
    }

    return lines;
}

std::runtime_error parse_error(const std::string &path, std::size_t row) {
    return std::runtime_error("table: couldn't parse row " + std::to_string(row + 1) + " of " +
                              path);
}

} // namespace

double lookup1(std::span<const double> x, std::span<const double> y, double xq) {
    if (x.size() == 1)
        return y[0];

    const Segment s = find_segment(x, xq);
    return (1.0 - s.t) * y[s.k] + s.t * y[s.k + 1];
}

double lookup2(std::span<const double> x, std::span<const double> y, std::span<const double> z,
               double xq, double yq) {
    const std::size_t nx = x.size(), ny = y.size();
    const Segment sx = find_segment(x, xq);
    const Segment sy = find_segment(y, yq);

    // with only one breakpoint in a direction there's nothing to blend along it
    const std::size_t i1 = nx > 1 ? sx.k + 1 : sx.k;
    const std::size_t j1 = ny > 1 ? sy.k + 1 : sy.k;
    auto at = [&](std::size_t i, std::size_t j) { return z[i + nx * j]; };

    // bilinear, along x at both y's then along y
    const double z0 = (1.0 - sx.t) * at(sx.k, sy.k) + sx.t * at(i1, sy.k);
    const double z1 = (1.0 - sx.t) * at(sx.k, j1) + sx.t * at(i1, j1);
    return (1.0 - sy.t) * z0 + sy.t * z1;
}

void check_table(std::span<const double> x, std::span<const double> y) {
    check_breakpoints(x, "x");

    if (y.size() != x.size()) {
        throw std::invalid_argument("table: y isn't the same length as x");
    }
}

void check_table(std::span<const double> x, std::span<const double> y, std::span<const double> z) {
    check_breakpoints(x, "x");
    check_breakpoints(y, "y");

    if (z.size() != x.size() * y.size()) {
        throw std::invalid_argument("table: z isn't nx by ny");
    }
}

Table1D::Table1D(std::vector<double> x_, std::vector<double> y_)
    : x(std::move(x_)), y(std::move(y_)) {
    check_table(x, y);
}

Table2D::Table2D(std::vector<double> x_, std::vector<double> y_, std::vector<double> z_)
    : x(std::move(x_)), y(std::move(y_)), z(std::move(z_)) {
    check_table(x, y, z);
}

Table1D load_table1d(const std::string &path) {
    const std::vector<std::string> lines = read_lines(path);
    std::vector<double> x, y, row;

    for (std::size_t r = 0; r < lines.size(); r++) {
        const bool numbers = parse_row(lines[r], row);

        if (!numbers && r == 0) {
            continue; // header
        }

        if (!numbers || row.size() != 2) {
            throw parse_error(path, r);
        }

        x.push_back(row[0]);
        y.push_back(row[1]);
    }

    return Table1D(std::move(x), std::move(y));
}

Table2D load_table2d(const std::string &path) {
    const std::vector<std::string> lines = read_lines(path);

    if (lines.size() < 2) {
        throw std::runtime_error("table: " + path + " needs a breakpoint row and at least one row");
    }

    // first row is a label then the y breakpoints,
    const std::string &head = lines[0];
    const auto comma = head.find(',');
    std::vector<double> y;

    if (comma == std::string::npos || !parse_row(head.substr(comma + 1), y)) {
        throw parse_error(path, 0);
    }

    // rows are read in x order, but z wants to be column major (x down the columns)
    std::vector<double> x;
    std::vector<std::vector<double>> rows;
    std::vector<double> row;

    for (std::size_t r = 1; r < lines.size(); r++) {

        if (!parse_row(lines[r], row) || row.size() != y.size() + 1) {
            throw parse_error(path, r);
        }

        x.push_back(row[0]);
        rows.emplace_back(row.begin() + 1, row.end());
    }

    const std::size_t nx = x.size(), ny = y.size();
    std::vector<double> z(nx * ny);

    for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 0; j < ny; j++) {
            z[i + nx * j] = rows[i][j];
        }
    }

    return Table2D(std::move(x), std::move(y), std::move(z));
}

} // namespace BSIM
