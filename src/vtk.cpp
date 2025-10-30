#include "vtk.hpp"
#include <fstream>
#include <iomanip>

namespace hc
{
    void writeVTK(const std::string &filename,
                  const GridGeo &geo,
                  const GridArray &arr,
                  const std::vector<std::string> &components)
    {
        std::ofstream file(filename);
        if (!file)
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        const Grid &grid = arr.grid();
        const int n_cells = grid.Nx() * grid.Ny() * grid.Nz();
        const int n_comp = arr.n_comp();

        const double x0 = geo.start(0);
        const double y0 = geo.start(1);
        const double z0 = geo.start(2);
        const double dx = geo.dx();
        const double dy = geo.dy();
        const double dz = geo.dz();

        file << "# vtk DataFile Version 3.0\n";
        file << "Cartesian cell-centered data\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << grid.Nx() << " " << grid.Ny() << " " << grid.Nz() << "\n";
        file << "ORIGIN " << x0 + 0.5 * dx << " " << y0 + 0.5 * dy << " " << z0 + 0.5 * dz << "\n";
        file << "SPACING " << dx << " " << dy << " " << dz << "\n";
        file << "POINT_DATA " << n_cells << "\n";

        // Write each component as a separate scalar field
        for (int n = 0; n < n_comp; ++n)
        {
            file << "SCALARS " << components[n] << " double\n ";
            file
                << "LOOKUP_TABLE default\n";
            for (int k = 0; k < grid.Nz(); ++k)
            {
                for (int j = 0; j < grid.Ny(); ++j)
                {
                    for (int i = 0; i < grid.Nx(); ++i)
                    {
                        file << std::setprecision(10) << arr(i, j, k, n) << "\n";
                    }
                }
            }
        }

        file.close();
    }
}