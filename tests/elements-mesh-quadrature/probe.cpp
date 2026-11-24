#include <precice/Participant.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
  precice::Participant probe("Probe", "precice-config.xml", 0, 1);
  const std::array<double, 6> region{-5., 13., -1., 3., -1., 4.};
  for (const auto mesh : {"Solid-Nodes", "Tets", "Hexes"}) {
    probe.setMeshAccessRegion(mesh, region);
  }
  probe.initialize();

  for (const std::string mesh : {"Tets", "Hexes"}) {
    const auto count = probe.getMeshVertexSize(mesh);
    std::vector<precice::VertexID> ids(count);
    std::vector<double> coordinates(3 * count);
    probe.getMeshVertexIDsAndCoordinates(mesh, ids, coordinates);

    std::ofstream output(mesh + ".txt");
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << std::setprecision(17);
    for (auto i = 0; i < count; ++i) {
      output << coordinates[3 * i] << ' '
             << coordinates[3 * i + 1] << ' '
             << coordinates[3 * i + 2] << '\n';
    }
  }

  std::array<precice::VertexID, 1> node;
  std::array<double, 3> coordinates;
  std::array<double, 3> displacement;
  const double dt = probe.getMaxTimeStepSize();
  probe.getMeshVertexIDsAndCoordinates("Solid-Nodes", node, coordinates);
  probe.readData("Solid-Nodes", "Displacement", node, dt, displacement);
  if (!std::isfinite(displacement[2]) || displacement[2] >= 0.0) {
    throw std::runtime_error("Expected a negative Z displacement at the loaded node");
  }
  probe.advance(dt);

  if (probe.isCouplingOngoing()) {
    throw std::runtime_error("Expected one completed coupling window");
  }
  probe.finalize();
}
