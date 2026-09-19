#include <TCanvas.h>
#include <TEllipse.h>
#include <TArc.h>
#include <TLine.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TPad.h>
#include <TH1F.h>

#include <cmath>

namespace {
void rotate_assembly_y(double &x, double &z, double angle)
{
  const double x0 = x;
  const double z0 = z;
  x = std::cos(angle) * x0 + std::sin(angle) * z0;
  z = -std::sin(angle) * x0 + std::cos(angle) * z0;
}
}

void
irt_geometry()
{
  constexpr double mm = 1.e-3;
  double mirror_x = 1150.08 * mm, mirror_z = 939.0 * mm;
  constexpr double mirror_r = 2203.01 * mm;
  double mirror_pivot_x = 327.9 * mm, mirror_pivot_z = 3110.6 * mm;
  double detector_sphere_x = 1834.0 * mm, detector_sphere_z = 1414.0 * mm;
  constexpr double detector_sphere_r = 1100.0 * mm;
  double detector_x = 1163.3 * mm, detector_z = 2282.6 * mm;
  constexpr double emission_z = 2534.0 * mm;
  constexpr double assembly_rotation_y = -0.081479300;

  rotate_assembly_y(mirror_x, mirror_z, assembly_rotation_y);
  rotate_assembly_y(mirror_pivot_x, mirror_pivot_z, assembly_rotation_y);
  rotate_assembly_y(detector_sphere_x, detector_sphere_z, assembly_rotation_y);
  rotate_assembly_y(detector_x, detector_z, assembly_rotation_y);

  auto canvas = new TCanvas("cIRTGeometry", "IRT nominal geometry", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  auto frame = canvas->DrawFrame(-1., -1., 4., 4.);
  frame->SetTitle("Nominal IRT geometry; z (m); x (m)");
  frame->GetXaxis()->SetTitleOffset(1.5);
  frame->GetYaxis()->SetTitleOffset(1.5);
  canvas->SetBit(TPad::kClipFrame);

  auto mirror = new TEllipse(mirror_z, mirror_x, mirror_r, mirror_r);
  mirror->SetFillStyle(0); mirror->SetLineColor(kBlack); mirror->SetLineStyle(2);
  mirror->SetLineWidth(1); mirror->Draw("same");
  // TArc uses the plotted coordinates (z, x).  The mirror surface spans
  // -3 to +7 degrees around the mirror-pivot direction.
  const double pivot_angle = std::atan2(mirror_pivot_x - mirror_x,
                                       mirror_pivot_z - mirror_z) * 180. / std::acos(-1.);
  auto mirror_surface = new TArc(mirror_z, mirror_x, mirror_r,
                                  pivot_angle - 3., pivot_angle + 7.);
  mirror_surface->SetNoEdges();
  mirror_surface->SetLineColor(kAzure - 3); mirror_surface->SetLineStyle(1);
  mirror_surface->SetLineWidth(5);
  mirror_surface->Draw("same");
  auto detector_sphere = new TEllipse(detector_sphere_z, detector_sphere_x,
                                      detector_sphere_r, detector_sphere_r);
  detector_sphere->SetFillStyle(0); detector_sphere->SetLineColor(kBlack);
  detector_sphere->SetLineStyle(2); detector_sphere->SetLineWidth(1);
  detector_sphere->Draw("same");

  // The detector plane is represented by a short tangent segment in this projection.
  const double plane_half = 100. * mm;
  // In this z-x projection the detector plane is tangent when its direction
  // is perpendicular to the radius from the sphere centre to the plane centre.
  const double radius_z = detector_z - detector_sphere_z;
  const double radius_x = detector_x - detector_sphere_x;
  const double radius_norm = std::hypot(radius_z, radius_x);
  const double tangent_z = -radius_x / radius_norm;
  const double tangent_x = radius_z / radius_norm;
  const double dz = plane_half * tangent_z;
  const double dx = plane_half * tangent_x;
  auto plane = new TLine(detector_z - dz, detector_x - dx,
                         detector_z + dz, detector_x + dx);
  plane->SetLineColor(kAzure - 3); plane->SetLineWidth(5); plane->Draw("same");

  auto track = new TArrow(-1., 0., 3500. * mm, 0., .015, ">");
  track->SetLineColor(kBlack); track->SetFillColor(kBlack);
  track->SetArrowSize(0.008);
  track->SetLineWidth(1); track->Draw(">");
  auto ip = new TMarker(0., 0., 20);
  ip->SetMarkerColor(kBlack); ip->SetMarkerSize(0.8); ip->Draw("same");
  auto emission = new TMarker(emission_z, 0., 20);
  emission->SetMarkerColor(kAzure - 3); emission->SetMarkerSize(0.8); emission->Draw("same");

  auto mirror_pivot = new TMarker(mirror_pivot_z, mirror_pivot_x, 24);
  mirror_pivot->SetMarkerColor(kAzure - 3); mirror_pivot->SetMarkerSize(0.8);
  mirror_pivot->Draw("same");

  // Representative photon path; the numerical ray is calculated by irt.
  const double mirror_arc_z = mirror_z + mirror_r * std::cos(pivot_angle * std::acos(-1.) / 180.);
  const double mirror_arc_x = mirror_x + mirror_r * std::sin(pivot_angle * std::acos(-1.) / 180.);
  auto ray_in = new TLine(emission_z, 0., mirror_arc_z, mirror_arc_x);
  ray_in->SetLineColor(kAzure - 3); ray_in->SetLineStyle(2);
  ray_in->SetLineWidth(2); ray_in->Draw("same");
  auto ray_out = new TLine(mirror_arc_z, mirror_arc_x, detector_z, detector_x);
  ray_out->SetLineColor(kAzure - 3); ray_out->SetLineStyle(2);
  ray_out->SetLineWidth(2); ray_out->Draw("same");

  canvas->Modified();
  canvas->Update();
}
