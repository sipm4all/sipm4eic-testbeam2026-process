#include "../lib/trigger_reader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

constexpr double ring_x0_mean = 6.17236;
constexpr double ring_x0_sigma = 0.793097;
constexpr double ring_y0_mean = -2.18484;
constexpr double ring_y0_sigma = 0.825066;

void
clusters(const char *input_filename = "triggered.root",
         const char *output_filename = "clusters.root",
         double adjacency_distance = 3.2,
         int large_cluster_size = 4)
{
  trigger_reader_t reader;
  if (!reader.open(input_filename))
    return;
  TFile output(output_filename, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << output_filename << std::endl;
    return;
  }
  auto size_hist = new TH1D("hClusterSize", "Cherenkov cluster size;cluster size;clusters", 256, -.5, 255.5);
  auto count_hist = new TH1D("hClustersPerFrame", "clusters per frame;clusters/frame;frames", 256, -.5, 255.5);
  auto largest_hist = new TH1D("hLargestCluster", "largest cluster;hits;frames", 256, -.5, 255.5);
  auto fraction_hist = new TH1D("hLargestClusterFraction", "largest cluster fraction;largest cluster / total hits;frames", 100, 0., 1.);
  auto count_vs_hits = new TH2D("hClustersVsHits", "clusters versus hits;hits/frame;clusters/frame", 512, -.5, 511.5, 256, -.5, 255.5);
  auto count_vs_largest = new TH2D("hClustersVsLargestCluster",
                                   "clusters versus largest cluster;clusters/frame;largest cluster size",
                                   256, -.5, 255.5, 256, -.5, 255.5);
  auto hits_vs_largest = new TH2D("hHitsVsLargestCluster",
                                  "hits versus largest cluster;largest cluster size;hits/frame",
                                  256, -.5, 255.5, 512, -.5, 511.5);
  auto large_count_hist = new TH1D("hLargeClusters",
                                   "large clusters per frame;clusters with size >= threshold;frames",
                                   128, -.5, 127.5);
  auto large_vs_hits = new TH2D("hLargeClustersVsHits",
                                "large clusters versus hits;hits/frame;large clusters/frame",
                                512, -.5, 511.5, 128, -.5, 127.5);
  auto large_vs_largest = new TH2D("hLargeClustersVsLargestCluster",
                                   "large clusters versus largest cluster;largest cluster size;large clusters/frame",
                                   256, -.5, 255.5, 128, -.5, 127.5);
  auto inliers_all = new TH1D("hRingInliersAll",
                              "ring inliers, exactly one ring;inliers;frames",
                              256, -.5, 255.5);
  auto inliers_clean = new TH1D("hRingInliersLargestCluster2",
                                "ring inliers, largest cluster <= 2;inliers;frames",
                                256, -.5, 255.5);
  auto inliers_clean3 = new TH1D("hRingInliersLargestCluster3",
                                 "ring inliers, largest cluster <= 3;inliers;frames",
                                 256, -.5, 255.5);
  auto inlier_clusters_all = new TH1D("hRingInlierClustersAll",
                                      "clusters among ring inliers, exactly one ring;clusters;frames",
                                      64, -.5, 63.5);
  auto inlier_clusters_clean2 = new TH1D("hRingInlierClustersLargestCluster2",
                                         "clusters among ring inliers, largest cluster <= 2;clusters;frames",
                                         64, -.5, 63.5);
  auto inlier_clusters_clean3 = new TH1D("hRingInlierClustersLargestCluster3",
                                         "clusters among ring inliers, largest cluster <= 3;clusters;frames",
                                         64, -.5, 63.5);
  long long frames = 0, total_clusters = 0;
  while (reader.next_spill()) while (reader.next_frame()) {
    if (reader.rings().size() != 1)
      continue;
    const auto &selected_ring = reader.rings().front();
    if (std::abs(selected_ring.x0 - ring_x0_mean) > 4. * ring_x0_sigma ||
        std::abs(selected_ring.y0 - ring_y0_mean) > 4. * ring_y0_sigma)
      continue;
    const auto &hits = reader.cherenkov_hits();
    std::vector<bool> visited(hits.size(), false);
    std::vector<int> stack, sizes;
    for (std::size_t seed = 0; seed < hits.size(); ++seed) {
      if (visited[seed]) continue;
      visited[seed] = true;
      stack.clear(); stack.push_back(static_cast<int>(seed));
      int size = 0;
      while (!stack.empty()) {
        const int current = stack.back(); stack.pop_back(); ++size;
        for (std::size_t other = 0; other < hits.size(); ++other) {
          if (visited[other]) continue;
          if (std::hypot(hits[current].x - hits[other].x,
                         hits[current].y - hits[other].y) <= adjacency_distance) {
            visited[other] = true;
            stack.push_back(static_cast<int>(other));
          }
        }
      }
      sizes.push_back(size); size_hist->Fill(size);
    }
    const int largest = sizes.empty() ? 0 : *std::max_element(sizes.begin(), sizes.end());
    count_hist->Fill(sizes.size()); largest_hist->Fill(largest);
    if (!hits.empty()) fraction_hist->Fill(static_cast<double>(largest) / hits.size());
    count_vs_hits->Fill(hits.size(), sizes.size());
    count_vs_largest->Fill(sizes.size(), largest);
    hits_vs_largest->Fill(largest, hits.size());
    const int nlarge = std::count_if(sizes.begin(), sizes.end(),
                                     [large_cluster_size](int size) {
                                       return size >= large_cluster_size;
                                     });
    large_count_hist->Fill(nlarge);
    large_vs_hits->Fill(hits.size(), nlarge);
    large_vs_largest->Fill(largest, nlarge);
    {
      const auto &ring = selected_ring;
      const int ninliers = ring.ninliers;
      inliers_all->Fill(ninliers);
      if (largest <= 2)
        inliers_clean->Fill(ninliers);
      if (largest <= 3)
        inliers_clean3->Fill(ninliers);

      // Reconstruct the ring-inlier set, since the ring tree stores only its count.
      std::vector<std::size_t> inlier_indices;
      for (std::size_t i = 0; i < hits.size(); ++i) {
        const double radial = std::hypot(hits[i].x - ring.x0,
                                         hits[i].y - ring.y0);
        if (std::abs(radial - ring.radius) <= 6. &&
            std::abs(hits[i].time - ring.time) <= 4.)
          inlier_indices.push_back(i);
      }
      std::vector<bool> inlier_visited(inlier_indices.size(), false);
      int n_inlier_clusters = 0;
      for (std::size_t seed = 0; seed < inlier_indices.size(); ++seed) {
        if (inlier_visited[seed]) continue;
        ++n_inlier_clusters;
        std::vector<std::size_t> pending{seed};
        inlier_visited[seed] = true;
        while (!pending.empty()) {
          const auto current = pending.back(); pending.pop_back();
          const auto current_hit = inlier_indices[current];
          for (std::size_t other = 0; other < inlier_indices.size(); ++other) {
            if (inlier_visited[other]) continue;
            const auto other_hit = inlier_indices[other];
            if (std::hypot(hits[current_hit].x - hits[other_hit].x,
                           hits[current_hit].y - hits[other_hit].y) < adjacency_distance) {
              inlier_visited[other] = true;
              pending.push_back(other);
            }
          }
        }
      }
      inlier_clusters_all->Fill(n_inlier_clusters);
      if (largest <= 2)
        inlier_clusters_clean2->Fill(n_inlier_clusters);
      if (largest <= 3)
        inlier_clusters_clean3->Fill(n_inlier_clusters);
    }
    total_clusters += sizes.size(); ++frames;
  }
  size_hist->Write(); count_hist->Write(); largest_hist->Write();
  fraction_hist->Write(); count_vs_hits->Write();
  count_vs_largest->Write();
  hits_vs_largest->Write();
  large_count_hist->Write();
  large_vs_hits->Write();
  large_vs_largest->Write();
  inliers_all->Write(); inliers_clean->Write(); inliers_clean3->Write();
  inlier_clusters_all->Write();
  inlier_clusters_clean2->Write();
  inlier_clusters_clean3->Write();
  output.Close();
  std::cout << "frames processed: " << frames << std::endl
            << "clusters found: " << total_clusters << std::endl
            << "adjacency distance [mm]: " << adjacency_distance << std::endl
            << "large cluster threshold: " << large_cluster_size << std::endl
            << "output: " << output_filename << std::endl;
}
