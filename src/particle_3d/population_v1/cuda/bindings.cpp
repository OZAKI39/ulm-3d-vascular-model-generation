#include <torch/extension.h>
#include <vector>
std::vector<torch::Tensor> fem_sample_cuda(std::vector<torch::Tensor> a, double padding);
torch::Tensor neighbor_count_cuda(std::vector<torch::Tensor> a, double cell_size);
std::vector<torch::Tensor> neighbor_fill_cuda(std::vector<torch::Tensor> a, double cell_size, torch::Tensor offsets, int64_t total);
std::vector<torch::Tensor> wall_query_cuda(std::vector<torch::Tensor> a, double roundoff);
torch::Tensor wall_path_cuda(std::vector<torch::Tensor> a, double roundoff);
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.def("fem_sample", &fem_sample_cuda, "Canonical float64 P1 FEM CUDA sample");
  m.def("neighbor_count", &neighbor_count_cuda);
  m.def("neighbor_fill", &neighbor_fill_cuda);
  m.def("wall_query", &wall_query_cuda);
  m.def("wall_path", &wall_path_cuda);
}
