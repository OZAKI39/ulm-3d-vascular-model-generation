#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>
namespace frozen {
using Vec=std::array<double,3>;
enum class Status { VALID=0, OUTSIDE=1, SOLID=2, MISSING=3, NONFINITE_POSITION=4 };
struct Query { Status status; Vec u; };
class FlowField {
public:
 explicit FlowField(const std::string& path);
 std::array<int64_t,3> dims; Vec origin; double dx;
 std::vector<uint64_t> indices,fluid_indices; std::vector<double> velocity;
 Status node(uint64_t id,size_t& slot) const;
};
class FlowFieldSampler {
public:
 explicit FlowFieldSampler(const FlowField& f):field(f){}
 Query query(Vec x) const;
private: const FlowField& field;
};
class TechnicalStokesDrag {
public: static Vec force(double mu,double diameter,Vec fluid,Vec particle);
};
}
