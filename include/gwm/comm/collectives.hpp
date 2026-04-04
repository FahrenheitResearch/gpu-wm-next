#pragma once

#include <cstdint>
#include <vector>

namespace gwm::comm {

void allreduce_sum_in_place(std::vector<double>& values);

void allreduce_sum_in_place(std::vector<std::uint64_t>& values);

}  // namespace gwm::comm
