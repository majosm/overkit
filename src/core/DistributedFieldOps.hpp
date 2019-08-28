// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_FIELD_OPS_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_FIELD_OPS_HPP_INCLUDED

#include <ovk/core/DistributedField.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Partition.hpp>

namespace ovk {
namespace core {

long long CountDistributedMask(const distributed_field<bool> &Mask);

}}

#endif
