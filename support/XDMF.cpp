// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/XDMF.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Handle.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Range.hpp>

#include "support/pugixml/pugixml.hpp"

#include <hdf5.h>
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <exception>
#include <string>
#include <sstream>
#include <utility>

namespace support {

xdmf::xdmf(std::string &&Path, int NumDims, ovk::comm &&Comm, ovk::array<xdmf_grid_meta> &&Grids,
  ovk::array<xdmf_attribute_meta> &&Attributes, ovk::core::handle<hid_t> &&HDF5File):
  Path_(std::move(Path)),
  NumDims_(NumDims),
  Comm_(std::move(Comm)),
  Grids_(std::move(Grids)),
  Attributes_(std::move(Attributes)),
  HDF5File_(std::move(HDF5File))
{}

xdmf xdmf::internal_Create(std::string &&Path, int NumDims, ovk::comm_view Comm_,
  ovk::array<xdmf_grid_meta> &&Grids, ovk::array<xdmf_attribute_meta> &&Attributes) {

  ovk::comm Comm = ovk::DuplicateComm(Comm_);

  std::size_t iNameStart = Path.find_last_of("/");
  if (iNameStart != std::string::npos) {
    iNameStart += 1;
  } else {
    iNameStart = 0;
  }
  std::size_t iNameEnd = Path.find_last_of(".");
  std::string Directory = Path.substr(0,iNameStart);
  std::string Name = Path.substr(iNameStart,iNameEnd-iNameStart);

  ovk::tuple<std::string> CoordNames = {"X", "Y", "Z"};

  if (Comm.Rank() == 0) {

    pugi::xml_document Document;

    pugi::xml_node DoctypeNode = Document.append_child(pugi::node_doctype);
    DoctypeNode.set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

    pugi::xml_node RootNode = Document.append_child("Xdmf");
    RootNode.append_attribute("Version") = "2.0";

    pugi::xml_node DomainNode = RootNode.append_child("Domain");

    for (auto &Grid : Grids) {

      pugi::xml_node GridNode = DomainNode.append_child("Grid");
      GridNode.append_attribute("Name") = Grid.Name().c_str();
      GridNode.append_attribute("GridType") = "Uniform";

      GridNode.append_child("Time").append_attribute("Value") = 0;

      std::string SizeString = std::to_string(Grid.Size()(NumDims-1));
      for (int iDim = NumDims-2; iDim >= 0; --iDim) {
        SizeString += " " + std::to_string(Grid.Size()(iDim));
      }

      pugi::xml_node TopologyNode = GridNode.append_child("Topology");
      TopologyNode.append_attribute("TopologyType") = "2DSMesh";
      TopologyNode.append_attribute("Dimensions") = SizeString.c_str();

      pugi::xml_node GeometryNode = GridNode.append_child("Geometry");
      std::string GeometryTypeString;
      GeometryTypeString += CoordNames(0);
      for (int iDim = 1; iDim < NumDims; ++iDim) {
        GeometryTypeString += "_" + CoordNames(iDim);
      }
      GeometryNode.append_attribute("GeometryType") = GeometryTypeString.c_str();

      for (int iDim = 0; iDim < NumDims; ++iDim) {
        pugi::xml_node DataItemNode = GeometryNode.append_child("DataItem");
        DataItemNode.append_attribute("Dimensions") = SizeString.c_str();
        DataItemNode.append_attribute("NumberType") = "Float";
        DataItemNode.append_attribute("Precision") = 8;
        DataItemNode.append_attribute("Format") = "HDF";
        std::string DataPath = Name + ".h5:/" + Grid.Name() + "/Geometry/" + CoordNames(iDim);
        DataItemNode.text().set(DataPath.c_str());
      }

      for (auto &Attribute : Attributes) {

        pugi::xml_node AttributeNode = GridNode.append_child("Attribute");
        AttributeNode.append_attribute("Name") = Attribute.Name().c_str();
        AttributeNode.append_attribute("AttributeType") = "Scalar";
        AttributeNode.append_attribute("Center") = "Node";

        pugi::xml_node DataItemNode = AttributeNode.append_child("DataItem");
        DataItemNode.append_attribute("Dimensions") = SizeString.c_str();
        switch (Attribute.Type()) {
        case xdmf_attribute_type::INT:
          DataItemNode.append_attribute("NumberType") = "Int";
          DataItemNode.append_attribute("Precision") = 4;
          break;
        case xdmf_attribute_type::LONG_LONG:
          DataItemNode.append_attribute("NumberType") = "Int";
          DataItemNode.append_attribute("Precision") = 8;
          break;
        case xdmf_attribute_type::DOUBLE:
          DataItemNode.append_attribute("NumberType") = "Float";
          DataItemNode.append_attribute("Precision") = 8;
          break;
        default:
          OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
          break;
        }
        DataItemNode.append_attribute("Format") = "HDF";
        std::string DataPath = Name + ".h5:/" + Grid.Name() + "/Attributes/" + Attribute.Name();
        DataItemNode.text().set(DataPath.c_str());

      }

    }

    Document.save_file(Path.c_str());

  }

  MPI_Barrier(Comm);

  std::string HDF5Path = Directory + Name + ".h5";

  auto FileProperties = ovk::core::MakeHandle(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
  H5Pset_fapl_mpio(FileProperties, Comm, MPI_INFO_NULL);
  hid_t FileRaw = H5Fcreate(HDF5Path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, FileProperties);
  if (FileRaw < 0) {
    if (Comm.Rank() == 0) {
      std::fprintf(stderr, "ERROR: Unable to create file '%s'.\n", HDF5Path.c_str());
      std::fflush(stderr);
    }
    throw xdmf_file_create_error(HDF5Path);
  }
  auto HDF5File = ovk::core::MakeHandle(FileRaw, H5Fclose);

  for (auto &Grid : Grids) {
    auto GridGroup = ovk::core::MakeHandle(H5Gcreate(HDF5File, Grid.Name().c_str(), H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    auto GeometryGroup = ovk::core::MakeHandle(H5Gcreate(GridGroup, "Geometry", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    ovk::tuple<hsize_t> DataspaceSize = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      DataspaceSize(iDim) = Grid.Size()(NumDims-1-iDim);
    }
    auto Dataspace = ovk::core::MakeHandle(H5Screate_simple(NumDims, DataspaceSize.Data(),
      nullptr), H5Sclose);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      ovk::core::MakeHandle(H5Dcreate(GeometryGroup, CoordNames(iDim).c_str(), H5T_NATIVE_DOUBLE,
        Dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    }
    auto AttributesGroup = ovk::core::MakeHandle(H5Gcreate(GridGroup, "Attributes", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    for (auto &Attribute : Attributes) {
      hid_t HDF5Datatype = 0;
      switch (Attribute.Type()) {
      case xdmf_attribute_type::INT:
        HDF5Datatype = H5T_NATIVE_INT;
        break;
      case xdmf_attribute_type::LONG_LONG:
        HDF5Datatype = H5T_NATIVE_LLONG;
        break;
      case xdmf_attribute_type::DOUBLE:
        HDF5Datatype = H5T_NATIVE_DOUBLE;
        break;
      default:
        OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
        break;
      }
      ovk::core::MakeHandle(H5Dcreate(AttributesGroup, Attribute.Name().c_str(), HDF5Datatype,
        Dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    }
  }

  return {std::move(Path), NumDims, std::move(Comm), std::move(Grids), std::move(Attributes),
    std::move(HDF5File)};

}

xdmf CreateXDMF(std::string Path, int NumDims, ovk::comm_view Comm, ovk::array<xdmf_grid_meta>
  Grids, ovk::array<xdmf_attribute_meta> Attributes) {
  return xdmf::internal_Create(std::move(Path), NumDims, Comm, std::move(Grids),
    std::move(Attributes));
}

ovk::optional<xdmf> CreateXDMF(std::string Path, int NumDims, ovk::comm_view Comm,
  ovk::array<xdmf_grid_meta> Grids, ovk::array<xdmf_attribute_meta> Attributes, captured_xdmf_error
  &Error) {

  Error.Reset();

  ovk::optional<xdmf> MaybeXDMF;

  try {
    MaybeXDMF = xdmf::internal_Create(std::move(Path), NumDims, Comm, std::move(Grids),
      std::move(Attributes));
  } catch (const xdmf_error &CreateError) {
    Error = CreateError.Capture();
  }

  return MaybeXDMF;

}

xdmf xdmf::internal_Open(std::string &&Path, ovk::comm_view Comm_) {

  ovk::comm Comm = ovk::DuplicateComm(Comm_);

  std::size_t iNameStart = Path.find_last_of("/");
  if (iNameStart != std::string::npos) {
    iNameStart += 1;
  } else {
    iNameStart = 0;
  }
  std::size_t iNameEnd = Path.find_last_of(".");
  std::string Directory = Path.substr(0,iNameStart);
  std::string Name = Path.substr(iNameStart,iNameEnd-iNameStart);

  ovk::tuple<std::string> CoordNames = {"X", "Y", "Z"};

  int NumDims;
  int NumGrids;
  int NumAttributes;

  ovk::array<xdmf_grid_meta> Grids;
  ovk::array<xdmf_attribute_meta> Attributes;

  if (Comm.Rank() == 0) {

    pugi::xml_document Document;

    Document.load_file(Path.c_str());

    pugi::xml_node RootNode = Document.child("Xdmf");
    pugi::xml_node DomainNode = RootNode.child("Domain");

    {
      pugi::xml_node GridNode = DomainNode.child("Grid");
      pugi::xml_node TopologyNode = GridNode.child("Topology");
      std::string SizeString = TopologyNode.attribute("Dimensions").value();
      std::stringstream SizeStringStream(SizeString);
      int DimSize;
      NumDims = 0;
      while (SizeStringStream >> DimSize) {
        ++NumDims;
      }
    }

    NumGrids = 0;
    for (pugi::xml_node GridNode = DomainNode.child("Grid"); GridNode; GridNode =
      GridNode.next_sibling("Grid")) {
      ++NumGrids;
    }

    Grids.Reserve(NumGrids);

    for (pugi::xml_node GridNode = DomainNode.child("Grid"); GridNode; GridNode =
      GridNode.next_sibling("Grid")) {

      std::string GridName = GridNode.attribute("Name").value();

      pugi::xml_node TopologyNode = GridNode.child("Topology");
      std::string SizeString = TopologyNode.attribute("Dimensions").value();
      std::stringstream SizeStringStream(SizeString);
      ovk::tuple<int> Size = ovk::MakeUniformTuple<int>(NumDims, 0, 1);
      int iDim = NumDims-1;
      while (SizeStringStream >> Size(iDim)) {
        --iDim;
      }

      Grids.Append(std::move(GridName), Size);

    }

    {

      pugi::xml_node GridNode = DomainNode.child("Grid");

      NumAttributes = 0;
      for (pugi::xml_node AttributeNode = GridNode.child("Attribute"); AttributeNode;
        AttributeNode = AttributeNode.next_sibling("Attribute")) {
        ++NumAttributes;
      }

      Attributes.Reserve(NumAttributes);

      for (pugi::xml_node AttributeNode = GridNode.child("Attribute"); AttributeNode;
        AttributeNode = AttributeNode.next_sibling("Attribute")) {

        std::string AttributeName = AttributeNode.attribute("Name").value();

        pugi::xml_node DataItemNode = AttributeNode.child("DataItem");

        xdmf_attribute_type Type;
        std::string NumberTypeString = DataItemNode.attribute("NumberType").value();
        int Precision = DataItemNode.attribute("NumberType").as_int();
        if (NumberTypeString == "Int") {
          Type = Precision == 8 ? xdmf_attribute_type::LONG_LONG : xdmf_attribute_type::INT;
        } else {
          Type = xdmf_attribute_type::DOUBLE;
        }

        Attributes.Append(std::move(AttributeName), Type);

      }

    }

  }

  MPI_Bcast(&NumDims, 1, MPI_INT, 0, Comm);
  MPI_Bcast(&NumGrids, 1, MPI_INT, 0, Comm);
  MPI_Bcast(&NumAttributes, 1, MPI_INT, 0, Comm);

  if (Comm.Rank() != 0) {
    Grids.Reserve(NumGrids);
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    std::string GridName;
    ovk::tuple<int> Size;
    if (Comm.Rank() == 0) {
      xdmf_grid_meta &Grid = Grids(iGrid);
      GridName = Grid.Name();
      Size = Grid.Size();
    }
    ovk::core::BroadcastString(GridName, 0, Comm);
    MPI_Bcast(Size.Data(), ovk::MAX_DIMS, MPI_INT, 0, Comm);
    if (Comm.Rank() != 0) {
      Grids.Append(std::move(GridName), Size);
    }
  }

  if (Comm.Rank() != 0) {
    Attributes.Reserve(NumAttributes);
  }

  for (int iAttribute = 0; iAttribute < NumAttributes; ++iAttribute) {
    std::string AttributeName;
    xdmf_attribute_type Type;
    if (Comm.Rank() == 0) {
      xdmf_attribute_meta &Attribute = Attributes(iAttribute);
      AttributeName = Attribute.Name();
      Type = Attribute.Type();
    }
    ovk::core::BroadcastString(AttributeName, 0, Comm);
    MPI_Bcast(&Type, 1, ovk::core::GetMPIDataType<xdmf_attribute_type>(), 0, Comm);
    if (Comm.Rank() != 0) {
      Attributes.Append(std::move(AttributeName), Type);
    }
  }

  MPI_Barrier(Comm);

  std::string HDF5Path = Directory + Name + ".h5";

  auto FileProperties = ovk::core::MakeHandle(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
  H5Pset_fapl_mpio(FileProperties, Comm, MPI_INFO_NULL);
  hid_t FileRaw = H5Fopen(HDF5Path.c_str(), H5F_ACC_RDWR, FileProperties);
  if (FileRaw < 0) {
    if (Comm.Rank() == 0) {
      std::fprintf(stderr, "ERROR: Unable to open file '%s'.\n", HDF5Path.c_str());
      std::fflush(stderr);
    }
    throw xdmf_file_open_error(HDF5Path);
  }
  auto HDF5File = ovk::core::MakeHandle(FileRaw, H5Fclose);

  return {std::move(Path), NumDims, std::move(Comm), std::move(Grids), std::move(Attributes),
    std::move(HDF5File)};

}

xdmf OpenXDMF(std::string Path, ovk::comm_view Comm) {
  return xdmf::internal_Open(std::move(Path), Comm);
}

ovk::optional<xdmf> OpenXDMF(std::string Path, ovk::comm_view Comm, captured_xdmf_error &Error) {

  Error.Reset();

  ovk::optional<xdmf> MaybeXDMF;

  try {
    MaybeXDMF = xdmf::internal_Open(std::move(Path), Comm);
  } catch (const xdmf_error &OpenError) {
    Error = OpenError.Capture();
  }

  return MaybeXDMF;

}

namespace {

template <typename T> void WriteHDF5Data(hid_t Dataset, int NumDims, hid_t Datatype,
  ovk::field_view<const T> Data, const ovk::range &WriteRange) {

  ovk::tuple<hsize_t> SourceBufferSize = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    SourceBufferSize(iDim) = Data.Extents().Size(NumDims-1-iDim);
  }

  ovk::tuple<hsize_t> SourceSlabStart = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
  ovk::tuple<hsize_t> SourceSlabSize = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    SourceSlabStart(iDim) = WriteRange.Begin(NumDims-1-iDim) - Data.Extents().Begin(NumDims-1-iDim);
    SourceSlabSize(iDim) = WriteRange.Size(NumDims-1-iDim);
  }

  auto SourceDataspace = ovk::core::MakeHandle(H5Screate_simple(NumDims, SourceBufferSize.Data(),
    nullptr), H5Sclose);

  H5Sselect_hyperslab(SourceDataspace, H5S_SELECT_SET, SourceSlabStart.Data(), nullptr,
    SourceSlabSize.Data(), nullptr);

  auto DestDataspace = ovk::core::MakeHandle(H5Dget_space(Dataset), H5Sclose);

  ovk::tuple<hsize_t> DestSlabStart = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
  ovk::tuple<hsize_t> DestSlabSize = ovk::MakeUniformTuple<hsize_t>(NumDims, 0, 1);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    DestSlabStart(iDim) = WriteRange.Begin(NumDims-1-iDim);
    DestSlabSize(iDim) = WriteRange.Size(NumDims-1-iDim);
  }

  H5Sselect_hyperslab(DestDataspace, H5S_SELECT_SET, DestSlabStart.Data(), nullptr,
    DestSlabSize.Data(), nullptr);

  auto TransferProperties = ovk::core::MakeHandle(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
  H5Pset_dxpl_mpio(TransferProperties, H5FD_MPIO_COLLECTIVE);
  
  H5Dwrite(Dataset, Datatype, SourceDataspace, DestDataspace, TransferProperties, Data.Data());

}

}

xdmf &xdmf::WriteGeometry(const std::string &GridName, int Dimension, ovk::field_view<const double>
  Data) {

  return WriteGeometry(GridName, Dimension, Data, Data.Extents());

  return *this;

}

xdmf &xdmf::WriteGeometry(const std::string &GridName, int Dimension, ovk::field_view<const double>
  Data, const ovk::range &WriteRange) {

  if (OVK_DEBUG) {
    auto Iter = std::find_if(Grids_.Begin(), Grids_.End(), [&GridName](const xdmf_grid_meta &Grid)
      -> bool {
      return Grid.Name() == GridName;
    });
    OVK_DEBUG_ASSERT(Iter != Grids_.End(), "Invalid grid name '%s'.", GridName);
  }

  auto GridGroup = ovk::core::MakeHandle(H5Gopen(HDF5File_, GridName.c_str(), H5P_DEFAULT),
    H5Gclose);
  auto GeometryGroup = ovk::core::MakeHandle(H5Gopen(GridGroup, "Geometry", H5P_DEFAULT), H5Gclose);

  ovk::tuple<std::string> CoordNames = {"X", "Y", "Z"};

  auto Dataset = ovk::core::MakeHandle(H5Dopen(GeometryGroup, CoordNames(Dimension).c_str(),
    H5P_DEFAULT), H5Dclose);

  WriteHDF5Data(Dataset, NumDims_, H5T_NATIVE_DOUBLE, Data, WriteRange);

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const int> Data) {

  return WriteAttribute(GridName, AttributeName, Data, Data.Extents());

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const int> Data, const ovk::range &WriteRange) {

  if (OVK_DEBUG) {
    auto GridsIter = std::find_if(Grids_.Begin(), Grids_.End(), [&GridName](const xdmf_grid_meta
      &Grid) -> bool {
      return Grid.Name() == GridName;
    });
    OVK_DEBUG_ASSERT(GridsIter != Grids_.End(), "Invalid grid name '%s'.", GridName);
    auto AttributesIter = std::find_if(Attributes_.Begin(), Attributes_.End(), [&AttributeName](
      const xdmf_attribute_meta &Attribute) -> bool {
      return Attribute.Name() == AttributeName;
    });
    OVK_DEBUG_ASSERT(AttributesIter != Attributes_.End(), "Invalid attribute name '%s'.",
      AttributeName);
  }

  auto GridGroup = ovk::core::MakeHandle(H5Gopen(HDF5File_, GridName.c_str(), H5P_DEFAULT),
    H5Gclose);
  auto AttributesGroup = ovk::core::MakeHandle(H5Gopen(GridGroup, "Attributes", H5P_DEFAULT),
    H5Gclose);

  auto Dataset = ovk::core::MakeHandle(H5Dopen(AttributesGroup, AttributeName.c_str(), H5P_DEFAULT),
    H5Dclose);

  WriteHDF5Data(Dataset, NumDims_, H5T_NATIVE_INT, Data, WriteRange);

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const long long> Data) {

  return WriteAttribute(GridName, AttributeName, Data, Data.Extents());

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const long long> Data, const ovk::range &WriteRange) {

  if (OVK_DEBUG) {
    auto GridsIter = std::find_if(Grids_.Begin(), Grids_.End(), [&GridName](const xdmf_grid_meta
      &Grid) -> bool {
      return Grid.Name() == GridName;
    });
    OVK_DEBUG_ASSERT(GridsIter != Grids_.End(), "Invalid grid name '%s'.", GridName);
    auto AttributesIter = std::find_if(Attributes_.Begin(), Attributes_.End(), [&AttributeName](
      const xdmf_attribute_meta &Attribute) -> bool {
      return Attribute.Name() == AttributeName;
    });
    OVK_DEBUG_ASSERT(AttributesIter != Attributes_.End(), "Invalid attribute name '%s'.",
      AttributeName);
  }

  auto GridGroup = ovk::core::MakeHandle(H5Gopen(HDF5File_, GridName.c_str(), H5P_DEFAULT),
    H5Gclose);
  auto AttributesGroup = ovk::core::MakeHandle(H5Gopen(GridGroup, "Attributes", H5P_DEFAULT),
    H5Gclose);

  auto Dataset = ovk::core::MakeHandle(H5Dopen(AttributesGroup, AttributeName.c_str(), H5P_DEFAULT),
    H5Dclose);

  WriteHDF5Data(Dataset, NumDims_, H5T_NATIVE_LLONG, Data, WriteRange);

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const double> Data) {

  return WriteAttribute(GridName, AttributeName, Data, Data.Extents());

  return *this;

}

xdmf &xdmf::WriteAttribute(const std::string &GridName, const std::string &AttributeName,
  ovk::field_view<const double> Data, const ovk::range &WriteRange) {

  if (OVK_DEBUG) {
    auto GridsIter = std::find_if(Grids_.Begin(), Grids_.End(), [&GridName](const xdmf_grid_meta
      &Grid) -> bool {
      return Grid.Name() == GridName;
    });
    OVK_DEBUG_ASSERT(GridsIter != Grids_.End(), "Invalid grid name '%s'.", GridName);
    auto AttributesIter = std::find_if(Attributes_.Begin(), Attributes_.End(), [&AttributeName](
      const xdmf_attribute_meta &Attribute) -> bool {
      return Attribute.Name() == AttributeName;
    });
    OVK_DEBUG_ASSERT(AttributesIter != Attributes_.End(), "Invalid attribute name '%s'.",
      AttributeName);
  }

  auto GridGroup = ovk::core::MakeHandle(H5Gopen(HDF5File_, GridName.c_str(), H5P_DEFAULT),
    H5Gclose);
  auto AttributesGroup = ovk::core::MakeHandle(H5Gopen(GridGroup, "Attributes", H5P_DEFAULT),
    H5Gclose);

  auto Dataset = ovk::core::MakeHandle(H5Dopen(AttributesGroup, AttributeName.c_str(), H5P_DEFAULT),
    H5Dclose);

  WriteHDF5Data(Dataset, NumDims_, H5T_NATIVE_DOUBLE, Data, WriteRange);

  return *this;

}

}
