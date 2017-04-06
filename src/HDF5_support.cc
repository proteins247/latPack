#include "HDF5_support.hh"


HDF5TrajWriter::HDF5TrajWriter(const char* filepath, OUT_MODE sim_out_mode_, size_t max_structure_length_, size_t chunk_size_)
        : chunk_size(chunk_size_), // default value is over 9000 (16384)
          sim_out_mode(sim_out_mode_),
          max_structure_length(max_structure_length_) // default value is 0
{
        file_id = H5Fcreate(filepath, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
        if (is_invalid())
                throw File_exists_error{};
        if (max_structure_length > 0) {
                str_memtype = H5Tcopy(H5T_C_S1);
                // str_filetype = H5Tcopy(H5T_FORTRAN_S1);
                str_filetype = H5Tcopy(H5T_C_S1);
                status = H5Tset_size(str_memtype, max_structure_length+1);
                status = H5Tset_size(str_filetype, max_structure_length+1);
                // allocate 2D array for structures
                structures = new char[chunk_size * (max_structure_length+1)];
                structure_count = 0;
                while (max_structure_length * chunk_size > 1000000)
                        chunk_size /= 2;
        }
        group_is_open = false;
}

HDF5TrajWriter::~HDF5TrajWriter()
{
        close_file();
}


bool
HDF5TrajWriter::is_invalid()
{
        if (file_id < 0)
                return true;
        return false;
}

void
HDF5TrajWriter::create_trajectory_group()
{
        // can't open a trajectory group if another is already open
        if (group_is_open)
                throw Group_is_open_error{};

        int traj_num = group_ids.size() + 1;
        std::string groupname = "traj" + std::to_string(traj_num);
        hid_t group = H5Gcreate(file_id, groupname.c_str(),
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        group_ids.push_back(group);
        
        // create dataspace with unlimited dimensions
        dims[0] = {0};          // (re)set dims of dataspace to 0
        hid_t dataspace = H5Screate_simple(1, dims, maxdims);

        // enable chunking, chunksize set by chunk_dims[0]
        hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_chunk(prop, 1, chunk_dims);

        // create datasets
        dataset_steps = H5Dcreate2(group, "steps", H5T_NATIVE_UINT, dataspace,
                                   H5P_DEFAULT, prop, H5P_DEFAULT);
        dataset_energies = H5Dcreate2(group, "energies", H5T_NATIVE_FLOAT, dataspace,
                                      H5P_DEFAULT, prop, H5P_DEFAULT);
        if (sim_out_mode == OUT_ES) {
                dataset_structures = H5Dcreate2(group, "structures", str_filetype, dataspace,
                                                H5P_DEFAULT, prop, H5P_DEFAULT);
        }
        // close resources
        H5Pclose(prop);
        H5Sclose(dataspace);

        group_is_open = true;

        // reset offset
        offset[0] = {0};
}

void
HDF5TrajWriter::write_attribute(const char* name, unsigned int attr)
{
        hsize_t attr_dims[1] = {1};
        // create dataspace
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // create attribute
        hid_t attribute = H5Acreate2(file_id, name, H5T_NATIVE_UINT, dataspace, 
                                     H5P_DEFAULT, H5P_DEFAULT);

        // write attribute
        status = H5Awrite(attribute, H5T_NATIVE_UINT, &attr);

        // close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
}

void
HDF5TrajWriter::write_attribute(const char* name, float attr)
{
        hsize_t attr_dims[1] = {1};
        // create dataspace
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // create attribute
        hid_t attribute = H5Acreate2(file_id, name, H5T_NATIVE_FLOAT, dataspace, 
                                     H5P_DEFAULT, H5P_DEFAULT);

        // write attribute
        status = H5Awrite(attribute, H5T_NATIVE_FLOAT, &attr);

        // close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
}

void
HDF5TrajWriter::write_attribute(const char* name, std::string attr)
{
        hsize_t attr_dims[1] = {1};
        size_t strlen = attr.size();

        // copy the string type
        hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
        status = H5Tset_size(filetype, strlen);
        hid_t memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, strlen+1);
        
        // create dataspace
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // create attribute
        hid_t attribute = H5Acreate2(file_id, name, filetype, dataspace, 
                                     H5P_DEFAULT, H5P_DEFAULT);

        // write attribute
        status = H5Awrite(attribute, memtype, attr.c_str());

        // close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
        status = H5Tclose(memtype);
        status = H5Tclose(filetype);
        
}

void
HDF5TrajWriter::write_buffered_traj(unsigned int step, float energy, std::string *structure)
{
        steps.push_back(step);
        energies.push_back(energy);
        if (structure) { // if not given as nullptr 
                structure->copy(&structures[structure_count*(max_structure_length+1)],
                                max_structure_length );
                structures[structure_count*(max_structure_length+1) + structure->size()] = '\0';
                structure_count++;
        }
        if (steps.size() == chunk_size)
                // time to write
                flush_buffers();
}

void
HDF5TrajWriter::close_trajectory_group()
{
        // perform final writes to file if the buffers still contain stuff
        if (steps.size() > 0)
                flush_buffers();

        status = H5Dclose(dataset_steps);
        status = H5Dclose(dataset_energies);
        if (sim_out_mode == OUT_ES)
                status = H5Dclose(dataset_structures);

        status = H5Gclose(group_ids.back());
        group_is_open = false;
}

void
HDF5TrajWriter::delete_last_group()
{
        unsigned int group_count = group_ids.size();
        std::string group_name = "traj" + std::to_string(group_count);
        // By removing the link from root to the group, which is the only link that exists,
        //   the group should be deleted.
        status = H5Ldelete(file_id, group_name.c_str(), H5P_DEFAULT);
        group_ids.pop_back();
}

void
HDF5TrajWriter::flush_buffers()
{
        hsize_t extend_dims[1] = {steps.size()};

        // grow the total size of dataspace
        dims[0] += steps.size(); 
        status = H5Dset_extent(dataset_steps, dims);
        status = H5Dset_extent(dataset_energies, dims);
        if (sim_out_mode == OUT_ES)
                status = H5Dset_extent(dataset_structures, dims);

        // get current filespace
        hid_t filespace = H5Dget_space(dataset_steps);

        // select hyperslab for each dataset
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     extend_dims, NULL);
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     extend_dims, NULL);
        if (sim_out_mode == OUT_ES)
                status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                             extend_dims, NULL);

        // define memory space
        hid_t memspace = H5Screate_simple(1, extend_dims, NULL); // RANK=1

        // write!
        status = H5Dwrite(dataset_steps, H5T_NATIVE_UINT, memspace, filespace,
                          H5P_DEFAULT, steps.data());
        status = H5Dwrite(dataset_energies, H5T_NATIVE_FLOAT, memspace, filespace,
                          H5P_DEFAULT, energies.data());
        if (sim_out_mode == OUT_ES) {
                status = H5Dwrite(dataset_structures, str_memtype, memspace, filespace,
                                  H5P_DEFAULT, structures);
        }

        // close resources
        status = H5Sclose(memspace);
        status = H5Sclose(filespace);
                
        // after writing, clear the containers:
        steps.clear();
        energies.clear();
        structure_count = 0;    // just reset the count

        // adjust offset
        offset[0] += extend_dims[0];
}

void
HDF5TrajWriter::close_file()
{
        if (group_is_open)
                close_trajectory_group();
        if (max_structure_length > 0) {
                status = H5Tclose(str_memtype);
                status = H5Tclose(str_filetype);
                delete[] structures;
        }
        status = H5Fclose(file_id);
}
