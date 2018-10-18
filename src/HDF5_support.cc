#include "HDF5_support.hh"

HDF5TrajWriter::HDF5TrajWriter(const char* filepath,
                               OUT_MODE sim_out_mode_,
                               size_t max_structure_length_,
                               size_t chunk_size_)
        : sim_out_mode(sim_out_mode_),
          max_structure_length(max_structure_length_), // default value is 0
          chunk_size(chunk_size_) // default value is over 9000 (16384)
{
        file_id = H5Fcreate(filepath, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

        if (is_invalid())
                throw File_exists_error{};

        if (max_structure_length > 0)
        {
                str_memtype = H5Tcopy(H5T_C_S1);
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
        {
                return true;
        }
        return false;
}

void
HDF5TrajWriter::create_trajectory_group()
{
        // Can't open a trajectory group if another is already open
        if (group_is_open)
        {
                throw Group_is_open_error{};
        }

        std::string groupname = "traj" + std::to_string(group_ids.size() + 1);
        hid_t group = H5Gcreate(file_id, groupname.c_str(),
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        group_ids.push_back(group);
        
        // (Re)set dims of dataspace to 0
        dims[0] = {0};          
        // Reset offset
        offset[0] = {0};

        // Create dataspace with unlimited dimensions:
        hid_t dataspace = H5Screate_simple(1, dims, maxdims);

        // Enable chunking, chunksize set by chunk_dims[0]
        hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_chunk(prop, 1, chunk_dims);

        // Enable compression level 6
        status = H5Pset_deflate (prop, 6);         

        // Create datasets
        dataset_steps = H5Dcreate2(group, "steps",
                                   H5T_NATIVE_UINT, dataspace,
                                   H5P_DEFAULT, prop, H5P_DEFAULT);
        dataset_energies = H5Dcreate2(group, "energies",
                                      H5T_NATIVE_FLOAT, dataspace,
                                      H5P_DEFAULT, prop, H5P_DEFAULT);
        if (sim_out_mode == OUT_ES)
        {
                dataset_structures = H5Dcreate2(group, "structures",
                                                str_filetype, dataspace,
                                                H5P_DEFAULT, prop, H5P_DEFAULT);
        }

        // Close resources
        H5Pclose(prop);
        H5Sclose(dataspace);

        group_is_open = true;
}

void
HDF5TrajWriter::write_attribute(const char* name, unsigned int attr)
{
        hsize_t attr_dims[1] = {1};
        // Create dataspace
        // Note that I should have used H5SCreate(H5S_SCALAR)
        //   to store a single value. This creates a 1-sized array.
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // Create attribute
        hid_t attribute = H5Acreate2(file_id, name, H5T_NATIVE_UINT, dataspace,
                                     H5P_DEFAULT, H5P_DEFAULT);

        // Write attribute
        status = H5Awrite(attribute, H5T_NATIVE_UINT, &attr);

        // Close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
}

void
HDF5TrajWriter::write_attribute(const char* name, float attr)
{
        hsize_t attr_dims[1] = {1};
        // Create dataspace
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // Create attribute
        hid_t attribute = H5Acreate2(file_id, name, H5T_NATIVE_FLOAT, dataspace,
                                     H5P_DEFAULT, H5P_DEFAULT);

        // Write attribute
        status = H5Awrite(attribute, H5T_NATIVE_FLOAT, &attr);

        // Close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
}

void
HDF5TrajWriter::write_attribute(const char* name, std::string attr)
{
        hsize_t attr_dims[1] = {1};
        size_t strlen = attr.size();

        // Copy the string type
        hid_t strtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(strtype, strlen+1);
        // hid_t memtype = H5Tcopy(H5T_C_S1);
        // status = H5Tset_size(memtype, strlen+1);
        
        // Create dataspace
        hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

        // Create attribute
        hid_t attribute = H5Acreate2(file_id, name, strtype, dataspace, 
                                     H5P_DEFAULT, H5P_DEFAULT);

        // Write attribute
        status = H5Awrite(attribute, strtype, attr.c_str());

        // Close resources
        status = H5Aclose(attribute);
        status = H5Sclose(dataspace);
        // status = H5Tclose(memtype);
        status = H5Tclose(strtype);
        
}

void
HDF5TrajWriter::write_buffered_traj(unsigned int step,
                                    float energy,
                                    std::string *structure)
{
        steps.push_back(step);
        energies.push_back(energy);
        current_step = step;
        current_energy = energy;
        if (structure)   // If not given as nullptr
        {
                structure->copy(
                        &structures[structure_count*(max_structure_length+1)],
                        max_structure_length );
                structures[structure_count*(max_structure_length+1) + structure->size()] = '\0';
                structure_count++;
        }
        if (steps.size() == chunk_size)
                // time to write
                flush_buffers();
}

void
HDF5TrajWriter::close_trajectory_group(bool successfulRunMinE,
                                       bool foundFinalStructure)
        // Default values of false for both function args.
{
        // Perform final writes to file if the buffers still contain stuff
        if (steps.size() > 0)
                flush_buffers();

        status = H5Dclose(dataset_steps);
        status = H5Dclose(dataset_energies);
        if (sim_out_mode == OUT_ES)
                status = H5Dclose(dataset_structures);
        
        // Write attributes on sim success
        {
                hid_t strtype = H5Tcopy(H5T_C_S1);
                status = H5Tset_size(strtype, 4);
                // hid_t memtype = H5Tcopy(H5T_C_S1);
                // status = H5Tset_size(memtype, 4);
        
                // create dataspace
                hsize_t attr_dims[1] = {1};
                hid_t dataspace = H5Screate_simple(1, attr_dims, NULL);

                // create attribute
                hid_t attribute1 = H5Acreate2(group_ids.back(), "Reached min E",
                                              strtype, dataspace, 
                                              H5P_DEFAULT, H5P_DEFAULT);
                hid_t attribute2 = H5Acreate2(group_ids.back(), "Found final struct",
                                              strtype, dataspace, 
                                              H5P_DEFAULT, H5P_DEFAULT);
                hid_t attribute3 = H5Acreate2(group_ids.back(), "Last E",
                                              H5T_NATIVE_FLOAT, dataspace,
                                              H5P_DEFAULT, H5P_DEFAULT);
                hid_t attribute4 = H5Acreate2(group_ids.back(), "Last step",
                                              H5T_NATIVE_UINT, dataspace,
                                              H5P_DEFAULT, H5P_DEFAULT);

                // write attributes
                if (successfulRunMinE)
                        status = H5Awrite(attribute1, strtype, "yes");
                else
                        status = H5Awrite(attribute1, strtype, "no ");
                if (foundFinalStructure)
                        status = H5Awrite(attribute2, strtype, "yes");
                else
                        status = H5Awrite(attribute2, strtype, "no ");

                status = H5Awrite(attribute3, H5T_NATIVE_FLOAT, &current_energy);
                status = H5Awrite(attribute4, H5T_NATIVE_UINT, &current_step);

                // close resources
                status = H5Aclose(attribute1);
                status = H5Aclose(attribute2);
                status = H5Sclose(dataspace);
                // status = H5Tclose(memtype);
                status = H5Tclose(strtype);
        } // end sim success attributes writing

        status = H5Gclose(group_ids.back());
        group_is_open = false; 
}

void
HDF5TrajWriter::delete_last_group()
{
        std::string group_name = "traj" + std::to_string(group_ids.size());
        // By removing the link from root to the group, which is the
        //   only link that exists, the group should be deleted.
        status = H5Ldelete(file_id, group_name.c_str(), H5P_DEFAULT);
        group_ids.pop_back();
}

//////////////////////////////////////////////////
// protected:
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


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


HDF5TrajAnalyzer::HDF5TrajAnalyzer(const char* filepath, std::vector<std::string> * dataset_names_,
                                   size_t chunk_size_)
        : dataset_names(dataset_names_), chunk_size(chunk_size_)
{
        // open file in read/write mode
        file_id = H5Fopen(filepath, H5F_ACC_RDWR, H5P_DEFAULT);
        
        // read relevant attributes
        H5G_info_t group_info;
        status = H5Gget_info(file_id, &group_info);
        group_count = group_info.nlinks;

}

HDF5TrajAnalyzer::~HDF5TrajAnalyzer()
{
        close_file();
}

size_t
HDF5TrajAnalyzer::get_group_count()
{
        return group_count;
}

void
HDF5TrajAnalyzer::open_trajectory_group(size_t index)
{
        if (group_is_open)
                throw Group_is_open_error{};

        // open the group
        std::string groupname = "traj" + std::to_string(index);
        group_id = H5Gopen2(file_id, groupname.c_str(),  H5P_DEFAULT);
        group_is_open = true;
        
        // open the structures dataspace
        dataset_structures = H5Dopen(group_id, "structures", H5P_DEFAULT);
        filespace = H5Dget_space(dataset_structures);
        H5Sget_simple_extent_dims(filespace, trajectory_size, NULL);

        // for reading and writing. simple memspace of size 1
        memspace = H5Screate_simple(1, read_write_dims, NULL); 

        // create memtype
        hid_t filetype = H5Dget_type(dataset_structures);
        size_t sdim = H5Tget_size(filetype);
        memtype = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (memtype, sdim);

        // allocation of space to read data
        structure = new char[sdim];

        // chunking for new datasets
        hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_chunk(prop, 1, chunk_dims);

        // enable compression level 6
        status = H5Pset_deflate (prop, 6);         

        // creation of datasets
        for (auto &name : * dataset_names) {
                hid_t dataset = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_FLOAT, filespace,
                                           H5P_DEFAULT, prop, H5P_DEFAULT);
                datasets.insert({name, dataset});
        }

        // reset offsets
        read_offset[0] = 0;
        write_offset[0] = 0;

        // close those resources
        H5Tclose(filetype);
        H5Pclose(prop);
}

void
HDF5TrajAnalyzer::close_trajectory_group()
{
        // close the datasets
        for (auto it=datasets.begin(); it!=datasets.end(); ++it)
                H5Dclose(it->second);

        delete[] structure;
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Tclose(memtype);
        H5Dclose(dataset_structures);
        H5Gclose(group_id);
        group_is_open = false;
}

ssize_t
HDF5TrajAnalyzer::read_structure_traj(std::string *structure_str)
{
        if (read_offset[0] == trajectory_size[0])
                return -1;

        // setup reading filespace
        // filespace is the dataspace of dataset_structures

        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, read_offset, NULL,
                                     read_write_dims, NULL);
        // perform read operation
        status = H5Dread(dataset_structures, memtype, memspace, filespace, H5P_DEFAULT, structure);
        *structure_str = structure;

        // return the index of current read (increase offset for next time)
        return read_offset[0]++;
}

ssize_t
HDF5TrajAnalyzer::write_analysis(std::unordered_map<std::string, float> &data)
{
        if (write_offset[0] == trajectory_size[0])
                return -1;

        // select filespace for writing
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, write_offset, NULL,
                                     read_write_dims, NULL);

        // write to datasets
        for (auto &dataset : datasets) {
                H5Dwrite(dataset.second, H5T_NATIVE_FLOAT, memspace, filespace,
                         H5P_DEFAULT, &data[dataset.first]);
        }
        
        // return writing position, incrementing it for next write as well
        return write_offset[0]++;
}

// protected:
void
HDF5TrajAnalyzer::close_file()
{
        if (group_is_open)
                close_trajectory_group();
        H5Fclose(file_id);
}
