#ifndef HDF5_SUPPORT_HH
#define HDF5_SUPPORT_HH

//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <unordered_map>
#include "hdf5.h"

// for debugging: (comment SC_MinE.hh and uncomment enum)
#include "SC_MinE.hh"           // defines enum OUT_MODE
// enum OUT_MODE { OUT_NO, OUT_E, OUT_ES };

/*
This file should be designed so that the other files don't need
to know the hdf5 types.

We end up setting up for writing hdf5 files in a specific format
- attributes can only be written for the entire system
at the root level of the hierarchy
- we write trajectory groups one at a time, and they have a fixed structure

in a sense, we are replicating the text output of latFold and latFoldVec,
but with writes to an HDF5 file. hence things proceed in sequence:
- simulation attributes to root group
- trajectory data to trajectory groups, one by one

disadvantage is that there's no flexibility

1. HDF5TrajWriter 
2. HDF5TrajAnalyzer

*/

class HDF5TrajWriter {
public:
        // C'tor. max_length has default value 0 for var max_structure_length
        //   provide max_structure_length > 0 for sim_out_mode == OUT_ES
        //   chunk_size default should allow data to fit in default 1 MB HDF5 cache
        HDF5TrajWriter(const char* filepath, OUT_MODE sim_out_mode,
                       size_t max_structure_length = 0, size_t chunk_size = 16384);

        // Destructor. Calls protected function close_file()
        ~HDF5TrajWriter();
        
        // true if HDF5TrajWriter and underlying HDF5 objects are in invalid state
        bool is_invalid();

        // Create an HDF5 group within the file to store trajectory information.
        //   Sets up datasets within the group.
        void create_trajectory_group();

        // Write attribute of label /name/ to root group of HDF5 file
        //   these could/should be templaticized eventually
        //   but need to know how to handle creating for different types, esp. string attributes
        void write_attribute(const char* name, unsigned int attr);
        void write_attribute(const char* name, float attr);
        void write_attribute(const char* name, std::string attr);

        // template <typename T>
        // void write_attribute(int groupIndex, const char* name, T attr);

        // Take data and add to write buffers. Write to HDF5 file if buffer is full
        //   structure should be nullptr for OUT_E mode
        // Note, it's unclear whether I should have set up buffers to store data or whether
        //   HDF5 library's own buffering is sufficient. But if I wrote data each
        //   time this is called, I'd have to do a hyperslab selection each time
        void write_buffered_traj(unsigned int step, float energy, std::string *structure);

        // Perform final write (flush buffers) and close group
        void close_trajectory_group();

        // Delete the last used HDF5 group
        // memory/space is freed so long as the file has not been closed
        void delete_last_group();
        

protected:
        herr_t status;

        // File data will be stored in chunks of this value
        // Also the size to extend datasets each time
        // Decision details:
        // - default cache size is 1 MB
        // - need a chunk size less than that
        // - imagine a 100-residue structure (extreme) = 100 bytes (chars)
        // - 10000 of 100 bytes would be 1 MB
        unsigned int chunk_size;

        // Identifier for HDF5 file
        hid_t file_id;
        // Container for HDF5 groups, each of which corresponds to a trajectory
        std::vector<hid_t> group_ids;
        // Var to track whether last group in /group_ids/ is open
        bool group_is_open;

        // Output mode for HDF5TrajWriter
        const OUT_MODE sim_out_mode;

        // For sizing the string arrays in hdf5 file
        size_t max_structure_length; 

        // Each group will contain these datasets (no /dataset_structures/ if not OUT_ES)
        hid_t dataset_steps, dataset_energies, dataset_structures;
        // Memory and file datatypes for writing structure strings
        hid_t str_memtype, str_filetype;
        
        // Size of array held in dataset
        //   Value is (re)set in createTrajectoryGroup()
        hsize_t dims[1];
                        
        hsize_t chunk_dims[1] = {chunk_size};         // Our chunk size
        hsize_t offset[1] = {0};              // Var to store offset to do hyperslab selections
        hsize_t maxdims[1] = {H5S_UNLIMITED}; // maximum extensible size

        // Buffers to hold data in before writing
        std::vector<unsigned int> steps;
        std::vector<float> energies;
        char* structures;       // we need to manually manage char* to get
                                //   a contiguous char array that we write with
        // A count of how many structures stored
        size_t structure_count;

protected:
        // Write data stored in buffers to file
        void flush_buffers();

        // Write any buffered data (if any) and close group (if open) and file
        //   (called by d'tor)
        void close_file();      
};


/*
This second class is for opening already written trajectory h5 files.

facilities are provided for opening trajectory groups, reading trajectory data,
and writing additional datasets.

 */
class HDF5TrajAnalyzer {
public:
        HDF5TrajAnalyzer(const char* filepath, std::vector<std::string> * dataset_names, size_t chunk_size = 16384);
        ~HDF5TrajAnalyzer();

        size_t get_group_count();

        void open_trajectory_group(size_t index);

        void close_trajectory_group();

        ssize_t read_structure_traj(std::string *structure);

        ssize_t write_analysis(std::unordered_map<std::string, float> &data);

protected:
        hid_t file_id;
        herr_t status;

        // the number of groups in the root of the file
        // (corresponds to the number of trajectories)
        size_t group_count;

        // group associated stuff
        hid_t group_id;
        bool group_is_open;

        // additional datasets to add to the hdf5 file
        std::vector<std::string> * dataset_names;
        std::unordered_map<std::string, hid_t> datasets;

        size_t chunk_size;

        hid_t dataset_structures;  // hid_t to the "structures" dataset in a group
        hid_t memspace, filespace;        // dataspaces
        hid_t memtype;          // hdf5 type used for reading from dataset_structures

        char * structure;       // storage place for read string

        hsize_t trajectory_size[1];
        hsize_t chunk_dims[1] = {chunk_size};         // Our chunk size
        hsize_t read_write_dims[1] = {1};             // we read one at a time
        hsize_t read_offset[1] = {0};              // Var to store offset to do hyperslab selection
        hsize_t write_offset[1] = {0};             // same but for writes

protected:        // functions:
        void close_file();

};

class File_exists_error {
};

class Group_is_open_error {
};

#endif // HDF5_SUPPORT_HH
