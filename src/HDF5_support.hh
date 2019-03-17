#ifndef HDF5_SUPPORT_HH
#define HDF5_SUPPORT_HH

//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include "hdf5.h"

// If debugging: comment out #include "SC_MinE.hh" and uncomment enum
#include "SC_MinE.hh"           // defines enum OUT_MODE
// enum OUT_MODE { OUT_NO, OUT_E, OUT_ES };

/*
 * This file is designed so that the other files don't need to know
 * the hdf5 types. It is a wrapper around the hdf5 C code.
 * 
 * We create hdf5 files with a specific format.
 * - Arbitrary attributes can only be written for the entire system
 *   at the root level of the hierarchy.
 * - We write trajectory groups one at a time, and they have a fixed structure.
 * 
 * In a sense, we are replicating the text output of latFold and latFoldVec,
 * but with writes to an HDF5 file. Hence things proceed in sequence:
 * - Simulation attributes to root group
 * - Trajectory data to trajectory groups, one by one
 * 
 * Disadvantage of this setup is that there's no flexibility.
 * 
 * Code contains two classes:
 * 1. HDF5TrajWriter   (for latFold, latFoldVec)
 * 2. HDF5TrajAnalyzer (for latMapTraj)
 * 
 * My first C++ classes, so they're ugly.
 *
 */

class HDF5TrajWriter {
public:
        // C'tor:
        //   Provide max_structure_length > 0 for sim_out_mode == OUT_ES
        //   chunk_size default should allow data to fit in default 1 MB HDF5 cache
        //   but it is adjusted automatically.
        HDF5TrajWriter(const char* filepath,
                       OUT_MODE sim_out_mode,
                       size_t max_structure_length=0,
                       size_t chunk_size=16384);

        // Destructor. Calls protected function close_file()
        ~HDF5TrajWriter();
        
        // true if HDF5TrajWriter and underlying HDF5 objects are in invalid state
        // This function could be further developed.
        bool is_invalid();

        // Create an HDF5 group to store trajectory information.
        //   Sets up datasets within the group.
        void create_trajectory_group();

        // Write attribute of label /name/ to root group of HDF5 file
        //   These could/should be templaticized eventually
        //   But need to know how to handle creating attributes
        //   for different types, esp. string attributes
        void write_attribute(const char* name, unsigned int attr);
        void write_attribute(const char* name, float attr);
        void write_attribute(const char* name, std::string attr);

        // // (template idea, but don't know how to implement)
        // template <typename T>
        // void write_attribute(const char* name, T attr);

        // Take data from a simulation step and add to write
        //   buffers. Writes to HDF5 file once buffer is
        //   full. structure should be nullptr for OUT_E mode. Note,
        //   it's unclear whether I should have set up buffers to
        //   store data or whether HDF5 library's own buffering is
        //   sufficient. But if I wrote data each time this is called,
        //   I'd have to do a hyperslab selection each time.
        void write_buffered_traj(uint64_t step,
                                 float energy,
                                 std::string *structure);

        // Perform final write (flush buffers) and close group
        void close_trajectory_group(bool successfulRunMinE = false,
                                    bool foundFinalStructure = false,
                                    double targetFraction = 0.0,
                                    double targetEnergy = 0.0,
				    size_t stepsToTarget = 0,
				    double survivalSumFraction = 0.0);

        // Delete the last used HDF5 group. The memory/space is freed
        //   so long as the file has not been closed since the group's
        //   creation.
        void delete_last_group();
        

protected:
        herr_t status;

        // Identifier for HDF5 file
        hid_t file_id;

        // Output mode for HDF5TrajWriter
        const OUT_MODE sim_out_mode;

        // For sizing the string arrays in hdf5 file
        size_t max_structure_length; 

        // File data will be stored in chunks of this value
        // Also the size to extend datasets each time
        // Decision details:
        // - default cache size is 1 MB
        // - need a chunk size less than that
        // - imagine a 100-residue structure (extreme) = 100 bytes (chars)
        // - 10000 of 100 bytes would be 1 MB
        unsigned int chunk_size;


        // Container for HDF5 groups, each of which corresponds to a trajectory
        std::vector<hid_t> group_ids;

        // Var to track whether last group in /group_ids/ is open
        bool group_is_open;

        // Each group will contain these datasets (no /dataset_structures/ if not OUT_ES)
        hid_t dataset_steps, dataset_energies, dataset_structures;

        // Stores the current step and current energy saved
        uint64_t current_step;
        float current_energy;

        // Memory and file datatypes for writing structure strings
        hid_t str_memtype, str_filetype;
        
        // Size of array held in dataset
        //   Value is (re)set in createTrajectoryGroup()
        hsize_t dims[1];
        hsize_t chunk_dims[1] = {chunk_size};         // Our chunk size
        hsize_t offset[1] = {0};              // Var to store offset to do hyperslab selections
        hsize_t maxdims[1] = {H5S_UNLIMITED}; // maximum extensible size

        // Buffers to hold data in before writing
        std::vector<uint64_t> steps;
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
 * This second class is for opening already written trajectory h5
 * files.
 *
 * Facilities are provided for opening trajectory groups, reading
 * trajectory data, and writing additional datasets.
 *
 * This analyzer class could be improved to handle hdf5 files which
 * the analyzer has previously examined
 *
 */
class HDF5TrajAnalyzer {
public:
        // filepath - name of file
        // dataset_names - pointer to vector of names of new datasets to create in each trajectory group
        // chunk_size - chunking for data in new datasets
        HDF5TrajAnalyzer(const char* filepath, std::vector<std::string> * dataset_names, size_t chunk_size = 16384);
        ~HDF5TrajAnalyzer();

        // true if HDF5TrajAnalyzer and underlying HDF5 objects are in invalid state
        // This function could be further developed.
        bool is_invalid();

        size_t get_group_count();

        // Open an HDF5 group corresponding to index. First traj is index = 1 (traj1).
        void open_trajectory_group(size_t index);

        // Close currently open trajectory group. (only one can be open at a time)
        void close_trajectory_group();

        // Read the trajectory
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

class File_opening_error {
};

// class File_exists_error {
// };

class Group_is_open_error {
};

class Group_creation_error {
};

class Group_opening_error {
};

class Data_write_error {
};

#endif // HDF5_SUPPORT_HH
