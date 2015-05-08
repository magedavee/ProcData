#include "HDF5Reader.hh"
#include "HDF5_StructInfo.hh"

#include <sstream>
#include <fstream>
#include <iostream>

SimIoniReader::SimIoniReader(const string& f_in):
ioni_reader("ScIoni", IoniCluster_offsets, IoniCluster_sizes, 1024),
prim_reader("Prim", ParticleVertex_offsets, ParticleVertex_sizes, 1024),
ncapt_reader("NCapt", NCapt_offsets, NCapt_sizes, 1024) {
    infile_id = H5Fopen(f_in.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(infile_id <= 0) { printf("Unable to open file '%s'!\n", f_in.c_str()); return; }
    assert(infile_id > 0);
    
    herr_t err = H5TBget_table_info(infile_id, "Prim", &nfields, &nrecords );
    assert(err >= 0);
    printf("Reading table 'Prim' with %llu records from '%s'.\n", nrecords, f_in.c_str());
    
    err = H5TBget_table_info(infile_id, "ScIoni", &nfields, &nrecords );
    assert(err >= 0);
    printf("Reading table 'ScIoni' with %llu records from '%s'.\n", nrecords, f_in.c_str());

    ioni_reader.setFile(infile_id);
    prim_reader.setFile(infile_id);
    ncapt_reader.setFile(infile_id);
}

SimIoniReader::~SimIoniReader() { 
    if(infile_id > 0) {
        herr_t err = H5Fclose(infile_id);
        assert(err >= 0);
    }
}

/// comparison function for time-sorting hits
bool compare_hit_times(const s_IoniCluster& a, const s_IoniCluster& b) { return a.t < b.t; }

bool SimIoniReader::loadMergedIoni() {
    if(!ioni_reader.loadEvent()) return false;
    
    map<Int_t, vector<s_IoniCluster> > volClusts;
    for(auto it = ioni_reader.event_read.begin(); it != ioni_reader.event_read.end(); it++) {
        if(!it->E) continue; // ignore 0-energy (neutrino tag) entries
        assert(it->t == it->t); // NaN check
        Int_t V = it->vol;
        if(!volClusts.count(V)) volClusts[V].push_back(*it);
        else {
            if(it->t > volClusts[V].back().t + dt_max)
                volClusts[V].push_back(*it);
            else
                volClusts[V].back() += *it;
        }
    }

    // merge and sort into single detector history
    merged.clear();
    for(auto it = volClusts.begin(); it != volClusts.end(); it++)
        merged.insert(merged.end(), it->second.begin(), it->second.end());
    std::sort(merged.begin(), merged.end(), compare_hit_times);
    
    return true;
}

double SimIoniReader::getAttr(const string& tbl, const string& attr) {
    double val = 0;
    herr_t err = H5LTget_attribute_double(infile_id, tbl.c_str(), attr.c_str(), &val);
    assert(err >= 0);
    return val;
}

/////////////////////////////
/////////////////////////////
/////////////////////////////

template<>
int64_t HDF5_Table_Cache<s_Event>::getIdentifier(const s_Event& val) { return val.N; }
template<>
void  HDF5_Table_Cache<s_Event>::setIdentifier(s_Event& val, int64_t newID) { val.N = newID; }

template<>
int64_t HDF5_Table_Cache<s_IoniCluster>::getIdentifier(const s_IoniCluster& val) { return val.evt; }
template<>
void  HDF5_Table_Cache<s_IoniCluster>::setIdentifier(s_IoniCluster& val, int64_t newID) { val.evt = newID; }

template<>
int64_t HDF5_Table_Cache<s_ParticleVertex>::getIdentifier(const s_ParticleVertex& val) { return val.evt; }
template<>
void  HDF5_Table_Cache<s_ParticleVertex>::setIdentifier(s_ParticleVertex& val, int64_t newID) { val.evt = newID; }

template<>
int64_t HDF5_Table_Cache<s_NCapt>::getIdentifier(const s_NCapt& val) { return val.evt; }
template<>
void  HDF5_Table_Cache<s_NCapt>::setIdentifier(s_NCapt& val, int64_t newID) { val.evt = newID; }



SimEventSelector::SimEventSelector():
evt_xfer("Evt", Event_offsets, Event_sizes, nchunk),
ioni_xfer("ScIoni", IoniCluster_offsets, IoniCluster_sizes, nchunk),
prim_xfer("Prim", ParticleVertex_offsets, ParticleVertex_sizes, nchunk),
ncapt_xfer("NCapt", NCapt_offsets, NCapt_sizes, nchunk) { }

void SimEventSelector::setOutfile(const string& filename) {
    assert(!outfile_id);
    outfile_id = H5Fcreate(filename.c_str(), // file name
                           H5F_ACC_TRUNC, // access_mode : overwrite old file with new data
                           H5P_DEFAULT,   // create_ID
                           H5P_DEFAULT    // access_ID
    );
    
    makeScIoniTable(outfile_id, nchunk, compress);
    makePrimTable(outfile_id, nchunk, compress);
    makeNCaptTable(outfile_id, nchunk, compress);
    makeEvtTable(outfile_id, nchunk, compress);
    
    evt_xfer.tableOut.setFile(outfile_id);
    ioni_xfer.tableOut.setFile(outfile_id);
    prim_xfer.tableOut.setFile(outfile_id);
    ncapt_xfer.tableOut.setFile(outfile_id);
}

void SimEventSelector::setInfile(const string& filename) {
    infile_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(infile_id > 0);
    
    evt_xfer.tableIn.setFile(infile_id);
    ioni_xfer.tableIn.setFile(infile_id);
    prim_xfer.tableIn.setFile(infile_id);
    ncapt_xfer.tableIn.setFile(infile_id);
}

void SimEventSelector::transfer(const vector<int64_t>& ids) {
    printf("Transferring %zu events...", ids.size()); fflush(stdout);
    
    evt_xfer.transferIDs(ids, nTransferred);
    ioni_xfer.transferIDs(ids, nTransferred);
    prim_xfer.transferIDs(ids, nTransferred);
    ncapt_xfer.transferIDs(ids, nTransferred);
    
    nTransferred += ids.size();

    printf(" Done.\n");
}

void SimEventSelector::transfer(const string& evtfile) {
    vector<int64_t> v;
    std::ifstream infl(evtfile.c_str());
    int64_t e = -1;
    while (infl) {
        if(e >= 0) v.push_back(e);
        infl >> e;
    }
    transfer(v);
}

void SimEventSelector::setTotalTime(double t) {
    herr_t err = H5LTset_attribute_double(outfile_id, "Evt", "runtime", &t, 1);
    assert(err >= 0);
    err = H5LTset_attribute_double(outfile_id, "ScIoni", "runtime", &t, 1);
    assert(err >= 0);
    err = H5LTset_attribute_double(outfile_id, "Prim", "runtime", &t, 1);
    assert(err >= 0);
    err = H5LTset_attribute_double(outfile_id, "NCapt", "runtime", &t, 1);
    assert(err >= 0);
}

