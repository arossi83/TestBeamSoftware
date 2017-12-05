#ifndef __Event_h
#define __Event_h
#include <vector>
#include <map>
#include <string>
#include <vector>
#include <stdint.h>
#include "TObject.h"
#include "Cluster.h"
#include "Hit.h"
#include "Cbc.h"
#include "Stub.h"
namespace tbeam {
  class Event :  public TObject {
  public:
    Event();
    virtual ~Event(){};
    void reset();
    //header part
    unsigned int run;//
    unsigned int lumiSection;//
    unsigned int event;//
    unsigned long long time;//
    unsigned long int unixtime;//
    unsigned long int dt;//is it used?
    uint32_t key;//key and value can be dropped
    uint32_t value;//key and value can ber dropped
    int condData;
    int glibStatus;
    
    uint8_t dataFormatVersion;//
    uint8_t debugMode;//
    uint8_t readoutMode;//
    uint8_t dataType;//
    uint64_t glibStatusCode;//
    uint16_t numberOfCBC;//
    //condition data
    unsigned int tdcPhase;//
    unsigned int HVsettings;//
    unsigned int DUTangle;//

    uint32_t window;
    uint32_t offset;
    uint32_t cwd;
    uint32_t tilt;
    uint32_t vcth;
    uint32_t stubLatency;
    uint32_t triggerLatency;

    std::map<uint32_t,uint32_t > conddatamap;//this is kept for now.If not used in offline analysis, it can be dropped
    std::map<std::string,std::vector<tbeam::cbc> > cbcs;
    std::map< std::string,std::vector<tbeam::hit> > dutHits;//
    std::map< std::string, std::vector<tbeam::cluster> > cbcClusters;//only in sparsified mode
    std::map<std::string, std::vector<tbeam::stub>> cbcStubs;//
    
    std::map< std::string, std::vector<tbeam::cluster> > offlineClusters;//only in unsparsified mode
    std::map< std::string, std::vector<tbeam::stub> > offlineStubs;//
    ClassDef(Event,1)
  };
}
#endif
