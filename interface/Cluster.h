#ifndef __Cluster__h
#define __Cluster__h
#include<stdint.h>
namespace tbeam {
class cluster : public TObject {
 public:
  cluster() :
   firstStrip(9999),
   firstRow(9999),
  edge(9999),
  column(999),
  size(999),
  threshold(0)
  {}

  cluster(unsigned int fStrip, unsigned int ftRow, unsigned int ed, unsigned int col, uint16_t s) :
  firstStrip(fStrip),
  firstRow(ftRow),
  edge(ed),
  column(col),
  size(s),
  threshold(0)
  {}

  cluster(unsigned int fStrip, unsigned int ftRow, unsigned int ed, unsigned int col, uint16_t s, uint16_t th) :
  firstStrip(fStrip),
  firstRow(ftRow),
  edge(ed),
  column(col),
  size(s),
  threshold(th)
  {}
  virtual ~cluster(){}
  float center() const { return float(firstStrip) + 0.5*(float(size) - 1.); }
  std::pair< float, float > barycenter() const { return std::make_pair(column, center()); }
  
  unsigned int firstStrip;
  unsigned int firstRow;
  unsigned int edge;
  unsigned int column;
  uint16_t size;
  uint16_t threshold;
  ClassDef(cluster,1)
  };
}
#endif
