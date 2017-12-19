#ifndef __Track__h
#define __Track__h
#include<stdint.h>
#include<iostream>
#include<iomanip>
namespace tbeam {
class Track : public TObject {
 public:
  Track():
    xPos_(0.),yPos_(0.),dxdz_(0.),dydz_(0.),chi2_(0.),ndof_(0.),xPosError_(0.),yPosError_(0.)
    {
      xPosPrevHit_ = 0.;
      yPosPrevHit_ = 0.;
      xPosErrsPrevHit_ = 0.;
      yPosErrsPrevHit_ = 0.;
      xPosNextHit_ = 0.;
      yPosNextHit_ = 0.;
      xPosErrNextHit_ = 0.;
      yPosErrNextHit_ = 0.;
    } 
  Track(const float x, const float y, const float dxdz, const float dydz, 
        const float chi2, const float ndof, const float xErr, const float yErr):
    xPos_(x),yPos_(y),dxdz_(dxdz),dydz_(dydz),chi2_(chi2),ndof_(ndof),xPosError_(xErr),yPosError_(yErr)
    {
      xPosPrevHit_ = 0.;
      yPosPrevHit_ = 0.;
      xPosErrsPrevHit_ = 0.;
      yPosErrsPrevHit_ = 0.;
      xPosNextHit_ = 0.;
      yPosNextHit_ = 0.;
      xPosErrNextHit_ = 0.;
      yPosErrNextHit_ = 0.;
    } 
  Track(const float x, const float y, const float dxdz, const float dydz, 
        const float chi2, const float ndof, const float xErr, const float yErr,
        const float xposPrevHit, const float yposPrevHit, const float xposErrsPrevHit, const float yposErrsPrevHit,
        const float xposNextHit, const float yposNextHit, const float xposErrNextHit, const float yposErrNextHit):
    xPos_(x),yPos_(y),dxdz_(dxdz),dydz_(dydz),chi2_(chi2),ndof_(ndof),xPosError_(xErr),yPosError_(yErr)
    {
      xPosPrevHit_ = xposPrevHit;
      yPosPrevHit_ = yposPrevHit;
      xPosErrsPrevHit_ = xposErrsPrevHit;
      yPosErrsPrevHit_ = yposErrsPrevHit;
      xPosNextHit_ = xposNextHit;
      yPosNextHit_ = yposNextHit;
      xPosErrNextHit_ = xposErrNextHit;
      yPosErrNextHit_ = yposErrNextHit;
    } 
    
    float xPos() const { return xPos_; }
    float yPos() const { return yPos_; }
    float dxdz() const { return dxdz_; }
    float dydz() const { return dydz_; }
    float chi2() const { return chi2_; }
    float ndof() const { return ndof_; }
    float xPosErr() const { return xPosError_;}
    float yPosErr() const { return yPosError_; }

    float xPosPrevHit() const     { return xPosPrevHit_; }
    float yPosPrevHit() const     { return yPosPrevHit_; }
    float xPosErrsPrevHit() const { return xPosErrsPrevHit_; }
    float yPosErrsPrevHit() const { return yPosErrsPrevHit_; }
    float xPosNextHit() const     { return xPosNextHit_; }
    float yPosNextHit() const     { return yPosNextHit_; }
    float xPosErrNextHit() const  { return xPosErrNextHit_; }
    float yPosErrNextHit() const  { return yPosErrNextHit_; }
    
    friend std::ostream &operator<<(std::ostream& os, const tbeam::Track& tk) {    
      os << "X="       << std::setw(6) << tk.xPos() 
	 << " Y="      << std::setw(6) << tk.yPos()
	 << " dxdz="   << std::setw(6) << tk.dxdz()
	 << " dydz="   << std::setw(6) << tk.dydz()
	 << " chi2="   << std::setw(6) << tk.chi2()
	 << " ndof="   << std::setw(6) << tk.ndof()
	 << " xError=" << std::setw(6) << tk.xPosErr()
	 << " yError=" << std::setw(6) << tk.yPosErr()
	 << std::endl; 
      return os;
    }
    
  private:
    float xPos_;//Track impact x-position at the reference plane
    float yPos_;//Track impact y-position at the reference plane
    float dxdz_;//x-Gradient along beam direction
    float dydz_;//y-Gradient along beam direction
    float chi2_;//Chi-square from track fitting
    float ndof_;//#degrees of freedom
    
    //optional requirement
    float xPosError_;//Track impact x-position error
    float yPosError_;//Track impact y-position error 
    //additional info available from FNAL beam test
    float xPosPrevHit_;
    float yPosPrevHit_;
    float xPosErrsPrevHit_;
    float yPosErrsPrevHit_;
    float xPosNextHit_;
    float yPosNextHit_;
    float xPosErrNextHit_;
    float yPosErrNextHit_;
    ClassDef(Track,1)
      };
}
#endif
