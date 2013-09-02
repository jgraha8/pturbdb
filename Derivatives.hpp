#ifndef DERIVATIVES_H
#define DERIVATIVES_H

namespace DERIVATIVES_H {
  
  // Finite difference class
  class FiniteDiff {
    
  public:
    int norder;

  public:
    // Constructor
    FiniteDiff(){};
    FiniteDiff( int _ndim, int *_dims );
    FiniteDiff( const FiniteDiff &g );
    // Deconstructor
    ~FiniteDiff();

  protected:
    void FiniteDiffInit( int _norder );

  };

}

#endif
