// -*- C++ -*-
#ifndef MCATNLO_HiggsAnalysis_H
#define MCATNLO_HiggsAnalysis_H
//
// This is the declaration of the HiggsAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

namespace MCatNLO {

using namespace ThePEG;
using namespace Herwig;

/**
 * Here is the documentation of the HiggsAnalysis class.
 *
 * @see \ref HiggsAnalysisInterfaces "The interfaces"
 * defined for HiggsAnalysis.
 */
class HiggsAnalysis: public AnalysisHandler {

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<HiggsAnalysis> initHiggsAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsAnalysis & operator=(const HiggsAnalysis &);

  void printPlots(bool=true);
  int modPrint;

private:

  HistogramPtr _hpt[6],_hptlog[2],_hm[2],_hy[2],_jeteta[2];

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiggsAnalysis. */
template <>
struct BaseClassTrait<MCatNLO::HiggsAnalysis,1> {
  /** Typedef of the first base class of HiggsAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiggsAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<MCatNLO::HiggsAnalysis>
  : public ClassTraitsBase<MCatNLO::HiggsAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "MCatNLO::HiggsAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiggsAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class HiggsAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiggsAnalysis.so"; }
};

/** @endcond */

}

#endif /* MCATNLO_HiggsAnalysis_H */
