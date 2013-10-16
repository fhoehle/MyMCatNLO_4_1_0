// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsAnalysis class.
//

#include "HiggsAnalysis.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;
using namespace MCatNLO;

namespace {

  // check particle doesn't have a shower t->t g branching
  bool isLastInShower(const Particle & p) {
    return p.children().size() > 1 
      && p.children()[0]->id() != p.id()
      && p.children()[1]->id() != p.id();
  }

  // selector struct
  struct Higgs {
    static bool AllCollisions() { return false; }
    static bool AllSteps() { return true; }
    // ===
    // pick the last instance from the shower
    static bool FinalState() { return false; }
    static bool Intermediate() { return true; }
    // ===
    static bool Check(const Particle & p) {
      return (abs(p.id()) == ParticleID::h0); 
    }
  };
}

void HiggsAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tParticleVector higgses;
  // find the t tbar at the end of the shower
  event->select(back_inserter(higgses), ThePEG::ParticleSelector<Higgs>());

  // define cuts as in mcatnlo_hwanhgg.f
  double ycut  = 2.0;
  double ptcut = 20.0;

  // only fill Histos if there is at least one Higgs
  if( higgses.size() > 0 ) {
    Lorentz5Momentum higgs=higgses[0]->momentum();
    // Higgs kinematics
    double pt = higgs.perp()/GeV;
    double yh = higgs.rapidity();
    double mh = higgs.mass()/GeV;
    // Jet kinematics
    //double jetpt = pt;
    double jeteta = yh;
    
    // all cut combinations (as in mcatnlo_hwanhgg.f)
    bool passcuts1 = ( pt > ptcut );            // cut on pthiggs
    bool passcuts2 = ( abs(yh) < ycut );         // cuts on |y|_higgs
    bool passcuts3=  ( passcuts2 && passcuts1 ); // pass both cuts
    
    double weight=event->weight();               // event weight
    
    //std::cout<<weight<<std::endl;
    

    // Fill all the histograms
    for(int ix=0;ix<passcuts1+1;++ix) {
      // y Higgs & eta Jet
      _hy[ix]      ->addWeighted(yh,weight);          // Higgs rapidity
      _jeteta[ix]  ->addWeighted(jeteta,weight);        // Jet pseudo-rapidity
    }
    for(int ix=0;ix<passcuts2+1;++ix) {
      // pt Higgs
      _hpt[ix]     ->addWeighted(pt,weight);          // pt in 0-200
      _hpt[2+ix]   ->addWeighted(pt,weight);          // pt in 200-400
      _hpt[4+ix]   ->addWeighted(pt,weight);          // pt in 0-500
      if(pt>0.) 
	_hptlog[ix]->addWeighted(log10(pt),weight);   // pt in log scale
    }
    for(int ix=0;ix<passcuts3+1;++ix)
      // Higgs mass
      _hm[ix]->addWeighted(mh/GeV,weight);
  }

  if(!(event->number()%modPrint)) printPlots(false);
  return;

}

IBPtr HiggsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr HiggsAnalysis::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<HiggsAnalysis> HiggsAnalysis::initHiggsAnalysis;
// Definition of the static class description member.

void HiggsAnalysis::Init() {

  static ClassDocumentation<HiggsAnalysis> documentation
    ("There is no documentation for the HiggsAnalysis class");

}

void HiggsAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  printPlots(true);

}

void HiggsAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();

  modPrint=10000;

  for(unsigned int ix=0;ix<2;++ix) {
    _hpt[ix  ]    = new_ptr(Histogram(0.,200.,100));
    _hptlog[ix]   = new_ptr(Histogram(.1,5.,98));
    _hy[ix]       = new_ptr(Histogram(-4.,4.,40));
    _jeteta[ix]   = new_ptr(Histogram(-4.,4.,40));
    _hm[ix]       = new_ptr(Histogram(300.,1000.,70));
    _hpt[ix+2]    = new_ptr(Histogram(200.,400.,100));
    _hpt[ix+4]    = new_ptr(Histogram(0.,500.,50));
  }
}


void HiggsAnalysis::printPlots(bool rescale) {

  
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;

  if(rescale) {
    for(unsigned int ix=0;ix<2;++ix) {
      _hpt[ix  ]->normaliseToCrossSection();
      _hpt[ix+2]->normaliseToCrossSection();
      _hpt[ix+4]->normaliseToCrossSection();
      _hptlog[ix]->normaliseToCrossSection();
      _hm[ix]->normaliseToCrossSection();
      _hy[ix]->normaliseToCrossSection();
      _jeteta[ix]->normaliseToCrossSection();
    }
  }
  
  for(unsigned int ix=0;ix<2;++ix) {
    string cc;
    if(ix==0)
      cc = string(" ");
    else
      cc = string(" cuts");
    
    _hpt[ix  ]->topdrawMCatNLO(outfile,Frame,"RED","Higgs pt"+cc);
    _hptlog[ix]->topdrawMCatNLO(outfile,Frame,"RED","Higgs log[pt]"+cc);
    _hy[ix]->topdrawMCatNLO(outfile,Frame,"RED","Higgs y"+cc);
    _jeteta[ix]->topdrawMCatNLO(outfile,Frame,"RED","jet y"+cc);
  }

  for(unsigned int ix=0;ix<2;++ix) {
    string cc;
    if(ix==0)
      cc = string(" ");
    else
      cc = string(" cuts");

    _hm[ix]->topdrawMCatNLO(outfile,Frame,"RED","mH"+cc);
    _hpt[ix+2]->topdrawMCatNLO(outfile,Frame,"RED","Higgs pt"+cc);
    _hpt[ix+4]->topdrawMCatNLO(outfile,Frame,"RED","Higgs pt"+cc);
  }
  
  outfile.close();
  return;

}
