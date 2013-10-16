// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopAnalysis class.
//

#include "TopAnalysis.h"
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
  struct TTBar {
    static bool AllCollisions() { return false; }
    static bool AllSteps() { return true; }
    // ===
    // pick the last instance from the shower
    static bool FinalState() { return false; }
    static bool Intermediate() { return true; }
    // ===
    static bool Check(const Particle & p) { 
      return abs(p.id()) == ParticleID::t && isLastInShower(p);
    }
  };

  // compute delta phi
  double deltaPhi(const Lorentz5Momentum & p1, const Lorentz5Momentum &p2) {
    Energy pt1 = p1.perp();
    Energy pt2 = p2.perp();
    if(pt1>0.*GeV&&pt2>0.*GeV) {
      double tmp = (p1.x()*p2.x()+p1.y()*p2.y())/(pt1*pt2);
      if(tmp>=1.)  tmp =  1.;
      if(tmp<=-1.) tmp = -1.;
      return acos(tmp);
    }
    else
      return 1e8;
  }
}

void TopAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tParticleVector particles;
  // find the t tbar at the end of the shower
  event->select(back_inserter(particles), ThePEG::ParticleSelector<TTBar>());
  // must be a t and tbar
  if(particles.size()!=2) {
    if(!(event->number()%modPrint)) printPlots(false);
    return;
  }

  // put the top first
  if(particles[0]->id()<0) swap(particles[0],particles[1]);
  // decide which cuts to use
  Energy ECMS = (event->incoming().first->momentum()+
		 event->incoming().second->momentum()).m();
  double ycut,ptcut;
  if(ECMS>5000.0*GeV) {
    ycut = 2.5;
    ptcut = 30.;
  }
  else {
    ycut = 1.0;
    ptcut = 15.;
  }
  // momentum of ttbar system
  Lorentz5Momentum ptt=particles[0]->momentum()+particles[1]->momentum();
  // pt of system
  double pt = ptt.perp()/GeV;
  // azimuth different of tops
  double azi     = deltaPhi(particles[0]->momentum(),particles[1]->momentum());
  double azinorm = (Constants::pi-azi)/Constants::pi;
  // rapidity of t and tbar
  double y1     = particles[0]->momentum().rapidity();
  double y2     = particles[1]->momentum().rapidity();
  // pseudorapidity
  double eta1   = particles[0]->momentum().eta();
  double eta2   = particles[1]->momentum().eta(); 
  // delta R
  double deltaR = sqrt(sqr(eta1-eta2)+sqr(azi)); 
  // pt of t and tbar
  double pt1=particles[0]->momentum().perp()/GeV;
  double pt2=particles[1]->momentum().perp()/GeV;

  bool siq1flag = ( pt1 > ptcut && abs(y1) < ycut ); 
  bool siq2flag = ( pt2 > ptcut && abs(y2) < ycut );
  bool ddflag = siq1flag && siq2flag ;

  for(int ix=0;ix<ddflag+1;++ix) {
    // properties of the ttbar system
    // pt
    _ttpt[ix]   ->addWeighted(pt,event->weight());
    _ttpt[2+ix]   ->addWeighted(pt,event->weight());
    if(pt>0.) _ttptlog[ix]->addWeighted(log10(pt),event->weight());
    // mass of ttbar system
    _ttinvm[ix]->addWeighted(ptt.m()/GeV,event->weight());
    // rapidity of ttbar system
    _raptt[ix]->addWeighted(ptt.rapidity(),event->weight());
    // azimuth difference of tops
    _ttazi[ix]->addWeighted(azi,event->weight());
    _ttaziB[ix]->addWeighted(azi,event->weight());
    if(azinorm>0.) _ttpiazi[ix]  ->addWeighted(log10(azinorm),
					       event->weight());
    // delta R
    _ttdeltaRB[ix]->addWeighted(deltaR,event->weight());
    _ttdeltaR[ix]->addWeighted(deltaR,event->weight());
    // delta eta
    _ttdeta[ix]->addWeighted(eta1-eta2,event->weight());
    // delta y 
    _deltay[ix]->addWeighted(y1-y2,event->weight());

    //Plot single-inclusive variables only when no corr cuts are applied
    if(ix==0) {
      // pT of top
      _toppt[ix]  ->addWeighted(pt1,event->weight());
      _toppt[2+ix]  ->addWeighted(pt1,event->weight());
      _topptlog[ix]->addWeighted(log10(pt1),event->weight());
      // pT of t bar
      _tbpt[ix]    ->addWeighted(pt2,event->weight());
      _tbpt[2+ix]  ->addWeighted(pt2,event->weight());
      _tbptlog[ix]->addWeighted(log10(pt2),event->weight());
      // rapidity of t and tbar
      _rapt[ix]     ->addWeighted(y1,event->weight());
      _raptb[ix]    ->addWeighted(y2,event->weight());
    }
  }

  //Plot single-inclusive variables with cuts
  if( abs(y2) < ycut ) {
    _tbpt[1]    ->addWeighted(pt2,event->weight());
    _tbpt[3]  ->addWeighted(pt2,event->weight());
    _tbptlog[1]->addWeighted(log10(pt2),event->weight());
  }
  if( pt2 > ptcut )
    _raptb[1]    ->addWeighted(y2,event->weight());

  if( abs(y1) < ycut ) {
    _toppt[1]  ->addWeighted(pt1,event->weight());
    _toppt[3]  ->addWeighted(pt1,event->weight());
    _topptlog[1]->addWeighted(log10(pt1),event->weight());
  }
  if( pt1 > ptcut )
    _rapt[1]     ->addWeighted(y1,event->weight());

  if(!(event->number()%modPrint)) printPlots(false);
  return;
  
}

IBPtr TopAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr TopAnalysis::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<TopAnalysis> TopAnalysis::initTopAnalysis;
// Definition of the static class description member.

void TopAnalysis::Init() {

  static ClassDocumentation<TopAnalysis> documentation
    ("There is no documentation for the TopAnalysis class");

}

void TopAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  this->printPlots(true);
}

void TopAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();

  // set the modPrint
  modPrint = 10000;


  for(unsigned int ix=0;ix<2;++ix) {
    _ttpt[ix  ]   = new_ptr(Histogram(0.,100.,50));
    _ttpt[ix+2]   = new_ptr(Histogram(80.,2000.,96));
    _ttptlog[ix]  = new_ptr(Histogram(0.1,5.,98));
    _ttinvm[ix]   = new_ptr(Histogram(300.,1000.,70));
    _ttazi[ix]    = new_ptr(Histogram(0.,Constants::pi,20));
    _ttdeltaR[ix] = new_ptr(Histogram(0.,3.*Constants::pi,60));
    _toppt[ix]    = new_ptr(Histogram(0.,500.,100));
    _toppt[ix+2]  = new_ptr(Histogram(400,2400,100));
    _topptlog[ix] = new_ptr(Histogram(0.1,5.,98));
    _tbpt[ix]     = new_ptr(Histogram(0.,500.,100));
    _tbpt[ix+2]   = new_ptr(Histogram(400,2400,100));
    _tbptlog[ix]  = new_ptr(Histogram(0.1,5.,98));
    _ttdeta[ix]   = new_ptr(Histogram(-4.,4.,40));
    _raptt[ix]    = new_ptr(Histogram(-4.,4.,80));
    _deltay[ix]   = new_ptr(Histogram(-4.,4.,40));
    _ttaziB[ix]   = new_ptr(Histogram(2.*Constants::pi/3.,Constants::pi,20));
    _ttdeltaRB[ix]= new_ptr(Histogram(2*Constants::pi/3,4*Constants::pi/3,40));
    _rapt[ix]     = new_ptr(Histogram(-4.,4.,80));
    _raptb[ix]    = new_ptr(Histogram(-4.,4.,80));
    _ttpiazi[ix]  = new_ptr(Histogram(-4.,0.1,78));
  }
}


void TopAnalysis::printPlots(bool rescale) {
  
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  
  for(unsigned int ix=0;ix<2;++ix) {
    if(rescale) {
      _ttpt[ix  ]   ->normaliseToCrossSection();
      _ttptlog[ix]  ->normaliseToCrossSection();
      _ttinvm[ix]   ->normaliseToCrossSection();
      _ttazi[ix]    ->normaliseToCrossSection();
      _ttdeltaR[ix] ->normaliseToCrossSection();
      _toppt[ix]    ->normaliseToCrossSection();
      _topptlog[ix] ->normaliseToCrossSection();
      _tbpt[ix]     ->normaliseToCrossSection();
      _tbptlog[ix]  ->normaliseToCrossSection();
      _ttdeta[ix]   ->normaliseToCrossSection();
      _raptt[ix]    ->normaliseToCrossSection();
      _deltay[ix]   ->normaliseToCrossSection();
      _ttaziB[ix]   ->normaliseToCrossSection();
      _ttdeltaRB[ix]->normaliseToCrossSection();
      _rapt[ix]     ->normaliseToCrossSection();
      _raptb[ix]    ->normaliseToCrossSection();
      _ttpiazi[ix]  ->normaliseToCrossSection();
    }
    
    string cc;
    if(ix==0) cc = string(" ");
    else cc = string(" cuts");
    
    //In order to have the same labels as in mcatnlo_hwan**.f, get
    //them from calls to mbook()
    _ttpt[ix]     ->topdrawMCatNLO(outfile,Frame,"RED","tt pt"+cc);
    _ttptlog[ix]  ->topdrawMCatNLO(outfile,Frame,"RED","tt log[pt]"+cc);
    _ttinvm[ix]   ->topdrawMCatNLO(outfile,Frame,"RED","tt inv m"+cc);
    _ttazi[ix]    ->topdrawMCatNLO(outfile,Frame,"RED","tt azimt"+cc);
    _ttdeltaR[ix] ->topdrawMCatNLO(outfile,Frame,"RED","tt del R"+cc);
    _tbpt[ix]     ->topdrawMCatNLO(outfile,Frame,"RED","tb pt"+cc);
    _tbptlog[ix]  ->topdrawMCatNLO(outfile,Frame,"RED","tb log[pt]"+cc);
    _toppt[ix]    ->topdrawMCatNLO(outfile,Frame,"RED","t pt"+cc);
    _topptlog[ix] ->topdrawMCatNLO(outfile,Frame,"RED","t log[pt]"+cc);
    _ttdeta[ix]   ->topdrawMCatNLO(outfile,Frame,"RED","tt delta eta"+cc);
    _raptt[ix]    ->topdrawMCatNLO(outfile,Frame,"RED","y_tt"+cc);
    _deltay[ix]   ->topdrawMCatNLO(outfile,Frame,"RED","delta y"+cc);
    _ttaziB[ix]   ->topdrawMCatNLO(outfile,Frame,"RED","tt azimt"+cc);
    _ttdeltaRB[ix]->topdrawMCatNLO(outfile,Frame,"RED","tt del R"+cc);  
    _raptb[ix]    ->topdrawMCatNLO(outfile,Frame,"RED","y_tb"+cc); 
    _rapt[ix]     ->topdrawMCatNLO(outfile,Frame,"RED","y_t"+cc);
    _ttpiazi[ix]  ->topdrawMCatNLO(outfile,Frame,"RED","tt log[pi-azimt]"+cc);
  }

  for(unsigned int ix=0;ix<2;++ix) {

    if(rescale) {
      _ttpt[ix+2]   ->normaliseToCrossSection();
      _toppt[ix+2]  ->normaliseToCrossSection();
      _tbpt[ix+2]   ->normaliseToCrossSection();
    }

    string cc = string(" ");
    if(ix != 0) cc = string(" cuts");
    
    _ttpt[ix+2]  ->topdrawMCatNLO(outfile,Frame,"RED","tt pt"+cc);
    _tbpt[ix+2]  ->topdrawMCatNLO(outfile,Frame,"RED","tb pt"+cc);
    _toppt[ix+2] ->topdrawMCatNLO(outfile,Frame,"RED","t pt"+cc);
  }

  outfile.close();
  return;
}
