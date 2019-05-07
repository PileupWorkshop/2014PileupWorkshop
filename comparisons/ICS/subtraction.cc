// Test of subtraction performance of different background mitigation techniques discussed in the ICS paper.
// This file was copied from comparisons/review/subtraction-tests.cc version 1.0.0 and further modified.
// The current implementation of the event-wide CS and ICS works only for CHS events

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "SimpleHist.hh"
#include "ProfileHist.hh"
#include "Matching.hh"
#include "helpers.hh"
#include "FunctionOfEvent.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include "fastjet/tools/Filter.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "puppiContainer.hh"

#include "TH1D.h"
#include "TFile.h"

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;
using namespace fastjet;

//========================================================================
// (statistical) analysis tools
//========================================================================

/// a specific class that handles Matteo's format
///
/// Usage: 
///   declare it with a name
///   for a given entry add
///     data.add_entry();   // for un-matched 
///     data.add_entry(hard_value, full_value); // for matched
class MatteoAverageAndError{
public:
  MatteoAverageAndError(){}
  
  AverageAndError match_eff, delta, njets, avg_hard;
  CorrelationCoefficient correl;
  
  void add_entry(){
    match_eff += 0.0;
  }

  void add_entry(double hard, double full){
    match_eff += 1.0;
    avg_hard += hard;
    delta += full-hard;
    correl.add_entry(hard, full);
  }

};

/// a class for output management
class OutputInfo{
public:
  OutputInfo(string outfile, unsigned int npu, double ptmin)
    : _outfile(outfile), _npu(npu), _ptmin(ptmin){}

  // what will get 
  string _outfile;
  ostringstream header;
  unsigned int _npu;
  double _ptmin;

  // the various infos
  map<string,SimpleHist> hists; //< normalised to the bin size
  map<string,ProfileHist> profiles;
  map<string,AverageAndError> averages;

  map<string,MatteoAverageAndError> matteos;

  /// output once in a while
  /// Arguments are:
  ///   iev    the number of events so far
  void write(const unsigned int iev){
    ofstream ostr(_outfile.c_str());

    cout << "Output of results for nev = " << iev << endl;

    ostr << header.str() << endl;
    ostr << "# nev = " << iev << endl;
  
    // write out averages
    for(map<string,AverageAndError>::iterator result = averages.begin(); 
        result != averages.end(); result++) {
      ostr << "# " << result->first << ": average = " << result->second.average()
           << " stddev = " << result->second.sd() << endl;
    }

    // Mattteo's output of the form 
    //   name  <delta>  sigma_delta   corr.coef  match.eff  match.eff_error
    ostr << "#----------------------------------------------------------------------" << endl;
    ostr << "# npu    ptmin        <delta>                sigma_delta        corr.coef        Njets(>20GeV)          match.eff              <hard>              name " << endl;
    for(map<string,MatteoAverageAndError>::iterator result = matteos.begin(); 
        result != matteos.end(); result++) {
      ostr << setprecision(4) 
           << setw(4) << _npu   << "    "
           << setw(4) << _ptmin << "    "
           << setw(8) << result->second.delta.average()     << " +- " << setw(8) << result->second.delta.error()       << "    "
           << setw(8) << result->second.delta.sd()          << " +- " << setw(8) << result->second.delta.error_on_sd() << "    "
           << setw(8) << result->second.correl.r()          << "    "
           << setw(8) << result->second.njets.average()     << " +- " << setw(8) << result->second.njets.error()       << "    "
           << setw(8) << result->second.match_eff.average() << " +- " << setw(8) << result->second.match_eff.error()   
           << setw(8) << result->second.avg_hard.average()  << " +- " << setw(8) << result->second.avg_hard.error()    << "    "
           << "   # " << setw(-20) << result->first << endl;
    }
    ostr << "#----------------------------------------------------------------------" << endl;

    // histograms (normalised to the number of hard jets) 
    for(map<string,SimpleHist>::iterator hist = hists.begin(); 
	hist != hists.end(); hist++) {
      ostr << "# " << hist->first << " [1/x-axis units] " << endl;
      output(hist->second, &ostr, 1.0/iev/hist->second.binsize());
      ostr << endl << endl;
    }

    // output profiles
    for(map<string,ProfileHist>::iterator hist = profiles.begin(); 
        hist != profiles.end(); hist++) {
      ostr << "# " << hist->first << " [average within each bin; stddev; error on avg; average of squares] " << endl;
      output_noNaN(hist->second, &ostr);
      ostr << endl << endl;
    }
  }
};

//----------------------------------------------------------------------
// declare histograms (they depend on boundaries so have to be declared)
void declare_output(OutputInfo &output, const string &tag, const double jetptmin){
  output.hists["hist_delta_pt_"   +tag].declare(-jetptmin, jetptmin, 200);
  output.hists["hist_delta_m_"    +tag].declare(-0.25*jetptmin, 0.25*jetptmin, 200);
  output.hists["hist_delta_eta_"   +tag].declare(-0.5, 0.5, 200);
  output.hists["hist_delta_width_"    +tag].declare(-0.25, 0.25, 200);

  output.hists["hist_pt_"   +tag].declare(0.0, 3*jetptmin, 300);
  output.hists["hist_m_"    +tag].declare(0.0, 0.5*jetptmin, 200);
  output.hists["hist_eta_"   +tag].declare(-5, 5, 200);
  output.hists["hist_width_"    +tag].declare(0, 1, 200);
}





//----------------------------------------------------------------------
// function to compute jet width
double width(const fastjet::PseudoJet &jet){
  const double pi=3.1415926;
  if (!jet.has_constituents()){
    return -0.1;
  }
  double width = 1e-8;
  double ptSum = 0;
  
  if (jet.constituents().size() < 2) return width;
  
  vector<fastjet::PseudoJet> constituents = jet.constituents();
  
  for (unsigned int iconst=0; iconst<constituents.size(); iconst++){
    double wdPhi = jet.phi() - constituents.at(iconst).phi();
    if(wdPhi>pi) wdPhi=2*pi-wdPhi;
    double wdEta = jet.eta() - constituents.at(iconst).eta();
    double dR = sqrt(wdPhi*wdPhi+wdEta*wdEta);
    double pt = constituents.at(iconst).pt();
    width += dR * pt;
    ptSum += pt;
  }
  
  return width/ptSum;
}



//----------------------------------------------------------------------
// record what has been going on
void record(const string &name,
            const vector<PseudoJet> &hard_jets, 
            const vector<PseudoJet> &subt_jets,
            Matching &matching,
            OutputInfo &output,
            double jet_rapmax,
            bool print){

  // write out jets in subtracted event
  if ( print ) {
    cerr << "Subtracted event with " << name << ":" << endl; 
    for (unsigned int i=0; i < 4U && i < subt_jets.size(); i++) {
      cerr << "  jet " << i << ": "  << subt_jets[i] << endl;
    }
  }

  // set up the set of full/subtracted jets from which to match
  matching.set_full_jets(subt_jets);
     
  // run over the hard jets
  for (unsigned int i=0; i < hard_jets.size(); i++) {
    // for each hard jet, find the corresponding full/subtracted jet that matches
    // (if any)
    const PseudoJet * match = matching.match(hard_jets[i]);
    if (match){
      output.matteos["pt_"+name].add_entry(hard_jets[i].pt(), match->pt());
      output.matteos["m_" +name].add_entry(hard_jets[i].m(),  match->m());
      output.matteos["eta_" +name].add_entry(hard_jets[i].pseudorapidity(),  match->pseudorapidity());
      output.matteos["width_" +name].add_entry(width(hard_jets[i]),  width(*match));

      output.hists["hist_delta_pt_"+name].add_entry(match->pt() - hard_jets[i].pt());
      output.hists["hist_delta_m_" +name].add_entry(match->m()  - hard_jets[i].m() );
      output.hists["hist_delta_eta_" +name].add_entry(match->pseudorapidity()  - hard_jets[i].pseudorapidity() );
      output.hists["hist_delta_width_" +name].add_entry(width(*match)  - width(hard_jets[i]) );
      
      output.hists["hist_pt_"+name].add_entry(match->pt());
      output.hists["hist_m_" +name].add_entry(match->m());
      output.hists["hist_eta_" +name].add_entry(match->pseudorapidity());
      output.hists["hist_width_" +name].add_entry(width(*match));
    } else {
      output.matteos["pt_"+name].add_entry();
      output.matteos["m_" +name].add_entry();
      output.matteos["eta_" +name].add_entry();
      output.matteos["width_" +name].add_entry();
    }    
  }
}





//========================================================================
// Jet transformers
//========================================================================


//------------------------------------------------------------------------
/// a quick wrapper for  transformer with an additional short name
///
/// It deletes the object at the end when going out of scope
class TransformerWithName : public Transformer{
public:
  /// ctor
  TransformerWithName(Transformer *base_transformer, const string name)
    : _base_transformer(base_transformer), _name(name){}

  /// description
  virtual string description() const{
    return _name + string(": ") + _base_transformer->description();
  }

  /// use the underlying transformer for doing the work
  virtual PseudoJet result(const PseudoJet &jet) const{
    return _base_transformer->result(jet);
  }

  /// retreive info
  Transformer* base() const{ return _base_transformer.get();}
  string name() const{ return _name;}

protected:
  SharedPtr<Transformer> _base_transformer;      ///< the object itself
  string _name; ///< its (short) name  
};


//------------------------------------------------------------------------
// identity transformer
class IdentityTransformer : public fastjet::Transformer{
public:
  IdentityTransformer(){}
  virtual std::string description() const{ return "identity transformer";}
  virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const{ return jet;}
};


//------------------------------------------------------------------------
// cleansing will require a transformer
class TJetCleanser : public fastjet::Transformer{
public:
  /// ctor
  ///  \param Rsub             radius for subjet clustering
  ///  \param cleansing_mode   JVF/linear cleansing or Gaussian cleamsing
  ///  \param input_mode       together or separate
  ///  \param gamma0           gamma0 parameter for all cleansing modes 
  ///  \param gamma1           gamma1 parameter for Gaussian cleansing
  ///  \param gamma0_width     gamma0 width for Gaussian cleansing
  ///  \param gamma1_width     gamma1 width for Gaussian cleansing
  ///  \param chs_rescaling    charged tracks from PU vertices will be 
  ///                          rescaled by that factor before calling
  ///                          cleansing (meant for use with CHS events)
  TJetCleanser(const double Rsub,
               const fastjet::contrib::JetCleanser::cleansing_mode cleansing_mode,
               const fastjet::contrib::JetCleanser::input_mode     input_mode,
               double ftrim, double ga0, 
               double ga1=1.0, 
               double ga0width=1.0, double ga1width=1.0,
               double chs_rescaling = 1.0) 
    : _cleanser(Rsub, cleansing_mode, input_mode), 
      _separate(input_mode==contrib::JetCleanser::input_nc_separate),
      _chs_rescaling(chs_rescaling){
    _cleanser.SetGaussianParameters(ga0, ga1, ga0width, ga1width);
    _cleanser.SetLinearParameters(ga0);
    _cleanser.SetTrimming(ftrim); 
  }

  virtual std::string description() const{
    ostringstream oss;
    oss << _cleanser.description() << " (with chs rescaling factor " << _chs_rescaling << ")";
    return oss.str();
  }

  virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const{
    
    // make sure no particles have 0 pt
    vector<PseudoJet> constits = (!SelectorIsPureGhost())(jet.constituents());
    if (constits.size()==0) return PseudoJet();
    
    // separate the charged from the neutrqls and then split the charged
    // according to their origin (leading vertex or PU vertices)
    vector<PseudoJet> neutrals, charged, charged_lv, charged_pu;
    SelectorCharged().sift(constits, charged, neutrals);
    SelectorHard().sift(charged, charged_lv, charged_pu);
    
    // rescale charged PU tracks
    for (unsigned int i=0; i< charged_pu.size(); i++)
      charged_pu[i] *= _chs_rescaling;

    return (_separate)
      ? _cleanser(neutrals, charged_lv, charged_pu)
      : _cleanser(jet,      charged_lv, charged_pu); 
  }

protected:
  fastjet::contrib::JetCleanser _cleanser;
  bool _separate;
  double _chs_rescaling;
};

//========================================================================
// Event transformers
//========================================================================

typedef FunctionOfEvent<vector<PseudoJet> > EventTransformer;

//------------------------------------------------------------------------
// process the event
//
// Note that this takes ownership of the base transformer
class EventTransformerWithName : public EventTransformer{
public:
  /// ctor
  EventTransformerWithName(EventTransformer *base_transformer, const string name)
  : _base_transformer(base_transformer), _name(name){}

  /// description
  virtual string description() const{
    return _name + string(": ") + _base_transformer->description();
  }
  
  /// use the underlying transformer for doing the work
  virtual vector<PseudoJet> result(const vector<PseudoJet> &event) const{
    return _base_transformer->result(event);
  }

  /// retreive info
  EventTransformer* base() const{ return _base_transformer.get();}
  string name() const{ return _name;}

protected:
  SharedPtr<EventTransformer> _base_transformer;      ///< the object itself
  string _name; ///< its (short) name  
};

//------------------------------------------------------------------------
// A wrapper around the SoftKiller to use it as an event transformer
class SoftKillerWrapper : public EventTransformer{
public:
  SoftKillerWrapper(double rapmax, double tile_size, const Selector &sel)
    : _sk(rapmax, tile_size, sel) {}
  virtual vector<PseudoJet> result(const vector<PseudoJet> & event) const{
    return _sk(event);
  }
  virtual string description() const{ return _sk.description(); }
protected:
  const contrib::SoftKiller _sk;
};


//------------------------------------------------------------------------
// SoftKiller+Zeroing to use it as an event transformer
class SoftKillerZeroing : public EventTransformer{
public:
  SoftKillerZeroing(double rapmax, double tile_size, double Rzero)
    : _sk(rapmax, tile_size, !SelectorCharged()), _Rzero(Rzero), pt2_protection(100.0){}
  virtual vector<PseudoJet> result(const vector<PseudoJet> & event) const{
    vector<PseudoJet> sk_event = _sk(event);
    vector<PseudoJet> zeroed_event;
    Selector sel_vicinity = SelectorCircle(_Rzero) * SelectorCharged() * SelectorHard();
    for (const PseudoJet &p : sk_event){
      sel_vicinity.set_reference(p);
        
      // see if there is a charged track from LV nearby or if large enough pt
      if ((p.pt2()>pt2_protection) ||
          (sel_vicinity.count(sk_event)>0))
        zeroed_event.push_back(p);          
    }

    return zeroed_event;
  }
  virtual string description() const{ return _sk.description()+string("+zeroing"); }
protected:
  const contrib::SoftKiller _sk;
  const double _Rzero;
  const double pt2_protection;
};



////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  cout << "# " << cmdline.command_line() << "\n#" << endl;

  // create mixer that will construct events by mixing hard and pileup
  // events read from files given from command line using 
  // -hard hard_events_file(.gz) -pileup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  

  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  

  // particle and jet selection & R
  double particle_rapmax = cmdline.value<double>("-particle-rapmax",4.);
  double jet_rapmax = cmdline.value<double>("-jet-rapmax",2.5);
  double jet_ptmin  = cmdline.value<double>("-jet-ptmin",20);
  double R = cmdline.value<double>("-R",0.4);

  string included_subs = cmdline.value<string>("-methods", "unsub,area,sk50,sk45z02,jetByJetCS,eventWideCS,ICS,CSSK,puppi");
  included_subs.insert(0,1,',');
  included_subs.push_back(',');

  // some containier for what we want to output
  OutputInfo output(cmdline.value<string>("-out"), mixer.npu(), jet_ptmin);
  int iev_periodic=10;

  output.header << "# " << cmdline.command_line() << "\n#" << endl;
  output.header << "# time now = " << cmdline.time_stamp() << endl;
  output.header << "# Mixer: "  << mixer.description() << endl;

  // matching condition
  double maxdeltaR = 0.75*R;
  Matching matching(maxdeltaR);
  cout << "# matching: " << matching.description() << endl;
  output.header << "# matching: " << matching.description() << endl;

  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);             // the jet definition
  AreaDefinition area_def(active_area_explicit_ghosts);  // the area definition
  cout << "# jet_def: "  << jet_def.description() << endl;            
  cout << "# area_def: "  << area_def.description() << endl;           
  output.header << "# jet_def: "  << jet_def.description() << endl;            
  output.header << "# area_def: "  << area_def.description() << endl;           

  // select 2 hardest with rap and pt cuts:
  Selector sel_jets = SelectorAbsRapMax(jet_rapmax)*SelectorPtMin(jet_ptmin)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;
  output.header << "# sel_jets: " << sel_jets.description() << endl;

  // take particles only within |y|<4
  Selector sel_particles = SelectorAbsRapMax(particle_rapmax);
  cout << "# sel_particles: " << sel_particles.description() << endl;
  output.header << "# sel_particles: " << sel_particles.description() << endl;

  AverageAndError npu;

  output.hists["hist_pt_hard"   ].declare(0.0,   3*jet_ptmin, 300);
  output.hists["hist_m_hard"    ].declare(0.0, 0.5*jet_ptmin, 200);
  output.hists["hist_eta_hard"    ].declare(-5, 5, 200);
  output.hists["hist_width_hard"    ].declare(0, 1, 200);

  //----------------------------------------------------------------------
  // declare all the jet-based methods we want to try
  // 
  vector<TransformerWithName> transformers;
  vector<EventTransformerWithName> event_transformers;

  bool is_chs = cmdline.present("-chs");
  assert(is_chs && "The code below assumes one uses CHS events");
  
  //......................................................................
  // [unsub] unsubtracted
  if (included_subs.find(string(",unsub,"))!=string::npos){
    transformers.push_back(TransformerWithName(new IdentityTransformer(), "unsub"));
    declare_output(output, "unsub", jet_ptmin);
  }

  //......................................................................
  // [area] [eventWideCS] [ICS] area--median subtraction and CS and ICS
  double grid_size = 0.6; // different grid size with respect to /comparisons/review/subtraction-tests.cc, version 1.0.0
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(particle_rapmax, grid_size);
  GridMedianBackgroundEstimator * gmbge_neutral = new GridMedianBackgroundEstimator(particle_rapmax, grid_size);
  BackgroundRescalingYPolynomial * rescaling = new BackgroundRescalingYPolynomial(0.9084412, 0, -0.0189777, 0, 4.05981e-05);
  //BackgroundRescalingYPolynomial * rescaling = new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  gmbge->set_rescaling_class(rescaling);
  gmbge_neutral->set_rescaling_class(rescaling);

  contrib::ConstituentSubtractor cs;
  cs.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  cs.set_max_distance(0.25); // free parameter for the maximal allowed distance between particle i and ghost k
  cs.set_alpha(1);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
  cs.set_ghost_area(0.005); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  cs.set_max_eta(particle_rapmax); // parameter for the maximal eta to distribute the ghosts. Using the maximal rapidity for selection of input particles here.
  cs.set_background_estimator(gmbge_neutral);
  cs.set_common_bge_for_rho_and_rhom(true);
  cs.initialize();

  contrib::IterativeConstituentSubtractor ics;
  ics.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  vector<double> max_distances={0.2,0.15};  // maximal allowed distances for the individual iterations
  vector<double> alphas={1,1};   // alpha values for the individual iterations
  ics.set_parameters(max_distances,alphas);
  ics.set_remove_remaining_proxies(false);  // if set to true, then the unsubtracted ghosts are discarded for the next iteration
  ics.set_ghost_area(0.005); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  ics.set_max_eta(particle_rapmax);   // parameter for the maximal eta to distribute the ghosts. Using the maximal rapidity for selection of input particles here.
  ics.set_background_estimator(gmbge_neutral);
  ics.set_common_bge_for_rho_and_rhom(true);
  ics.initialize();

  contrib::SoftKiller sk(particle_rapmax,0.6);



  if (included_subs.find(string(",eventWideCS,"))!=string::npos){
    declare_output(output, "eventWideCS", jet_ptmin);
  }
  if (included_subs.find(string(",CSSK,"))!=string::npos){
    declare_output(output, "CSSK", jet_ptmin);
  }
  if (included_subs.find(string(",ICS,"))!=string::npos){
    declare_output(output, "ICS", jet_ptmin);
  }

  
  Subtractor *area_subtractor = new Subtractor(gmbge_neutral);
  area_subtractor->set_safe_mass();
  if (!mixer.massless()) area_subtractor->set_use_rho_m();

  if (included_subs.find(string(",area,"))!=string::npos){
    transformers.push_back(TransformerWithName(area_subtractor, "area"));

    output.header << "# Parameters for area subtraction" << endl;
    output.header << "#   rapidity rescaling for rho included" << endl;
    output.header << "#   grid size for GMBGE        = " << grid_size << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;

    declare_output(output, "area", jet_ptmin);
  }

  //......................................................................
  //   [sk50]    SoftKiller (0.5)
  if (included_subs.find(string(",sk50,"))!=string::npos){
    double sk_grid = 0.5;
    Selector sel_soft_killer = !SelectorCharged();
    event_transformers.push_back(EventTransformerWithName(new SoftKillerWrapper(particle_rapmax, sk_grid, sel_soft_killer), "sk50"));
    output.header << "# Parameters for the SoftKiller" << endl;
    output.header << "#   grid size for SoftKiller = " << sk_grid << endl;
    output.header << "#   description     = " << event_transformers.back().description();
    
    declare_output(output, "sk50", jet_ptmin);
  }

  //......................................................................                      
  //   [sk45z02] SoftKiller(0.45) + zeroing(0.2)                                                                                                                       
  if (included_subs.find(string(",sk45z02,"))!=string::npos){
    double sk_grid = 0.45;
    event_transformers.push_back(EventTransformerWithName(new SoftKillerZeroing(particle_rapmax, sk_grid, 0.2), "sk45z02"));
    output.header << "# Parameters for the SoftKiller + eroing" << endl;
    output.header << "#   grid size for SoftKiller = " << sk_grid << endl;
    output.header << "#   Zeroing radius = " << 0.2 << endl;
    output.header << "#   description     = " << event_transformers.back().description();

    declare_output(output, "sk45z02", jet_ptmin);
  }

  //......................................................................
  //   [jetByJetCS] ConstituentSubtractor [external]
  if (included_subs.find(string(",jetByJetCS,"))!=string::npos){
    double csub_alpha =  0.0;
    double csub_maxDR = -1.0;
    transformers.push_back(TransformerWithName(new contrib::ConstituentSubtractor(gmbge_neutral, 0, csub_alpha, csub_maxDR), "jetByJetCS"));
    output.header << "# Parameters for the JetByJetCSuentSubtractor" << endl;
    output.header << "#   rapidity rescaling for rho" << endl;
    output.header << "#   grid size for GMBGE        = " << grid_size << endl;
    output.header << "#   alpha                      = " << csub_alpha << endl;
    output.header << "#   maxDeltaR                  = " << csub_maxDR << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;

    declare_output(output, "jetByJetCS", jet_ptmin);
  }
  
  //......................................................................
  //   [clns]    cleansing [external]
  if (included_subs.find(string(",clns,"))!=string::npos){
    if (is_chs){
      double gamma0 = 0.612;
      double Rsub  = 0.2;
      double ftrim = 0.0;
      transformers.push_back(TransformerWithName(new TJetCleanser(Rsub,
                                                                  contrib::JetCleanser::linear_cleansing, 
                                                                  contrib::JetCleanser::input_nc_separate, 
                                                                  ftrim, 
                                                                  gamma0, 0, 0, 0,
                                                                  1.0/mixer.chs_rescaling_factor()),
                                                 "clns"));
      output.header << "# Parameters for linear cleansing [Linear mode, input = neutral and charged separate]" << endl;
      output.header << "#   gamma0 = " << gamma0 << endl;
      output.header << "#   Rsub   = " << Rsub  << endl;
      output.header << "#   ftrim  = " << ftrim << endl;

      declare_output(output, "clns", jet_ptmin);
    }
  }

  //......................................................................
  //   [puppi]   PUPPI [external]
  // copied on 2014-05-17, 18:53.
  //   git rev-list --count HEAD
  //   158
  // Modified in 2 ways:
  //  - output PJ are done as weight * particle instead of (weight*px,
  //    weight*py,...) in order to preserve the UserInfo
  if (included_subs.find(string(",puppi,"))!=string::npos){
    declare_output(output, "puppi", jet_ptmin);
  }

  
  // make sure there are no unused command-line arguments
  cmdline.assert_all_options_used();
  
  // loop over events ----------------------------------------------------------------------
  int iev = 0;
  while ( mixer.next_event() && iev < nev ) {
    // increment event number    
    iev++;
     
    // extract particles from event
    vector<PseudoJet> full_event = sel_particles(mixer.particles());
    if ( iev <= maxprintout ) { cerr << "\nEvent " << iev << endl; }
    if ( iev <= maxprintout ) { cerr << "nPU = " << mixer.npu() << endl; }
    npu.add_entry(mixer.npu());
     
    // extract hard event
    vector<PseudoJet> hard_event, pileup_event;
    SelectorIsHard().sift(full_event, hard_event, pileup_event); 
    vector<PseudoJet> pure_chs_event = ((!SelectorCharged()) || (SelectorHard()))(full_event);
    vector<PseudoJet> neutral_event = (!SelectorCharged())(full_event);
    vector<PseudoJet> neutral_pileup_event = (!SelectorCharged())(pileup_event);
    vector<PseudoJet> charged_hard_event = (SelectorCharged() && SelectorHard())(full_event);

    /*for (unsigned int i=0;i<neutral_pileup_event.size();i++){
      //cout << i << "  "  << neutral_pileup_event[i].m() << "   " << neutral_pileup_event[i].eta() << "  " << neutral_pileup_event[i].rap() << endl;
      rap_dependence->Fill(neutral_pileup_event[i].rap(),neutral_pileup_event[i].pt());
      }*/
    
    // cluster hard event only (Note: area may not be needed here)                                      
    ClusterSequence cs_hard(hard_event,jet_def);
    // cluster full event (hard + pileup)
    ClusterSequenceArea cs_full(full_event,jet_def,area_def);
          
    // get all jets in full event
    vector<PseudoJet> full_jets = sorted_by_pt(cs_full.inclusive_jets());
    if ( iev <= maxprintout ) { cerr << "Full event" << endl; }
    for (unsigned int i=0; i < 4U && i < full_jets.size(); i++) {
      if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i] << endl; }
    }
             
    // select two hardest jets in hard event
    vector<PseudoJet> hard_jets = sel_jets(sorted_by_pt(cs_hard.inclusive_jets()));
    if ( iev <= maxprintout ) { cerr << "Hard event" << endl; }
    for (unsigned int i=0; i < hard_jets.size(); i++) {
      if ( iev <= maxprintout ) { cerr << "  jet " << i << ": " << hard_jets[i] << endl; }
      output.hists["hist_pt_hard"].add_entry(hard_jets[i].pt());
      output.hists["hist_m_hard" ].add_entry(hard_jets[i].m());
      output.hists["hist_eta_hard" ].add_entry(hard_jets[i].pseudorapidity());
      output.hists["hist_width_hard" ].add_entry(width(hard_jets[i]));
    }
     

    // generic needs for subtraction: feed particles to background estimator
    gmbge->set_particles(full_event);
    gmbge_neutral->set_particles(neutral_event);
    
    //----------------------------------------------------------------------
    // Jet transformers
    //----------------------------------------------------------------------
    for (unsigned int itrans=0; itrans<transformers.size(); itrans++){
      const TransformerWithName &sub = transformers[itrans];
      // cout << "------------------------------------------------------------------------" << endl;
      // cout << "DEALING WITH: " << sub.name() << endl;
      vector<PseudoJet> subtracted_jets = sub(full_jets); 
      record(sub.name(), hard_jets, subtracted_jets, matching, output, jet_rapmax, iev <= maxprintout);
      // cout << "------------------------------------------------------------------------" << endl;
    }

    //----------------------------------------------------------------------
    // Event Transformers
    //----------------------------------------------------------------------
    for (unsigned int itrans=0; itrans<event_transformers.size(); itrans++){
      const EventTransformerWithName &sub = event_transformers[itrans];
      // cout << "------------------------------------------------------------------------" << endl;
      // cout << "DEALING WITH: " << sub.name() << endl;
      vector<PseudoJet> subtracted_jets = jet_def(sub(pure_chs_event)); 
      record(sub.name(), hard_jets, subtracted_jets, matching, output, jet_rapmax, iev <= maxprintout);
      // cout << "------------------------------------------------------------------------" << endl;
    }


    // CS whole event
    if ((included_subs.find(string(",eventWideCS,"))!=string::npos) && is_chs){
      vector<PseudoJet> subtracted=cs.subtract_event(neutral_event);
      subtracted.insert(subtracted.end(), charged_hard_event.begin(), charged_hard_event.end());
      vector<PseudoJet> subtracted_jets = jet_def(subtracted); 
      record("eventWideCS", hard_jets, subtracted_jets, matching, output, jet_rapmax, iev <= maxprintout);
    }
    // CS+SK as in ATLAS
    if ((included_subs.find(string(",CSSK,"))!=string::npos) && is_chs){
      vector<PseudoJet> subtracted_CS=cs.subtract_event(neutral_event);
      vector<PseudoJet> subtracted;
      double pt_min;
      sk.apply(subtracted_CS,subtracted,pt_min);
      subtracted.insert(subtracted.end(), charged_hard_event.begin(), charged_hard_event.end());
      vector<PseudoJet> subtracted_jets = jet_def(subtracted); 
      record("CSSK", hard_jets, subtracted_jets, matching, output, jet_rapmax, iev <= maxprintout);
    }
    // ICS whole event
    if ((included_subs.find(string(",ICS,"))!=string::npos) && is_chs){
      vector<PseudoJet> subtracted=ics.subtract_event(neutral_event);
      subtracted.insert(subtracted.end(), charged_hard_event.begin(), charged_hard_event.end());
      vector<PseudoJet> subtracted_jets = jet_def(subtracted); 
      record("ICS", hard_jets, subtracted_jets, matching, output, jet_rapmax, iev <= maxprintout);
    }
    
    //----------------------------------------------------------------------
    // PUPPI 
    //
    // They assume CHS events, so we do not run it for the full and
    // for the CHS we scale chg=PU particles back up
    // ----------------------------------------------------------------------
    if ((included_subs.find(string(",puppi,"))!=string::npos) && is_chs){
      // rescale up chg-PU
      vector<PseudoJet> rescaled_pileup_event;
      for (unsigned int i=0; i<pileup_event.size(); i++){
        double scale = SelectorCharged().pass(pileup_event[i]) ? 1.0/mixer.chs_rescaling_factor() : 1.0;
        rescaled_pileup_event.push_back(scale * pileup_event[i]);
      }

      // apply PUPPI
      puppiContainer curEvent(hard_event, rescaled_pileup_event);
      vector<PseudoJet> puppi_event = curEvent.puppiFetch(mixer.npu());

      // cluster it
      ClusterSequence cs_puppi(puppi_event,jet_def);
      
      // get all jets in full event
      vector<PseudoJet> puppi_jets = sorted_by_pt(cs_puppi.inclusive_jets());
      record("puppi", hard_jets, puppi_jets, matching, output, jet_rapmax, iev <= maxprintout);
    }      

    // cout << "------------------------------------------------------------------------" << endl;
    // cout << "All done" << endl;

    // output from time to time
    if (iev % iev_periodic == 0){
       if (iev == 15*iev_periodic) iev_periodic*=10;
       output.write(iev);
    }
    // cout << "Moving to next event" << endl;
    // cout << "------------------------------------------------------------------------" << endl;
  }  // end loop over events

  output.write(iev);

  //rap_dependence->Write();
  //file->Close();
  
  // free allocated memory
  delete rescaling;
  delete gmbge;
  delete gmbge_neutral;
}

