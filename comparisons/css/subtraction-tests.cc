///////////////////////////////////////////////////////////////////////////////////////////////////
//
// test of subtraction performance done by the FJ authors
//
// We include the following methods:
//
// TODO:
//   - add hisrograms of the observable #DONE but not tested
//   - add histograms of the delta obs  #DONE but not tested
//   - switch cleansing to Rsub=0.2
//   - add area+trimming
//   - add angularity
//   - plot efficiencies
// 
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "SimpleHist.hh"
#include "ProfileHist.hh"
#include "Matching.hh"
#include "helpers.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

#include "fastjet/tools/Filter.hh"

#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/contrib/VertexJets.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include "puppiContainer.hh"

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;
using namespace fastjet;

bool do_width;

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
    // separate the charged from the neutrqls and then split the charged
    // according to their origin (leading vertex or PU vertices)
    vector<PseudoJet> neutrals, charged, charged_lv, charged_pu;
    SelectorCharged().sift(jet.constituents(), charged, neutrals);
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

//----------------------------------------------------------------------
// declare histograms (they depend on boundaries so have to be declared)
void declare_output(OutputInfo &output, const string &tag, const double jetptmin){
  output.hists["hist_delta_pt_"   +tag].declare(-jetptmin, jetptmin, 200);
  output.hists["hist_delta_m_"    +tag].declare(-0.25*jetptmin, 0.25*jetptmin, 200);
  output.hists["hist_delta_width_"+tag].declare(-0.2, 0.2, 200.0);

  output.hists["hist_pt_"   +tag].declare(0.0, 3*jetptmin, 300);
  output.hists["hist_m_"    +tag].declare(0.0, 0.5*jetptmin, 200);
  output.hists["hist_width_"+tag].declare(0.0, 0.4, 200.0);
}

//----------------------------------------------------------------------
// record what has been going on
void record(const string &name,
            const vector<PseudoJet> &hard_jets, 
            const vector<PseudoJet> &subt_jets,
            Matching &matching,
            OutputInfo &output,
            double jet_rapmax,
            const contrib::GenericSubtractor &gen_sub,
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

      output.hists["hist_delta_pt_"+name].add_entry(match->pt() - hard_jets[i].pt());
      output.hists["hist_delta_m_" +name].add_entry(match->m()  - hard_jets[i].m() );
      
      output.hists["hist_pt_"+name].add_entry(match->pt());
      output.hists["hist_m_" +name].add_entry(match->m() );

      // handle angularity
      // 
      // This has to be done differently for areasub and safeareasub, in
      // which case we have to use a generic subtractor
      if (do_width){
        Width width;
        double width_hard = width(hard_jets[i]);
        double width_subt = ((name == "areasub") || (name == "safeareasub") || (name == "trimmed_safeareasub"))
          ? max(gen_sub(width, *match), 0.0)
          : width(*match);
      
        output.matteos["width_"+name].add_entry(width_hard, width_subt);
        output.hists["hist_width_"+name].add_entry(width_subt);
        output.hists["hist_delta_width_"+name].add_entry(width_subt - width_hard);
      }
    } else {
      output.matteos["pt_"+name].add_entry();
      output.matteos["m_" +name].add_entry();
      if (do_width)
        output.matteos["width_" +name].add_entry();
    }    
  }
  
  // keep track of number of jets above 20 GeV within the jet rapidity window
  unsigned int njets = ( SelectorPtMin(20.0) && SelectorAbsRapMax(jet_rapmax) ).count(subt_jets);
  output.matteos["pt_"+name].njets += njets;
  output.matteos["m_"+name].njets += njets;
  if (do_width)
    output.matteos["width_"+name].njets += njets;
}

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

  string included_subs = cmdline.value<string>("-methods", "unsub,area,safearea,safenpc,cleansing,constit,puppi,killer");
  included_subs.insert(0,1,',');
  included_subs.push_back(',');

  // some containier for what we want to output
  OutputInfo output(cmdline.value<string>("-out"), mixer.npu(), jet_ptmin);
  int iev_periodic=10;

  output.header << "# " << cmdline.command_line() << "\n#" << endl;
  output.header << "# time now = " << cmdline.time_stamp() << endl;
  output.header << "# Mixer: "  << mixer.description() << endl;

  do_width = cmdline.present("-width");
  output.header << "# do_width: " << do_width << endl;

  // matching condition
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.3);
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

  // select 2 hardest > 20, 100, 500 in hard event, and then |y| < 2.5
  Selector sel_jets = SelectorAbsRapMax(jet_rapmax)*SelectorPtMin(jet_ptmin)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;
  output.header << "# sel_jets: " << sel_jets.description() << endl;

  // take particles only within |y|<4
  Selector sel_particles = SelectorAbsRapMax(particle_rapmax);
  cout << "# sel_particles: " << sel_particles.description() << endl;
  output.header << "# sel_particles: " << sel_particles.description() << endl;

  AverageAndError npu;
  //ProfileHist offset_v_rapidity(-jet_rapmax,jet_rapmax,0.50);

  output.hists["hist_pt_hard"   ].declare(0.0,   3*jet_ptmin, 300);
  output.hists["hist_m_hard"    ].declare(0.0, 0.5*jet_ptmin, 200);
  if (do_width)
    output.hists["hist_width_hard"].declare(0.0, 0.4, 200);

  //----------------------------------------------------------------------
  // declare all the jet-based methods we want to try
  // 
  vector<TransformerWithName> transformers;
  bool is_chs = cmdline.present("-chs");

  //......................................................................
  // unsubtracted
  if (included_subs.find(string(",unsub,"))!=string::npos){
    transformers.push_back(TransformerWithName(new IdentityTransformer(), "unsubtracted"));
    declare_output(output, "unsubtracted", jet_ptmin);
  }

  //......................................................................
  // default FJ subtraction
  bool rescale = ! cmdline.present("-norescale");
  double grid_size = cmdline.value<double>("-gridsize", 0.55);
  
  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(particle_rapmax, grid_size);
  // define (and then apply) function for rapidity rescaling of rho.
  // NB These parameters have actually been determined at 13 TeV, but they vary slowly
  // and should therefore also be aproproate for 14 TeV
  FunctionOfPseudoJet<double> * rescaling = new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  if ( rescale) { gmbge -> set_rescaling_class(rescaling); }

  // declare a generic subtractor for shapes using areasub
  contrib::GenericSubtractor gen_sub(gmbge);
  output.header << "# GenericSubtractor: " << gen_sub.description() << endl;

  if (included_subs.find(string(",area,"))!=string::npos){
    transformers.push_back(TransformerWithName(new Subtractor(gmbge), "areasub"));

    output.header << "# Parameters for area subtraction" << endl;
    output.header << "#   rapidity rescaling for rho = " << rescale << endl;
    output.header << "#   grid size for GMBGE        = " << grid_size << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;
    cout << "# Parameters for area subtraction" << endl;
    cout << "#   rapidity rescaling for rho = " << rescale << endl;
    cout << "#   grid size for GMBGE        = " << grid_size << endl;

    declare_output(output, "areasub", jet_ptmin);
  }

  //......................................................................
  // safe area subtraction  
  if (included_subs.find(string(",safearea,"))!=string::npos){
    if (is_chs){
      transformers.push_back(TransformerWithName(new contrib::SafeAreaSubtractor(gmbge, 0, SelectorCharged(), SelectorHard()), "safeareasub"));
    } else {
      transformers.push_back(TransformerWithName(new contrib::SafeAreaSubtractor(gmbge), "safeareasub"));
    }
    output.header << "# Parameters for safe area subtraction" << endl;
    output.header << "#   rapidity rescaling for rho = " << rescale << endl;
    output.header << "#   grid size for GMBGE        = " << grid_size << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;
    cout << "# Parameters for safe area subtraction" << endl;
    cout << "#   rapidity rescaling for rho = " << rescale << endl;
    cout << "#   grid size for GMBGE        = " << grid_size << endl;

    declare_output(output, "safeareasub", jet_ptmin);
  }

  // NOTHING TO BE INSERTED HERE (we're assuming areasub+trimming follows areasub)

  double Rsub  = cmdline.value("-Rsub",  0.2);
  double ftrim = cmdline.value("-ftrim", 0.0);
  
  //......................................................................
  // safe area subtraction + trimming (automatically included from safearea
  if (included_subs.find(string(",safearea,"))!=string::npos){
    // get the subtractor from the last addition
    Transformer *safeareasub = transformers.back().base();
    assert(dynamic_cast<contrib::SafeAreaSubtractor*>(safeareasub) != 0);

    Filter *trimmer = new Filter(JetDefinition(kt_algorithm, Rsub), SelectorPtFractionMin(ftrim));
    trimmer->set_subtractor(safeareasub);
    transformers.push_back(TransformerWithName(trimmer, "trimmed_safeareasub"));

    output.header << "# Parameters for safe area subtraction + trimming" << endl;
    output.header << "#   Rsub   = " << Rsub  << endl;
    output.header << "#   ftrim  = " << ftrim << endl;
    output.header << "#   others = cf. safe ares suntraction" << endl;
    output.header << "#   Note: unsubtracted jet given as reference" << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;
    cout << "# Parameters for safe area subtraction + trimming" << endl;
    cout << "#   Rsub   = " << Rsub  << endl;
    cout << "#   ftrim  = " << ftrim << endl;
    cout << "#   others = cf. safe ares suntraction" << endl;
    cout << "#   Note: unsubtracted jet given as reference" << endl;

    declare_output(output, "trimmed_safeareasub", jet_ptmin);
  }

  //......................................................................
  // NpC (CHS-only)
  double gamma0 = cmdline.value<double>("-gamma0", 0.612);
  double gamma_with_resc = (gamma0/(1-gamma0)) * mixer.chs_rescaling_factor(); 

  if (included_subs.find(string(",safenpc,"))!=string::npos){
    if (is_chs){
      transformers.push_back(TransformerWithName(new contrib::SafeNpCSubtractor(gamma_with_resc, SelectorIsCharged(), SelectorHard()), "safenpcsub"));

      output.header << "# Parameters for neutral-proportional-to-charged (CHS-only)" << endl;
      output.header << "#   gamma0      = " << gamma0 
                    << "  [rescaled to " << gamma_with_resc << "]" << endl;
      output.header << "#   description = " << transformers.back().description() << endl;
      cout << "# Parameters for neutral-proportional-to-charged (CHS-only)" << endl;
      cout << "#   gamma0 = " << gamma0 
           << "  [rescaled to " << gamma_with_resc << "]" << endl;

      declare_output(output, "safenpcsub", jet_ptmin);
    }
  }
  
  //......................................................................
  // Constituent subtractor
  double csub_alpha = cmdline.value("constit:alpha",  0.0);
  double csub_maxDR = cmdline.value("constit:maxDR", -1.0);

  if (included_subs.find(string(",constit,"))!=string::npos){
    transformers.push_back(TransformerWithName(new contrib::ConstituentSubtractor(gmbge, 0, csub_alpha, csub_maxDR), "constitsub"));
    output.header << "# Parameters for the ConstituentSubtractor" << endl;
    output.header << "#   rapidity rescaling for rho = " << rescale << endl;
    output.header << "#   grid size for GMBGE        = " << grid_size << endl;
    output.header << "#   alpha                      = " << csub_alpha << endl;
    output.header << "#   maxDeltaR                  = " << csub_maxDR << endl;
    output.header << "#   description                = " << transformers.back().description() << endl;
    cout << "# Parameters for the ConstituentSubtractor" << endl;
    cout << "#   rapidity rescaling for rho = " << rescale << endl;
    cout << "#   grid size for GMBGE        = " << grid_size << endl;
    cout << "#   alpha                      = " << csub_alpha << endl;
    cout << "#   maxDeltaR                  = " << csub_maxDR << endl;

    declare_output(output, "constitsub", jet_ptmin);
  }

  //......................................................................
  // Cleansing (CHS-only)
  if (included_subs.find(string(",cleansing,"))!=string::npos){
    if (is_chs){
      transformers.push_back(TransformerWithName(new TJetCleanser(Rsub,
                                                                  contrib::JetCleanser::linear_cleansing, 
                                                                  contrib::JetCleanser::input_nc_separate, 
                                                                  ftrim, 
                                                                  gamma0, 0, 0, 0,
                                                                  1.0/mixer.chs_rescaling_factor()),
                                                 "linear_cleansing"));
      output.header << "# Parameters for linear cleansing [Linear mode, input = neutral and charged separate]" << endl;
      output.header << "#   gamma0 = " << gamma0 << endl;
      output.header << "#   Rsub   = " << Rsub  << endl;
      output.header << "#   ftrim  = " << ftrim << endl;
      // disabled cleansing description output because it covers multiple lines
      //output.header << "#   description = " << transformers.back().description() << endl;
      cout << "# Parameters for linear cleansing [Linear mode, input = neutral and charged separate]" << endl;
      cout << "#   gamma0      = " << gamma0 << endl;
      cout << "#   Rsub        = " << Rsub  << endl;
      cout << "#   ftrim       = " << ftrim << endl;

      declare_output(output, "linear_cleansing", jet_ptmin);
    }
  }

  // //......................................................................
  // // corrJVF: check with Pascal what to run (CHS-only)
  // //
  // // Does not really subtract the jet: it returns a jet with extra
  // // user info. If a cut is set on corrJVF (as done below), it returns
  // // an empty-PJ if the cut is not satisfied (i.e. it acts as a tagger
  // // without modifying the original jet)
  // //
  // // disable it so far???
  // double corrjvf_scale = cmdline.value("-corrjvf:scale", 0.01);
  // double corrjvf_cut   = cmdline.value("-corrjvf:cut",   0.6);
  // Selector sel_chg_pu = SelectorCharged() * (!SelectorHard());
  // Selector sel_chg_lv = SelectorCharged() * SelectorHard();
  // output.header << "# Parameters for corrJVF" << endl;
  // output.header << "#   scale factor = " << corrjvf_scale << endl;
  // output.header << "#   cut          = " << corrjvf_cut   << endl;
  // cout << "# Parameters for corrJVF" << endl;
  // cout << "#   scale factor = " << corrjvf_scale << endl;
  // cout << "#   cut          = " << corrjvf_cut   << endl;
  // contrib::VertexJets *vertexjet = 0;
  // if (is_chs){
  //   vertexjet = new contrib::VertexJets(sel_chg_pu, sel_chg_lv);
  //   transformers.push_back(TransformerWithName(vertexjet, "corrJVF"));
  //   vertexjet->set_corrJVF_scale_factor (0.01);
  //   vertexjet->set_corrJVF_cut          (0.6);
  //   // disabled cleansing description output because it covers multiple lines
  //   //output.header << "   description  = " << vertexjet->description() << endl;
  // }


  //......................................................................
  // SoftKiller
  double sk_gridsize = cmdline.value("-sk:gridsize", (is_chs) ? 0.5 : 0.4);
  output.header << "# Parameters for the SoftKiller" << endl;
  output.header << "#   grid size for SoftKiller = " << sk_gridsize << endl;
  cout << "# Parameters for the SoftKiller" << endl;
  cout << "#   grid size for SoftKiller = " << grid_size << endl;
  Selector sel_soft_killer = is_chs ? (!SelectorCharged()) : SelectorIdentity();
  contrib::SoftKiller soft_killer(particle_rapmax, sk_gridsize, sel_soft_killer);
  output.header << "#   description     = " << soft_killer.description();

  if (included_subs.find(string(",killer,"))!=string::npos){
    declare_output(output, "soft_killer", jet_ptmin);
  }
  

  //......................................................................
  // PUPPI --- copied on 2014-05-17, 18:53.
  //   git rev-list --count HEAD
  //   105
  // Modified in 2 ways:
  //  - commented out the 'cout << "pass..."' lines
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

    // cluster hard event only (Note: area may not be needed here)                                      
    ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
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
      if (do_width)
        output.hists["hist_width_hard"].add_entry(Width()(hard_jets[i]));
    }
     

    // generic needs for subtraction
    // 1. feed particles to background estimator
    gmbge->set_particles(full_event);

    // // 2. set the number of charged PU tracks fot corrJVF
    // if (vertexjet){
    //   int n_pileup_tracks = sel_chg_pu(full_event).size(); 
    //   vertexjet->set_tot_n_pu_tracks(n_pileup_tracks);
    // }

    // loop over the subtraction methods
    for (unsigned int itrans=0; itrans<transformers.size(); itrans++){
      const TransformerWithName &sub = transformers[itrans];
      
      // subtract
      vector<PseudoJet> subtracted_jets = sub(full_jets); 
      record(sub.name(), hard_jets, subtracted_jets, matching, output, jet_rapmax, gen_sub, iev <= maxprintout);
    }

    //----------------------------------------------------------------------
    // SoftKiller
    //----------------------------------------------------------------------

    if (included_subs.find(string(",killer,"))!=string::npos){
      // apply the SoftKiller
      vector<PseudoJet> sk_event = soft_killer(full_event);
      
      // cluster it
      ClusterSequence cs_sk(sk_event,jet_def);
      
      // get all jets in full event
      vector<PseudoJet> sk_jets = sorted_by_pt(cs_sk.inclusive_jets());
      record("soft_killer", hard_jets, sk_jets, matching, output, jet_rapmax, gen_sub, iev <= maxprintout);
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

      // // check if there is any chg-PU particle coming out
      // cout << "#PUPPI event size = " << puppi_event.size() << endl;
      // cout << "#PUPPI Chg-PU count = " << (SelectorCharged()*(!SelectorHard())).count(puppi_event) << endl;
      // cout << "#PUPPI PU size = " << (!SelectorHard()).count(puppi_event) << endl;
      // if (puppi_event[0].has_user_info<PU14>())
      //   cout << "PUPPI has info" << endl;
      // else
      //   cout << "PUPPI does not have info" << endl;

      // cluster it
      ClusterSequence cs_puppi(puppi_event,jet_def);
          
      // get all jets in full event
      vector<PseudoJet> puppi_jets = sorted_by_pt(cs_puppi.inclusive_jets());
      record("puppi", hard_jets, puppi_jets, matching, output, jet_rapmax, gen_sub, iev <= maxprintout);
    }      

    // output from time to time
    if (iev % iev_periodic == 0){
       if (iev == 15*iev_periodic) iev_periodic*=10;
       output.write(iev);
    }
  }  // end loop over events

  output.write(iev);

  // free allocated memory
  delete rescaling;
  delete gmbge;
}

