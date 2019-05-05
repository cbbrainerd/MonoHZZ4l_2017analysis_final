#ifndef HZZ4LeptonsCommonRootTree_h
#define HZZ4LeptonsCommonRootTree_h

/** \class  HZZ4LeptonsCommonRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *          
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/FileBlock.h"

// MCTruth
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//@
//To retrieve LHE info
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
//To retrieve JEC Uncertainty
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//To compute Pt calibration
//#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

// Data format
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"    
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
//#include "DataFormats/METReco/interface/PFMET.h"
//#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"


// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// user include files
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//PU Jet ID
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

//Full Error
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

#include <TMatrixD.h>
#include <TLorentzVector.h>

//for muon Rochester correction

#include "roccor_Run2_v2/RoccoR.h"
#include "TRandom3.h"

//chisquare
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//ArrayVector hack
#include "HiggsAnalysis/HiggsToZZ4Leptons/interface/ArrayVector.h"

class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

// Namespaces
using namespace reco;
using namespace std;
using namespace pat;




class HZZ4LeptonsCommonRootTree : public edm::EDAnalyzer {
  
 public:
  
  HZZ4LeptonsCommonRootTree(const edm::ParameterSet& pset);
  
  ~HZZ4LeptonsCommonRootTree();
  
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void beginJob();
  void endJob();

  std::unique_ptr<RoccoR> calibrator; //muon calibrator 
 
  void respondToOpenInputFile(edm::FileBlock const& fb) {
    inputfileName = fb.fileName();
    cout << "Input Filename is=" << inputfileName.c_str() << endl;
  }
	
  void ReadParameters(const edm::ParameterSet& pset){ 
    std::cout << "Reading parameters from cfg" << std::endl;
    
    typedef std::vector<edm::InputTag> vtag;
    // Get the various input parameters
    decaychannel              = pset.getParameter<std::string>("decaychannel");
    rootFileName              = pset.getUntrackedParameter<std::string>("rootFileName");
    useRECOformat             = pset.getUntrackedParameter<bool>("useRECOformat");

    // search param in cfg
    module_to_search=pset.getUntrackedParameter<vector<std::string> >("module_to_search");
    par_to_search= pset.getUntrackedParameter<std::string>("par_to_search");
    

    // Get PU simulation info
    fillPUinfo                = pset.getUntrackedParameter<bool>("fillPUinfo");
    PileupSrc_                = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("PileupSrc"));

    // Generator
    generator_                = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator"));

    //@ for LHE informations for Jets
    LHE_                      = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag> ("LHEProduct"));

    // Get HLT flags
    fillHLTinfo               = pset.getUntrackedParameter<bool>("fillHLTinfo");
    HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
    HLTAnalysisinst           = pset.getParameter<string>("HLTAnalysisinst");
    flagHLTnames              = pset.getParameter<vtag>("flagHLTnames");
    // Get HLT matching
//AOD    triggerEvent              = consumes<trigger::TriggerEvent >(pset.getParameter<edm::InputTag>("triggerEvent"));
    triggerObjects_              = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("triggerobjects"));

    triggerBits_              = consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerbits"));

    //Reham 
    noiseFilterTag_          =consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("noiseFilterTag"));
    HBHENoiseFilter_Selector_ = pset.getParameter<string>("HBHENoiseFilter_Selector_");
    EEBadScNoiseFilter_Selector_ = pset.getParameter<string>("EEBadScNoiseFilter_Selector_");
    GoodVtxNoiseFilter_Selector_ =  pset.getParameter<string>("GoodVtxNoiseFilter_Selector_");
    GlobalSuperTightHalo2016NoiseFilter_Selector_  =  pset.getParameter<string>("GlobalSuperTightHalo2016NoiseFilter_Selector_");
    HBHENoiseIsoFilter_Selector_ = pset.getParameter<string>("HBHENoiseIsoFilter_Selector_");
    EcalDeadCellTriggerPrimitiveNoiseFilter_Selector_ = pset.getParameter<string>("EcalDeadCellTriggerPrimitiveNoiseFilter_Selector_");
    BadPFMuonFilter_Selector_ = pset.getParameter<string>("BadPFMuonFilter_Selector_");
    BadChargedCandidateFilter_Selector_ = pset.getParameter<string>("BadChargedCandidateFilter_Selector_");                          
    EcalBadCalibFilter_Selector_ = pset.getParameter<string>("EcalBadCalibFilter_Selector_");
   
    triggerFilter             = pset.getParameter<std::string>("triggerFilter");
    triggerMatchObject        = consumes<edm::Association<std::vector<pat::TriggerObjectStandAlone> > >(pset.getParameter<edm::InputTag>("triggerMatchObject"));
    triggerMatchObjectEle     = pset.getParameter<edm::InputTag>("triggerMatchObjectEle");
    triggerHLTcollection      = pset.getParameter<std::string>("triggerHLTcollection");

    // Get flags
    flaginst                  = pset.getParameter<std::string>("flaginst");
    flagtags                  = pset.getParameter<std::vector<std::string> >("flagtags");
    // MCtruth tags
    fillMCTruth               = pset.getUntrackedParameter<bool>("fillMCTruth");
    MCcollName                = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("MCcollName"));
    genParticles_             = consumes<vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genParticles"));
    fourgenleptons_           = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("fourgenleptons"));
    digenZ_                   = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("digenZ"));

    
    // RECO tags
    RECOcollNameBest2e2mu     = pset.getParameter<vtag>("RECOcollNameBest2e2mu");
    RECOcollNameBest4mu       = pset.getParameter<vtag>("RECOcollNameBest4mu");
    RECOcollNameBest4e        = pset.getParameter<vtag>("RECOcollNameBest4e");
 
    // RECO additional tags
    useAdditionalRECO         = pset.getUntrackedParameter<bool>("useAdditionalRECO");
    RECOcollNameZ             = pset.getParameter<vtag>("RECOcollNameZ");
    RECOcollNameZss           = pset.getParameter<vtag>("RECOcollNameZss");
    RECOcollNameDiLep_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameDiLep"));

    RECOcollNameDiLep         = pset.getParameter<edm::InputTag>("RECOcollNameDiLep");
    RECOcollNameMMMM_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameMMMM"));
    RECOcollNameEEEE          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEEE")); 
    RECOcollNameEEMM          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEMM"));
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
    RECOcollNameLLLL         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameLLLL"));   

//    RECOcollNameMMMM          = pset.getParameter<vtag>("RECOcollNameMMMM");
//    RECOcollNameEEEE          = pset.getParameter<vtag>("RECOcollNameEEEE");
//    RECOcollNameEEMM          = pset.getParameter<vtag>("RECOcollNameEEMM");
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
//    RECOcollNameLLLL          = pset.getParameter<edm::InputTag>("RECOcollNameLLLL");

    // electrons and muons tags
    use2011EA                 = pset.getUntrackedParameter<bool>("use2011EA");
    muonTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
    muonCorrPtErrorMapTag_    = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsCorrPtErrorMapLabel"));
    
    //@
    isData                    = pset.getParameter<bool>("isData");
     slimmedMuonsTag_          = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("SlimmedMuonsLabel"));
     
     muonPFTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("PFMuonsLabel"));
     
     
     clusterCollectionTag_    = pset.getParameter<edm::InputTag>("SuperClustersLabel");
     gsftrackCollection_      = pset.getParameter<edm::InputTag>("GsfTracksElectronsLabel");
     
     
     electronEgmTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("ElectronsEgmLabel"));
     mvaElectronTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("mvaElectronTag"));
    ///@@@/// TEST Electron ID  REHAM
    //@//    mvaTrigV0MapTag_         = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaTrigV0MapTag"));
    //@// mvaNonTrigV0MapTag_      = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaNonTrigV0MapTag"));
    electronsMiniAODToken_    = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronsMiniAOD")); //REHAM
    eleIDToken_ = consumes<edm::ValueMap<bool>>(pset.getParameter<edm::InputTag>("elecID")); //REHAM 
    //   mvaValuesMapToken_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaValuesMap"));
    //  mvaCategoriesMapToken_= consumes<edm::ValueMap<int> >(pset.getParameter<edm::InputTag>("mvaCategoriesMap"));   
    
    // PF photons
//    pfphotonsTag_                 = consumes<edm::View<reco::PFCandidate>>(pset.getParameter<edm::InputTag>("PFPhotonsLabel"));
    pfTag_                 = consumes<edm::View<pat::PackedCandidate>>(pset.getParameter<edm::InputTag>("pfCands"));

    // vertexing 
    // 3D w.r.t primary vertex DA
    muonTag_Vert           = pset.getParameter<edm::InputTag>("MuonsLabelVert");
    muonMapTag_Vert        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVert"));
    muonMapTag_VertValue   = consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsMapLabelVertValue"));
    muonMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertError"));
    // KF
    muonMapTag_VertKF        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertKF"));
    muonMapTag_VertValueKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertValueKF"));
    muonMapTag_VertErrorKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertErrorKF"));
    
    
    // GD, Std and Kin vertex
    muonMapTag_VertGD      = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGD");
    muonMapTag_VertGDMMMM  = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGDMMMM");
    muonMapTag_VertStd     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStd");
    muonMapTag_VertStdMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStdMMMM");
    muonMapTag_VertKin     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKin");
    muonMapTag_VertKinMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKinMMMM");
    
    // STIP SLIP
    muonSTIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVert"));
    muonSLIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVert"));
    
    muonSTIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertValue"));
    muonSLIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertValue"));
    muonSTIPMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertError"));
    muonSLIPMapTag_VertError   =consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertError"));

    //electronID
    eleIDTag_                 = pset.getParameter<vtag>("eleIDLabel");
    
    // Other objets
    photonsTag_              = consumes<edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("PhotonsLabel"));
    //    tracksTag_               = consumes<vector<reco::Track> >(pset.getParameter<edm::InputTag>("TracksLabel"));
    jetsTag_                 = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsLabel")); 
    //   jetsDataTag_             = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsDataLabel"));
    jetsMVATag_              = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsMVALabel"));
    jetsQGTag_               = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsLabel"));//REHAM
    jetQGMapTag_             = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("jetQGMapTag"));//REHAM
    jetQGMapTag_axis2_             = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("jetQGMapTag_axis2"));//REHAM
    jetQGMapTag_ptd_             = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("jetQGMapTag_ptd"));//REHAM
    jetQGMapTag_mult_             = consumes<edm::ValueMap<int> >(pset.getParameter<edm::InputTag>("jetQGMapTag_mult"));//REHAM
    
    
    rhojetsTag_              = consumes<double>(pset.getParameter<edm::InputTag>("RhoJetsLabel"));
    verticesTag_             = consumes<vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VerticesLabel"));
    
    // GenJet 
    genjetTag_              = consumes<vector<reco::GenJet> >(pset.getParameter<edm::InputTag>("GenJetLabel")); // GenJet
    
    // PF MET
    pfmetTag_               = consumes<vector<pat::MET> >(pset.getParameter<edm::InputTag>("PfMETLabel"));
    // ecalBadCalibFilterUpdate_token   = consumes<bool>(edm::InputTag("ecalBadCalibReducedMINIAODFilter")); new met filter to be used (under test)
    // BadChCandFilterToken_            = consumes<bool>(edm::InputTag("BadChargedCandidateFilter"));
    // BadPFMuonFilterToken_            = consumes<bool>(edm::InputTag("BadPFMuonFilter"));
      
    
    // btagging
    //   bDiscriminators_         =pset.getParameter<std::vector<std::string> >("bDiscriminators");
    
    tCHighEff_bTag_         =pset.getParameter<std::string>("tCHighEff_bTagLabel");
    tCHighPur_bTag_         =pset.getParameter<std::string>("tCHighPur_bTagLabel");
    cSV_bTag1_         =pset.getParameter<std::string>("cSV_bTagLabel1");
    cSV_bTag2_         =pset.getParameter<std::string>("cSV_bTagLabel2");

    
    // Matching
    // /* to add matching informations Reham

    goodElectronMCMatch_  = consumes<edm::Association<std::vector<reco::GenParticle> >>(pset.getParameter<edm::InputTag>("goodElectronMCMatch"));
    myElectrons_          = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myElectrons"));
    goodMuonMCMatch_      = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodMuonMCMatch")); 
    myMuons_              = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myMuons"));
    goodGammaMCMatch_     = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodGammaMCMatch"));
    myGammas_             = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myGammas"));
    //  */

    // Beam Spot
    offlineBeamSpot_      = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("offlineBeamSpot"));
  }
  
  
  void DefineBranches(TTree *Tree_){
    // Run event
    Tree_->Branch("Run",&irun,"irun/i");
    Tree_->Branch("Event",&ievt,"ievt/i");
    Tree_->Branch("LumiSection",&ils,"ils/i");
    Tree_->Branch("Avginstlumi",&Avginstlumi,"Avginstlumi/F");
    
    // MC
    Tree_->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
    Tree_->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
    Tree_->Branch("MC_weighting",&MC_weighting,"MC_weighting/F");

    // HLT 
    Tree_->Branch("RECO_nMuHLTMatch",&RECO_nMuHLTMatch,"RECO_nMuHLTMatch/I");
    Tree_->Branch("RECOMU_PT_MuHLTMatch",(std::vector<float >*)(&RECOMU_PT_MuHLTMatch));//"std::vector<float>"); TEST
    Tree_->Branch("RECOMU_ETA_MuHLTMatch",(std::vector<float >*)(&RECOMU_ETA_MuHLTMatch));
    Tree_->Branch("RECOMU_PHI_MuHLTMatch",(std::vector<float >*)(&RECOMU_PHI_MuHLTMatch));

    Tree_->Branch("RECO_nEleHLTMatch",&RECO_nEleHLTMatch,"RECO_nEleHLTMatch/I");
    Tree_->Branch("RECOELE_PT_EleHLTMatch",(std::vector<float >*)(&RECOELE_PT_EleHLTMatch));
    Tree_->Branch("RECOELE_ETA_EleHLTMatch",(std::vector<float >*)(&RECOELE_ETA_EleHLTMatch));
    Tree_->Branch("RECOELE_PHI_EleHLTMatch",(std::vector<float >*)(&RECOELE_PHI_EleHLTMatch));

    Tree_->Branch("HLTPathsFired",HLTPathsFired,"HLTPathsFired/C");

   
    // MC block 
    Tree_->Branch("MC_E",(std::vector<float >*)(&MC_E)); 
    Tree_->Branch("MC_PT",(std::vector<float >*)(&MC_PT)); 
    Tree_->Branch("MC_ETA",(std::vector<float >*)(&MC_ETA)); 
    Tree_->Branch("MC_THETA",(std::vector<float >*)(&MC_THETA));
    Tree_->Branch("MC_PHI",(std::vector<float >*)(&MC_PHI));
    Tree_->Branch("MC_MASS",(std::vector<float >*)(&MC_MASS));
    Tree_->Branch("MC_PDGID",(std::vector<float >*)(&MC_PDGID));
   

    //  Leptons stable ordered in PT
    Tree_->Branch("MC_LEPT_PT",(std::vector<float >*)(&MC_LEPT_PT)); 
    Tree_->Branch("MC_LEPT_ETA",(std::vector<float >*)(&MC_LEPT_ETA)); 
    Tree_->Branch("MC_LEPT_PHI",(std::vector<float >*)(&MC_LEPT_PHI)); 
    Tree_->Branch("MC_LEPT_THETA",(std::vector<float >*)(&MC_LEPT_THETA)); 
    Tree_->Branch("MC_LEPT_PDGID",(std::vector<float >*)(&MC_LEPT_PDGID));

    // MC Z1, Z2 (first index) and daughter leptons and photon (second index, order: L1, L2, P1, P2 with PT ordering of the leptons and P1 associated to L1)
    Tree_->Branch("MC_Z_PT",(std::vector<std::vector<float > >*)(&MC_Z_PT));
    Tree_->Branch("MC_Z_ETA",(std::vector<std::vector<float > >*)(&MC_Z_ETA));
    Tree_->Branch("MC_Z_PHI",(std::vector<std::vector<float > >*)(&MC_Z_PHI));
    Tree_->Branch("MC_Z_THETA",(std::vector<std::vector<float > >*)(&MC_Z_THETA));
    Tree_->Branch("MC_Z_MASS",(std::vector<std::vector<float > >*)(&MC_Z_MASS));
    Tree_->Branch("MC_Z_PDGID",(std::vector<std::vector<float > >*)(&MC_Z_PDGID));

    // 4l from stable particles
    Tree_->Branch("MC_fourl_MASS",(std::vector<std::vector<float > >*)(&MC_fourl_MASS));
    Tree_->Branch("MC_fourl_PT",(std::vector<std::vector<float > >*)(&MC_fourl_PT));
    Tree_->Branch("MC_fourl_PDGID",(std::vector<std::vector<float > >*)(&MC_fourl_PDGID));

    // diZ
    Tree_->Branch("MC_ZZ_MASS",(std::vector<std::vector<float > >*)(&MC_ZZ_MASS));
    Tree_->Branch("MC_ZZ_PT",(std::vector<std::vector<float > >*)(&MC_ZZ_PT));
    Tree_->Branch("MC_ZZ_ETA",(std::vector<std::vector<float > >*)(&MC_ZZ_ETA));
    Tree_->Branch("MC_ZZ_PHI",(std::vector<std::vector<float > >*)(&MC_ZZ_PHI));
    Tree_->Branch("MC_ZZ_THETA",(std::vector<std::vector<float > >*)(&MC_ZZ_THETA));
    Tree_->Branch("MC_ZZ_PDGID",(std::vector<std::vector<float > >*)(&MC_ZZ_PDGID));


    // GenJet
    Tree_->Branch( "MC_GENJET_PT",(std::vector<float >*)(&MC_GENJET_PT));
    Tree_->Branch( "MC_GENJET_ETA",(std::vector<float >*)(&MC_GENJET_ETA));
    Tree_->Branch( "MC_GENJET_PHI",(std::vector<float >*)(&MC_GENJET_PHI));

    // GenMET  
    Tree_->Branch("MC_GENMET", &genmet, "MC_GENMET/F");
    
  
     // RECORF block 2e2mu
    
    Tree_->Branch("RECORF_2e2mu_cosTheta1_spin",(std::vector<double >*)(&RECORF_2e2mu_cosTheta1_spin));
    Tree_->Branch("RECORF_2e2mu_cosTheta2_spin",(std::vector<double >*)(&RECORF_2e2mu_cosTheta2_spin));
    Tree_->Branch("RECORF_2e2mu_cosThetaStar_spin",(std::vector<double >*)(&RECORF_2e2mu_cosThetaStar_spin));
    Tree_->Branch("RECORF_2e2mu_Phi_spin",(std::vector<double >*)(&RECORF_2e2mu_Phi_spin));
    Tree_->Branch("RECORF_2e2mu_Phi1_spin",(std::vector<double >*)(&RECORF_2e2mu_Phi1_spin));
    Tree_->Branch("RECORF_2e2mu_Phi2_spin",(std::vector<double >*)(&RECORF_2e2mu_Phi2_spin));
    Tree_->Branch("RECORF_2e2mu_phi1RF_spin",(std::vector<double >*)(&RECORF_2e2mu_phi1RF_spin));
    Tree_->Branch("RECORF_2e2mu_phi2RF_spin",(std::vector<double >*)(&RECORF_2e2mu_phi2RF_spin));
    Tree_->Branch("RECORF_2e2mu_MELA",(std::vector<double >*)(&RECORF_2e2mu_MELA));

   
    Tree_->Branch("RECORF_4e_cosTheta1_spin",(std::vector<double >*)(&RECORF_4e_cosTheta1_spin));
    Tree_->Branch("RECORF_4e_cosTheta2_spin",(std::vector<double >*)(&RECORF_4e_cosTheta2_spin));
    Tree_->Branch("RECORF_4e_cosThetaStar_spin",(std::vector<double >*)(&RECORF_4e_cosThetaStar_spin));
    Tree_->Branch("RECORF_4e_Phi_spin",(std::vector<double >*)(&RECORF_4e_Phi_spin));
    Tree_->Branch("RECORF_4e_Phi1_spin",(std::vector<double >*)(&RECORF_4e_Phi1_spin));
    Tree_->Branch("RECORF_4e_Phi2_spin",(std::vector<double >*)(&RECORF_4e_Phi2_spin));
    Tree_->Branch("RECORF_4e_phi1RF_spin",(std::vector<double >*)(&RECORF_4e_phi1RF_spin));
    Tree_->Branch("RECORF_4e_phi2RF_spin",(std::vector<double >*)(&RECORF_4e_phi2RF_spin));
    Tree_->Branch("RECORF_4e_MELA",(std::vector<double >*)(&RECORF_4e_MELA));
     
    Tree_->Branch("RECORF_4mu_cosTheta1_spin",(std::vector<double >*)(&RECORF_4mu_cosTheta1_spin));
    Tree_->Branch("RECORF_4mu_cosTheta2_spin",(std::vector<double >*)(&RECORF_4mu_cosTheta2_spin));
    Tree_->Branch("RECORF_4mu_cosThetaStar_spin",(std::vector<double >*)(&RECORF_4mu_cosThetaStar_spin));
    Tree_->Branch("RECORF_4mu_Phi_spin",(std::vector<double >*)(&RECORF_4mu_Phi_spin));
    Tree_->Branch("RECORF_4mu_Phi1_spin",(std::vector<double >*)(&RECORF_4mu_Phi1_spin));
    Tree_->Branch("RECORF_4mu_Phi2_spin",(std::vector<double >*)(&RECORF_4mu_Phi2_spin));
    Tree_->Branch("RECORF_4mu_phi1RF_spin",(std::vector<double >*)(&RECORF_4mu_phi1RF_spin));
    Tree_->Branch("RECORF_4mu_phi2RF_spin",(std::vector<double >*)(&RECORF_4mu_phi2RF_spin));
    Tree_->Branch("RECORF_4mu_MELA",(std::vector<double >*)(&RECORF_4mu_MELA));



    // RECO additional block for reconstructed higgs, Z and their daughters
    Tree_->Branch("RECO_ZMM_MASS",(std::vector<float >*)(&RECO_ZMM_MASS));
    Tree_->Branch("RECO_ZEE_MASS",(std::vector<float >*)(&RECO_ZEE_MASS));
    Tree_->Branch("RECO_DiLep_MASS",(std::vector<float >*)(&RECO_DiLep_MASS));
    Tree_->Branch("RECO_ZMM_PT",(std::vector<std::vector<float > >*)(&RECO_ZMM_PT));
    Tree_->Branch("RECO_ZEE_PT",(std::vector<std::vector<float > >*)(&RECO_ZEE_PT));  
    Tree_->Branch("RECO_DiLep_PT",(std::vector<std::vector<float > >*)(&RECO_DiLep_PT));  
    Tree_->Branch("RECO_ZMM_ETA",(std::vector<std::vector<float > >*)(&RECO_ZMM_ETA));
    Tree_->Branch("RECO_ZEE_ETA",(std::vector<std::vector<float > >*)(&RECO_ZEE_ETA));
    Tree_->Branch("RECO_DiLep_ETA",(std::vector<std::vector<float > >*)(&RECO_DiLep_ETA));  
    Tree_->Branch("RECO_ZMM_PHI",(std::vector<std::vector<float > >*)(&RECO_ZMM_PHI));
    Tree_->Branch("RECO_ZEE_PHI",(std::vector<std::vector<float > >*)(&RECO_ZEE_PHI));
    Tree_->Branch("RECO_DiLep_PHI",(std::vector<std::vector<float > >*)(&RECO_DiLep_PHI));  
                              
    Tree_->Branch("RECO_ZMMss_MASS",(std::vector<float >*)(&RECO_ZMMss_MASS));
    Tree_->Branch("RECO_ZEEss_MASS",(std::vector<float >*)(&RECO_ZEEss_MASS));
    Tree_->Branch("RECO_ZEM_MASS",(std::vector<float >*)(&RECO_ZEM_MASS));
    Tree_->Branch("RECO_ZMMss_PT",(std::vector<std::vector<float > >*)(&RECO_ZMMss_PT));
    Tree_->Branch("RECO_ZEEss_PT",(std::vector<std::vector<float > >*)(&RECO_ZEEss_PT));
    Tree_->Branch("RECO_ZEM_PT",(std::vector<std::vector<float > >*)(&RECO_ZEM_PT));
    Tree_->Branch("RECO_ZMMss_ETA",(std::vector<std::vector<float > >*)(&RECO_ZMMss_ETA));
    Tree_->Branch("RECO_ZEEss_ETA",(std::vector<std::vector<float > >*)(&RECO_ZEEss_ETA));
    Tree_->Branch("RECO_ZEM_ETA",(std::vector<std::vector<float > >*)(&RECO_ZEM_ETA));
    Tree_->Branch("RECO_ZMMss_PHI",(std::vector<std::vector<float > >*)(&RECO_ZMMss_PHI));
    Tree_->Branch("RECO_ZEEss_PHI",(std::vector<std::vector<float > >*)(&RECO_ZEEss_PHI));
    Tree_->Branch("RECO_ZEM_PHI",(std::vector<std::vector<float > >*)(&RECO_ZEM_PHI));

    
    Tree_->Branch("RECO_MMMM_MASS",(std::vector<std::vector<float > >*)(&RECO_MMMM_MASS));
    Tree_->Branch("RECO_MMMM_PT",(std::vector<std::vector<float > >*)(&RECO_MMMM_PT));
    Tree_->Branch("RECO_MMMM_ETA",(std::vector<std::vector<float > >*)(&RECO_MMMM_ETA));
    Tree_->Branch("RECO_MMMM_PHI",(std::vector<std::vector<float > >*)(&RECO_MMMM_PHI));
    Tree_->Branch("RECO_MMMM_MASS_REFIT",(std::vector<float >*)(&RECO_MMMM_MASS_REFIT));

    Tree_->Branch("RECO_EEEE_MASS",(std::vector<std::vector<float > >*)(&RECO_EEEE_MASS));
    Tree_->Branch("RECO_EEEE_PT",(std::vector<std::vector<float > >*)(&RECO_EEEE_PT)); 
    Tree_->Branch("RECO_EEEE_ETA",(std::vector<std::vector<float > >*)(&RECO_EEEE_ETA));
    Tree_->Branch("RECO_EEEE_PHI",(std::vector<std::vector<float > >*)(&RECO_EEEE_PHI));
    Tree_->Branch("RECO_EEEE_MASS_REFIT",(std::vector<float >*)(&RECO_EEEE_MASS_REFIT));

    Tree_->Branch("RECO_EEMM_MASS",(std::vector<std::vector<float > >*)(&RECO_EEMM_MASS));
    Tree_->Branch("RECO_EEMM_PT",(std::vector<std::vector<float > >*)(&RECO_EEMM_PT));
    Tree_->Branch("RECO_EEMM_ETA",(std::vector<std::vector<float > >*)(&RECO_EEMM_ETA));
    Tree_->Branch("RECO_EEMM_PHI",(std::vector<std::vector<float > >*)(&RECO_EEMM_PHI));
    Tree_->Branch("RECO_EEMM_MASS_REFIT",(std::vector<float >*)(&RECO_EEMM_MASS_REFIT));

    Tree_->Branch("RECO_LLL0_MASS",(std::vector<float >*)(&RECO_LLL0_MASS)); 
    Tree_->Branch("RECO_LLL1_MASS",(std::vector<float >*)(&RECO_LLL1_MASS)); 
    Tree_->Branch("RECO_LLL2_MASS",(std::vector<float >*)(&RECO_LLL2_MASS)); 
    Tree_->Branch("RECO_LLL3_MASS",(std::vector<float >*)(&RECO_LLL3_MASS)); 
    Tree_->Branch("RECO_LLL0_PT",(std::vector<std::vector<float > >*)(&RECO_LLL0_PT)); 
    Tree_->Branch("RECO_LLL1_PT",(std::vector<std::vector<float > >*)(&RECO_LLL1_PT)); 
    Tree_->Branch("RECO_LLL2_PT",(std::vector<std::vector<float > >*)(&RECO_LLL2_PT)); 
    Tree_->Branch("RECO_LLL3_PT",(std::vector<std::vector<float > >*)(&RECO_LLL3_PT)); 

    Tree_->Branch("RECO_LLLl0_MASS",(std::vector<float >*)(&RECO_LLLl0_MASS)); 
    Tree_->Branch("RECO_LLLl1_MASS",(std::vector<float >*)(&RECO_LLLl1_MASS)); 
    Tree_->Branch("RECO_LLLl0_PT",(std::vector<std::vector<float > >*)(&RECO_LLLl0_PT)); 
    Tree_->Branch("RECO_LLLl1_PT",(std::vector<std::vector<float > >*)(&RECO_LLLl1_PT)); 

    Tree_->Branch("RECO_LLLL0ss_MASS",(std::vector<float >*)(&RECO_LLLL0ss_MASS)); 
    Tree_->Branch("RECO_LLLL0ss_PT",(std::vector<std::vector<float > >*)(&RECO_LLLL0ss_PT)); 
    Tree_->Branch("RECO_LLLL1ss_MASS",(std::vector<float >*)(&RECO_LLLL1ss_MASS)); 
    Tree_->Branch("RECO_LLLL1ss_PT",(std::vector<std::vector<float > >*)(&RECO_LLLL1ss_PT)); 
    Tree_->Branch("RECO_LLLL2ss_MASS",(std::vector<float >*)(&RECO_LLLL2ss_MASS)); 
    Tree_->Branch("RECO_LLLL2ss_PT",(std::vector<std::vector<float > >*)(&RECO_LLLL2ss_PT)); 
       
    //Tree_->Branch("RECOcollNameLLLLssos_MASS",(std::vector<float >*)(&RECOcollNameLLLLssos_MASS));
    //Tree_->Branch("RECOcollNameLLLLssos_PT",(std::vector<std::vector<float > >*)(&RECOcollNameLLLLssos_PT));

    Tree_->Branch("RECO_LLLL_MASS",(std::vector<std::vector<float > >*)(&RECO_LLLL_MASS));
    Tree_->Branch("RECO_LLLL_PT",(std::vector<std::vector<float > >*)(&RECO_LLLL_PT));
    Tree_->Branch("RECO_LLLL_ETA",(std::vector<std::vector<float > >*)(&RECO_LLLL_ETA));
    Tree_->Branch("RECO_LLLL_PHI",(std::vector<std::vector<float > >*)(&RECO_LLLL_PHI));
 
   
    // Electron block
    Tree_->Branch("RECOELE_E",(std::vector<float >*)(&RECOELE_E)); 
    Tree_->Branch("RECOELE_PT",(std::vector<float >*)(&RECOELE_PT));
    Tree_->Branch("RECOELE_PTError",(std::vector<float >*)(&RECOELE_PTError));
    Tree_->Branch("RECOELE_P",(std::vector<float >*)(&RECOELE_P));
    Tree_->Branch("RECOELE_ETA",(std::vector<float >*)(&RECOELE_ETA)); 
    Tree_->Branch("RECOELE_THETA",(std::vector<float >*)(&RECOELE_THETA)); 
    Tree_->Branch("RECOELE_PHI",(std::vector<float >*)(&RECOELE_PHI)); 
    Tree_->Branch("RECOELE_MASS",(std::vector<float >*)(&RECOELE_MASS)); 
    Tree_->Branch("RECOELE_CHARGE",(std::vector<float >*)(&RECOELE_CHARGE)); 
  
    //Reham
    Tree_->Branch("RECOELE_ID",(std::vector<float >*)(&RECOELE_ID));
    Tree_->Branch("RECOELE_PT_uncorr",(std::vector<float >*)(&RECOELE_PT_uncorr));

    // Core attributes
    Tree_->Branch("RECOELE_isEcalDriven",(std::vector<int >*)(&RECOELE_isEcalDriven));   
    Tree_->Branch("RECOELE_isTrackerDriven",(std::vector<int >*)(&RECOELE_isTrackerDriven));  
    Tree_->Branch("RECOELE_gsftrack_NPixHits",(std::vector<float >*)(&RECOELE_gsftrack_NPixHits));
    Tree_->Branch("RECOELE_gsftrack_NStripHits",(std::vector<float >*)(&RECOELE_gsftrack_NStripHits));
    Tree_->Branch("RECOELE_gsftrack_chi2",(std::vector<float >*)(&RECOELE_gsftrack_chi2));
    Tree_->Branch("RECOELE_gsftrack_dxyB",(std::vector<float >*)(&RECOELE_gsftrack_dxyB));
    Tree_->Branch("RECOELE_gsftrack_dxy",(std::vector<float >*)(&RECOELE_gsftrack_dxy));
    Tree_->Branch("RECOELE_gsftrack_dxyError",(std::vector<float >*)(&RECOELE_gsftrack_dxyError));
    Tree_->Branch("RECOELE_gsftrack_dzB",(std::vector<float >*)(&RECOELE_gsftrack_dzB));
    Tree_->Branch("RECOELE_gsftrack_dz",(std::vector<float >*)(&RECOELE_gsftrack_dz));
    Tree_->Branch("RECOELE_gsftrack_dzError",(std::vector<float >*)(&RECOELE_gsftrack_dzError));
    Tree_->Branch("RECOELE_gsftrack_losthits",(std::vector<int >*)(&RECOELE_gsftrack_losthits));
    Tree_->Branch("RECOELE_gsftrack_validhits",(std::vector<int >*)(&RECOELE_gsftrack_validhits));
    Tree_->Branch("RECOELE_gsftrack_expected_inner_hits",(std::vector<int >*)(&RECOELE_gsftrack_expected_inner_hits)) ; 
    Tree_->Branch("RECOELE_scl_E",(std::vector<float >*)(&RECOELE_scl_E));
    Tree_->Branch("RECOELE_scl_Et",(std::vector<float >*)(&RECOELE_scl_Et));
    Tree_->Branch("RECOELE_scl_Eta",(std::vector<float >*)(&RECOELE_scl_Eta));
    Tree_->Branch("RECOELE_scl_Phi",(std::vector<float >*)(&RECOELE_scl_Phi));    
    Tree_->Branch("RECOELE_ecalEnergy",(std::vector<float >*)(&RECOELE_ecalEnergy));

    // Track-Cluster Matching    
    Tree_->Branch("RECOELE_ep",(std::vector<float >*)(&RECOELE_ep));
    Tree_->Branch("RECOELE_eSeedp",(std::vector<float >*)(&RECOELE_eSeedp));
    Tree_->Branch("RECOELE_eSeedpout",(std::vector<float >*)(&RECOELE_eSeedpout));
    Tree_->Branch("RECOELE_eElepout",(std::vector<float >*)(&RECOELE_eElepout));
    Tree_->Branch("RECOELE_deltaEtaIn",(std::vector<float >*)(&RECOELE_deltaEtaIn));
    Tree_->Branch("RECOELE_deltaEtaSeed",(std::vector<float >*)(&RECOELE_deltaEtaSeed));
    Tree_->Branch("RECOELE_deltaEtaEle",(std::vector<float >*)(&RECOELE_deltaEtaEle));
    Tree_->Branch("RECOELE_deltaPhiIn",(std::vector<float >*)(&RECOELE_deltaPhiIn));
    Tree_->Branch("RECOELE_deltaPhiSeed",(std::vector<float >*)(&RECOELE_deltaPhiSeed));
    Tree_->Branch("RECOELE_deltaPhiEle",(std::vector<float >*)(&RECOELE_deltaPhiEle));
    // Fiducial flags
    Tree_->Branch("RECOELE_isbarrel",(std::vector<int >*)(&RECOELE_isbarrel));   
    Tree_->Branch("RECOELE_isendcap",(std::vector<int >*)(&RECOELE_isendcap));   
    Tree_->Branch("RECOELE_isGap",(std::vector<int >*)(&RECOELE_isGap));
    Tree_->Branch("RECOELE_isEBetaGap",(std::vector<int >*)(&RECOELE_isEBetaGap));   
    Tree_->Branch("RECOELE_isEBphiGap",(std::vector<int >*)(&RECOELE_isEBphiGap));   
    Tree_->Branch("RECOELE_isEEdeeGap",(std::vector<int >*)(&RECOELE_isEEdeeGap));   
    Tree_->Branch("RECOELE_isEEringGap",(std::vector<int >*)(&RECOELE_isEEringGap));   
    // Shower shape
    Tree_->Branch("RECOELE_sigmaIetaIeta",(std::vector<float >*)(&RECOELE_sigmaIetaIeta));   
    Tree_->Branch("RECOELE_sigmaEtaEta",(std::vector<float >*)(&RECOELE_sigmaEtaEta));   
    Tree_->Branch("RECOELE_e15",(std::vector<float >*)(&RECOELE_e15));   
    Tree_->Branch("RECOELE_e25max",(std::vector<float >*)(&RECOELE_e25max));   
    Tree_->Branch("RECOELE_e55",(std::vector<float >*)(&RECOELE_e55));   
    Tree_->Branch("RECOELE_he",(std::vector<float >*)(&RECOELE_he));   
    Tree_->Branch("RECOELE_r9",(std::vector<float >*)(&RECOELE_r9));   
    // Particle flow
    Tree_->Branch("RECOELE_mva",(std::vector<float >*)(&RECOELE_mva));   
    // Brem & Classifaction
    Tree_->Branch("RECOELE_fbrem",(std::vector<float >*)(&RECOELE_fbrem));   
    Tree_->Branch("RECOELE_nbrems",(std::vector<int >*)(&RECOELE_nbrems));   
    //  golden/bigbrem/(narrow)/showering/crack
    Tree_->Branch("RECOELE_Class",(std::vector<int >*)(&RECOELE_Class));  
    //fBrem addition
    Tree_->Branch("RECOELE_fbrem_mode",(std::vector<double >*)(&RECOELE_fbrem_mode));
    Tree_->Branch("RECOELE_fbrem_mean",(std::vector<double >*)(&RECOELE_fbrem_mean));
    
    // Isolation
    Tree_->Branch("RECOELE_EGMTRACKISO",(std::vector<float >*)(&RECOELE_EGMTRACKISO));  
    Tree_->Branch("RECOELE_EGMHCALISO",(std::vector<float >*)(&RECOELE_EGMHCALISO));  
    Tree_->Branch("RECOELE_EGMECALISO",(std::vector<float >*)(&RECOELE_EGMECALISO)); 
    Tree_->Branch("RECOELE_EGMX",(std::vector<float >*)(&RECOELE_EGMX)); 

    // PF isolation
    Tree_->Branch("RECOELE_PFchAllPart",(std::vector<double >*)(&RECOELE_PFchAllPart));
    Tree_->Branch("RECOELE_PFchHad",(std::vector<double >*)(&RECOELE_PFchHad));
    Tree_->Branch("RECOELE_PFneuHad",(std::vector<double >*)(&RECOELE_PFneuHad));
    Tree_->Branch("RECOELE_PFphoton",(std::vector<double >*)(&RECOELE_PFphoton));
    Tree_->Branch("RECOELE_PFPUchAllPart",(std::vector<double >*)(&RECOELE_PFPUchAllPart));
    Tree_->Branch("RECOELE_PFX_dB",(std::vector<double >*)(&RECOELE_PFX_dB));
    Tree_->Branch("RECOELE_PFX_rho",(std::vector<double >*)(&RECOELE_PFX_rho));

    // Electron Regression
    Tree_->Branch("RECOELE_regEnergy",(std::vector<double >*)(&RECOELE_regEnergy));
    Tree_->Branch("RECOELE_regEnergyError",(std::vector<double >*)(&RECOELE_regEnergyError));
    
    // Vertexing DA and KF
    Tree_->Branch("RECOELE_SIP",(std::vector<float >*)(&RECOELE_SIP)); 
    Tree_->Branch("RECOELE_IP",(std::vector<float >*)(&RECOELE_IP)); 
    Tree_->Branch("RECOELE_IPERROR",(std::vector<float >*)(&RECOELE_IPERROR)); 
    Tree_->Branch("RECOELE_SIP_KF",(std::vector<float >*)(&RECOELE_SIP_KF)); 
    Tree_->Branch("RECOELE_IP_KF",(std::vector<float >*)(&RECOELE_IP_KF)); 
    Tree_->Branch("RECOELE_IPERROR_KF",(std::vector<float >*)(&RECOELE_IPERROR_KF)); 

     // GD vertex
    Tree_->Branch("RECOELE_SIP_GD",(std::vector<float >*)(&RECOELE_SIP_GD)); //2e2mu
    Tree_->Branch("RECOELE_SIP_GDEEEE",(std::vector<float >*)(&RECOELE_SIP_GDEEEE));  //4e
    // Std vertex
    Tree_->Branch("RECOELE_SIP_Std",(std::vector<float >*)(&RECOELE_SIP_Std)); //2e2mu
    Tree_->Branch("RECOELE_SIP_StdEEEE",(std::vector<float >*)(&RECOELE_SIP_StdEEEE));  //4e
    // Kin vertex
    Tree_->Branch("RECOELE_SIP_Kin",(std::vector<float >*)(&RECOELE_SIP_Kin)); //2e2mu
    Tree_->Branch("RECOELE_SIP_KinEEEE",(std::vector<float >*)(&RECOELE_SIP_KinEEEE));  //4e


    Tree_->Branch("RECOELE_STIP",(std::vector<float >*)(&RECOELE_STIP)); 
    Tree_->Branch("RECOELE_SLIP",(std::vector<float >*)(&RECOELE_SLIP)); 
    Tree_->Branch("RECOELE_TIP",(std::vector<float >*)(&RECOELE_TIP)); 
    Tree_->Branch("RECOELE_LIP",(std::vector<float >*)(&RECOELE_LIP)); 
    Tree_->Branch("RECOELE_TIPERROR",(std::vector<float >*)(&RECOELE_TIPERROR)); 
    Tree_->Branch("RECOELE_LIPERROR",(std::vector<float >*)(&RECOELE_LIPERROR)); 

    Tree_->Branch("RECOELE_sclRawE",(std::vector<double >*)(&ele_sclRawE));
    Tree_->Branch("RECOELE_sclX",(std::vector<double >*)(&ele_sclX)); 
    Tree_->Branch("RECOELE_sclY",(std::vector<double >*)(&ele_sclY)); 
    Tree_->Branch("RECOELE_sclZ",(std::vector<double >*)(&ele_sclZ));
    Tree_->Branch("RECOELE_seedSubdet1",(std::vector<double >*)(&ele_seedSubdet1));
    Tree_->Branch("RECOELE_seedDphi1",(std::vector<double >*)(&ele_seedDphi1)); 
    Tree_->Branch("RECOELE_seedDrz1",(std::vector<double >*)(&ele_seedDrz1));
    Tree_->Branch("RECOELE_seedSubdet2",(std::vector<double >*)(&ele_seedSubdet2));
    Tree_->Branch("RECOELE_seedDphi2",(std::vector<double >*)(&ele_seedDphi2)); 
    Tree_->Branch("RECOELE_seedDrz2",(std::vector<double >*)(&ele_seedDrz2));
    Tree_->Branch("RECOELE_eidVeryLoose",(std::vector<double >*)(&ele_eidVeryLoose)); 
    Tree_->Branch("RECOELE_eidLoose",(std::vector<double >*)(&ele_eidLoose)); 
    Tree_->Branch("RECOELE_eidMedium",(std::vector<double >*)(&ele_eidMedium)); 
    Tree_->Branch("RECOELE_eidTight",(std::vector<double >*)(&ele_eidTight)); 
    Tree_->Branch("RECOELE_eidHZZVeryLoose",(std::vector<double >*)(&ele_eidHZZVeryLoose)); 
    Tree_->Branch("RECOELE_eidHZZLoose",(std::vector<double >*)(&ele_eidHZZLoose)); 
    Tree_->Branch("RECOELE_eidHZZMedium",(std::vector<double >*)(&ele_eidHZZMedium)); 
    Tree_->Branch("RECOELE_eidHZZTight",(std::vector<double >*)(&ele_eidHZZTight)); 
    Tree_->Branch("RECOELE_mvaTrigV0",(std::vector<double >*)(&RECOELE_mvaTrigV0));     
    Tree_->Branch("RECOELE_mvaNonTrigV0",(std::vector<double >*)(&RECOELE_mvaNonTrigV0)); 
    Tree_->Branch("RECOELE_COV",(std::vector<std::vector<std::vector<double > > >*)(&RECOELE_COV)); 

    Tree_->Branch("RECOELE_TLE_ParentSC_X",(std::vector<float >*)(&RECOELE_TLE_ParentSC_X));
    Tree_->Branch("RECOELE_TLE_ParentSC_Y",(std::vector<float >*)(&RECOELE_TLE_ParentSC_Y));
    Tree_->Branch("RECOELE_TLE_ParentSC_Z",(std::vector<float >*)(&RECOELE_TLE_ParentSC_Z));


    //Reham Electron systematic variables

    Tree_->Branch("RECOELE_ecalTrkEnergyPreCorr",(std::vector<float >*)(&RECOELE_ecalTrkEnergyPreCorr));
    Tree_->Branch("RECOELE_ecalTrkEnergyErrPreCorr",(std::vector<float >*)(&RECOELE_ecalTrkEnergyErrPreCorr));
    Tree_->Branch("RECOELE_ecalTrkEnergyErrPostCorr",(std::vector<float >*)(&RECOELE_ecalTrkEnergyErrPostCorr));
    Tree_->Branch("RECOELE_energyScaleValue",(std::vector<float >*)(&RECOELE_energyScaleValue));       
    Tree_->Branch("RECOELE_energySigmaValue",(std::vector<float >*)(&RECOELE_energySigmaValue));
    Tree_->Branch("RECOELE_energyScaleUp",(std::vector<float >*)(&RECOELE_energyScaleUp));     
    Tree_->Branch("RECOELE_energyScaleDown",(std::vector<float >*)(&RECOELE_energyScaleDown));       
    Tree_->Branch("RECOELE_energyScaleStatUp",(std::vector<float >*)(&RECOELE_energyScaleStatUp));       
    Tree_->Branch("RECOELE_energyScaleStatDown",(std::vector<float >*)(&RECOELE_energyScaleStatDown));        
    Tree_->Branch("RECOELE_energyScaleSystUp",(std::vector<float >*)(&RECOELE_energyScaleSystUp));        
    Tree_->Branch("RECOELE_energyScaleSystDown",(std::vector<float >*)(&RECOELE_energyScaleSystDown));        
    Tree_->Branch("RECOELE_energyScaleGainUp",(std::vector<float >*)(&RECOELE_energyScaleGainUp));        
    Tree_->Branch("RECOELE_energyScaleGainDown",(std::vector<float >*)(&RECOELE_energyScaleGainDown));      
    Tree_->Branch("RECOELE_energyScaleEtUp",(std::vector<float >*)(&RECOELE_energyScaleEtUp));       
    Tree_->Branch("RECOELE_energyScaleEtDown",(std::vector<float >*)(&RECOELE_energyScaleEtDown));       
    Tree_->Branch("RECOELE_energySigmaUp",(std::vector<float >*)(&RECOELE_energySigmaUp));         
    Tree_->Branch("RECOELE_energySigmaDown",(std::vector<float >*)(&RECOELE_energySigmaDown));       
    Tree_->Branch("RECOELE_energySigmaPhiUp",(std::vector<float >*)(&RECOELE_energySigmaPhiUp));        
    Tree_->Branch("RECOELE_energySigmaPhiDown",(std::vector<float >*)(&RECOELE_energySigmaPhiDown));     
    Tree_->Branch("RECOELE_energySigmaRhoUp",(std::vector<float >*)(&RECOELE_energySigmaRhoUp));        
    Tree_->Branch("RECOELE_energySigmaRhoDown",(std::vector<float >*)(&RECOELE_energySigmaRhoDown));      

    // Muon block
    Tree_->Branch("RECOMU_isPFMu",(std::vector<int >*)(&RECOMU_isPFMu));
    Tree_->Branch("RECOMU_isGlobalMu",(std::vector<int >*)(&RECOMU_isGlobalMu));
    Tree_->Branch("RECOMU_isStandAloneMu",(std::vector<int >*)(&RECOMU_isStandAloneMu));
    Tree_->Branch("RECOMU_isTrackerMu",(std::vector<int >*)(&RECOMU_isTrackerMu));
    Tree_->Branch("RECOMU_isCaloMu",(std::vector<int >*)(&RECOMU_isCaloMu));
    Tree_->Branch("RECOMU_isTrackerHighPtMu",(std::vector<int >*)(&RECOMU_isTrackerHighPtMu));
    Tree_->Branch("RECOMU_isME0Muon",(std::vector<int >*)(&RECOMU_isME0Muon));
    Tree_->Branch("RECOMU_E",(std::vector<float >*)(&RECOMU_E)); 
    Tree_->Branch("RECOMU_PT",(std::vector<float >*)(&RECOMU_PT));
    Tree_->Branch("RECOMU_P",(std::vector<float >*)(&RECOMU_P)); 
    Tree_->Branch("RECOMU_ETA",(std::vector<float >*)(&RECOMU_ETA)); 
    Tree_->Branch("RECOMU_THETA",(std::vector<float >*)(&RECOMU_THETA)); 
    Tree_->Branch("RECOMU_PHI",(std::vector<float >*)(&RECOMU_PHI)); 
    Tree_->Branch("RECOMU_MASS",(std::vector<float >*)(&RECOMU_MASS)); 
    Tree_->Branch("RECOMU_CHARGE",(std::vector<float >*)(&RECOMU_CHARGE));
    // Tree_->Branch("RECOMU_Roch_calib_error",(std::vector<float >*)(&RECOMU_Roch_calib_error));
    Tree_->Branch("RECOMU_PT_uncorr",(std::vector<float >*)(&RECOMU_PT_uncorr));

    Tree_->Branch("RECOMU_COV",(std::vector<std::vector<std::vector<double > > >*)(&RECOMU_COV)); 

    Tree_->Branch("RECOMU_TRACKISO",(std::vector<float >*)(&RECOMU_TRACKISO));  
    Tree_->Branch("RECOMU_TRACKISO_SUMPT",(std::vector<float >*)(&RECOMU_TRACKISO_SUMPT));  
    Tree_->Branch("RECOMU_HCALISO",(std::vector<float >*)(&RECOMU_HCALISO));  
    Tree_->Branch("RECOMU_ECALISO",(std::vector<float >*)(&RECOMU_ECALISO)); 
    Tree_->Branch("RECOMU_X",(std::vector<float >*)(&RECOMU_X));

    Tree_->Branch("RECOMU_PFchHad",(std::vector<double >*)(&RECOMU_PFchHad));
    Tree_->Branch("RECOMU_PFneuHad",(std::vector<double >*)(&RECOMU_PFneuHad));
    Tree_->Branch("RECOMU_PFphoton",(std::vector<double >*)(&RECOMU_PFphoton));
    Tree_->Branch("RECOMU_PFPUchAllPart",(std::vector<double >*)(&RECOMU_PFPUchAllPart));
    Tree_->Branch("RECOMU_PFX_dB",(std::vector<double >*)(&RECOMU_PFX_dB));
    Tree_->Branch("RECOMU_PFX_rho",(std::vector<double >*)(&RECOMU_PFX_rho));


    // photon
    Tree_->Branch("RECOPFPHOT_PFchHad",(std::vector<double >*)(&RECOPFPHOT_PFchHad));
    Tree_->Branch("RECOPFPHOT_PFneuHad",(std::vector<double >*)(&RECOPFPHOT_PFneuHad));
    Tree_->Branch("RECOPFPHOT_PFphoton",(std::vector<double >*)(&RECOPFPHOT_PFphoton));
    Tree_->Branch("RECOPFPHOT_PFPUchAllPart",(std::vector<double >*)(&RECOPFPHOT_PFPUchAllPart));
    Tree_->Branch("RECOPFPHOT_PFX_rho",(std::vector<double >*)(&RECOPFPHOT_PFX_rho));

    // vertexing DA and KF
    Tree_->Branch("RECOMU_SIP",(std::vector<float >*)(&RECOMU_SIP)); 
    Tree_->Branch("RECOMU_IP",(std::vector<float >*)(&RECOMU_IP)); 
    Tree_->Branch("RECOMU_IPERROR",(std::vector<float >*)(&RECOMU_IPERROR)); 
    Tree_->Branch("RECOMU_SIP_KF",(std::vector<float >*)(&RECOMU_SIP_KF)); 
    Tree_->Branch("RECOMU_IP_KF",(std::vector<float >*)(&RECOMU_IP_KF)); 
    Tree_->Branch("RECOMU_IPERROR_KF",(std::vector<float >*)(&RECOMU_IPERROR_KF)); 

    // GD vertex
    Tree_->Branch("RECOMU_SIP_GD",(std::vector<float >*)(&RECOMU_SIP_GD)); //2e2mu
    Tree_->Branch("RECOMU_SIP_GDMMMM",(std::vector<float >*)(&RECOMU_SIP_GDMMMM));  //4mu
    // Std vertex
    Tree_->Branch("RECOMU_SIP_Std",(std::vector<float >*)(&RECOMU_SIP_Std)); //2e2mu
    Tree_->Branch("RECOMU_SIP_StdMMMM",(std::vector<float >*)(&RECOMU_SIP_StdMMMM));  //4mu
    // Kin vertex
    Tree_->Branch("RECOMU_SIP_Kin",(std::vector<float >*)(&RECOMU_SIP_Kin)); //2e2mu
    Tree_->Branch("RECOMU_SIP_KinMMMM",(std::vector<float >*)(&RECOMU_SIP_KinMMMM));  //4mu



    Tree_->Branch("RECOMU_STIP",(std::vector<float >*)(&RECOMU_STIP)); 
    Tree_->Branch("RECOMU_SLIP",(std::vector<float >*)(&RECOMU_SLIP)); 
    Tree_->Branch("RECOMU_TIP",(std::vector<float >*)(&RECOMU_TIP)); 
    Tree_->Branch("RECOMU_LIP",(std::vector<float >*)(&RECOMU_LIP)); 
    Tree_->Branch("RECOMU_TIPERROR",(std::vector<float >*)(&RECOMU_TIPERROR)); 
    Tree_->Branch("RECOMU_LIPERROR",(std::vector<float >*)(&RECOMU_LIPERROR)); 
    
 

    Tree_->Branch("RECOMU_caloCompatibility",(std::vector<float >*)(&RECOMU_caloCompatibility));
    Tree_->Branch("RECOMU_segmentCompatibility",(std::vector<float >*)(&RECOMU_segmentCompatibility)); 
    Tree_->Branch("RECOMU_numberOfMatches",(std::vector<int >*)(&RECOMU_numberOfMatches));
    Tree_->Branch("RECOMU_numberOfMatchedStations",(std::vector<int >*)(&RECOMU_numberOfMatchedStations));
    Tree_->Branch("RECOMU_glbmuPromptTight",(std::vector<int >*)(&RECOMU_glbmuPromptTight));
 
    // track variables from muons:
    Tree_->Branch( "RECOMU_trkmuArbitration",(std::vector<int >*)(&RECOMU_trkmuArbitration));
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityLoose",(std::vector<int >*)(&RECOMU_trkmu2DCompatibilityLoose));
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityTight",(std::vector<int >*)(&RECOMU_trkmu2DCompatibilityTight));
    Tree_->Branch( "RECOMU_trkmuOneStationLoose",(std::vector<int >*)(&RECOMU_trkmuOneStationLoose));
    Tree_->Branch( "RECOMU_trkmuOneStationTight",(std::vector<int >*)(&RECOMU_trkmuOneStationTight));
    Tree_->Branch( "RECOMU_trkmuLastStationLoose",(std::vector<int >*)(&RECOMU_trkmuLastStationLoose));
    Tree_->Branch( "RECOMU_trkmuLastStationTight",(std::vector<int >*)(&RECOMU_trkmuLastStationTight));
    Tree_->Branch( "RECOMU_trkmuOneStationAngLoose",(std::vector<int >*)(&RECOMU_trkmuOneStationAngLoose));
    Tree_->Branch( "RECOMU_trkmuOneStationAngTight",(std::vector<int >*)(&RECOMU_trkmuOneStationAngTight));
    Tree_->Branch( "RECOMU_trkmuLastStationAngLoose",(std::vector<int >*)(&RECOMU_trkmuLastStationAngLoose));
    Tree_->Branch( "RECOMU_trkmuLastStationAngTight",(std::vector<int >*)(&RECOMU_trkmuLastStationAngTight));
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtLoose",(std::vector<int >*)(&RECOMU_trkmuLastStationOptimizedLowPtLoose));
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtTight",(std::vector<int >*)(&RECOMU_trkmuLastStationOptimizedLowPtTight));

    Tree_->Branch( "RECOMU_mutrkPT",(std::vector<float >*)(&RECOMU_mutrkPT));
    Tree_->Branch( "RECOMU_mutrkPTError",(std::vector<float >*)(&RECOMU_mutrkPTError));
    Tree_->Branch( "RECOMU_mutrkDxy",(std::vector<float >*)(&RECOMU_mutrkDxy));
    Tree_->Branch( "RECOMU_mutrkDxyError",(std::vector<float >*)(&RECOMU_mutrkDxyError));
    Tree_->Branch( "RECOMU_mutrkDxyB",(std::vector<float >*)(&RECOMU_mutrkDxyB));
    Tree_->Branch( "RECOMU_mutrkDz",(std::vector<float >*)(&RECOMU_mutrkDz));
    Tree_->Branch( "RECOMU_mutrkDzError",(std::vector<float >*)(&RECOMU_mutrkDzError));
    Tree_->Branch( "RECOMU_mutrkDzB",(std::vector<float >*)(&RECOMU_mutrkDzB));
    Tree_->Branch( "RECOMU_mutrkChi2PerNdof",(std::vector<float >*)(&RECOMU_mutrkChi2PerNdof));
    Tree_->Branch( "RECOMU_mutrkCharge",(std::vector<float >*)(&RECOMU_mutrkCharge));
    Tree_->Branch( "RECOMU_mutrkNHits",(std::vector<float >*)(&RECOMU_mutrkNHits));
    Tree_->Branch( "RECOMU_mutrkNStripHits",(std::vector<float >*)(&RECOMU_mutrkNStripHits));
    Tree_->Branch( "RECOMU_mutrkNPixHits",(std::vector<float >*)(&RECOMU_mutrkNPixHits));
    Tree_->Branch( "RECOMU_mutrkNMuonHits",(std::vector<float >*)(&RECOMU_mutrkNMuonHits));
    Tree_->Branch( "RECOMU_mutrktrackerLayersWithMeasurement",(std::vector<float >*)(&RECOMU_mutrktrackerLayersWithMeasurement));
    
    Tree_->Branch( "RECOMU_muInnertrkDxy",(std::vector<float >*)(&RECOMU_muInnertrkDxy));
    Tree_->Branch( "RECOMU_muInnertrkDxyError",(std::vector<float >*)(&RECOMU_muInnertrkDxyError));
    Tree_->Branch( "RECOMU_muInnertrkDxyB",(std::vector<float >*)(&RECOMU_muInnertrkDxyB));
    Tree_->Branch( "RECOMU_muInnertrkDz",(std::vector<float >*)(&RECOMU_muInnertrkDz));
    Tree_->Branch( "RECOMU_muInnertrkDzError",(std::vector<float >*)(&RECOMU_muInnertrkDzError));
    Tree_->Branch( "RECOMU_muInnertrkDzB",(std::vector<float >*)(&RECOMU_muInnertrkDzB));
    Tree_->Branch( "RECOMU_muInnertrkChi2PerNdof",(std::vector<float >*)(&RECOMU_muInnertrkChi2PerNdof));
    Tree_->Branch( "RECOMU_muInnertrktrackerLayersWithMeasurement",(std::vector<float >*)(&RECOMU_muInnertrktrackerLayersWithMeasurement));
    Tree_->Branch( "RECOMU_muInnertrkPT",(std::vector<float >*)(&RECOMU_muInnertrkPT));
    Tree_->Branch( "RECOMU_muInnertrkPTError",(std::vector<float >*)(&RECOMU_muInnertrkPTError));
    Tree_->Branch( "RECOMU_muInnertrkCharge",(std::vector<float >*)(&RECOMU_muInnertrkCharge));
    Tree_->Branch( "RECOMU_muInnertrkNHits",(std::vector<float >*)(&RECOMU_muInnertrkNHits));
    Tree_->Branch( "RECOMU_muInnertrkNStripHits",(std::vector<float >*)(&RECOMU_muInnertrkNStripHits));
    Tree_->Branch( "RECOMU_muInnertrkNPixHits",(std::vector<float >*)(&RECOMU_muInnertrkNPixHits));
    // best tracks for 13 TeV analysis
    Tree_->Branch( "RECOMU_mubesttrkType",(std::vector<int >*)(&RECOMU_mubesttrkType));
    Tree_->Branch( "RECOMU_mubesttrkDxy",(std::vector<float >*)(&RECOMU_mubesttrkDxy));
    Tree_->Branch( "RECOMU_mubesttrkDxyError",(std::vector<float >*)(&RECOMU_mubesttrkDxyError));
    Tree_->Branch( "RECOMU_mubesttrkDz",(std::vector<float >*)(&RECOMU_mubesttrkDz));
    Tree_->Branch( "RECOMU_mubesttrkDzError",(std::vector<float >*)(&RECOMU_mubesttrkDzError));
    Tree_->Branch( "RECOMU_mubesttrkPTError",(std::vector<float >*)(&RECOMU_mubesttrkPTError));
    Tree_->Branch( "RECOMU_Rochester_Error",(std::vector<float >*)(&RECOMU_Rochester_Error));


    // Geom. Discri.
    Tree_->Branch("ftsigma",(std::vector<double >*)(&ftsigma));
    Tree_->Branch("gdX",(std::vector<double >*)(&gdX));
    Tree_->Branch("gdY",(std::vector<double >*)(&gdY));
    Tree_->Branch("gdZ",(std::vector<double >*)(&gdZ));
    Tree_->Branch("ftsigmalag",(std::vector<double >*)(&ftsigmalag));
    Tree_->Branch("gdlagX",(std::vector<double >*)(&gdlagX));
    Tree_->Branch("gdlagY",(std::vector<double >*)(&gdlagY));
    Tree_->Branch("gdlagZ",(std::vector<double >*)(&gdlagZ));
    Tree_->Branch("gdlagProb",(std::vector<double >*)(&gdlagProb));
    Tree_->Branch("gdlagNdof",(std::vector<double >*)(&gdlagNdof));
    Tree_->Branch("ftsigmaMMMM",(std::vector<double >*)(&ftsigmaMMMM));
    Tree_->Branch("gdXMMMM",(std::vector<double >*)(&gdXMMMM));
    Tree_->Branch("gdYMMMM",(std::vector<double >*)(&gdYMMMM));
    Tree_->Branch("gdZMMMM",(std::vector<double >*)(&gdZMMMM));
    Tree_->Branch("ftsigmalagMMMM",(std::vector<double >*)(&ftsigmalagMMMM));
    Tree_->Branch("gdlagXMMMM",(std::vector<double >*)(&gdlagXMMMM));
    Tree_->Branch("gdlagYMMMM",(std::vector<double >*)(&gdlagYMMMM));
    Tree_->Branch("gdlagZMMMM",(std::vector<double >*)(&gdlagZMMMM));
    Tree_->Branch("gdlagProbMMMM",(std::vector<double >*)(&gdlagProbMMMM));
    Tree_->Branch("gdlagNdofMMMM",(std::vector<double >*)(&gdlagNdofMMMM));
    Tree_->Branch("ftsigmaEEEE",(std::vector<double >*)(&ftsigmaEEEE));
    Tree_->Branch("gdXEEEE",(std::vector<double >*)(&gdXEEEE));
    Tree_->Branch("gdYEEEE",(std::vector<double >*)(&gdYEEEE));
    Tree_->Branch("gdZEEEE",(std::vector<double >*)(&gdZEEEE));
    Tree_->Branch("ftsigmalagEEEE",(std::vector<double >*)(&ftsigmalagEEEE));
    Tree_->Branch("gdlagXEEEE",(std::vector<double >*)(&gdlagXEEEE));
    Tree_->Branch("gdlagYEEEE",(std::vector<double >*)(&gdlagYEEEE));
    Tree_->Branch("gdlagZEEEE",(std::vector<double >*)(&gdlagZEEEE));
    Tree_->Branch("gdlagProbEEEE",(std::vector<double >*)(&gdlagProbEEEE));
    Tree_->Branch("gdlagNdofEEEE",(std::vector<double >*)(&gdlagNdofEEEE));
    
    // ConstraintFit 4l
    Tree_->Branch("StdFitVertexX",(std::vector<double >*)(&StdFitVertexX));
    Tree_->Branch("StdFitVertexY",(std::vector<double >*)(&StdFitVertexY));
    Tree_->Branch("StdFitVertexZ",(std::vector<double >*)(&StdFitVertexZ));
    Tree_->Branch("StdFitVertexChi2r",(std::vector<double >*)(&StdFitVertexChi2r));
    Tree_->Branch("StdFitVertexProb",(std::vector<double >*)(&StdFitVertexProb));
    Tree_->Branch("StdFitVertexTrack_PT",(std::vector<std::vector<float > >*)(&StdFitVertexTrack_PT));
    Tree_->Branch("StdFitVertexTrack_ETA",(std::vector<std::vector<float > >*)(&StdFitVertexTrack_ETA));
    Tree_->Branch("StdFitVertexTrack_PHI",(std::vector<std::vector<float > >*)(&StdFitVertexTrack_PHI));
    Tree_->Branch("KinFitVertexX",(std::vector<double >*)(&KinFitVertexX));
    Tree_->Branch("KinFitVertexY",(std::vector<double >*)(&KinFitVertexY));
    Tree_->Branch("KinFitVertexZ",(std::vector<double >*)(&KinFitVertexZ));
    Tree_->Branch("KinFitVertexChi2r",(std::vector<double >*)(&KinFitVertexChi2r));
    Tree_->Branch("KinFitVertexProb",(std::vector<double >*)(&KinFitVertexProb));

    Tree_->Branch("StdFitVertexXMMMM",(std::vector<double >*)(&StdFitVertexXMMMM));
    Tree_->Branch("StdFitVertexYMMMM",(std::vector<double >*)(&StdFitVertexYMMMM));
    Tree_->Branch("StdFitVertexZMMMM",(std::vector<double >*)(&StdFitVertexZMMMM));
    Tree_->Branch("StdFitVertexChi2rMMMM",(std::vector<double >*)(&StdFitVertexChi2rMMMM));
    Tree_->Branch("StdFitVertexProbMMMM",(std::vector<double >*)(&StdFitVertexProbMMMM));
    Tree_->Branch("StdFitVertexTrackMMMM_PT",(std::vector<std::vector<float > >*)(&StdFitVertexTrackMMMM_PT));
    Tree_->Branch("StdFitVertexTrackMMMM_ETA",(std::vector<std::vector<float > >*)(&StdFitVertexTrackMMMM_ETA));
    Tree_->Branch("StdFitVertexTrackMMMM_PHI",(std::vector<std::vector<float > >*)(&StdFitVertexTrackMMMM_PHI));
    Tree_->Branch("KinFitVertexXMMMM",(std::vector<double >*)(&KinFitVertexXMMMM));
    Tree_->Branch("KinFitVertexYMMMM",(std::vector<double >*)(&KinFitVertexYMMMM));
    Tree_->Branch("KinFitVertexZMMMM",(std::vector<double >*)(&KinFitVertexZMMMM));
    Tree_->Branch("KinFitVertexChi2rMMMM",(std::vector<double >*)(&KinFitVertexChi2rMMMM));
    Tree_->Branch("KinFitVertexProbMMMM",(std::vector<double >*)(&KinFitVertexProbMMMM));
    

    Tree_->Branch("StdFitVertexXEEEE",(std::vector<double >*)(&StdFitVertexXEEEE));
    Tree_->Branch("StdFitVertexYEEEE",(std::vector<double >*)(&StdFitVertexYEEEE));
    Tree_->Branch("StdFitVertexZEEEE",(std::vector<double >*)(&StdFitVertexZEEEE));
    Tree_->Branch("StdFitVertexChi2rEEEE",(std::vector<double >*)(&StdFitVertexChi2rEEEE));
    Tree_->Branch("StdFitVertexProbEEEE",(std::vector<double >*)(&StdFitVertexProbEEEE));
    Tree_->Branch("StdFitVertexTrackEEEE_PT",(std::vector<std::vector<float > >*)(&StdFitVertexTrackEEEE_PT));
    Tree_->Branch("StdFitVertexTrackEEEE_ETA",(std::vector<std::vector<float > >*)(&StdFitVertexTrackEEEE_ETA));
    Tree_->Branch("StdFitVertexTrackEEEE_PHI",(std::vector<std::vector<float > >*)(&StdFitVertexTrackEEEE_PHI));
    Tree_->Branch("KinFitVertexXEEEE",(std::vector<double >*)(&KinFitVertexXEEEE));
    Tree_->Branch("KinFitVertexYEEEE",(std::vector<double >*)(&KinFitVertexYEEEE));
    Tree_->Branch("KinFitVertexZEEEE",(std::vector<double >*)(&KinFitVertexZEEEE));
    Tree_->Branch("KinFitVertexChi2rEEEE",(std::vector<double >*)(&KinFitVertexChi2rEEEE));
    Tree_->Branch("KinFitVertexProbEEEE",(std::vector<double >*)(&KinFitVertexProbEEEE));

    // constrintFit 3l
    Tree_->Branch("StdFitVertexChi2rMMM",(std::vector<double >*)(&StdFitVertexChi2rMMM));
    Tree_->Branch("StdFitVertexProbMMM",(std::vector<double >*)(&StdFitVertexProbMMM));
    Tree_->Branch("StdFitVertexChi2rMME",(std::vector<double >*)(&StdFitVertexChi2rMME));
    Tree_->Branch("StdFitVertexProbMME",(std::vector<double >*)(&StdFitVertexProbMME));
    Tree_->Branch("StdFitVertexChi2rEEE",(std::vector<double >*)(&StdFitVertexChi2rEEE));
    Tree_->Branch("StdFitVertexProbEEE",(std::vector<double >*)(&StdFitVertexProbEEE));
    Tree_->Branch("StdFitVertexChi2rMEE",(std::vector<double >*)(&StdFitVertexChi2rMEE));
    Tree_->Branch("StdFitVertexProbMEE",(std::vector<double >*)(&StdFitVertexProbMEE));


     // constrintFit Dileptons
    Tree_->Branch("StdFitVertexChi2rDiLep",(std::vector<double >*)(&StdFitVertexChi2rDiLep));
    Tree_->Branch("StdFitVertexProbDiLep",(std::vector<double >*)(&StdFitVertexProbDiLep));

    // Conversions
    Tree_->Branch("ConvMapDist",(std::vector<float >*)(&ConvMapDist));
    Tree_->Branch("ConvMapDcot",(std::vector<float >*)(&ConvMapDcot));



    //MatchingMC:
    //Muons:
    Tree_->Branch("RECOMU_MatchingMCTruth",(std::vector<int >*)(&RECOMU_MatchingMCTruth));
    Tree_->Branch("RECOMU_MatchingMCpT",(std::vector<float >*)(&RECOMU_MatchingMCpT));
    Tree_->Branch("RECOMU_MatchingMCEta",(std::vector<float >*)(&RECOMU_MatchingMCEta));
    Tree_->Branch("RECOMU_MatchingMCPhi",(std::vector<float >*)(&RECOMU_MatchingMCPhi));

    //Electrons:
    Tree_->Branch("RECOELE_MatchingMCTruth",(std::vector<int >*)(&RECOELE_MatchingMCTruth));
    Tree_->Branch("RECOELE_MatchingMCpT",(std::vector<float >*)(&RECOELE_MatchingMCpT));
    Tree_->Branch("RECOELE_MatchingMCEta",(std::vector<float >*)(&RECOELE_MatchingMCEta));
    Tree_->Branch("RECOELE_MatchingMCPhi",(std::vector<float >*)(&RECOELE_MatchingMCPhi));

    //Gamma:
    Tree_->Branch("RECOPHOT_MatchingMCTruth",(std::vector<int >*)(&RECOPHOT_MatchingMCTruth));
    Tree_->Branch("RECOPHOT_MatchingMCpT",(std::vector<float >*)(&RECOPHOT_MatchingMCpT));
    Tree_->Branch("RECOPHOT_MatchingMCEta",(std::vector<float >*)(&RECOPHOT_MatchingMCEta));
    Tree_->Branch("RECOPHOT_MatchingMCPhi",(std::vector<float >*)(&RECOPHOT_MatchingMCPhi));

    //ZtoMuMu:
    Tree_->Branch("RECOzMuMu_MatchingMCTruth",(std::vector<int >*)(&RECOzMuMu_MatchingMCTruth));
    Tree_->Branch("RECOzMuMu_MatchingMCpT",(std::vector<float >*)(&RECOzMuMu_MatchingMCpT));
    Tree_->Branch("RECOzMuMu_MatchingMCmass",(std::vector<float >*)(&RECOzMuMu_MatchingMCmass));
    Tree_->Branch("RECOzMuMu_MatchingMCEta",(std::vector<float >*)(&RECOzMuMu_MatchingMCEta));
    Tree_->Branch("RECOzMuMu_MatchingMCPhi",(std::vector<float >*)(&RECOzMuMu_MatchingMCPhi));

    //ZtoEE:
    Tree_->Branch("RECOzEE_MatchingMCTruth",(std::vector<int >*)(&RECOzEE_MatchingMCTruth));
    Tree_->Branch("RECOzEE_MatchingMCpT",(std::vector<float >*)(&RECOzEE_MatchingMCpT));
    Tree_->Branch("RECOzEE_MatchingMCmass",(std::vector<float >*)(&RECOzEE_MatchingMCmass));
    Tree_->Branch("RECOzEE_MatchingMCEta",(std::vector<float >*)(&RECOzEE_MatchingMCEta));
    Tree_->Branch("RECOzEE_MatchingMCPhi",(std::vector<float >*)(&RECOzEE_MatchingMCPhi));

    //HtoZtoEEEE:
    Tree_->Branch("RECOHzzEEEE_MatchingMCTruth",(std::vector<int >*)(&RECOHzzEEEE_MatchingMCTruth));
    Tree_->Branch("RECOHzzEEEE_MatchingMCpT",(std::vector<float >*)(&RECOHzzEEEE_MatchingMCpT));
    Tree_->Branch("RECOHzzEEEE_MatchingMCmass",(std::vector<float >*)(&RECOHzzEEEE_MatchingMCmass));
    Tree_->Branch("RECOHzzEEEE_MatchingMCEta",(std::vector<float >*)(&RECOHzzEEEE_MatchingMCEta));
    Tree_->Branch("RECOHzzEEEE_MatchingMCPhi",(std::vector<float >*)(&RECOHzzEEEE_MatchingMCPhi));

    //HtoZtoEEMM:
    Tree_->Branch("RECOHzzEEMM_MatchingMCTruth",(std::vector<int >*)(&RECOHzzEEMM_MatchingMCTruth));
    Tree_->Branch("RECOHzzEEMM_MatchingMCpT",(std::vector<float >*)(&RECOHzzEEMM_MatchingMCpT));
    Tree_->Branch("RECOHzzEEMM_MatchingMCmass",(std::vector<float >*)(&RECOHzzEEMM_MatchingMCmass));
    Tree_->Branch("RECOHzzEEMM_MatchingMCEta",(std::vector<float >*)(&RECOHzzEEMM_MatchingMCEta));
    Tree_->Branch("RECOHzzEEMM_MatchingMCPhi",(std::vector<float >*)(&RECOHzzEEMM_MatchingMCPhi));

    //HtoZtoMMMM:
    Tree_->Branch("RECOHzzMMMM_MatchingMCTruth",(std::vector<int >*)(&RECOHzzMMMM_MatchingMCTruth));
    Tree_->Branch("RECOHzzMMMM_MatchingMCpT",(std::vector<float >*)(&RECOHzzMMMM_MatchingMCpT));
    Tree_->Branch("RECOHzzMMMM_MatchingMCmass",(std::vector<float >*)(&RECOHzzMMMM_MatchingMCmass));
    Tree_->Branch("RECOHzzMMMM_MatchingMCEta",(std::vector<float >*)(&RECOHzzMMMM_MatchingMCEta));
    Tree_->Branch("RECOHzzMMMM_MatchingMCPhi",(std::vector<float >*)(&RECOHzzMMMM_MatchingMCPhi));



    //Global Event 
    Tree_->Branch( "RECO_NMU", &RECO_NMU, "RECO_NMU/I"); 
    Tree_->Branch( "RECO_NELE", &RECO_NELE, "RECO_NELE/I"); 
    
    // Tracks
    Tree_->Branch( "RECO_NTRACK", &RECO_NTRACK, "RECO_NTRACK/I");
    Tree_->Branch( "RECO_TRACK_PT",(std::vector<float >*)(&RECO_TRACK_PT));
    Tree_->Branch( "RECO_TRACK_ETA",(std::vector<float >*)(&RECO_TRACK_ETA));
    Tree_->Branch( "RECO_TRACK_PHI",(std::vector<float >*)(&RECO_TRACK_PHI));
    Tree_->Branch( "RECO_TRACK_CHI2",(std::vector<float >*)(&RECO_TRACK_CHI2));
    Tree_->Branch( "RECO_TRACK_CHI2RED",(std::vector<float >*)(&RECO_TRACK_CHI2RED));
    Tree_->Branch( "RECO_TRACK_CHI2PROB",(std::vector<float >*)(&RECO_TRACK_CHI2PROB));
    Tree_->Branch( "RECO_TRACK_NHITS",(std::vector<int >*)(&RECO_TRACK_NHITS));
    Tree_->Branch( "RECO_TRACK_DXY",(std::vector<float >*)(&RECO_TRACK_DXY));
    Tree_->Branch( "RECO_TRACK_DXYERR",(std::vector<float >*)(&RECO_TRACK_DXYERR));
    Tree_->Branch( "RECO_TRACK_DZ",(std::vector<float >*)(&RECO_TRACK_DZ));
    Tree_->Branch( "RECO_TRACK_DZERR",(std::vector<float >*)(&RECO_TRACK_DZERR));
    
    // Photons
    Tree_->Branch("RECO_NPHOT", &RECO_NPHOT, "RECO_NPHOT/I");
    Tree_->Branch("RECOPHOT_PT",(std::vector<float >*)(&RECOPHOT_PT)); 
    Tree_->Branch("RECOPHOT_ETA",(std::vector<float >*)(&RECOPHOT_ETA)); 
    Tree_->Branch("RECOPHOT_PHI",(std::vector<float >*)(&RECOPHOT_PHI)); 
    Tree_->Branch("RECOPHOT_THETA",(std::vector<float >*)(&RECOPHOT_THETA)); 
    Tree_->Branch("RECOPHOT_TLE_ParentSC_X",(std::vector<float >*)(&RECOPHOT_TLE_ParentSC_X));
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Y",(std::vector<float >*)(&RECOPHOT_TLE_ParentSC_Y));
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Z",(std::vector<float >*)(&RECOPHOT_TLE_ParentSC_Z));

    Tree_->Branch("RECO_NPFPHOT", &RECO_NPFPHOT, "RECO_NPFPHOT/I");
    Tree_->Branch("RECOPFPHOT_PT",(std::vector<float >*)(&RECOPFPHOT_PT)); 
    Tree_->Branch("RECOPFPHOT_PTError",(std::vector<float >*)(&RECOPFPHOT_PTError));  
    Tree_->Branch("RECOPFPHOT_ETA",(std::vector<float >*)(&RECOPFPHOT_ETA)); 
    Tree_->Branch("RECOPFPHOT_PHI",(std::vector<float >*)(&RECOPFPHOT_PHI)); 
    Tree_->Branch("RECOPFPHOT_THETA",(std::vector<float >*)(&RECOPFPHOT_THETA)); 
    
    //Reham

     Tree_->Branch("RECOPFPHOT_PT_uncorr",(std::vector<float >*)(&RECOPFPHOT_PT_uncorr));

    //Reham photn systematic variables

    Tree_->Branch("RECOPFPHOT_ecalEnergyPreCorr",(std::vector<float >*)(&RECOPFPHOT_ecalEnergyPreCorr));
    Tree_->Branch("RECOPFPHOT_ecalEnergyErrPreCorr",(std::vector<float >*)(&RECOPFPHOT_ecalEnergyErrPreCorr));   
    Tree_->Branch("RECOPFPHOT_ecalEnergyErrPostCorr",(std::vector<float >*)(&RECOPFPHOT_ecalEnergyErrPostCorr));
    Tree_->Branch("RECOPFPHOT_energyScaleValue",(std::vector<float >*)(&RECOPFPHOT_energyScaleValue));       
    Tree_->Branch("RECOPFPHOT_energySigmaValue",(std::vector<float >*)(&RECOPFPHOT_energySigmaValue));
    Tree_->Branch("RECOPFPHOT_energyScaleUp",(std::vector<float >*)(&RECOPFPHOT_energyScaleUp));     
    Tree_->Branch("RECOPFPHOT_energyScaleDown",(std::vector<float >*)(&RECOPFPHOT_energyScaleDown));       
    Tree_->Branch("RECOPFPHOT_energyScaleStatUp",(std::vector<float >*)(&RECOPFPHOT_energyScaleStatUp));       
    Tree_->Branch("RECOPFPHOT_energyScaleStatDown",(std::vector<float >*)(&RECOPFPHOT_energyScaleStatDown));        
    Tree_->Branch("RECOPFPHOT_energyScaleSystUp",(std::vector<float >*)(&RECOPFPHOT_energyScaleSystUp));        
    Tree_->Branch("RECOPFPHOT_energyScaleSystDown",(std::vector<float >*)(&RECOPFPHOT_energyScaleSystDown));        
    Tree_->Branch("RECOPFPHOT_energyScaleGainUp",(std::vector<float >*)(&RECOPFPHOT_energyScaleGainUp));        
    Tree_->Branch("RECOPFPHOT_energyScaleGainDown",(std::vector<float >*)(&RECOPFPHOT_energyScaleGainDown));      
    Tree_->Branch("RECOPFPHOT_energyScaleEtUp",(std::vector<float >*)(&RECOPFPHOT_energyScaleEtUp));       
    Tree_->Branch("RECOPFPHOT_energyScaleEtDown",(std::vector<float >*)(&RECOPFPHOT_energyScaleEtDown));       
    Tree_->Branch("RECOPFPHOT_energySigmaUp",(std::vector<float >*)(&RECOPFPHOT_energySigmaUp));         
    Tree_->Branch("RECOPFPHOT_energySigmaDown",(std::vector<float >*)(&RECOPFPHOT_energySigmaDown));       
    Tree_->Branch("RECOPFPHOT_energySigmaPhiUp",(std::vector<float >*)(&RECOPFPHOT_energySigmaPhiUp));        
    Tree_->Branch("RECOPFPHOT_energySigmaPhiDown",(std::vector<float >*)(&RECOPFPHOT_energySigmaPhiDown));     
    Tree_->Branch("RECOPFPHOT_energySigmaRhoUp",(std::vector<float >*)(&RECOPFPHOT_energySigmaRhoUp));        
    Tree_->Branch("RECOPFPHOT_energySigmaRhoDown",(std::vector<float >*)(&RECOPFPHOT_energySigmaRhoDown)); 
    
    //Beam Spot position
    Tree_->Branch("BeamSpot_X",&BeamSpot_X,"BeamSpot_X/D");
    Tree_->Branch("BeamSpot_Y",&BeamSpot_Y,"BeamSpot_Y/D");
    Tree_->Branch("BeamSpot_Z",&BeamSpot_Z,"BeamSpot_Z/D");
    // Vertices
    Tree_->Branch( "RECO_NVTX", &RECO_NVTX, "RECO_NVTX/I");
    Tree_->Branch( "RECO_VERTEX_x",(std::vector<float >*)(&RECO_VERTEX_x));
    Tree_->Branch( "RECO_VERTEX_y",(std::vector<float >*)(&RECO_VERTEX_y));
    Tree_->Branch( "RECO_VERTEX_z",(std::vector<float >*)(&RECO_VERTEX_z));
    Tree_->Branch( "RECO_VERTEX_ndof",(std::vector<float >*)(&RECO_VERTEX_ndof));
    Tree_->Branch( "RECO_VERTEX_chi2",(std::vector<float >*)(&RECO_VERTEX_chi2));
    Tree_->Branch( "RECO_VERTEX_ntracks",(std::vector<int >*)(&RECO_VERTEX_ntracks));
    Tree_->Branch( "RECO_VERTEXPROB",(std::vector<float >*)(&RECO_VERTEXPROB));
    Tree_->Branch( "RECO_VERTEX_isValid",(std::vector<int >*)(&RECO_VERTEX_isValid));
    Tree_->Branch( "RECO_VERTEX_TRACK_PT",(std::vector<std::vector<float > >*)(&RECO_VERTEX_TRACK_PT));
    
    // PFJets
    Tree_->Branch( "RECO_PFJET_N",   &RECO_PFJET_N,   "RECO_PFJET_N/I");
    Tree_->Branch( "RECO_PFJET_CHARGE",(std::vector<int >*)(&RECO_PFJET_CHARGE));
    Tree_->Branch( "RECO_PFJET_ET",(std::vector<float >*)(&RECO_PFJET_ET));
    Tree_->Branch( "RECO_PFJET_PT",(std::vector<float >*)(&RECO_PFJET_PT));
    Tree_->Branch( "RECO_PFJET_ETA",(std::vector<float >*)(&RECO_PFJET_ETA));
    Tree_->Branch( "RECO_PFJET_PHI",(std::vector<float >*)(&RECO_PFJET_PHI));
    Tree_->Branch( "RECO_PFJET_PUID_loose",(std::vector<int >*)(&RECO_PFJET_PUID_loose));
    Tree_->Branch( "RECO_PFJET_PUID_medium",(std::vector<int >*)(&RECO_PFJET_PUID_medium));
    Tree_->Branch( "RECO_PFJET_PUID",(std::vector<int >*)(&RECO_PFJET_PUID));
    Tree_->Branch( "RECO_PFJET_PUID_MVA",(std::vector<float >*)(&RECO_PFJET_PUID_MVA));
    Tree_->Branch( "RECO_PFJET_QG_Likelihood",(std::vector<float >*)(&RECO_PFJET_QG_Likelihood));//REHAM QG tagger
    Tree_->Branch( "RECO_PFJET_QG_axis2",(std::vector<float >*)(&RECO_PFJET_QG_axis2));//REHAM QG tagger
    Tree_->Branch( "RECO_PFJET_QG_ptd",(std::vector<float >*)(&RECO_PFJET_QG_ptd));//REHAM QG tagger
    Tree_->Branch( "RECO_PFJET_QG_mult",(std::vector<int >*)(&RECO_PFJET_QG_mult));//REHAM QG tagger
    Tree_->Branch( "RHO_ele", &RHO_ele, "RHO_ele/D");
    Tree_->Branch( "RHO_mu", &RHO_mu, "RHO_mu/D");

 //@
    Tree_->Branch("LHE_PARTON_N", &LHE_PARTON_N, "LHE_PARTON_N/I");
    Tree_->Branch("LHE_PARTON_CLEAR",(std::vector<int >*)(&LHE_PARTON_CLEAR));
    Tree_->Branch("LHE_PARTON_PDGID",(std::vector<int >*)(&LHE_PARTON_PDGID));
    Tree_->Branch("LHE_PARTON_PT",(std::vector<float >*)(&LHE_PARTON_PT));
    Tree_->Branch("LHE_PARTON_ETA",(std::vector<float >*)(&LHE_PARTON_ETA));
    Tree_->Branch("LHE_PARTON_PHI",(std::vector<float >*)(&LHE_PARTON_PHI));
    Tree_->Branch("LHE_PARTON_E",(std::vector<float >*)(&LHE_PARTON_E));
    Tree_->Branch("RECO_PFJET_PT_UncUp",(std::vector<float >*)(&RECO_PFJET_PT_UncUp));
    Tree_->Branch("RECO_PFJET_PT_UncDn",(std::vector<float >*)(&RECO_PFJET_PT_UncDn));
    Tree_->Branch("RECO_PFJET_AREA",(std::vector<float >*)(&RECO_PFJET_AREA));
    Tree_->Branch("RECO_PFJET_PTD",(std::vector<float >*)(&RECO_PFJET_PTD));
    Tree_->Branch("RECO_PFJET_CHARGED_HADRON_ENERGY",(std::vector<float >*)(&RECO_PFJET_CHARGED_HADRON_ENERGY));
    Tree_->Branch("RECO_PFJET_NEUTRAL_HADRON_ENERGY",(std::vector<float >*)(&RECO_PFJET_NEUTRAL_HADRON_ENERGY));
    Tree_->Branch("RECO_PFJET_PHOTON_ENERGY",(std::vector<float >*)(&RECO_PFJET_PHOTON_ENERGY));
    Tree_->Branch("RECO_PFJET_ELECTRON_ENERGY",(std::vector<float >*)(&RECO_PFJET_ELECTRON_ENERGY));
    Tree_->Branch("RECO_PFJET_MUON_ENERGY",(std::vector<float >*)(&RECO_PFJET_MUON_ENERGY));
    Tree_->Branch("RECO_PFJET_HF_HADRON_ENERGY",(std::vector<float >*)(&RECO_PFJET_HF_HADRON_ENERGY));
    Tree_->Branch("RECO_PFJET_HF_EM_ENERGY",(std::vector<float >*)(&RECO_PFJET_HF_EM_ENERGY));
    Tree_->Branch("RECO_PFJET_CHARGED_EM_ENERGY",(std::vector<float >*)(&RECO_PFJET_CHARGED_EM_ENERGY));
    Tree_->Branch("RECO_PFJET_CHARGED_MU_ENERGY",(std::vector<float >*)(&RECO_PFJET_CHARGED_MU_ENERGY));
    Tree_->Branch("RECO_PFJET_NEUTRAL_EM_ENERGY",(std::vector<float >*)(&RECO_PFJET_NEUTRAL_EM_ENERGY));
    Tree_->Branch("RECO_PFJET_CHARGED_HADRON_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_CHARGED_HADRON_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_PHOTON_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_PHOTON_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_ELECTRON_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_ELECTRON_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_MUON_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_MUON_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_HF_HADRON_MULTIPLICTY",(std::vector<int >*)(&RECO_PFJET_HF_HADRON_MULTIPLICTY));
    Tree_->Branch("RECO_PFJET_HF_EM_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_HF_EM_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_CHARGED_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_CHARGED_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_NEUTRAL_MULTIPLICITY",(std::vector<int >*)(&RECO_PFJET_NEUTRAL_MULTIPLICITY));
    Tree_->Branch("RECO_PFJET_NCOMPONENTS",(std::vector<int >*)(&RECO_PFJET_NCOMPONENTS));
    Tree_->Branch("RECO_PFJET_COMPONENT_PDGID",(std::vector<std::vector<int > >*)(&RECO_PFJET_COMPONENT_PDGID));
    Tree_->Branch("RECO_PFJET_COMPONENT_PT",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_PT));
    Tree_->Branch("RECO_PFJET_COMPONENT_ETA",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_ETA));
    Tree_->Branch("RECO_PFJET_COMPONENT_PHI",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_PHI));
    Tree_->Branch("RECO_PFJET_COMPONENT_E",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_E));
    Tree_->Branch("RECO_PFJET_COMPONENT_CHARGE",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_CHARGE));
    Tree_->Branch("RECO_PFJET_COMPONENT_TRANSVERSE_MASS",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_TRANSVERSE_MASS));
    Tree_->Branch("RECO_PFJET_COMPONENT_XVERTEX",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_XVERTEX));
    Tree_->Branch("RECO_PFJET_COMPONENT_YVERTEX",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_YVERTEX));
    Tree_->Branch("RECO_PFJET_COMPONENT_ZVERTEX",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_ZVERTEX));
    Tree_->Branch("RECO_PFJET_COMPONENT_VERTEX_CHI2",(std::vector<std::vector<float > >*)(&RECO_PFJET_COMPONENT_VERTEX_CHI2));
    
    
    //CaloMET
    Tree_->Branch( "RECO_CALOMET",          &calomet,          "RECO_CALOMET/F");
 /*    Tree_->Branch( "RECO_CALOMETHO",        &calometho,        "RECO_CALOMETHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHFHO",    &calometnohfho,    "RECO_CALOMETNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHF",      &calometnohf,      "RECO_CALOMETNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTHO",     &calometoptho,     "RECO_CALOMETOPTHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHFHO", &calometoptnohfho, "RECO_CALOMETOPTNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHF",   &calometoptnohf,   "RECO_CALOMETOPTNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPT",       &calometopt,       "RECO_CALOMETOPT/F"); */
 
   //Particle Flow MET
    //Type1 correction //Reham

    Tree_->Branch( "RECO_PFMET", &pfmet, "RECO_PFMET/F");
    Tree_->Branch( "RECO_PFMET_X", &pfmet_x, "RECO_PFMET_X/F");
    Tree_->Branch( "RECO_PFMET_Y", &pfmet_y, "RECO_PFMET_Y/F");
    Tree_->Branch( "RECO_PFMET_PHI", &pfmet_phi, "RECO_PFMET_PHI/F");
    Tree_->Branch( "RECO_PFMET_THETA", &pfmet_theta, "RECO_PFMET_THETA/F");

    //Uncorrected //Reham
    Tree_->Branch( "RECO_PFMET_uncorr", &pfmet_uncorr, "RECO_PFMET_uncorr/F");
    Tree_->Branch( "RECO_PFMET_X_uncorr", &pfmet_x_uncorr, "RECO_PFMET_X_uncorr/F");
    Tree_->Branch( "RECO_PFMET_Y_uncorr", &pfmet_y_uncorr, "RECO_PFMET_Y_uncorr/F");
    Tree_->Branch( "RECO_PFMET_PHI_uncorr", &pfmet_phi_uncorr, "RECO_PFMET_PHI_uncorr/F");
    Tree_->Branch( "RECO_PFMET_THETA_uncorr", &pfmet_theta_uncorr, "RECO_PFMET_THETA_uncorr/F");

    //MET uncertinities

    Tree_->Branch( "RECO_PFMET_JetEnUp", &pfmet_JetEnUp, "RECO_PFMET_JetEnUp/F");
    Tree_->Branch( "RECO_PFMET_JetEnDn", &pfmet_JetEnDn, "RECO_PFMET_JetEnDn/F");
    Tree_->Branch( "RECO_PFMET_ElectronEnUp", &pfmet_ElectronEnUp, "RECO_PFMET_ElectronEnUp/F");
    Tree_->Branch( "RECO_PFMET_ElectronEnDn", &pfmet_ElectronEnDn, "RECO_PFMET_ElectronEnDn/F");
    Tree_->Branch( "RECO_PFMET_MuonEnUp", &pfmet_MuonEnUp, "RECO_PFMET_MuonEnUp/F");
    Tree_->Branch( "RECO_PFMET_MuonEnDn", &pfmet_MuonEnDn, "RECO_PFMET_MuonEnDn/F");
    Tree_->Branch( "RECO_PFMET_JetResUp", &pfmet_JetResUp, "RECO_PFMET_JetResUp/F");
    Tree_->Branch( "RECO_PFMET_JetResDn", &pfmet_JetResDn, "RECO_PFMET_JetResDn/F");
    Tree_->Branch( "RECO_PFMET_UnclusteredEnUp", &pfmet_UnclusteredEnUp, "RECO_PFMET_UnclusteredEnUp/F");
    Tree_->Branch( "RECO_PFMET_UnclusteredEnDn", &pfmet_UnclusteredEnDn, "RECO_PFMET_UnclusteredEnDn/F");
    Tree_->Branch( "RECO_PFMET_PhotonEnUp", &pfmet_PhotonEnUp, "RECO_PFMET_PhotonEnUp/F");
    Tree_->Branch( "RECO_PFMET_PhotonEnDn", &pfmet_PhotonEnDn, "RECO_PFMET_PhotonEnDn/F");
    Tree_->Branch( "RECO_PFMET_TauEnUp ", &pfmet_TauEnUp , "RECO_PFMET_TauEnUp /F");
    Tree_->Branch( "RECO_PFMET_TauEnDown", &pfmet_TauEnDown, "RECO_PFMET_TauEnDown/F");

    // decision of MET filter
    //Tree_->Branch( "RECO_PFMET_passecalBadCalibFilterUpdate", &PassecalBadCalibFilterUpdated, "RECO_PFMET_passecalBadCalibFilterUpdate/I");//new
    // Tree_->Branch( "RECO_PFMET_filterbadChCandidate", &filterbadChCandidate, "RECO_PFMET_filterbadChCandidate/I");
    // Tree_->Branch( "RECO_PFMET_filterbadPFMuon", &filterbadPFMuon, "RECO_PFMET_filterbadPFMuon/I");
   
    Tree_->Branch( "RECO_PFMET_GoodVtxNoiseFilter",&passFilterGoodVtxNoise, "RECO_PFMET_GoodVtxNoiseFilter/I");
    Tree_->Branch( "RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter",&passFilterGlobalSuperTightHalo2016NoiseFilter, "RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter/I");
    Tree_->Branch( "RECO_PFMET_HBHENoiseFilter",&passFilterHBHENoise, "RECO_PFMET_HBHENoiseFilter/I");
    Tree_->Branch( "RECO_PFMET_HBHENoiseIsoFilter",&passFilterHBHENoiseIso, "RECO_PFMET_HBHENoiseIsoFilter/I");
    Tree_->Branch( "RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter",&passFilterEcalDeadCellTriggerPrimitiveNoise, "RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter/I");
    Tree_->Branch( "RECO_PFMET_BadPFMuonFilter",& passFilterBadPFMuon, "RECO_PFMET_BadPFMuonFilter/I");
    Tree_->Branch( "RECO_PFMET_BadChargedCandidateFilter",&passFilterBadChargedCandidate, "RECO_PFMET_BadChargedCandidateFilter/I");
    Tree_->Branch( "RECO_PFMET_EEBadScNoiseFilter",&passFilterEEBadScNoise, "RECO_PFMET_EEBadScNoiseFilter/I");
    Tree_->Branch( "RECO_PFMET_EcalBadCalibFilter",&passFilterEcalBadCalib, "RECO_PFMET_EcalBadCalibFilter/I");
    
 
    //Track Corrected MET
    Tree_->Branch( "RECO_TCMET", &tcmet, "RECO_TCMET/F");
    //Type I correction MET
    Tree_->Branch( "RECO_CORMETMUONS",  &cormetmuons,  "RECO_CORMETMUONS/F");
   

    // Btagging jets and discriminators
    Tree_->Branch("tCHighEff_BTagJet_PT",(std::vector<float >*)(&tCHighEff_BTagJet_PT));
    Tree_->Branch("tCHighPur_BTagJet_PT",(std::vector<float >*)(&tCHighPur_BTagJet_PT));
    Tree_->Branch("cSV_BTagJet_PT",(std::vector<float >*)(&cSV_BTagJet_PT));
    Tree_->Branch("tCHighEff_BTagJet_ETA",(std::vector<float >*)(&tCHighEff_BTagJet_ETA));
    Tree_->Branch("tCHighPur_BTagJet_ETA",(std::vector<float >*)(&tCHighPur_BTagJet_ETA));
    Tree_->Branch("cSV_BTagJet_ETA",(std::vector<float >*)(&cSV_BTagJet_ETA));
    Tree_->Branch("tCHighEff_BTagJet_PHI",(std::vector<float >*)(&tCHighEff_BTagJet_PHI));
    Tree_->Branch("tCHighPur_BTagJet_PHI",(std::vector<float >*)(&tCHighPur_BTagJet_PHI));
    Tree_->Branch("cSV_BTagJet_PHI",(std::vector<float >*)(&cSV_BTagJet_PHI));
    Tree_->Branch("cSV_BTagJet_ET",(std::vector<float >*)(&cSV_BTagJet_ET));
    Tree_->Branch("tCHighEff_BTagJet_DISCR",(std::vector<float >*)(&tCHighEff_BTagJet_DISCR));
    Tree_->Branch("tCHighPur_BTagJet_DISCR",(std::vector<float >*)(&tCHighPur_BTagJet_DISCR));
    Tree_->Branch("cSV_BTagJet_DISCR",(std::vector<float >*)(&cSV_BTagJet_DISCR));     
  }
  
  
  void Initialize() {}
  
  void fillPU(const edm::Event& iEvent){
      edm::Handle<vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(PileupSrc_, PupInfo);

      if(!PupInfo.isValid()) return;
      for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) { 
	//@//	std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
	if (cand->getBunchCrossing() == 0) num_PU_vertices=cand->getTrueNumInteractions();;
	// num_PU_vertices=cand->getPU_NumInteractions(); in-time,out-of-time pileup
	PU_BunchCrossing=cand->getBunchCrossing();
      }	
  }

  void EventsMCReWeighting(const edm::Event& iEvent){

    // get the weight                                                                                                                                                     
    MC_weighting=0.;
    float EventWeight = 1.0;
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    iEvent.getByToken(generator_, gen_ev_info);
    if(!gen_ev_info.isValid()) return;
    EventWeight = gen_ev_info->weight();
    //std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;
                                                                                                                                                                        
    float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
    //std::cout<<"mc_weight = "<< mc_weight <<std::endl;                                                                                                                  
    MC_weighting=mc_weight;

  }


  void fillHLTFired(const edm::Event& iEvent){
    edm::Handle<vector<std::string> > HLTfired_;
    iEvent.getByLabel(HLTInfoFired,HLTfired_);

    vector<string> HLTimported;
    string tmpstring="";

    for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
      unsigned int i=cand-HLTfired_->begin();
      HLTimported.push_back(cand->c_str());
      string newstr=HLTimported.at(i) + ":" + tmpstring;
      tmpstring=newstr;
    }

    //@//  cout << "HLTFiredString= " << tmpstring.c_str() << endl;
    if (!tmpstring.empty()) sprintf(HLTPathsFired,tmpstring.c_str());
        

  }


  void triggermatching(const edm::Event& iEvent){

    //@//  cout << "Start Trigger matching for muons" << endl;
    // check HLTrigger/Configuration/python/HLT_GRun_cff.py

//AOD    edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects );
    //    const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
    size_t nMuHLT =0, nEleHLT=0;
    
    
    //Reham
    ////
     iEvent.getByToken(triggerBits_, triggerBits);//Reham
     iEvent.getByToken(triggerObjects_, triggerObjects);//Reham
    
     const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits); //Reham   
    
    
    // cout<<"number 1"<<endl;
    
    std::vector<pat::TriggerObjectStandAlone>  HLTMuMatched,HLTEleMatched;
    std::vector<string> HLTMuMatchedNames,HLTEleMatchedNames;
    
    //MiniAOD
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {

       obj.unpackPathNames(names);//Reham
       obj.unpackFilterLabels(iEvent,*triggerBits);//Reham

      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
	std::string fullname = obj.filterLabels()[h];
	std::string name;
	size_t p = fullname.find_first_of(':');
	if ( p != std::string::npos) {
	  name = fullname.substr(0, p);
	}
	else {
	  name = fullname;
	} 
	
	//cout<<"number 2"<<endl;
	
	if (name == triggerFilter.c_str()) {
	  HLTMuMatched.push_back(obj);
	  HLTMuMatchedNames.push_back(name);
	  //@//  cout << "Matching " << triggerFilter.c_str()  << endl;
	  nMuHLT++;
	}
	if (name == triggerEleFilter.c_str()) {
	  HLTEleMatched.push_back(obj);
	  HLTEleMatchedNames.push_back(name);
	  //@// cout << "Matching " << triggerEleFilter.c_str()  << endl;
	  nEleHLT++;
	}
      }
    }
    
    
    
    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);
    float maxDeltaR_=0.2;
    float maxDPtRel_=1.0;
    int nMuHLTMatch=0;
    for (edm::View<pat::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){
      unsigned int i=iCand-MuCandidates->begin();
      //@// cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;
      if (IsMuMatchedToHLTMu(*iCand,  HLTMuMatched , HLTMuMatchedNames, maxDeltaR_, maxDPtRel_)==true){
	nMuHLTMatch++;
	//@//	cout << "Muon HLT Matched with pT= " << iCand->pt() << endl;
	RECOMU_PT_MuHLTMatch[i] =iCand->pt();
	RECOMU_ETA_MuHLTMatch[i]=iCand->eta();
	RECOMU_PHI_MuHLTMatch[i]=iCand->phi();
      }
    }

    //@// cout << "N. Muons HLT Matched= " << nMuHLTMatch << " FiredString:" << HLTPathsFired << endl;
    RECO_nMuHLTMatch    = nMuHLTMatch;


    //@// cout << "Start Trigger matching for electron" << endl;

    int nEleHLTMatch=0;
    edm::Handle<edm::View<pat::Electron> > EleCandidates;
    iEvent.getByToken(electronEgmTag_, EleCandidates);

    for (edm::View<pat::Electron>::const_iterator iCand = EleCandidates->begin(); iCand != EleCandidates->end(); ++iCand){

      unsigned int i=iCand-EleCandidates->begin();
      //@// cout << "Electron with pt= " << iCand->pt() << ": check trigger matching" << endl;
      if (IsEleMatchedToHLTEle(*iCand,  HLTEleMatched , HLTEleMatchedNames, maxDeltaR_, maxDPtRel_)==true){
	    cout << "Electron HLT Matched with pT= " << iCand->pt() << endl;
	    nEleHLTMatch++;
	    RECOELE_PT_EleHLTMatch[i]=iCand->pt();
	    RECOELE_ETA_EleHLTMatch[i]=iCand->eta();
	    RECOELE_PHI_EleHLTMatch[i]=iCand->phi();
      }
    }

    //@// cout << "N. Electrons HLT Matched= " << nEleHLTMatch << endl;

    RECO_nEleHLTMatch = nEleHLTMatch;


  }

  bool IsMuMatchedToHLTMu ( const pat::Muon &mu, std::vector<pat::TriggerObjectStandAlone> HLTMu , std::vector<string> HLTMuNames, double DR, double DPtRel ) {
    size_t dim =  HLTMu.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT mu filter is= " << HLTMuNames[k].c_str() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
      if (  (deltaR(HLTMu[k], mu) < DR)   && (fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt()<DPtRel)){ 
	//@//	cout << "HLT mu filter is= " << HLTMuNames[k].c_str() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool IsEleMatchedToHLTEle ( const pat::Electron &ele, std::vector<pat::TriggerObjectStandAlone> HLTEle , std::vector<string> HLTEleNames, double DR, double DPtRel ) {
    size_t dim =  HLTEle.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
      if (  (deltaR(HLTEle[k], ele) < DR)   && (fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt()<DPtRel)){ 
	//@//	cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool isTrackerHighPtMu (const reco::Muon &mu, math::XYZPoint primaryVertex){
    return( mu.numberOfMatchedStations() > 1 &&
            (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3
            && std::abs(mu.muonBestTrack()->dxy(primaryVertex)) < 0.2
            && std::abs(mu.muonBestTrack()->dz(primaryVertex)) < 0.5
            && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
            && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
            );
  }

  
  // PDT
  std::string getParticleName(int id) const{
    const ParticleData * pd = pdt_->particle( id );
    if (!pd) {
      std::ostringstream ss;
      ss << "P" << id;
      return ss.str();
    } else
      return pd->name();
  }

  
  void fillgenjets(const edm::Event& iEvent){
    edm::Handle<reco::GenJetCollection> genjetHandle;
    iEvent.getByToken(genjetTag_,genjetHandle);
    int i=0;
    for ( GenJetCollection::const_iterator igen=genjetHandle->begin(); igen!=genjetHandle->end(); igen++) {
      if (i>99) break;
      MC_GENJET_PT[i]=igen->pt();
      MC_GENJET_ETA[i]=igen->eta();
      MC_GENJET_PHI[i]=igen->phi();
      i++;
    }
  }

  // GenParticles  
  void fillgenparticles(const edm::Event& iEvent, const edm::EventSetup &es){    
    
    // get gen particle candidates 
    edm::Handle<reco::GenParticleCollection> genCandidates;           
    iEvent.getByToken(genParticles_, genCandidates);

    es.getData( pdt_ );
    vector<float> leptonpt,leptoneta,leptonphi,leptontheta,leptonpdgid;
    
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // lepton stable      
   
      if ( ( abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13 || abs(mcIter->pdgId())==15  ) && mcIter->status()==1 ){	
	
	//@//	std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ; 
	leptonpt.push_back(mcIter->pt());
	leptoneta.push_back(mcIter->eta());
	leptonphi.push_back(mcIter->phi());
	leptontheta.push_back(mcIter->theta());
        leptonpdgid.push_back(float(mcIter->pdgId()));

	if (mcIter->numberOfMothers()>0){ 
	  //@// std::cout << getParticleName(int(mcIter->mother(0)->pdgId())) << " << ";
	  if (mcIter->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	    //@//  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->pdgId())) << " << ";
	    if (mcIter->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      //@//  std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		//@//	std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  //@// std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && 
mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		    //@//  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << " << std::endl;
		    if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  }
		}
	      }
	    }
	  }
	}

    
      }
    }

    sort(leptonpt.rbegin(),leptonpt.rend());
    
    if (leptonpt.size()>0) { 
      //@// cout << "\n The 4 highest pt leptons in MCtruth are= ";
      for (unsigned int i=0;i<leptonpt.size();i++){
	if (i>3) continue;
	cout << leptonpt.at(i) << " ";
	MC_LEPT_PT[i]=leptonpt.at(i);
	MC_LEPT_ETA[i]=leptoneta.at(i);
	MC_LEPT_PHI[i]=leptonphi.at(i);
	MC_LEPT_THETA[i]=leptontheta.at(i);
        MC_LEPT_PDGID[i]=leptonpdgid.at(i);
      }
      cout << endl;
    }

    // Check if Z bosons are generated (specially for Z + jet MC)
    int i=0;
      //j=0;
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // Select the Z decay: status 3 
      if ( abs(mcIter->pdgId())==23 && mcIter->status()==3 ) {
	int j=0;
	//@//	std::cout << "\n Found generated Z with PDGId = " << mcIter->pdgId() << " and status = "<< mcIter->status() << " and " << mcIter->numberOfDaughters() << " daughts |";
	//@//	std::cout << " pt = "<<std::setw(12)<<mcIter->pt()<<" GeV/c | eta = "<<std::setw(12)<<mcIter->eta()<<" | phi = "<<std::setw(12)<<mcIter->phi();
	//@//	std::cout << "   ===> Filled Z Bo at ["<<i<<"]["<<j<<"]";
	MC_Z_PT[i][0]=mcIter->pt();       MC_Z_ETA[i][0]=mcIter->eta();   MC_Z_PHI[i][0]=mcIter->phi();
        MC_Z_THETA[i][0]=mcIter->theta(); MC_Z_MASS[i][0]=mcIter->mass(); MC_Z_PDGID[i][0]=mcIter->pdgId();  
	
	// Check daughters of Z decay: status 3
	reco::GenParticle::daughters dst3 = mcIter->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_dst3 = dst3.begin(), est3 = dst3.end(); it_dst3 != est3; ++it_dst3) {
	  // Select status 3 electron, muon or tau
	  if (((abs((**it_dst3).pdgId()) == 11) || (abs((**it_dst3).pdgId()) == 13) || (abs((**it_dst3).pdgId()) == 15)) && (**it_dst3).status() == 3) {
	    //@//  std::cout<<"\n |--> Z daughter Particle  "<<std::setw(18)<<"| id = "<<std::setw(5)<<(**it_dst3).pdgId()<<" | st = "<<std::setw(5)<<(**it_dst3).status()<<" | pt = ";
	    std::cout<<std::setw(12)<<(**it_dst3).pt()<<" GeV/c | eta = "<<std::setw(12)<<(**it_dst3).eta()<<" | phi = "<<std::setw(12)<<(**it_dst3).phi();
	    // daughters of the status 3 electron, muon or tau
	    reco::GenParticle::daughters dst1 = (**it_dst3).daughterRefVector();
	    bool flag_emu_found = false; 
	    // Ask for status 1 daughters as long as no status 1 daughters are found 
	    while (flag_emu_found == false) {
	      reco::GenParticle::daughters::const_iterator it_dst1 = dst1.begin(), e_dst1 = dst1.end();
	      for (it_dst1 = dst1.begin(); it_dst1 != e_dst1; ++it_dst1) {
		//@//	std::cout<<"\n |--> Z daughter Particle  "<<std::setw(18)<<"| id = "<<std::setw(5)<<(**it_dst1).pdgId()<<" | st = "<<std::setw(5)<<(**it_dst1).status()<<" | pt = ";
		//@//	std::cout<<std::setw(12)<<(**it_dst1).pt()<<" GeV/c | eta = "<<std::setw(12)<<(**it_dst1).eta()<<" | phi = "<<std::setw(12)<<(**it_dst1).phi();
		if (((abs((**it_dst1).pdgId()) == 11) || (abs((**it_dst1).pdgId()) == 13)) && (**it_dst1).status() == 1) {
		  flag_emu_found = true;
		  ++j;
		  //@// std::cout<<"   ===> Filled E/MU at ["<<i<<"]["<<j<<"]";
		  MC_Z_PT[i][j]=(**it_dst1).pt();          MC_Z_ETA[i][j]=(**it_dst1).eta();      MC_Z_PHI[i][j]=(**it_dst1).phi();
		  MC_Z_THETA[i][j]=(**it_dst1).theta();    MC_Z_MASS[i][j]=(**it_dst1).mass();    MC_Z_PDGID[i][j]=(**it_dst1).pdgId();
		}
		else if (((abs((**it_dst1).pdgId()) == 11) || (abs((**it_dst1).pdgId()) == 13) || (abs((**it_dst1).pdgId()) == 15) || 
			  (abs((**it_dst1).pdgId()) == 24)) && (**it_dst1).status() == 2) {
		  // std::cout<<"   ===> Not Status 1 || Provided daughterRefVector";
		  dst1 = (**it_dst1).daughterRefVector();
		}
		else if (abs((**it_dst1).pdgId()) > 22 ) { // exotic or hadronic tau detected (tau --> nu W)
		  // std::cout<<"   ===> Hadronic tau || Quit Loop";
		  flag_emu_found = true; // quit loop || no e mu candidate found
		}
		else if ((abs((**it_dst1).pdgId()) == 22) && ((**it_dst1).status() == 1)) { // FSR Photon
		  int ph_j = 0; if (flag_emu_found == true) ph_j = j+2; else ph_j = j+3;
		  if ( (**it_dst1).pt() > MC_Z_PT[i][ph_j] ) { // Save only highest pt photon
		    //@// std::cout<<"   ===> Filled PHOT at ["<<i<<"]["<<ph_j<<"]";
		    MC_Z_PT[i][ph_j]=(**it_dst1).pt();          MC_Z_ETA[i][ph_j]=(**it_dst1).eta();      MC_Z_PHI[i][ph_j]=(**it_dst1).phi();
		    MC_Z_THETA[i][ph_j]=(**it_dst1).theta();    MC_Z_MASS[i][ph_j]=(**it_dst1).mass();    MC_Z_PDGID[i][ph_j]=(**it_dst1).pdgId();
		  }
		  else {
		    //@//  std::cout<<"   ===> NOT Filled, pt = "<<MC_Z_PT[i][ph_j]<<" GeV/c";
		  }
		}
		else {} // other particle than e or mu or photon
	      }
	    }
	  }
	}
	++i;
	//@//	std::cout<<"\n-----------------------------------------------------------------"<<std::endl;
      }
    }
    //@//  std::cout<<"\n"<<std::endl;




    // get 4l candidates
    i=0; 
    //j=0;
    edm::Handle<edm::View<Candidate> >  fourlCandidates;
    iEvent.getByToken(fourgenleptons_, fourlCandidates);
    for (edm::View<Candidate>::const_iterator mcIter=fourlCandidates->begin(); mcIter!=fourlCandidates->end(); ++mcIter ) {
      if (i>49 ) continue;
      //@// cout << " MC 4l Mass= " << mcIter->mass()
      //@//   << " Charge= " 
      //@//   << mcIter->daughter(0)->daughter(0)->charge() << " " 
      //@//  << mcIter->daughter(0)->daughter(1)->charge() << " "
      //@//   << mcIter->daughter(0)->daughter(2)->charge() << " " 
      //@//   << mcIter->daughter(1)->charge() << " " 
      //@//   << endl;
      MC_fourl_MASS[i][0]=mcIter->p4().mass();
      MC_fourl_PT[i][0]=mcIter->p4().pt();
      MC_fourl_PDGID[i][0]=mcIter->pdgId();

      int ii=0; // l=0;
      for (unsigned j = 0; j < mcIter->numberOfDaughters(); ++j ) {
	//cout << "j= " << j << " " << abs(mcIter->daughter(j)->pdgId()) << endl;
	if (j==0){
	  for (unsigned k = 0; k < mcIter->daughter(j)->numberOfDaughters(); ++k ) {
	    //cout << abs(mcIter->daughter(j)->daughter(k)->pdgId()) << endl;
	    if ( abs(mcIter->daughter(j)->daughter(k)->pdgId())==13 || 
		 abs(mcIter->daughter(j)->daughter(k)->pdgId())==15 || 
		 abs(mcIter->daughter(j)->daughter(k)->pdgId())==11){
	      
	      //cout << "k+1= " << k+1 << endl;
	      MC_fourl_MASS[i][k+1]=mcIter->daughter(j)->daughter(k)->p4().mass();
	      MC_fourl_PT[i][k+1]=mcIter->daughter(j)->daughter(k)->p4().pt();
	      MC_fourl_PDGID[i][k+1]=mcIter->daughter(j)->daughter(k)->pdgId();
	      
	    }
	  }
	}
	
	if (j==1) {
	  //cout << "k+1= " << mcIter->daughter(0)->numberOfDaughters()+1 << endl;
	  MC_fourl_MASS[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->p4().mass();
	  MC_fourl_PT[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->p4().pt();
	  MC_fourl_PDGID[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->pdgId();
	}

	ii++;
	
      }

      i++;
    }
    
    // get ZZ candidates
    i =0;
    edm::Handle<edm::View<Candidate> >  ZZCandidates;
    iEvent.getByToken(digenZ_, ZZCandidates);
    
    for (edm::View<Candidate>::const_iterator mcIterZZ=ZZCandidates->begin(); mcIterZZ!=ZZCandidates->end(); ++mcIterZZ ) {
      if (i>3 ) continue;
      //@// cout << "MC ZZ Mass= " << mcIterZZ->p4().mass() 
      //@//  << " and pT= " << mcIterZZ->p4().pt()  
      //@//  << endl;
      
      
      MC_ZZ_MASS[i][0]   = mcIterZZ->p4().mass();
      MC_ZZ_PT[i][0]     = mcIterZZ->p4().pt();
      MC_ZZ_ETA[i][0]    = mcIterZZ->p4().eta();
      MC_ZZ_PHI[i][0]    = mcIterZZ->p4().phi();
      MC_ZZ_THETA[i][0]  = mcIterZZ->p4().theta();
      MC_ZZ_PDGID[i][0]  = mcIterZZ->pdgId();
      
      int ii=0,l=0;
      for (unsigned j = 0; j < mcIterZZ->numberOfDaughters(); ++j ) {
	// cout << " j= " << j << " " << abs(mcIterZZ->daughter(j)->pdgId()) << endl;
	
	MC_ZZ_MASS[i][j+1] = mcIterZZ->daughter(j)->p4().mass();
	MC_ZZ_PT[i][j+1]   = mcIterZZ->daughter(j)->p4().pt();
	MC_ZZ_ETA[i][j+1]  = mcIterZZ->daughter(j)->p4().eta();
	MC_ZZ_PHI[i][j+1]  = mcIterZZ->daughter(j)->p4().phi();
	MC_ZZ_THETA[i][j+1]= mcIterZZ->daughter(j)->p4().theta();
	MC_ZZ_PDGID[i][j+1]= mcIterZZ->daughter(j)->pdgId();
	
	//cout << mcIterZZ->daughter(j)->numberOfDaughters()<< endl;
	
	int kk=0;
	for (unsigned k = 0; k < mcIterZZ->daughter(j)->numberOfDaughters(); ++k ) {
	  // cout << " k= " << k << abs(mcIterZZ->daughter(j)->daughter(k)->pdgId()) << endl;
	  if ( abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==13 || 
	       abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==15 || 
	       abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==11){
	    l=ii+j+kk+3; 	      
	    
	    //cout << " l= " << l << " " << abs(mcIterZZ->daughter(j)->daughter(k)->pdgId()) << endl;
	    MC_ZZ_MASS[i][l] = mcIterZZ->daughter(j)->daughter(k)->p4().mass();
	    MC_ZZ_PT[i][l]   = mcIterZZ->daughter(j)->daughter(k)->p4().pt();
	    MC_ZZ_ETA[i][l]  = mcIterZZ->daughter(j)->daughter(k)->p4().eta();
	    MC_ZZ_PHI[i][l]  = mcIterZZ->daughter(j)->daughter(k)->p4().phi();
	    MC_ZZ_THETA[i][l]= mcIterZZ->daughter(j)->daughter(k)->p4().theta();
	    MC_ZZ_PDGID[i][l]= mcIterZZ->daughter(j)->daughter(k)->pdgId();
	    kk++;
	  }
	}
	ii++;	    	  
      }
    }
      
    
  }
  

  // MC Higgs 
  void fillmc(const edm::Event& iEvent){
    edm::Handle<edm::View<Candidate> > Candidates;
    iEvent.getByToken(MCcollName, Candidates);
    //@//   cout << "running fillmc" << endl;
    //@// cout << "size= " << Candidates->end()-Candidates->begin() << endl;
    for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
      //@//  cout << "Filling MC variables" << endl;
      const reco::Candidate& theParticle = (*cand);
      SetMCValues(theParticle,0);   
      int i=0,l=0;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	if (cand->daughter(j)->pdgId()==23){	  
	  const reco::Candidate& theParticle = (*cand->daughter(j));
	  SetMCValues(theParticle,j+1);
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if ( abs(cand->daughter(j)->daughter(k)->pdgId())==13 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==15 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==11){
	      l=i+j+k+3; 	      
	      const reco::Candidate& theParticle = (*cand->daughter(j)->daughter(k));
	      SetMCValues(theParticle,l);
	    }
	  }
	  i++;
	}
      }
    }    
    //fillMCCP(iEvent);
  }

  
  
  struct SortCandByDecreasingPt {
    bool operator()( const reco::Candidate &c1, const reco::Candidate &c2) const {
      return c1.pt() > c2.pt();
    }
  };
  

  
  void fillAdditionalRECO(const edm::Event& iEvent){
   
    //Matching ZtoMuMu:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZMuMu;
    iEvent.getByToken(goodZtoMuMuMCMatch_, GenParticlesMatchZMuMu);
    if (GenParticlesMatchZMuMu.isValid() ){
      //@//  cout << "The matched map collection has size= " <<  GenParticlesMatchZMuMu->size() << endl;
    }
    
    //Matching ZtoEE:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZEE;
    iEvent.getByToken(goodZtoEEMCMatch_, GenParticlesMatchZEE);
    if (GenParticlesMatchZEE.isValid() ){
      //@//  cout << "The matched map collection has size= " <<  GenParticlesMatchZEE->size() << endl;
    }

    // di-leptons OS
    leptonscands_Z0->clear();
    leptonscands_Z1->clear();
   
    
    //@// cout << "RECOcollNameZ size " << RECOcollNameZ.size() << endl;
    //@// cout << "RECOcollNameZ string" << RECOcollNameZ.at(0) << endl;
    for (unsigned int i=0; i<RECOcollNameZ.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZ;
      if(i==0) iEvent.getByToken(zToMuMu, CandidatesZ); 
      else if(i==1) iEvent.getByToken(zToEE, CandidatesZ);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZ->begin();cand != CandidatesZ->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_ZMM_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_ZMM_PT[0][kk]=cand->p4().pt();
	if (i==0) RECO_ZMM_ETA[0][kk]=cand->p4().eta();
	if (i==0) RECO_ZMM_PHI[0][kk]=cand->p4().phi();
	if (i==1) RECO_ZEE_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_ZEE_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_ZEE_ETA[0][kk]=cand->p4().eta();
	if (i==1) RECO_ZEE_PHI[0][kk]=cand->p4().phi();
	//@//	cout << "di-lepton candidate of type=" << RECOcollNameZ.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  //@//  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEE_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEE_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEE_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z1->push_back( cand->daughter(j)->clone());
	  }
	}
	
	if (i==0 && fillMCTruth==true){
	  // Matching ZtoMuMu
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZMuMu = (*GenParticlesMatchZMuMu)[Ref]; 
	  if (!genrefZMuMu.isNull()){
	    //@//  cout << "Gen Z with pT= " << genrefZMuMu->p4().pt() << " and mass="<< genrefZMuMu->p4().mass() << endl;	  
	    RECOzMuMu_MatchingMCTruth[kk]= true;
	    RECOzMuMu_MatchingMCpT[kk]= genrefZMuMu->p4().pt();
	    RECOzMuMu_MatchingMCmass[kk]= genrefZMuMu->p4().mass();
	    RECOzMuMu_MatchingMCEta[kk]= genrefZMuMu->p4().eta();
	    RECOzMuMu_MatchingMCPhi[kk]= genrefZMuMu->p4().phi();
	  }
	}
	if (i==1 && fillMCTruth==true){
	  // Matching ZtoEE
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZEE = (*GenParticlesMatchZEE)[Ref]; 
	  if (!genrefZEE.isNull()){
	    //@//  cout << "Gen Z with pT= " << genrefZEE->p4().pt() << " and mass="<< genrefZEE->p4().mass() << endl;	  
	    RECOzEE_MatchingMCTruth[kk]= true;
	    RECOzEE_MatchingMCpT[kk]= genrefZEE->p4().pt();
	    RECOzEE_MatchingMCmass[kk]= genrefZEE->p4().mass();
	    RECOzEE_MatchingMCEta[kk]= genrefZEE->p4().eta();
	    RECOzEE_MatchingMCPhi[kk]= genrefZEE->p4().phi();
	  }
	}


	kk++;
      }
    }

    // di-leptons SS and cross-leptons
    leptonscands_Zss0->clear();
    leptonscands_Zss1->clear();
    leptonscands_Zcross->clear();


    //@//  cout << "RECOcollNameZss size " << RECOcollNameZss.size() << endl; 
    for (unsigned int i=0; i<RECOcollNameZss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZss;
      if(i==0)   iEvent.getByToken(zToMuMussmerge, CandidatesZss);
      if(i==1)   iEvent.getByToken(zToEEssmerge, CandidatesZss);
      if(i==2)   iEvent.getByToken(zToCrossLeptons, CandidatesZss); 
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZss->begin();cand != CandidatesZss->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0)  RECO_ZMMss_MASS[kk]=cand->p4().mass();
	if (i==0)  RECO_ZMMss_PT[0][kk]=cand->p4().pt();
	if (i==0)  RECO_ZMMss_ETA[0][kk]=cand->p4().eta();
	if (i==0)  RECO_ZMMss_PHI[0][kk]=cand->p4().phi();

	if (i==1)  RECO_ZEEss_MASS[kk]=cand->p4().mass();
	if (i==1)  RECO_ZEEss_PT[0][kk]=cand->p4().pt();
	if (i==1)  RECO_ZEEss_ETA[0][kk]=cand->p4().eta();
	if (i==1)  RECO_ZEEss_PHI[0][kk]=cand->p4().phi();

	if (i==2)  RECO_ZEM_MASS[kk]=cand->p4().mass();
	if (i==2)  RECO_ZEM_PT[0][kk]=cand->p4().pt();
	if (i==2)  RECO_ZEM_ETA[0][kk]=cand->p4().eta();
	if (i==2)  RECO_ZEM_PHI[0][kk]=cand->p4().phi();

	//@//cout << "di-lepton candidate of type=" << RECOcollNameZss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  //@//  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMMss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMMss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMMss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEEss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEEss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEEss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_ZEM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zcross->push_back( cand->daughter(j)->clone());	  
	  }
	}
	kk++;
      }
    }


    // di-leptons ALL
    leptonscands_DiLep->clear();

    edm::Handle<edm::View<Candidate> > CandidatesDiLep;
    iEvent.getByToken(RECOcollNameDiLep_, CandidatesDiLep); 
    int kkk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesDiLep->begin();cand != CandidatesDiLep->end(); ++ cand ) { 
      if (kkk>49) continue;
      RECO_DiLep_MASS[kkk]=cand->p4().mass();
      RECO_DiLep_PT[0][kkk]=cand->p4().pt();
      RECO_DiLep_ETA[0][kkk]=cand->p4().eta();
      RECO_DiLep_PHI[0][kkk]=cand->p4().phi();
      
      //@//  cout << "di-lepton candidate of type=" << RECOcollNameDiLep.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	//@//	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl;
	RECO_DiLep_PT[j+1][kkk]=cand->daughter(j)->p4().pt();
	RECO_DiLep_ETA[j+1][kkk]=cand->daughter(j)->p4().eta();
	RECO_DiLep_PHI[j+1][kkk]=cand->daughter(j)->p4().phi();
	leptonscands_DiLep->push_back( cand->daughter(j)->clone());
      }
      kkk++;
    }



    // MuMuMuMu
    int i=1; 
    leptonscands_MMMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByToken(RECOcollNameMMMM_, CandidatesMMMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHMMMM;
    iEvent.getByToken(goodHiggsTozzToMMMMMCMatch_, GenParticlesMatchHMMMM);
    if (GenParticlesMatchHMMMM.isValid()){
      //@//  cout << "The matched map collection has size= " <<  GenParticlesMatchHMMMM->size() << endl;
    }

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) { 
      if (kk>99) break;
      RECO_MMMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMMM_PT[i-1][kk]=cand->p4().pt();
      RECO_MMMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_MMMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	} 
	l++;
      }


      // Matching HiggsToMMMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesMMMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzMMMM = (*GenParticlesMatchHMMMM)[Ref]; 
	if (!genrefHzzMMMM.isNull()){
	  //@// cout << "Gen Z with pT= " << genrefHzzMMMM->p4().pt() << " and mass="<< genrefHzzMMMM->p4().mass() << endl;	  
	  RECOHzzMMMM_MatchingMCTruth[kk]= true;
	  RECOHzzMMMM_MatchingMCpT[kk]= genrefHzzMMMM->p4().pt();
	  RECOHzzMMMM_MatchingMCmass[kk]= genrefHzzMMMM->p4().mass();
	  RECOHzzMMMM_MatchingMCEta[kk]= genrefHzzMMMM->p4().eta();
	  RECOHzzMMMM_MatchingMCPhi[kk]= genrefHzzMMMM->p4().phi();
	}
      }

      kk++;
    }

    
    // EEEE
    i=1;
    leptonscands_EEEE->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByToken(RECOcollNameEEEE, CandidatesEEEE);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEEE;
    iEvent.getByToken(goodHiggsTozzToEEEEMCMatch_, GenParticlesMatchHEEEE);
    if (GenParticlesMatchHEEEE.isValid()){
      //@// cout << "The matched map collection has size= " <<  GenParticlesMatchHEEEE->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEEE_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEEE_PT[i-1][kk]=cand->p4().pt();
      RECO_EEEE_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEEE_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEEE_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_EEEE->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEEE
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEEE,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEEE = (*GenParticlesMatchHEEEE)[Ref]; 
	if (!genrefHzzEEEE.isNull()){
	  //@// cout << "Gen Z with pT= " << genrefHzzEEEE->p4().pt() << " and mass="<< genrefHzzEEEE->p4().mass() << endl;	  
	  RECOHzzEEEE_MatchingMCTruth[kk]= true;
	  RECOHzzEEEE_MatchingMCpT[kk]= genrefHzzEEEE->p4().pt();
	  RECOHzzEEEE_MatchingMCmass[kk]= genrefHzzEEEE->p4().mass();
	  RECOHzzEEEE_MatchingMCEta[kk]= genrefHzzEEEE->p4().eta();
	  RECOHzzEEEE_MatchingMCPhi[kk]= genrefHzzEEEE->p4().phi();
	}
      }


      kk++;
    }


    // EEMM
    i=1;
    leptonscands_EEMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByToken(RECOcollNameEEMM, CandidatesEEMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEMM;
    iEvent.getByToken(goodHiggsTozzToEEMMMCMatch_, GenParticlesMatchHEEMM);
    if (GenParticlesMatchHEEMM.isValid()){
      //@//  cout << "The matched map collection has size= " <<  GenParticlesMatchHEEMM->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEMM_PT[i-1][kk]=cand->p4().pt();
      RECO_EEMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();

	  leptonscands_EEMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEMM = (*GenParticlesMatchHEEMM)[Ref]; 
	if (!genrefHzzEEMM.isNull()){
	  //@// cout << "Gen Z with pT= " << genrefHzzEEMM->p4().pt() << " and mass="<< genrefHzzEEMM->p4().mass() << endl;	  
	  RECOHzzEEMM_MatchingMCTruth[kk]= true;
	  RECOHzzEEMM_MatchingMCpT[kk]= genrefHzzEEMM->p4().pt();
	  RECOHzzEEMM_MatchingMCmass[kk]= genrefHzzEEMM->p4().mass();
	  RECOHzzEEMM_MatchingMCEta[kk]= genrefHzzEEMM->p4().eta();
	  RECOHzzEEMM_MatchingMCPhi[kk]= genrefHzzEEMM->p4().phi();
	}
      }

      kk++;
    }
  }


  void SetMCValues(const reco::Candidate& cand, int nMC){
    MC_E[nMC]     = cand.p4().energy();
    MC_PT[nMC]    = cand.p4().pt();
    MC_ETA[nMC]   = cand.p4().eta();
    MC_THETA[nMC] = cand.p4().theta();
    MC_PHI[nMC]   = cand.p4().phi();
    MC_MASS[nMC]  = cand.p4().mass();
    MC_PDGID[nMC] = cand.pdgId();
  }

 
 
  
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll){
    
    bool found=false;
    for( reco::CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
      if ((fabs(pp->p4().mass()-mass)  <0.001 ) &&
	  (fabs(pp->p4().pt()  -pt)    <0.001 ) &&
	  (fabs(pp->charge()   -charge)<0.001 )  ){
	found=true;
	//std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
      }
    }
    return found;
  }
  
  bool matchParticle(double mass, double pt, int charge, const reco::Candidate *c1){
    
    bool found=false;
    if ((fabs(c1->p4().mass()-mass)  <0.001 ) &&
	(fabs(c1->p4().pt()  -pt)    <0.001 ) &&
	(fabs(c1->charge()   -charge)<0.001 )  ){
      found=true;
      //std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
    }    
    return found;
  }
  
  
  void fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){
   
 
    // Supercluster collection
    edm::Handle<reco::SuperClusterCollection> clusters;
    iEvent.getByLabel(clusterCollectionTag_,clusters);
    

    // EG isolation
    //edm::Handle<reco::GsfElectronRefVector> EleRefs;
    //iEvent.getByLabel(electronEgmTag_, EleRefs);
    edm::Handle<edm::View<pat::Electron> > EleRefs;
    iEvent.getByToken(electronEgmTag_, EleRefs);
    
    
    
    //Electron ID MVA Trig and Non Trig
    //edm::Handle<edm::ValueMap<float> >  mvaTrigV0map;
    //iEvent.getByToken(mvaTrigV0MapTag_, mvaTrigV0map);
    
    // edm::Handle<edm::ValueMap<float> >  mvaNonTrigV0map;
    //iEvent.getByToken(mvaNonTrigV0MapTag_, mvaNonTrigV0map);


    // primary vertex
   

 
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    //@// cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  //@//
	  /*cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;*/
	  break;
	}
      }
    } 
    
    int index=0;
    RECO_NELE=EleRefs->size();
   
    ///* to add matching informations
    // Matching
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchEle;
    iEvent.getByToken(goodElectronMCMatch_, GenParticlesMatchEle);
    edm::Handle<reco::CandidateCollection > CollEle;
    iEvent.getByToken(myElectrons_, CollEle);

    if (GenParticlesMatchEle.isValid()){
      //@// cout << endl<< "Electrons:"<<endl<<"The reco collection to be matched has size= " <<  CollEle->size() << endl;
      //@// cout << "The matched map collection has size= " <<  GenParticlesMatchEle->size() << endl;
    }//*/
    //

    for (edm::View<pat::Electron>::const_iterator cand = EleRefs->begin(); cand != EleRefs->end(); ++cand) {
      
      if(index>99) break;
      
      edm::Ref<edm::View<pat::Electron> > eletrackref(EleRefs,index);
      // edm::Ref<edm::View<pat::Electron> > eletrackrefv(VertEleCandidates,index);
     
    
      /* // kinematic //Reham comment to check the new corrected pt
	 RECOELE_E[index]       = cand->p4().energy();
	 RECOELE_PT[index]      = cand->p4().pt();
	 
	 
	 RECOELE_P[index]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
	 RECOELE_ETA[index]     = cand->p4().eta();
	 RECOELE_THETA[index]   = cand->p4().theta();
	 RECOELE_PHI[index]     = cand->p4().phi();
	 RECOELE_MASS[index]    = cand->p4().mass();
	 RECOELE_CHARGE[index]  = cand->charge();*/
      
      ///// For the elctron energy correction
      //Reham
      
      float corrEt = cand->et() * cand->userFloat("ecalTrkEnergyPostCorr")/cand->energy();
      auto corrP4  = cand->p4() * cand->userFloat("ecalTrkEnergyPostCorr")/cand->energy();
      cout<<"corrEt= "<<corrEt<<"corrP4 = "<<corrP4<<endl;
     
      cout<<"TESTTTT Kinematics corrPT = "<<corrP4.pt()<<" corrE = "<<corrP4.energy()<<"corrEta = "<<corrP4.eta()<<"corr phi = "<<corrP4.phi()<<"corrTheta = "<<corrP4.theta()<<"corrcharge = "<<cand->charge()<<endl;

      cout<<"TESTTTT Kinematics BEFORE corr. PT = "<< cand->p4().pt()<<" E = "<<cand->p4().energy()<<"Eta = "<<cand->p4().eta()<<" phi = "<<cand->p4().phi()<<"Theta = "<<cand->p4().theta()<<"charge = "<<cand->charge()<<endl;
 

      if (corrP4.pt() <7 )continue;
      if (corrP4.eta() >2.5)continue;
     
      RECOELE_E[index]       = corrP4.energy();
      RECOELE_PT[index]      = corrP4.pt();   
      RECOELE_P[index]=sqrt(corrP4.px()*corrP4.px()+corrP4.py()*corrP4.py()+corrP4.pz()*corrP4.pz());
      RECOELE_ETA[index]     = corrP4.eta();
      RECOELE_THETA[index]   = corrP4.theta();
      RECOELE_PHI[index]     = corrP4.phi();
      RECOELE_MASS[index]    = corrP4.mass();
      RECOELE_CHARGE[index]  = cand->charge();
      RECOELE_PT_uncorr[index]   = cand->pt();

	//Reham error in PT and some systematic variables 
      
      // cout<<"ecalTrkEnergyErrPostCorr = "<<cand->userFloat("ecalTrkEnergyErrPostCorr")<<endl;
      // cout<<"ecalTrkEnergyPreCorr"<<cand->userFloat("ecalTrkEnergyPreCorr")<<" ecalTrkEnergyErrPreCorr = "<<cand->userFloat("ecalTrkEnergyErrPreCorr");
   
      RECOELE_ecalTrkEnergyPreCorr[index]  = cand->userFloat("ecalTrkEnergyPreCorr");
      RECOELE_ecalTrkEnergyErrPreCorr[index]  = cand->userFloat("ecalTrkEnergyErrPreCorr");
      RECOELE_ecalTrkEnergyErrPostCorr[index]  = cand->userFloat("ecalTrkEnergyErrPostCorr");
      RECOELE_energyScaleValue[index]          = cand->userFloat("energyScaleValue");
      RECOELE_energySigmaValue[index]          = cand->userFloat("energySigmaValue");
      RECOELE_energyScaleUp[index]          = cand->userFloat("energyScaleUp");
      RECOELE_energyScaleDown[index]          = cand->userFloat("energyScaleDown");
      RECOELE_energyScaleStatUp[index]          = cand->userFloat("energyScaleStatUp");
      RECOELE_energyScaleStatDown[index]          = cand->userFloat("energyScaleStatDown");
      RECOELE_energyScaleSystUp[index]          = cand->userFloat("energyScaleSystUp");
      RECOELE_energyScaleSystDown[index]          = cand->userFloat("energyScaleSystDown");
      RECOELE_energyScaleGainUp[index]          = cand->userFloat("energyScaleGainUp");
      RECOELE_energyScaleGainDown[index]          = cand->userFloat("energyScaleGainDown");
      RECOELE_energyScaleEtUp[index]          = cand->userFloat("energyScaleEtUp");
      RECOELE_energyScaleEtDown[index]          = cand->userFloat("energyScaleEtDown");
      RECOELE_energySigmaUp[index]               = cand->userFloat("energySigmaUp");
      RECOELE_energySigmaDown[index]          = cand->userFloat("energySigmaDown");
      RECOELE_energySigmaPhiUp[index]          = cand->userFloat("energySigmaPhiUp");
      RECOELE_energySigmaPhiDown[index]          = cand->userFloat("energySigmaPhiDown");
      RECOELE_energySigmaRhoUp[index]          = cand->userFloat("energySigmaRhoUp");
      RECOELE_energySigmaRhoDown[index]          = cand->userFloat("energySigmaRhoDown");

      ///// */    

      
        /* std::cout << "--kinematic:" 		 */
	/* 	<< "  pT="     << RECOELE_PT[index]   */
	/* 	<< "  E="      << RECOELE_E[index]   */
	/* 	<< "  p="      << RECOELE_P[index]  */
	/* 	<< "  eta="    << RECOELE_ETA[index]  */
	/* 	<< "  theta="  << RECOELE_THETA[index]  */
	/* 	<< "  phi="    << RECOELE_PHI[index]  */
	/* 	<< "  mass="   << RECOELE_MASS[index]  */
	/* 	<< "  charge=" << RECOELE_CHARGE[index] */
	/* 	<< "  pT_uncorr ="     << RECOELE_PT_uncorr[index]  */
	/* 	<< std::endl; */

	/* std::cout << "Electron syst ==" */
	/* 	  <<"\n RECOELE_ecalTrkEnergyPreCorr[index] = "<<RECOELE_ecalTrkEnergyPreCorr[index] */
	/* 	  <<"\n RECOELE_ecalTrkEnergyErrPreCorr[index] = "<<RECOELE_ecalTrkEnergyErrPreCorr[index] */
	/* 	  <<"\n RECOELE_ecalTrkEnergyErrPostCorr[index]"<< RECOELE_ecalTrkEnergyErrPostCorr[index] */
	/* 	  <<"\n RECOELE_energyScaleValue[index]"<<RECOELE_energyScaleValue[index] */
	/* 	  <<"\n RECOELE_energySigmaValue[index]"<<RECOELE_energySigmaValue[index] */
	/* 	  <<"\n RECOELE_energyScaleUp[index]"<<RECOELE_energyScaleUp[index] */
	/* 	  <<"\n RECOELE_energyScaleDown[index]"<<RECOELE_energyScaleDown[index] */
	/* 	  <<"\n RECOELE_energyScaleStatUp[index]"<< RECOELE_energyScaleStatUp[index] */
	/* 	  <<"\n RECOELE_energyScaleStatDown[index]"<< RECOELE_energyScaleStatDown[index]      */
	/* 	  <<"\n RECOELE_energyScaleSystUp[index]"<<RECOELE_energyScaleSystUp[index]         */
	/* 	  <<"\n RECOELE_energyScaleSystDown[index]"<<RECOELE_energyScaleSystDown[index]        */
	/* 	  <<"\n RECOELE_energyScaleGainUp[index]"<<RECOELE_energyScaleGainUp[index]          */
	/* 	  <<"\n RECOELE_energyScaleGainDown[index]"<<RECOELE_energyScaleGainDown[index]        */
	/* 	  <<"\n RECOELE_energyScaleEtUp[index]"<< RECOELE_energyScaleEtUp[index]        */
	/* 	  <<"\n RECOELE_energyScaleEtDown[index]"<< RECOELE_energyScaleDown[index] */
	/* 	  <<"\n RECOELE_energySigmaUp[index]"<< RECOELE_energySigmaUp[index]        */
	/* 	  <<"\n RECOELE_energySigmaDown[index]"<< RECOELE_energySigmaDown[index]        */
	/* 	  <<"\n RECOELE_energySigmaPhiUp[index]"<<RECOELE_energySigmaPhiUp[index]          */
	/* 	  <<"\n RECOELE_energySigmaPhiDown[index]"<<RECOELE_energySigmaPhiDown[index]         */
	/* 	  <<"\n RECOELE_energySigmaRhoUp[index]"<<RECOELE_energySigmaRhoUp[index]        */
	/* 	  <<"\n RECOELE_energySigmaRhoDown[index]"<<RECOELE_energySigmaRhoDown[index]        */
	/* 	  <<std::endl; */
	
	
	cout<<"#########Electron step 1 "<<endl;

     // Global variables
      RECOELE_isEcalDriven[index]    = cand->ecalDrivenSeed();
      RECOELE_isTrackerDriven[index] = cand->trackerDrivenSeed();

     //@//
       std::cout << "\n Electron in the event: "
	 << "  isEcalDriven="    << RECOELE_isEcalDriven[index]  
	 << "  isTrackerDriven=" << RECOELE_isTrackerDriven[index] 
	 << std::endl;

      ///@@@/// Electron ID  REHAM
      

     cout << cand->electronID("mvaEleID-Fall17-iso-V2-wpHZZ") << " " << std::abs(cand->dB(pat::Electron::PV3D))/cand->edB(pat::Electron::PV3D) << endl;
     cout<<">>>>>MVA value " <<cand->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values")<<endl;

     RECOELE_ID[index] = cand->electronID("mvaEleID-Fall17-iso-V2-wpHZZ");
     RECOELE_mvaNonTrigV0[index] =  cand->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");

     cout<<"RECOELE_ID[index] = "<<RECOELE_ID[index]<<" RECOELE_mvaNonTrigV0[index] = "<< RECOELE_mvaNonTrigV0[index]<<endl;

     //@@@///  Electron SIP Reham

     RECOELE_SIP[index]= std::abs(cand->dB(pat::Electron::PV3D))/cand->edB(pat::Electron::PV3D);
     RECOELE_IP[index]= std::abs(cand->dB(pat::Electron::PV3D));
     RECOELE_IPERROR[index]= cand->edB(pat::Electron::PV3D);

     // cout<<"ID = "<<RECOELE_ID[index]<<" SIP = "<<RECOELE_SIP[index]<<endl;
      
      // SuperCluster
      reco::SuperClusterRef sclRef = cand->superCluster();
      math::XYZPoint sclPos = cand->superClusterPosition();

      //if(!cand->ecalDrivenSeed() && cand->trackerDrivenSeed())
      //clRef = cand->pflowSuperCluster();
      double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
      double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
      ele_sclRawE[index]   = sclRef->rawEnergy() ;
      RECOELE_scl_E[index]   = sclRef->energy() ;
      RECOELE_scl_Et[index]  = sclRef->energy()*(Rt/R) ;
      RECOELE_scl_Eta[index] = sclRef->eta() ;
      RECOELE_scl_Phi[index] = sclRef->phi() ;
      ele_sclX[index]  = sclPos.X();
      ele_sclY[index] =  sclPos.Y();
      ele_sclZ[index] =  sclPos.Z();
			
      //@//
      /* std::cout << "--supercluster properties:"
	<< "  E_sc="    << RECOELE_scl_E[index]
	<< "  Et_sc="   << RECOELE_scl_Et[index]
	<< "  Eeta_sc=" << RECOELE_scl_Eta[index]
	<< "  phi_sc="  << RECOELE_scl_Phi[index]
	<< std::endl;*/

      //Covariance matrix
      TMatrixDSym bigCov;
      double dp = 0.;
      if (cand->ecalDriven()) {
	dp = cand->p4Error(pat::Electron::P4_COMBINATION);	
	RECOELE_PTError[index]=dp;
      }
      
      else {
	// Parametrization from Claude Charlot, 
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZAnalysis/AnalysisStep/src/ZZMassErrors.cc?revision=1.2&view=markup
#if CMSSW_VERSION<500
	double ecalEnergy = cand->ecalEnergy() ;
#else
	double ecalEnergy = cand->correctedEcalEnergy() ;
#endif
	double err2 = 0.0;
	if (cand->isEB()) {
	  err2 += (5.24e-02*5.24e-02)/ecalEnergy;  
	  err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.00e-02*1.00e-02;
	} else if (cand->isEE()) {
	  err2 += (1.46e-01*1.46e-01)/ecalEnergy;  
	  err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.94e-03*1.94e-03;
	}
	dp = ecalEnergy * sqrt(err2);
	RECOELE_PTError[index]=dp;
      }
      // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
      //            [ Px/P  ]
      //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
      //            [ Pz/P  ]
      AlgebraicMatrix31 ptop3;
      ptop3(0,0) = cand->px()/cand->p();
      ptop3(1,0) = cand->py()/cand->p();
      ptop3(2,0) = cand->pz()/cand->p();
      AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
      for (int i = 0; i < 3; ++i) { 
	for (int j = 0; j < 3; ++j) {
	  RECOELE_COV[index][i][j]=mat[i][j];
	}
      }
      
      cout<<"#########Electron step 2 "<<endl;
      cout<<"RECOELE_PTError[index]= "<<RECOELE_PTError[index]<<endl;
      cout<<"cand ecal energy = "<<cand->correctedEcalEnergy()<<endl;

      // Egamma isolation
      //RECOELE_EGMTRACKISO[index]=(*egmisoTkelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMECALISO[index]=(*egmisoEcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMHCALISO[index]=(*egmisoHcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMX[index]=( (*egmisoTkelemap)[eletrackref] + (*egmisoEcalelemap)[eletrackref] + (*egmisoHcalelemap)[eletrackref])/eletrackref->pt();

      // temporary solution for the problem of rechits
      RECOELE_EGMECALISO[index]=(eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt();
      RECOELE_EGMHCALISO[index]=(eletrackref->dr03HcalTowerSumEt())/eletrackref->pt();
      RECOELE_EGMX[index]=RECOELE_EGMTRACKISO[index]+RECOELE_EGMECALISO[index]+RECOELE_EGMHCALISO[index];
   
      //@//   
      /* std::cout << "--EG isolation: electron"
		<< "  X="     << RECOELE_EGMX[index]
		<< "  Track=" << RECOELE_EGMTRACKISO[index]
		<< "  Ecal= " << RECOELE_EGMECALISO[index]
		<< "  Hcal=" << RECOELE_EGMHCALISO[index]
		<< std::endl;*/


      // PF isolation  (R unknown, TODO)
      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt);
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt);
      RECOELE_PFPUchAllPart[index]= (cand->pfIsolationVariables().sumPUPt);
      RECOELE_PFchAllPart[index]= (cand->pfIsolationVariables().sumChargedParticlePt);

      //RECOELE_PFchAllPart[index]      =  (*isoPFChargedAllelemap)[eletrackref];
      //RECOELE_PFchHad[index]          =  (*isoPFChargedelemap)[eletrackref];
      //RECOELE_PFneuHad[index]         =  (*isoPFNeutralelemap)[eletrackref];
      //RECOELE_PFphoton[index]         =  (*isoPFGammaelemap)[eletrackref];
     // RECOELE_PFPUchAllPart[index]    =  (*isoPFPUelemap)[eletrackref];

      RECOELE_PFX_dB[index]           =  (RECOELE_PFchHad[index] + max(RECOELE_PFphoton[index]+RECOELE_PFneuHad[index]-0.5*RECOELE_PFPUchAllPart[index],0.))/eletrackref->pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                       EffectiveArea = 0.29;
      }
      else { 
	// 7_4_X use eta 
	//if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
	//if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
	//if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
	//if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
	//if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;

	// 7_6_X use eta supercluster   
//	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
        // 8_0_X  
        if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
      }

      cout<<"#########Electron step 3 "<<endl;

      //RECOELE_PFX_rho[index]=(RECOELE_PFchHad[index]+
      //		      max( (RECOELE_PFneuHad[index]+RECOELE_PFphoton[index]-max(RHO_ele,0.0)*EffectiveArea),0.0) )/cand->p4().pt(); 

      RECOELE_PFX_rho[index]= (cand->pfIsolationVariables().sumChargedHadronPt+
			       std::max(
					cand->pfIsolationVariables().sumPhotonEt+
					cand->pfIsolationVariables().sumNeutralHadronEt-
					max(RHO_ele,0.0)*EffectiveArea,
					0.0)
			       )/cand->p4().pt();

      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt); 
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt); 
      //RECOELE_PFsumPUPt[index]= (cand->pfIsolationVariables().sumPUPt); 
      RECOELE_PFPUchAllPart[index]=(cand->pfIsolationVariables().sumPUPt);   




      cout<<"#########Electron step 4 "<<endl;

      // GsfTrack
      RECOELE_gsftrack_NPixHits[index]  = cand->gsfTrack()->hitPattern().numberOfValidPixelHits();
      RECOELE_gsftrack_NStripHits[index]= cand->gsfTrack()->hitPattern().numberOfValidStripHits();
      RECOELE_gsftrack_chi2[index]      = cand->gsfTrack()->normalizedChi2();
      RECOELE_gsftrack_dxyB[index]      = cand->gsfTrack()->dxy(bs.position()) ;
      RECOELE_gsftrack_dxy[index]       = cand->gsfTrack()->dxy(pVertex);
      RECOELE_gsftrack_dxyError[index]  = cand->gsfTrack()->dxyError();
      RECOELE_gsftrack_dzB[index]       = cand->gsfTrack()->dz(bs.position());
      RECOELE_gsftrack_dz[index]        = cand->gsfTrack()->dz(pVertex);
      RECOELE_gsftrack_dzError[index]   = cand->gsfTrack()->dzError();

      //Conversion variables
      RECOELE_gsftrack_losthits[index]  = cand->gsfTrack()->numberOfLostHits();
      RECOELE_gsftrack_validhits[index] = cand->gsfTrack()->numberOfValidHits();
      // RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);//original
       RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);//Reham
      // RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->hitPattern().numberOfValidHits();//Reham
      // RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->trackerExpectedHitsInner().numberOfHits();//Nicola Not work
     
       //@//
       /*  std::cout << "--gfstrack properties: " 
	<< "  nPixhits=" << RECOELE_gsftrack_NPixHits[index]
	<< "  nStriphits=" << RECOELE_gsftrack_NStripHits[index]
	<< "  chi2="     << RECOELE_gsftrack_chi2[index] 
	<< "  dxy="      << RECOELE_gsftrack_dxy[index] 
	<< "  dxyB="     << RECOELE_gsftrack_dxyB[index] 
	<< "  dxyError=" << RECOELE_gsftrack_dxyError[index] 
	<< "  dz="       << RECOELE_gsftrack_dz[index]
	<< "  dzB="      << RECOELE_gsftrack_dzB[index]
        << "  dzError="  << RECOELE_gsftrack_dzError[index] 
	<< "  losthits=" << RECOELE_gsftrack_losthits[index] 
	<< "  validhits="<< RECOELE_gsftrack_validhits[index] 
        << "  innerhits="<< RECOELE_gsftrack_expected_inner_hits[index] 
	<< std::endl;*/

      // Track-Cluster matching attributes
      RECOELE_ep[index]             = cand->eSuperClusterOverP(); 
      RECOELE_eSeedp[index]         = cand->eSeedClusterOverP();     
      RECOELE_eSeedpout[index]      = cand->eSeedClusterOverPout(); 
      RECOELE_eElepout[index]       = cand->eEleClusterOverPout();        
      
      RECOELE_deltaEtaIn [index]  = cand->deltaEtaSuperClusterTrackAtVtx(); 
      RECOELE_deltaEtaSeed[index] = cand->deltaEtaSeedClusterTrackAtCalo(); 
      RECOELE_deltaEtaEle[index]  = cand->deltaEtaEleClusterTrackAtCalo();  
      RECOELE_deltaPhiIn[index]   = cand->deltaPhiSuperClusterTrackAtVtx();
      RECOELE_deltaPhiSeed[index] = cand->deltaPhiSeedClusterTrackAtCalo(); 
      RECOELE_deltaPhiEle[index]  = cand->deltaPhiEleClusterTrackAtCalo() ;  

      cout<<"#########Electron step 5 "<<endl;

      //@// 
     /* std::cout << "--track-cluster matching: " 
	<< "  eSC_p="        << RECOELE_ep[index] 
	<< "  eSeed_p="      << RECOELE_eSeedp[index] 
	<< "  eSeed_pout="   << RECOELE_eSeedpout[index] 
	<< "  eEle_pout="    << RECOELE_eElepout[index] 
	<< "  deltaEtaIn="   << RECOELE_deltaEtaIn[index] 
        << "  deltaEtaSeed=" << RECOELE_deltaEtaSeed[index] 
        << "  deltaEtaEle="  << RECOELE_deltaEtaEle[index] 
        << "  deltaPhiIn="   << RECOELE_deltaPhiIn[index] 
        << "  deltaPhiSeed=" << RECOELE_deltaPhiSeed[index] 
	<< "  deltaPhiEle="  << RECOELE_deltaPhiEle[index] 
	<< std::endl;*/
      
       // Fiducial flags
      if (cand->isEB()) RECOELE_isbarrel[index] = 1 ; 
      else  RECOELE_isbarrel[index] = 0 ;
      if (cand->isEE()) RECOELE_isendcap[index] = 1 ; 
      else  RECOELE_isendcap[index] = 0 ;
      if (cand->isEBEtaGap())  RECOELE_isEBetaGap[index]  = 1 ;  
      if (cand->isEBPhiGap())  RECOELE_isEBphiGap[index]  = 1 ;  
      if (cand->isEEDeeGap())  RECOELE_isEEdeeGap[index]  = 1 ;  
      if (cand->isEERingGap()) RECOELE_isEEringGap[index] = 1 ;
      if (cand->isGap())       RECOELE_isGap[index] = 1 ;

      //@//
      /*  std::cout << "--fiducial flags: " 
	<< "  isEB="        << RECOELE_isbarrel[index] 
	<< "  isEBetagap="  << RECOELE_isEBetaGap[index] 
	<< "  isEBphigap="  << RECOELE_isEBphiGap[index] 
	<< "  isEE="        << RECOELE_isendcap[index] 
        << "  isEEdeegap="  << RECOELE_isEEdeeGap[index] 
	<< "  isEEringgap=" << RECOELE_isEEringGap[index] 
        << "  isGap="       << RECOELE_isGap[index]
	<< std::endl;*/
      

      // Shower shape
      RECOELE_sigmaIetaIeta[index] = cand->sigmaIetaIeta() ; 
      RECOELE_sigmaEtaEta[index]   = cand->sigmaEtaEta() ;
      RECOELE_e15[index]           = cand->e1x5() ;
      RECOELE_e25max[index]        = cand->e2x5Max() ;
      RECOELE_e55[index]           = cand->e5x5() ;
      RECOELE_he[index]            = cand->hadronicOverEm() ;
      // RECOELE_r9[index]            = cand->r9() ;

      //@//
      /*  std::cout << "--shower shape:" 
	<< "  sigmaIetaIeta=" << RECOELE_sigmaIetaIeta[index] 
	<< "  sigmaEtaEta=" << RECOELE_sigmaEtaEta[index]  
	<< "  e15=" << RECOELE_e15[index] 
	<< "  e25max=" << RECOELE_e25max[index] 
	<< "  e55=" << RECOELE_e55[index]  
	<< "  he=" << RECOELE_he[index] 
	<< "  r9=" << RECOELE_r9[index] 
	<< std::endl;*/
      
      // also variables on hcal
      // Isolation
      // Particle flow
      //RECOELE_mva[index] = cand->mva();
      //@//  std::cout << "--PF: mva = " << RECOELE_mva[index] << std::endl;

      // Brem & Classifaction
      RECOELE_fbrem[index]         = cand->fbrem() ;
      RECOELE_nbrems[index]        = cand->numberOfBrems();
      RECOELE_Class[index]         = cand->classification() ;
      if (useRECOformat) RECOELE_fbrem_mean[index]=1. - cand->gsfTrack()->outerMomentum().R()/cand->gsfTrack()->innerMomentum().R();
      RECOELE_fbrem_mode[index]=cand->fbrem();

      //@//
      /* std::cout << "--brem & classification: fbrem/nbrems/Class/fbrem_mean/fbrem_mode = "
		<< "  fbrem="      << RECOELE_fbrem[index]
		<< "  nbrems="     << RECOELE_nbrems[index]
        	<< "  class="      << RECOELE_Class[index]
		<< "  fbrem_mean=" << RECOELE_fbrem_mean[index]
		<< "  fbrem_mode=" << RECOELE_fbrem_mode[index]
		<< std::endl;*/

      // Corrections
      /*        RECOELE_isEcalEnergyCorrected[index] = cand->isEcalEnergyCorrected(); */
              RECOELE_ecalEnergy[index]            = cand->ecalEnergy(); 
      /*        RECOELE_ecalEnergyError[index]       = cand->ecalEnergyError(); */
      /*        RECOELE_isMomentumCorrected[index]   = cand->isMomentumCorrected(); */
      /*        RECOELE_trackMomentumError[index]    = cand->trackMomentumError(); */
      /*        RECOELE_electronMomentumError[index] = cand->electronMomentumError(); */
      
	      cout<<"#########Electron step 6 "<<endl;
    
      // Seed Collection
      if (useRECOformat) {
	edm::RefToBase<TrajectorySeed> seed = cand->gsfTrack()->extra()->seedRef();
	reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
	ele_seedSubdet2[index] = int(MyS->subDet2());
	if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2[index] = double(MyS->dPhi2());
	if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2[index]  = double(MyS->dRz2());
	
	ele_seedSubdet1[index] = int(MyS->subDet1());
	if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1[index] = double(MyS->dPhi1());
	if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1[index]  = double(MyS->dRz1());
      }
      //

      cout<<"#########Electron step 6a "<<endl;


      //add Gen Level info for electrons 

      ////////

      // Matching
      if (fillMCTruth==true){
      	//	int i=0;
      	//	for ( reco::CandidateCollection::const_iterator hIter=CollEle->begin(); hIter!= CollEle->end(); ++hIter ){
      	  //cout << "Reco Electron with pT= " << hIter->pt() << " " << RECOELE_PT[index] << " " << fabs(hIter->pt()-RECOELE_PT[index]) << " and mass="<< hIter->mass()<< endl;
      	// if (fabs(hIter->pt()-RECOELE_PT[index])<0.01){
      	//  i=hIter-(CollEle->begin());
      	// CandidateRef Ref( CollEle, i );
      	    edm::Ref<std::vector<reco::GenParticle> > genrefEle = (*GenParticlesMatchEle)[eletrackref];
      	    if (!genrefEle.isNull()){
      	      cout << "GenElectron with pT= " << genrefEle->p4().pt() << " and mass="<< genrefEle->p4().mass()<< endl;
      	      RECOELE_MatchingMCTruth[index]= true;
      	      RECOELE_MatchingMCpT[index]= genrefEle->p4().pt();
      	      RECOELE_MatchingMCEta[index]= genrefEle->p4().eta();
      	      RECOELE_MatchingMCPhi[index]= genrefEle->p4().phi();
      	    }
      	    else {
      	      cout << "There is no a reference to a genElectron" << endl;
      	    }
      	    // }
      	  //	}
      }//end match

      cout<<"RECOELE_MatchingMCTruth[index] = "<<RECOELE_MatchingMCTruth[index]<<endl;
      cout<<" RECOELE_MatchingMCpT[index] = "<< RECOELE_MatchingMCpT[index]<<endl;
      cout<<"RECOELE_MatchingMCEta[index]= "<<RECOELE_MatchingMCEta[index]<<endl;
      cout<<" RECOELE_MatchingMCPhi[index]= "<< RECOELE_MatchingMCPhi[index]<<endl;
      ///////

      index ++;
    }  
  } 

  void fillMuons(const edm::Event& iEvent,const edm::EventSetup& iSetup){

    //@//  cout<<"filling muons"<<endl;

    // Get the B-field
    //edm::ESHandle<MagneticField> magfield_;
    iSetup.get<IdealMagneticFieldRecord>().get( magfield_ );        


    // Muons
    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);

    edm::Handle<edm::ValueMap<float> > corrpterrormumap;
    iEvent.getByToken(muonCorrPtErrorMapTag_,corrpterrormumap);

    ////////////////////////////////////////////
    //Kalman Muon Calibrator for 2016 data

    /*   //@ Miquias
    edm::Handle<edm::View<pat::Muon> > SlimmedMuons;
    iEvent.getByToken(slimmedMuonsTag_, SlimmedMuons);
    vector<double> vcorrPt, vcorrPtError;
    double corrPt=0.,corrPtError=0.;
    double smearedPt=0., smearedPtError=0.;
    TLorentzVector p4;
    for(edm::View<pat::Muon>::const_iterator SlimmedMuonsIter = SlimmedMuons->begin(); SlimmedMuonsIter != SlimmedMuons->end(); ++SlimmedMuonsIter){
      smearedPt=SlimmedMuonsIter->pt();
      smearedPtError=SlimmedMuonsIter->muonBestTrack()->ptError();
      if (SlimmedMuonsIter->muonBestTrackType() == 1 && SlimmedMuonsIter->pt()<=200.){
    	if (isData){
    	  if (SlimmedMuonsIter->pt()>2.0 && fabs(SlimmedMuonsIter->eta())<2.4){
    	    KalmanMuonCalibrator calibrator("DATA_80X_13TeV");
    	    corrPt = calibrator.getCorrectedPt(SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(), SlimmedMuonsIter->charge());
    	    corrPtError = corrPt * calibrator.getCorrectedError(corrPt, SlimmedMuonsIter->eta(), SlimmedMuonsIter->bestTrack()->ptError()/corrPt );
    	    smearedPt=corrPt; // no smearing on data
    	    smearedPtError=corrPtError; // no smearing on data
    	  }
    	}
    	else { // isMC - calibration from data + smearing
    	  KalmanMuonCalibrator calibrator("MC_80X_13TeV");
    	  corrPt = calibrator.getCorrectedPt(SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(), SlimmedMuonsIter->charge());
    	  corrPtError = corrPt * calibrator.getCorrectedError(corrPt, SlimmedMuonsIter->eta(), SlimmedMuonsIter->bestTrack()->ptError()/corrPt );
    	  // smearedPt = calibrator.smearForSync(corrPt, SlimmedMuonsIter->eta()); // for synchronization
    	  smearedPt = calibrator.smear(corrPt, SlimmedMuonsIter->eta());
    	  smearedPtError = corrPtError;
    	  //smearedPt * calibrator.getCorrectedErrorAfterSmearing(smearedPt, SlimmedMuonsIter->eta(), corrPtError /smearedPt );
    	}
      }
      vcorrPt.push_back(smearedPt);
      vcorrPtError.push_back(smearedPtError);
    }
    */  
  /////////////////////////////////////////////////////////////////
   
     //Reham To get the error in muon PT

    
      // To add matching informations Reham
    // Matching

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchMu;
    iEvent.getByToken(goodMuonMCMatch_, GenParticlesMatchMu);
    edm::Handle<reco::CandidateCollection > CollMu;
    iEvent.getByToken(myMuons_, CollMu);

    if (GenParticlesMatchMu.isValid()){
      //@//  cout << endl<< "#Muons:"<<endl<<"#The reco collection to be matched has size= " <<  CollMu->size() << endl;
      //@//  cout << "#The matched map collection has size= " <<  GenParticlesMatchMu->size() << endl;
    }
    
    
    //using Rochester muon correction for 2017 data

    edm::Handle<edm::View<pat::Muon> > SlimmedMuons;
    iEvent.getByToken(slimmedMuonsTag_, SlimmedMuons);


    vector<double> vcorrPt, vcorrPtError,vuncorrpt, GenMuon_pt, GenMuon_eta, GenMuon_phi;
    vector<bool>Muon_match;
    double oldPt=0., oldPtError=0.;
    double  smearedPt=0., smearedPtError=0.;
    double correction=1, correction_error=0;
    TLorentzVector p4;
    int nl;

    cout<<"######## //ROOT TREE Muon correction// ########"<<endl;
   
    //calibrator->init(edm::FileInPath("RoccoR2017.txt").fullPath()); 

    edm::FileInPath corrPath("roccor_Run2_v2/data/RoccoR2017.txt");
    calibrator = std::unique_ptr<RoccoR>(new RoccoR(corrPath.fullPath()));

    cout<<"#ROOT TREE open the txt file for muon corrections"<<endl;
   
    int ii=0;

    for(edm::View<pat::Muon>::const_iterator SlimmedMuonsIter = SlimmedMuons->begin(); SlimmedMuonsIter != SlimmedMuons->end(); ++SlimmedMuonsIter){

      edm::Ref<edm::View<pat::Muon> > muref(SlimmedMuons,ii);

    double u = rgen_->Rndm();
    std::cout << "Random number: " << u << std::endl;

    bool Muon_Match = false;
    double Gen_Mu_pt= 0., Gen_Mu_eta= 0., Gen_Mu_phi= 0.;


      oldPt=SlimmedMuonsIter->pt();
      oldPtError=SlimmedMuonsIter->muonBestTrack()->ptError();
      
      cout<<"Slimmed Muon Pt = "<<oldPt<<"old pt error from muon best track= "<<oldPtError<<endl;
      ////
      
      /* // To add matching informations Reham */
      /* // Matching */
        if (fillMCTruth==true){ 
      /*  	int i=0; */
      /*  	for ( reco::CandidateCollection::const_iterator hIter=CollMu->begin(); hIter!= CollMu->end(); ++hIter ){ */
      /*  	  cout << "ROOT TREE Reco Muon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl; */
      /*  	  if (fabs(hIter->pt()- SlimmedMuonsIter->pt()) <0.01){ */
      /*  	    cout<<"Match found "<<endl; */
      /*  	      i=hIter-(CollMu->begin()); */
      /*  	      cout<<"i = "<<i<<endl; */
	  //        	      CandidateRef Ref(SlimmedMuons, ii ); 
        	      edm::Ref<std::vector<reco::GenParticle> > genrefMu = (*GenParticlesMatchMu)[muref]; 
        	      if (!genrefMu.isNull()){ 
        		cout<<"#ROOT TREE Matched found "<<endl; 
        		cout << "#ROOT TREE GenMuon with pT= " << genrefMu->p4().pt() << " and mass="<< genrefMu->p4().mass()<< endl; 
        		 Muon_Match = true; 
        		 Gen_Mu_pt = genrefMu->p4().pt(); 
			 Gen_Mu_eta = genrefMu->p4().eta();
			 Gen_Mu_phi = genrefMu->p4().phi();
			 Muon_match.push_back(Muon_Match);
			 GenMuon_pt.push_back(Gen_Mu_pt);
			 GenMuon_eta.push_back(Gen_Mu_eta);
			 GenMuon_phi.push_back(Gen_Mu_phi);
        	      } 
        	      else { 
        		cout << "ROOT TREE There is no reference to a genMuon" << endl; 
        		Muon_Match = false;
			Muon_match.push_back(Muon_Match);
			GenMuon_pt.push_back(Gen_Mu_pt); 
			GenMuon_eta.push_back(Gen_Mu_eta);
			GenMuon_phi.push_back(Gen_Mu_phi);
        	      } 
      /*  	  } */
      /*  	} //end of matching */
        } 
      ////
      
	if(Muon_Match == true){cout<<"Gen muon pt ="<<Gen_Mu_pt<<"Gen muon eta = "<< Gen_Mu_eta<<"Gen muon phi = "<< Gen_Mu_phi<<endl; } 
	

     cout<<"#ROOT TREE slimmed muon pt = "<<SlimmedMuonsIter->pt()<<" and best track type = "<<SlimmedMuonsIter->muonBestTrackType()<<endl;
     
      if (SlimmedMuonsIter->muonBestTrackType() == 1 && SlimmedMuonsIter->pt()<=200.){

	nl =  SlimmedMuonsIter->track()->hitPattern().trackerLayersWithMeasurement();

    	if (isData && nl>5 ){ //on data correction oly

	  cout<<"#ROOT TREE Data Muon correction "<<endl;

    	  if (SlimmedMuonsIter->pt()>2.0 && fabs(SlimmedMuonsIter->eta())<2.4){

	    //RoccoR Calibrator;
	    correction = calibrator->kScaleDT( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(), 0 ,0); 
	     correction_error = calibrator->kScaleDTerror( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi()); }//end eta > 2.4
	    else {// with pt<2 or not in the acceptance
	       correction =1;
	       correction_error =0;}//end else
	}//end data
   
    	else if(!isData && nl>5 ){ // isMC - calibration from data + smearing

	  cout<<"#ROOT TREE This is MC Muon correction"<<endl;

	  if ( Muon_Match == true){// when matched gen muon available
	    
	    cout<<"ROOT TREE MC and matched gen muon "<<endl;	    
	    cout<<"ROOT TREE gen muon pt = "<<Gen_Mu_pt<<endl;
	    
	    correction = calibrator->kSpreadMC( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(),  Gen_Mu_pt , 0 ,0 ); 
	    correction_error = calibrator->kSpreadMCerror( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(),  Gen_Mu_pt);
	    cout<<"TREE Corrrection = "<<correction<<endl;}
	  
	  else { //when matched gen muon not available
	    
	    cout<<"#rOOT TREE MC and No matched gen muon "<<endl;		 
	    
	    correction = calibrator->kSmearMC( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(), nl , u, 0 ,0 ); 
	    correction_error = calibrator->kSmearMCerror( SlimmedMuonsIter->charge(), SlimmedMuonsIter->pt(), SlimmedMuonsIter->eta(), SlimmedMuonsIter->phi(), nl ,u );
 cout<<"TREE Corrrection = "<<correction<<endl;}
	  
	}//ends if !isdata
      }//end if tracker , pt<200
   
    smearedPt = oldPt*correction;
    smearedPtError = correction_error;
    cout<<"ROOT TREE #Old Muon Pt = "<<oldPt<<" , correction = "<<correction<<" smeared muon pt = "<<smearedPt<<endl;
    cout<<"ROOT TREE #Old Muon Pt error = "<<oldPtError<<" , correction = "<<correction_error<<" smeared muon pt Error = "<<smearedPtError<<endl;
    vuncorrpt.push_back(oldPt); 
    vcorrPt.push_back(smearedPt);
    vcorrPtError.push_back(smearedPtError);
    
    ii++;

 }//end of for loop

  
    for (unsigned int k=0; k < Muon_match.size(); k++){
      cout<<"muon match vector = "<<Muon_match.at(k)<<"gen muon vector = "<< GenMuon_pt.at(k)<<"gen muon eta = "<<GenMuon_eta.at(k)<<"gen muon phi = "<<GenMuon_phi.at(k); }
    ///////////////////////////////////////////////////////////
 

   // primary vertex
 
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    //@// cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  //@//
	  /*  cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;*/
	  break;
	}
      }
    }
    
    //pVertex = math::XYZPoint(                                                                                                                                                    
    //			     PV.at(2).position().x(),                                                                                                                       
    //			     PV.at(2).position().y(),                                                                                                                               
    //			     PV.at(2).position().z());                            

    int indexbis=0;
    RECO_NMU=MuCandidates->size();


    for (edm::View<pat::Muon>::const_iterator cand = MuCandidates->begin(); cand != MuCandidates->end(); ++cand) {
            
      if(indexbis>99) break;

      edm::Ref<edm::View<pat::Muon> > muonref(MuCandidates,indexbis);

           /////Reham need to check this part 
      //Reham comment for now 
      
     //@ finds the pt error
      float calibratorPtError = 0, mu_uncorr_pt = 0, gen_muon_pt=0, gen_muon_eta=0, gen_muon_phi=0;
      bool MUMATCH=false;
      unsigned int NslimmedMuons = vcorrPt.size();
      cout << ">>>>> NslimmedMuons = " << NslimmedMuons <<" NMuCandidates = " << MuCandidates->size() << endl;
      for(unsigned int iref=0; iref<NslimmedMuons; ++iref){
	cout<<"Vcorrpt = "<<vcorrPt.at(iref)<<" cand pt = "<<cand->pt()<<endl;
	float candpt = cand->pt();
	float corrpt = vcorrPt.at(iref);
	if( candpt == corrpt ){
	//if( (vcorrPt.at(iref)) == (cand->pt()) ){
	  calibratorPtError = vcorrPtError.at(iref);
	  mu_uncorr_pt =vuncorrpt.at(iref);
	  cout<< "MuCandPt = " << cand->pt() << ", calibratorPt = " << vcorrPt.at(iref) << ", calibratorPtError = " << vcorrPtError.at(iref) << endl;
	  cout<<"mu uncorr pt = "<<mu_uncorr_pt<<endl;
	 
	  if (fillMCTruth==true){
	  if(Muon_match.at(iref)== true) {
	    MUMATCH = Muon_match.at(iref);
	    gen_muon_pt = GenMuon_pt.at(iref);
	    gen_muon_eta = GenMuon_eta.at(iref);
	    gen_muon_phi = GenMuon_phi.at(iref);
	    cout<<"MUMATCH = "<<MUMATCH<<"gen_muon = "<<gen_muon_pt<<"eta = "<<gen_muon_eta<<" phi = "<<gen_muon_phi<<endl;
	  }
	  }
	} 
      }
      
      RECOMU_MatchingMCTruth[indexbis]=MUMATCH;
      RECOMU_MatchingMCpT[indexbis]= gen_muon_pt;
      RECOMU_MatchingMCEta[indexbis]= gen_muon_eta;
      RECOMU_MatchingMCPhi[indexbis]= gen_muon_phi;
      
     cout<<" RECOMU_MatchingMCTruth[indexbis] = "<< RECOMU_MatchingMCTruth[indexbis]<<endl;
      cout<<" RECOMU_MatchingMCpT[indexbis] = "<< RECOMU_MatchingMCpT[indexbis]<<endl;
      cout<<" RECOMU_MatchingMCEta[indexbis] = "<< RECOMU_MatchingMCEta[indexbis]<<endl;
      cout<<" RECOMU_MatchingMCPhi[indexbis] = "<< RECOMU_MatchingMCPhi[indexbis]<<endl;
      
      edm::Ref<edm::View<pat::Muon> > mutrackref(MuCandidates,indexbis); 
      //@// edm::Ref<edm::View<pat::Muon> > mutrackrefv(VertMuCandidates,indexbis); 

    
      RECOMU_isPFMu[indexbis]=cand->isPFMuon();     
      RECOMU_isGlobalMu[indexbis]=cand->isGlobalMuon();
      RECOMU_isStandAloneMu[indexbis]=cand->isStandAloneMuon();
      RECOMU_isTrackerMu[indexbis]=cand->isTrackerMuon();
      RECOMU_isCaloMu[indexbis]=cand->isCaloMuon();

      if(cand->track().isNonnull()){//Reham
      // if(cand->innerTrack().isNonnull()){//Reham
      RECOMU_isTrackerHighPtMu[indexbis]=isTrackerHighPtMu(*cand,pVertex); //Reham error here 
       }//Reham
      //@//   
   /* std::cout << "\n Muon in the event: "
	        <<   "  isPF=" << RECOMU_isPFMu[indexbis]
		<<   "  isGB=" << RECOMU_isGlobalMu[indexbis]
		<< "   isSTA=" << RECOMU_isStandAloneMu[indexbis]
		<<   "  isTM=" << RECOMU_isTrackerMu[indexbis]
		<< "  isCalo=" << RECOMU_isCaloMu[indexbis]
		<< std::endl;*/

      // Kinematic of muon
      RECOMU_E[indexbis]=cand->p4().energy();
      RECOMU_PT[indexbis]=cand->p4().pt();
      RECOMU_P[indexbis]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOMU_ETA[indexbis]=cand->p4().eta();
      RECOMU_THETA[indexbis]=cand->p4().theta();
      RECOMU_PHI[indexbis]=cand->p4().phi();
      RECOMU_MASS[indexbis]=cand->p4().mass();
      RECOMU_CHARGE[indexbis]=cand->charge();
      RECOMU_PT_uncorr[indexbis]=mu_uncorr_pt;

      //@//    
       std::cout << "--kinematic:"
		<< "  pT="     << RECOMU_PT[indexbis]
	        << "  E="      << RECOMU_E[indexbis]
		<< "  p="      << RECOMU_P[indexbis]
		<< "  eta="    << RECOMU_ETA[indexbis]
		<< "  theta="  << RECOMU_THETA[indexbis]
		<< "  phi="    << RECOMU_PHI[indexbis]
		<< "  mass="   << RECOMU_MASS[indexbis]
		<< "  charge=" << RECOMU_CHARGE[indexbis]
		 <<"uncorr pt = "<< RECOMU_PT_uncorr[indexbis]
		<< std::endl;
	
      // Covariance matrix
      
      if(cand->innerTrack().isAvailable()){

	GlobalTrajectoryParameters gp(GlobalPoint(cand->innerTrack()->vx(), cand->innerTrack()->vy(),  cand->innerTrack()->vz()),
				      GlobalVector(cand->innerTrack()->px(),cand->innerTrack()->py(),cand->innerTrack()->pz()),
				      cand->innerTrack()->charge(),
				      magfield_.product());
	JacobianCurvilinearToCartesian curv2cart(gp);
	CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), cand->innerTrack()->covariance());
	const AlgebraicSymMatrix66 mat = cartErr.matrix();
	for (int i = 0; i < 3; ++i) { 
	  for (int j = 0; j < 3; ++j) { 
	    //bigCov(i,j) = mat[i+3][j+3]; 
	    //	  cout << "index= " << indexbis << "i= " << i << "j= " << j <<endl;
	    RECOMU_COV[indexbis][i][j]=mat[i+3][j+3]; 
	  }  
	}  
	
      }

      // Isolation
      //      RECOMU_TRACKISO[indexbis]=(*isoTkmumap)[mutrackref]/cand->p4().pt();
      RECOMU_TRACKISO_SUMPT[indexbis]=(cand->isolationR03().sumPt)/cand->p4().pt();
      //RECOMU_ECALISO[indexbis]=(*isoEcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_HCALISO[indexbis]=(*isoHcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_X[indexbis]=[(*isoTkmumap)[mutrackref]+(*isoEcalmumap)[mutrackref]+(*isoHcalmumap)[mutrackref]]/cand->p4().pt();
      // temporary solution for reducedrechit problem
      RECOMU_ECALISO[indexbis]=(cand->isolationR03().emEt)/cand->p4().pt();
      RECOMU_HCALISO[indexbis]=(cand->isolationR03().hadEt)/cand->p4().pt();
      RECOMU_X[indexbis]=RECOMU_TRACKISO[indexbis]+RECOMU_ECALISO[indexbis]+RECOMU_HCALISO[indexbis];

      RECOMU_PFchHad[indexbis]  = (cand->pfIsolationR03().sumChargedHadronPt);
      RECOMU_PFneuHad[indexbis] = (cand->pfIsolationR03().sumNeutralHadronEt);
      RECOMU_PFphoton[indexbis] = (cand->pfIsolationR03().sumPhotonEt);
      RECOMU_PFPUchAllPart[indexbis]= (cand->pfIsolationR03().sumPUPt);

//      RECOMU_PFchHad[indexbis]          =  (*isoPFChargedmumap)[mutrackref];
//      RECOMU_PFneuHad[indexbis]         =  (*isoPFNeutralmumap)[mutrackref];
//      RECOMU_PFphoton[indexbis]         =  (*isoPFGammamumap)[mutrackref];
//      RECOMU_PFPUchAllPart[indexbis]    =  (*isoPFPUmumap)[mutrackref];

      RECOMU_PFX_dB[indexbis]=(RECOMU_PFchHad[indexbis]+max(0.,RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-0.5*RECOMU_PFPUchAllPart[indexbis]))/cand->p4().pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.132;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.120;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.114;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.139;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.168;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.189;
      }
      else {       
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.674;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.565;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.442;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.515;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.821;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.660;
      }
       
      RECOMU_PFX_rho[indexbis]=(RECOMU_PFchHad[indexbis]+max( (RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-max(RHO_mu,0.0)*(EffectiveArea)),0.0) 
)/double(cand->p4().pt());

 

      //Reham
      RECOMU_SIP[indexbis]= std::abs(cand->dB(pat::Muon::PV3D))/cand->edB(pat::Muon::PV3D);
      RECOMU_IP[indexbis]= std::abs(cand->dB(pat::Muon::PV3D));
      RECOMU_IPERROR[indexbis]= cand->edB(pat::Muon::PV3D);

      cout<<"SIP = "<<RECOMU_SIP[indexbis]<<"IP , error = "<<RECOMU_IP[indexbis]<<", "<<RECOMU_IPERROR[indexbis]<<"std::abs(cand->dB(pat::Muon::PV3D))/cand->edB(pat::Muon::PV3D) = "<<std::abs(cand->dB(pat::Muon::PV3D))/cand->edB(pat::Muon::PV3D)<<"cand->dB(pat::Muon::PV3D) ="<<cand->dB(pat::Muon::PV3D)<<"cand->edB(pat::Muon::PV3D)= "<<cand->edB(pat::Muon::PV3D)<<endl;


      
      // Other properties
      RECOMU_numberOfMatches[indexbis]=cand->numberOfMatches();
      RECOMU_numberOfMatchedStations[indexbis]=cand->numberOfMatchedStations();
      RECOMU_caloCompatibility[indexbis]=cand->caloCompatibility();
      RECOMU_segmentCompatibility[indexbis]=(muon::segmentCompatibility( (*cand)));
      RECOMU_glbmuPromptTight[indexbis]=(muon::isGoodMuon( (*cand),muon::GlobalMuonPromptTight));

      //@//
      /*   std::cout	<< "--other properties:"
		<< "  n.matches="   << RECOMU_numberOfMatches[indexbis]
		<< "  caloComp="    << RECOMU_caloCompatibility[indexbis]
		<< "  segmentComp=" << RECOMU_segmentCompatibility[indexbis]
	        << "  glbmuPromptTight=" << RECOMU_glbmuPromptTight[indexbis]
		<< std::endl;*/

  
      // Track properties
      if(cand->muonBestTrack().isAvailable()){
	RECOMU_mubesttrkType[indexbis]=cand->muonBestTrackType();
	RECOMU_mubesttrkDxy[indexbis]=cand->muonBestTrack()->dxy(pVertex);
	RECOMU_mubesttrkDxyB[indexbis]=cand->muonBestTrack()->dxy(bs.position());
	RECOMU_mubesttrkDxyError[indexbis]=cand->muonBestTrack()->dxyError();
	RECOMU_mubesttrkDz[indexbis]=cand->muonBestTrack()->dz(pVertex);
	RECOMU_mubesttrkDzB[indexbis]=cand->muonBestTrack()->dz(bs.position());
	RECOMU_mubesttrkDzError[indexbis]=cand->muonBestTrack()->dzError();
	//RECOMU_mubesttrkPTError[indexbis]=(*corrpterrormumap)[mutrackref];;
      }

      //@ Reham comment for now

       RECOMU_mubesttrkPTError[indexbis]=cand->muonBestTrack()->ptError();
       RECOMU_Rochester_Error[indexbis]=calibratorPtError;
       cout<<">>>> check muon best track pt error "<< RECOMU_mubesttrkPTError[indexbis]<<" Muon rochester error ="<< RECOMU_Rochester_Error[indexbis]<<endl;

      if(cand->globalTrack().isAvailable()){
	RECOMU_mutrkPT[indexbis]=cand->globalTrack()->pt();
	RECOMU_mutrkPTError[indexbis]=cand->globalTrack()->ptError();
	//RECOMU_mutrkPTError[indexbis]=cand->bestTrack()->ptError(); // bestTrack
	RECOMU_mutrkDxy[indexbis]=cand->globalTrack()->dxy(pVertex);
	RECOMU_mutrkDxyError[indexbis]=cand->globalTrack()->dxyError();
	RECOMU_mutrkDxyB[indexbis]=cand->globalTrack()->dxy(bs.position()) ;
	RECOMU_mutrkDz[indexbis]=cand->globalTrack()->dz(pVertex);
	RECOMU_mutrkDzError[indexbis]=cand->globalTrack()->dzError();
	RECOMU_mutrkDzB[indexbis]=cand->globalTrack()->dz(bs.position());
	RECOMU_mutrkChi2PerNdof[indexbis]=cand->globalTrack()->normalizedChi2();
	RECOMU_mutrkCharge[indexbis]=cand->globalTrack()->charge();
 	RECOMU_mutrkNHits[indexbis]=cand->globalTrack()->numberOfValidHits(); 
	RECOMU_mutrkNPixHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_mutrkNStripHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidStripHits();
	RECOMU_mutrkNMuonHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidMuonHits(); 
	RECOMU_mutrktrackerLayersWithMeasurement[indexbis]=cand->globalTrack()->hitPattern().trackerLayersWithMeasurement(); 

	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack

	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();

      }
      else if(cand->innerTrack().isAvailable()){
	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack
	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();
      }

      if(cand->globalTrack().isAvailable() || cand->innerTrack().isAvailable() ){

	//@//
	/*	std::cout << "--muon track properties: "
	          << "  pt="        << RECOMU_mutrkPT[indexbis] 
	          << "  bestTrackType= " << RECOMU_mubesttrkType[indexbis]
		  << "  dxy="       << RECOMU_mubesttrkDxy[indexbis]
		  << "  dxyError="  << RECOMU_mubesttrkDxyError[indexbis]
		  << "  dxyB="      << RECOMU_mubesttrkDxyB[indexbis] 
		  << "  dz="        << RECOMU_mubesttrkDz[indexbis] 
		  << "  dzError="   << RECOMU_mubesttrkDzError[indexbis]
		  << "  dzB="       << RECOMU_mubesttrkDzB[indexbis]
		  << "  PtError="   << RECOMU_mubesttrkPTError[indexbis]
		  << "  chi2_nodf=" << RECOMU_mutrkChi2PerNdof[indexbis]
		  << "  charge="    << RECOMU_mutrkCharge[indexbis]
		  << "  nhits="     << RECOMU_mutrkNHits[indexbis] 
		  << "  nPixhits="   << RECOMU_mutrkNPixHits[indexbis]
	          << "  nStriphits=" << RECOMU_mutrkNStripHits[indexbis]
	          << "  nMuonhits="  << RECOMU_mutrkNMuonHits[indexbis]
		  << std::endl;*/
	
	
	// Tracker muon properties
	// RECOMU_trkmuArbitration[indexbis]=(muon::isGoodMuon( (*cand),muon::TrackerMuonArbitrated));
	RECOMU_trkmuArbitration[indexbis]=(muon::segmentCompatibility((*cand),pat::Muon::SegmentAndTrackArbitration));
	RECOMU_trkmu2DCompatibilityLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityLoose));
	RECOMU_trkmu2DCompatibilityTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityTight));
	RECOMU_trkmuOneStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationLoose));
	RECOMU_trkmuOneStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationTight));
	RECOMU_trkmuLastStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationLoose));
	RECOMU_trkmuLastStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationTight));
	RECOMU_trkmuOneStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngLoose));
	RECOMU_trkmuOneStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngTight));
	RECOMU_trkmuLastStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngLoose));
	RECOMU_trkmuLastStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngTight));
	RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtLoose));
	RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtTight));
	

	//@//
	/*	std::cout << "--tracker muon properties:"
		  << "  arbitration="              << RECOMU_trkmuArbitration[indexbis]
		  << "  2DCompLoose="              << RECOMU_trkmu2DCompatibilityLoose[indexbis]
		  << "  2DCompTight="              << RECOMU_trkmu2DCompatibilityTight[indexbis]
		  << "  1StationLoose="            << RECOMU_trkmuOneStationLoose[indexbis]
		  << "  1StationTight="            << RECOMU_trkmuOneStationTight[indexbis]
		  << "  LastStationLoose="         << RECOMU_trkmuLastStationLoose[indexbis]
		  << "  LastStationTight="         << RECOMU_trkmuLastStationTight[indexbis]
		  << "  1StationAngLoose="         << RECOMU_trkmuOneStationAngLoose[indexbis]
		  << "  1StationAngTight="         << RECOMU_trkmuOneStationAngTight[indexbis]
		  << "  LastStationAngLoose="      << RECOMU_trkmuLastStationAngLoose[indexbis]
		  << "  LastStationAngTight="      << RECOMU_trkmuLastStationAngTight[indexbis]
		  << "  LastStationOptLowptLoose=" << RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]
		  << "  LastStationOptLowptTight=" << RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]
		  << std::endl;*/
      }
  
         //    To add matching informations Reham */
      /* // Matching */
      /* if (fillMCTruth==true){  */
      /* 	int i=0;  */
      /* 	for ( reco::CandidateCollection::const_iterator hIter=CollMu->begin(); hIter!= CollMu->end(); ++hIter ){  */
      /* 	  cout << "@@@@@@@@@Reco Muon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl; */
      /* 	  cout<<"RECOMU_PT[indexbis] = "<<RECOMU_PT[indexbis]<<endl; */
      /* 	  if (fabs(hIter->pt()-RECOMU_PT[indexbis])<0.01){ */
      /* 	    i=hIter-(CollMu->begin()); */
      /* 	    CandidateRef Ref( CollMu, i ); */
      /* 	    edm::Ref<std::vector<reco::GenParticle> > genrefMu = (*GenParticlesMatchMu)[Ref];  */
      /* 	    if (!genrefMu.isNull()){  */
      /* 	      cout << "@@@@@@GenMuon with pT= " << genrefMu->p4().pt() << " and mass="<< genrefMu->p4().mass()<< endl;  */
      /* 	      RECOMU_MatchingMCTruth[i]= true; */
      /* 	      RECOMU_MatchingMCpT[i]= genrefMu->p4().pt(); */
      /* 	      RECOMU_MatchingMCEta[i]= genrefMu->p4().eta(); */
      /* 	      RECOMU_MatchingMCPhi[i]= genrefMu->p4().phi(); */
      /* 	      //cout<<"GenMuon pT check = "<< RECOMU_MatchingMCpT[indexbis]<<endl; */
      /* 	      //cout<<"GenMuon eta check = "<< RECOMU_MatchingMCEta[indexbis]<<endl; */
      /* 	      // cout<<"GenMuon phi check = "<< RECOMU_MatchingMCPhi[indexbis]<<endl; */
      /* 	    }  */
      /* 	    else{  */
      /* 	      cout << "There is no reference to a genMuon" << endl;  */
      /* 		}  */
      /* 	  } */
      /* 	} */
      /* } //end of matching  */
      
      indexbis++;
    }
  }
  


  void fillPhotons(const edm::Event& iEvent){
    // Photons
    //edm::Handle<edm::View<reco::Candidate> > photons;
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonsTag_, photons);
    RECO_NPHOT=photons->size();

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchPhot;
    iEvent.getByToken(goodGammaMCMatch_, GenParticlesMatchPhot);
    edm::Handle<reco::CandidateCollection > CollPhot;
    iEvent.getByToken(myGammas_, CollPhot);
    bool ismyGammas=false;
    if (CollPhot.isValid()) ismyGammas=true;
    
    int iphot=0;
    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {
      if (iphot>19) break;

      
      RECOPHOT_PT[iphot]=cand->pt();
      RECOPHOT_ETA[iphot]=cand->eta();
      RECOPHOT_PHI[iphot]=cand->phi();
      RECOPHOT_THETA[iphot]=cand->theta();

      /////

      // Matching
      int i=0;
      if (ismyGammas){
	for ( reco::CandidateCollection::const_iterator hIter=CollPhot->begin(); hIter!= CollPhot->end(); ++hIter ){
	  //cout << "Reco Photon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	  if (fabs(hIter->pt()-RECOPHOT_PT[iphot])<0.01){
	    i=hIter-(CollPhot->begin());
	    CandidateRef Ref( CollPhot, i );
	    edm::Ref<std::vector<reco::GenParticle> > genrefPhot = (*GenParticlesMatchPhot)[Ref];
	    if (!genrefPhot.isNull()){
	      cout << "GenPhot with pT= " << genrefPhot->p4().pt() << " and mass="<< genrefPhot->p4().mass()<< endl;
	      RECOPHOT_MatchingMCTruth[i]= true;
	      RECOPHOT_MatchingMCpT[i]= genrefPhot->p4().pt();
	      RECOPHOT_MatchingMCEta[i]= genrefPhot->p4().eta();
	      RECOPHOT_MatchingMCPhi[i]= genrefPhot->p4().phi();	    
	    } 
	  }   
	}
      }

      //////

      iphot++;

       cout<<"photon pt = "<<RECOPHOT_PT[iphot]<<"eta= "<<RECOPHOT_ETA[iphot]<<endl;

       cout<<" RECOPHOT_MatchingMCTruth[i] = "<< RECOPHOT_MatchingMCTruth[i]<<endl;
       cout<<"RECOPHOT_MatchingMCpT[i] = "<<RECOPHOT_MatchingMCpT[i]<<endl;
       cout<<"RECOPHOT_MatchingMCEta[i] = "<<RECOPHOT_MatchingMCEta[i]<<endl;
       cout<<"RECOPHOT_MatchingMCPhi[i] = "<<RECOPHOT_MatchingMCPhi[i]<<endl;
      }
    // PF ISR photon

    iphot=0;
    RECO_NPFPHOT=photons->size();

    //@// cout << "There are " << photons->size() << " photons" << endl;

    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {

      if (iphot>19) break;

      edm::Ref<edm::View<pat::Photon> > phtrackref(photons,iphot); 

       RECOPFPHOT_PT_uncorr[iphot]=cand->pt();

      cout<<"PF photon uncorr pt = "<< RECOPFPHOT_PT_uncorr[iphot]<<"eta= "<<cand->eta()<<endl;

      //Reham correction of photon energy 

      auto corrP4  = cand->p4() * cand->userFloat("ecalEnergyPostCorr")/cand->energy();
     
      RECOPFPHOT_PT[iphot]= corrP4.pt();

      // cout<<"corrP4"<<corrP4<<endl;
      cout<<"Corrected PF photon pt = "<<corrP4.pt()<<" eta = "<<corrP4.eta()<<" phi = "<<corrP4.phi()<<" energy E = "<<corrP4.E()<<" energy Energy= "<<corrP4.energy()<<" energy Et= "<<corrP4.Et()<<endl;


      // error on pt of photon taken from #include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
      double C;
      double S;
      double N;
      if(TMath::Abs(cand->eta())<1.48){C=0.35/100; S=5.51/100; N=98./1000.;}
      else{C=0; S=12.8/100; N=440./1000.;} 
      double  perr = TMath::Sqrt(C*C*cand->p4().e()*cand->p4().e() + S*S*cand->p4().e() + N*N);
      double pterr = perr*cand->pt()/cand->p(); 
      RECOPFPHOT_PTError[iphot]=float(pterr);
      
        cout << "PF Photon : pT= " << RECOPFPHOT_PT[iphot] << " pTerr= " << RECOPFPHOT_PTError[iphot]<< endl;

      RECOPFPHOT_ETA[iphot]=cand->eta();
      RECOPFPHOT_PHI[iphot]=cand->phi();
      RECOPFPHOT_THETA[iphot]=cand->theta();
/*      
      RECOPFPHOT_PFchAllPart[iphot]      =  (*isoPFChargedAllphmap)[phtrackref]; 
      RECOPFPHOT_PFchHad[iphot]          =  (*isoPFChargedphmap)[phtrackref]; 
      RECOPFPHOT_PFneuHad[iphot]         =  (*isoPFNeutralphmap)[phtrackref]; 
      RECOPFPHOT_PFphoton[iphot]         =  (*isoPFGammaphmap)[phtrackref]; 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (*isoPFPUphmap)[phtrackref]; 
*/
//TODO      RECOPFPHOT_PFchAllPart[iphot]      =  (cand->PflowIsolationVariables().sumChargedParticlePt); 
      RECOPFPHOT_PFchHad[iphot]          =  (cand->chargedHadronIso());
      RECOPFPHOT_PFneuHad[iphot]         =  (cand->neutralHadronIso()); 
      RECOPFPHOT_PFphoton[iphot]         =  (cand->photonIso()); 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (cand->puChargedHadronIso());
    
      RECOPFPHOT_PFX_rho[iphot]=( 
 				 RECOPFPHOT_PFchHad[iphot]+ 
 				 RECOPFPHOT_PFneuHad[iphot]+ 
 				 RECOPFPHOT_PFphoton[iphot]+ 
 				 RECOPFPHOT_PFPUchAllPart[iphot] 
				  )/double(cand->p4().pt()); 

     	//Reham error in PT and some systematic variables 
      
      cout<<"ecalEnergyErrPostCorr = "<<cand->userFloat("ecalEnergyErrPostCorr")<<endl;
      cout<<"ecalEnergyPreCorr = "<<cand->userFloat("ecalEnergyPreCorr")<<endl;
      cout<<"ecalEnergyErrPreCorr = "<<cand->userFloat("ecalEnergyErrPreCorr")<<endl;

      RECOPFPHOT_ecalEnergyPreCorr[iphot]  = cand->userFloat("ecalEnergyPreCorr");
      RECOPFPHOT_ecalEnergyErrPreCorr[iphot]  = cand->userFloat("ecalEnergyErrPreCorr");
      RECOPFPHOT_ecalEnergyErrPostCorr[iphot]  = cand->userFloat("ecalEnergyErrPostCorr");
      RECOPFPHOT_energyScaleValue[iphot]          = cand->userFloat("energyScaleValue");
      RECOPFPHOT_energySigmaValue[iphot]          = cand->userFloat("energySigmaValue");
      RECOPFPHOT_energyScaleUp[iphot]          = cand->userFloat("energyScaleUp");
      RECOPFPHOT_energyScaleDown[iphot]          = cand->userFloat("energyScaleDown");
      RECOPFPHOT_energyScaleStatUp[iphot]          = cand->userFloat("energyScaleStatUp");
      RECOPFPHOT_energyScaleStatDown[iphot]          = cand->userFloat("energyScaleStatDown");
      RECOPFPHOT_energyScaleSystUp[iphot]          = cand->userFloat("energyScaleSystUp");
      RECOPFPHOT_energyScaleSystDown[iphot]          = cand->userFloat("energyScaleSystDown");
      RECOPFPHOT_energyScaleGainUp[iphot]          = cand->userFloat("energyScaleGainUp");
      RECOPFPHOT_energyScaleGainDown[iphot]          = cand->userFloat("energyScaleGainDown");
      RECOPFPHOT_energyScaleEtUp[iphot]          = cand->userFloat("energyScaleEtUp");
      RECOPFPHOT_energyScaleEtDown[iphot]          = cand->userFloat("energyScaleEtDown");
      RECOPFPHOT_energySigmaUp[iphot]               = cand->userFloat("energySigmaUp");
      RECOPFPHOT_energySigmaDown[iphot]          = cand->userFloat("energySigmaDown");
      RECOPFPHOT_energySigmaPhiUp[iphot]          = cand->userFloat("energySigmaPhiUp");
      RECOPFPHOT_energySigmaPhiDown[iphot]          = cand->userFloat("energySigmaPhiDown");
      RECOPFPHOT_energySigmaRhoUp[iphot]          = cand->userFloat("energySigmaRhoUp");
      RECOPFPHOT_energySigmaRhoDown[iphot]          = cand->userFloat("energySigmaRhoDown");

	/* std::cout << "Photon syst ==" */
	/* 	  <<"\n RECOPFPHOT_ecalEnergyPreCorr[iphot]= "<< RECOPFPHOT_ecalEnergyPreCorr[iphot] */
	/* 	  <<"\n RECOPFPHOT_ecalEnergyErrPreCorr[iphot]= "<< RECOPFPHOT_ecalEnergyErrPreCorr[iphot] */
	/* 	  <<"\n RECOPFPHOT_ecalEnergyErrPostCorr[iphot]= "<< RECOPFPHOT_ecalEnergyErrPostCorr[iphot] */
	/* 	  <<"\n RECOPFPHOT_energyScaleValue[iphot]= "<<RECOPFPHOT_energyScaleValue[iphot] */
	/* 	  <<"\n RECOPFPHOT_energySigmaValue[iphot]= "<<RECOPFPHOT_energySigmaValue[iphot] */
	/* 	  <<"\n RECOPFPHOT_energyScaleUp[iphot]= "<<RECOPFPHOT_energyScaleUp[iphot] */
	/* 	  <<"\n RECOPFPHOT_energyScaleDown[iphot]= "<<RECOPFPHOT_energyScaleDown[iphot] */
	/* 	  <<"\n RECOPFPHOT_energyScaleStatUp[iphot]= "<< RECOPFPHOT_energyScaleStatUp[iphot] */
	/* 	  <<"\n RECOPFPHOT_energyScaleStatDown[iphot]= "<< RECOPFPHOT_energyScaleStatDown[iphot]      */
	/* 	  <<"\n RECOPFPHOT_energyScaleSystUp[iphot]= "<<RECOPFPHOT_energyScaleSystUp[iphot]         */
	/* 	  <<"\n RECOPFPHOT_energyScaleSystDown[iphot]= "<<RECOPFPHOT_energyScaleSystDown[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energyScaleGainUp[iphot]= "<<RECOPFPHOT_energyScaleGainUp[iphot]          */
	/* 	  <<"\n RECOPFPHOT_energyScaleGainDown[iphot]= "<<RECOPFPHOT_energyScaleGainDown[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energyScaleEtUp[iphot]= "<< RECOPFPHOT_energyScaleEtUp[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energyScaleEtDown[iphot]= "<< RECOPFPHOT_energyScaleDown[iphot] */
	/* 	  <<"\n RECOPFPHOT_energySigmaUp[iphot]= "<< RECOPFPHOT_energySigmaUp[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energySigmaDown[iphot]= "<< RECOPFPHOT_energySigmaDown[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energySigmaPhiUp[iphot]= "<<RECOPFPHOT_energySigmaPhiUp[iphot]          */
	/* 	  <<"\n RECOPFPHOT_energySigmaPhiDown[iphot]= "<<RECOPFPHOT_energySigmaPhiDown[iphot]         */
	/* 	  <<"\n RECOPFPHOT_energySigmaRhoUp[iphot]= "<<RECOPFPHOT_energySigmaRhoUp[iphot]        */
	/* 	  <<"\n RECOPFPHOT_energySigmaRhoDown[iphot]= "<<RECOPFPHOT_energySigmaRhoDown[iphot]        */
	/* 	  <<std::endl;  */
     
      iphot++;
      
    }
    
  }

void fillTracks(const edm::Event& iEvent){
    // Tracks
  using namespace edm; using namespace std; using namespace reco;
  edm::Handle<edm::View<pat::PackedCandidate>> cands;
  iEvent.getByToken(pfTag_,cands);
 
  int countk=0;
  RECO_NTRACK=0;
  
//  for(edm::View<pat::PackedCandidate>::const_iterator i=cands->begin(); i!=cands->end(); i++){
  for(unsigned int i=0;i<cands->size();i++){
     if (countk>199) break;
     const pat::PackedCandidate & c = (*cands)[i];
     if(!(c.charge() != 0 && c.numberOfHits()> 0)) continue;
     RECO_NTRACK++;
    RECO_TRACK_PT[countk]=c.pseudoTrack().pt();
    RECO_TRACK_ETA[countk]=c.pseudoTrack().eta();
    RECO_TRACK_PHI[countk]=c.pseudoTrack().phi();
    RECO_TRACK_CHI2[countk]=c.pseudoTrack().chi2();
    RECO_TRACK_CHI2RED[countk]=c.pseudoTrack().normalizedChi2();
    //RECO_TRACK_CHI2PROB=TMath::Prob(i->chi2(),i->ndof());
    RECO_TRACK_CHI2PROB[countk]=ChiSquaredProbability(c.pseudoTrack().chi2(),c.pseudoTrack().ndof());
    RECO_TRACK_NHITS[countk]=c.pseudoTrack().numberOfValidHits();
    RECO_TRACK_DXY[countk]=c.pseudoTrack().dxy();
    RECO_TRACK_DXYERR[countk]=c.pseudoTrack().dxyError();
    RECO_TRACK_DZ[countk]=c.pseudoTrack().dz();
    RECO_TRACK_DZERR[countk]=c.pseudoTrack().dzError();
    countk++;
  }
  //@// cout << "Number of Tracks in the event= " << RECO_NTRACK << endl;
}
        


  void fillMET(const edm::Event& iEvent){

    //Type I Muon Correction MET
    //cormetmuons = SetCaloMET(iEvent,cormetMuTag_);
    
    //PFMET
    //@// std::cout<<"met uncer"<< pat::MET::METUncertainty(12) <<endl;
    edm::Handle<pat::METCollection> pfmetHandle;
    iEvent.getByToken(pfmetTag_,pfmetHandle);
    for ( pat::METCollection::const_iterator i=pfmetHandle->begin(); i!=pfmetHandle->end(); i++) {
      
      //Read the decision of MET Filter Reham

      //new filter to be use (under test)

      /* edm::Handle <bool> passecalBadCalibFilterUpdate ; */
      /* iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate); */
      /* PassecalBadCalibFilterUpdated =  (*passecalBadCalibFilterUpdate ); */
      /* cout<<" PassecalBadCalibFilterUpdated = "<< PassecalBadCalibFilterUpdated<<endl; */ 

      //old filters

      /* edm::Handle<bool> ifilterbadChCand; */
      /* iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand); */
      /* filterbadChCandidate = (*ifilterbadChCand);  */

      /* cout<<"filterbadChCandidate= "<<filterbadChCandidate<<endl; */

      /* edm::Handle<bool> ifilterbadPFMuon; */
      /* iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon); */
      /* filterbadPFMuon = (*ifilterbadPFMuon); */

      /* cout<<"filterbadPFMuon = "<<filterbadPFMuon<<endl; */
      
      cormetmuons = i->pt();
      
      //Type1 correction Reham
      pfmet       = i->corPt();     
      pfmet_x     = i->corPx();
      pfmet_y     = i->corPy();
      pfmet_phi   = i->corPhi();
      pfmet_theta = i->corP3().theta();
      
      //uncorrected met Reham
      pfmet_uncorr       = i->uncorPt();     
      pfmet_x_uncorr     = i->uncorPx();
      pfmet_y_uncorr     = i->uncorPy();
      pfmet_phi_uncorr   = i->uncorPhi();
      pfmet_theta_uncorr = i->uncorP3().theta();


      
      //std::cout<<"met uncer"<< i->MET::METUncertainty(0) <<endl;

      //MET uncertinities 

      pfmet_JetEnUp         = i->shiftedPt(pat::MET::JetEnUp, pat::MET::Type1);
      pfmet_JetEnDn         = i->shiftedPt(pat::MET::JetEnDown, pat::MET::Type1);
      pfmet_ElectronEnUp    = i->shiftedPt(pat::MET::ElectronEnUp, pat::MET::Type1);
      pfmet_ElectronEnDn    = i->shiftedPt(pat::MET::ElectronEnDown, pat::MET::Type1);
      pfmet_MuonEnUp        = i->shiftedPt(pat::MET::MuonEnUp, pat::MET::Type1);
      pfmet_MuonEnDn        = i->shiftedPt(pat::MET::MuonEnDown, pat::MET::Type1);
      pfmet_JetResUp        = i->shiftedPt(pat::MET::JetResUp, pat::MET::Type1);
      pfmet_JetResDn        = i->shiftedPt(pat::MET::JetResDown, pat::MET::Type1);
      pfmet_UnclusteredEnUp = i->shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Type1);
      pfmet_UnclusteredEnDn = i->shiftedPt(pat::MET::UnclusteredEnDown, pat::MET::Type1);
      pfmet_PhotonEnUp      = i->shiftedPt(pat::MET::PhotonEnUp, pat::MET::Type1);
      pfmet_PhotonEnDn      = i->shiftedPt(pat::MET::PhotonEnDown, pat::MET::Type1);
      pfmet_TauEnUp         = i->shiftedPt(pat::MET::TauEnUp , pat::MET::Type1);
      pfmet_TauEnDown       = i->shiftedPt(pat::MET::TauEnDown, pat::MET::Type1); 

    }
     
  }


  void fillMETFilters(const edm::Event& iEvent){
    
  //filters decision from trigger results

   edm::Handle<edm::TriggerResults> noiseFilterBits_; //our trigger result object
   iEvent.getByToken(noiseFilterTag_,noiseFilterBits_);
   const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);

   for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
     // cout <<"&&&&&&&"<< names.triggerName(i) << endl;
     

     if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_ ){
       cout<<"&&&&& pass GoodVtxNoiseFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterGoodVtxNoise = int(noiseFilterBits_->accept(i));};

     if (names.triggerName(i) == GlobalSuperTightHalo2016NoiseFilter_Selector_ ){
       cout<<"&&&&& pass GlobalSuperTightHalo2016NoiseFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterGlobalSuperTightHalo2016NoiseFilter = int(noiseFilterBits_->accept(i));};

     if (names.triggerName(i) == HBHENoiseFilter_Selector_){
       cout<<"&&&& Pass HBHENoiseFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterHBHENoise = int(noiseFilterBits_->accept(i));};
     
     if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_ ){
       cout<<"&&&&& pass HBHENoiseIsoFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterHBHENoiseIso = int(noiseFilterBits_->accept(i));};

     if (names.triggerName(i) == EcalDeadCellTriggerPrimitiveNoiseFilter_Selector_ ){
       cout<<"&&&&& pass EcalDeadCellTriggerPrimitiveNoiseFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterEcalDeadCellTriggerPrimitiveNoise = int(noiseFilterBits_->accept(i));};

    if (names.triggerName(i) == BadPFMuonFilter_Selector_ ){
       cout<<"&&&&& pass EcalBadPFMuonFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterBadPFMuon = int(noiseFilterBits_->accept(i));};

    if (names.triggerName(i) == BadChargedCandidateFilter_Selector_ ){
       cout<<"&&&&& pass BadChargedCandidateFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterBadChargedCandidate = int(noiseFilterBits_->accept(i));};

     if (names.triggerName(i) == EEBadScNoiseFilter_Selector_ ){
       cout<<"&&&&& pass EEBadScNoiseFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterEEBadScNoise = int(noiseFilterBits_->accept(i));};

    if (names.triggerName(i) == EcalBadCalibFilter_Selector_ ){
       cout<<"&&&&& pass EcalBadCalibFilter = "<< noiseFilterBits_->accept(i)<<endl;
       passFilterEcalBadCalib = int(noiseFilterBits_->accept(i));};
      
   };
   // cout<<"passFilterHBHE ="<< passFilterHBHE <<", "
   //  <<"passFilterEEBadSC ="<< passFilterEEBadSC 
   //  << endl;
  }

  void fillGD2e2mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD;
    iEvent.getByLabel(ftsigma_Vert, GeomD);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD->begin(); cand!=GeomD->end(); ++cand){
      ftsigma[jjj]=(*cand);
      jjj++;
    }

   
    edm::Handle<vector<double> > GeomDlag;
    iEvent.getByLabel(ftsigmalag_Vert, GeomDlag);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag->begin(); cand!=GeomDlag->end(); ++cand){
      ftsigmalag[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdX_;
    iEvent.getByLabel(gdX_Vert, gdX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdX_->begin(); cand!=gdX_->end(); ++cand){
      gdX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagX_;
    iEvent.getByLabel(gdlagX_Vert, gdlagX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagX_->begin(); cand!=gdlagX_->end(); ++cand){
      gdlagX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdY_;
    iEvent.getByLabel(gdY_Vert, gdY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdY_->begin(); cand!=gdY_->end(); ++cand){
      gdY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagY_;
    iEvent.getByLabel(gdlagY_Vert, gdlagY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagY_->begin(); cand!=gdlagY_->end(); ++cand){
      gdlagY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZ_;
    iEvent.getByLabel(gdZ_Vert, gdZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZ_->begin(); cand!=gdZ_->end(); ++cand){
      gdZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZ_;
    iEvent.getByLabel(gdlagZ_Vert, gdlagZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZ_->begin(); cand!=gdlagZ_->end(); ++cand){
      gdlagZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProb_;
    iEvent.getByLabel(gdlagProb_Vert, gdlagProb_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProb_->begin(); cand!=gdlagProb_->end(); ++cand){
      gdlagProb[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdof_;
    iEvent.getByLabel(gdlagNdof_Vert, gdlagNdof_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdof_->begin(); cand!=gdlagNdof_->end(); ++cand){
      gdlagNdof[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_MMMM;
    iEvent.getByLabel(ftsigma_VertMMMM, GeomD_MMMM);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_MMMM->begin(); cand!=GeomD_MMMM->end(); ++cand){
      ftsigmaMMMM[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_MMMM;
    iEvent.getByLabel(ftsigmalag_VertMMMM, GeomDlag_MMMM);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_MMMM->begin(); cand!=GeomDlag_MMMM->end(); ++cand){
      ftsigmalagMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdXMMMM_;
    iEvent.getByLabel(gdX_VertMMMM, gdXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXMMMM_->begin(); cand!=gdXMMMM_->end(); ++cand){
      gdXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXMMMM_;
    iEvent.getByLabel(gdlagX_VertMMMM, gdlagXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXMMMM_->begin(); cand!=gdlagXMMMM_->end(); ++cand){
      gdlagXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYMMMM_;
    iEvent.getByLabel(gdY_VertMMMM, gdYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYMMMM_->begin(); cand!=gdYMMMM_->end(); ++cand){
      gdYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYMMMM_;
    iEvent.getByLabel(gdlagY_VertMMMM, gdlagYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYMMMM_->begin(); cand!=gdlagYMMMM_->end(); ++cand){
      gdlagYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZMMMM_;
    iEvent.getByLabel(gdZ_VertMMMM, gdZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZMMMM_->begin(); cand!=gdZMMMM_->end(); ++cand){
      gdZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZMMMM_;
    iEvent.getByLabel(gdlagZ_VertMMMM, gdlagZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZMMMM_->begin(); cand!=gdlagZMMMM_->end(); ++cand){
      gdlagZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbMMMM_;
    iEvent.getByLabel(gdlagProb_VertMMMM, gdlagProbMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbMMMM_->begin(); cand!=gdlagProbMMMM_->end(); ++cand){
      gdlagProbMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofMMMM_;
    iEvent.getByLabel(gdlagNdof_VertMMMM, gdlagNdofMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofMMMM_->begin(); cand!=gdlagNdofMMMM_->end(); ++cand){
      gdlagNdofMMMM[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4e(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_EEEE;
    iEvent.getByLabel(ftsigma_VertEEEE, GeomD_EEEE);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_EEEE->begin(); cand!=GeomD_EEEE->end(); ++cand){
      ftsigmaEEEE[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_EEEE;
    iEvent.getByLabel(ftsigmalag_VertEEEE, GeomDlag_EEEE);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_EEEE->begin(); cand!=GeomDlag_EEEE->end(); ++cand){
      ftsigmalagEEEE[jjj]=(*cand);
      jjj++;
    }


    edm::Handle<vector<double> > gdXEEEE_;
    iEvent.getByLabel(gdX_VertEEEE, gdXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXEEEE_->begin(); cand!=gdXEEEE_->end(); ++cand){
      gdXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXEEEE_;
    iEvent.getByLabel(gdlagX_VertEEEE, gdlagXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXEEEE_->begin(); cand!=gdlagXEEEE_->end(); ++cand){
      gdlagXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYEEEE_;
    iEvent.getByLabel(gdY_VertEEEE, gdYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYEEEE_->begin(); cand!=gdYEEEE_->end(); ++cand){
      gdYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYEEEE_;
    iEvent.getByLabel(gdlagY_VertEEEE, gdlagYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYEEEE_->begin(); cand!=gdlagYEEEE_->end(); ++cand){
      gdlagYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZEEEE_;
    iEvent.getByLabel(gdZ_VertEEEE, gdZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZEEEE_->begin(); cand!=gdZEEEE_->end(); ++cand){
      gdZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZEEEE_;
    iEvent.getByLabel(gdlagZ_VertEEEE, gdlagZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZEEEE_->begin(); cand!=gdlagZEEEE_->end(); ++cand){
      gdlagZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbEEEE_;
    iEvent.getByLabel(gdlagProb_VertEEEE, gdlagProbEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbEEEE_->begin(); cand!=gdlagProbEEEE_->end(); ++cand){
      gdlagProbEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofEEEE_;
    iEvent.getByLabel(gdlagNdof_VertEEEE, gdlagNdofEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofEEEE_->begin(); cand!=gdlagNdofEEEE_->end(); ++cand){
      gdlagNdofEEEE[jjj]=(*cand);
      jjj++;
    }

  }

  void fillConstraintVtx2e2mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertex, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexX[jjj]=cand->position().x();
      StdFitVertexY[jjj]=cand->position().y();
      StdFitVertexZ[jjj]=cand->position().z();
      StdFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      //@//   
      /*  cout << "Std Fit: " 
	   <<  StdFitVertexX[jjj] << " " 
	   <<  StdFitVertexY[jjj] << " " 
	   <<  StdFitVertexZ[jjj] << " " 
	   <<  StdFitVertexChi2r[jjj] << " " 
	   <<  StdFitVertexProb[jjj] << endl;*/

      
      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){	 
	  if (i > 3) break;
	  //@//  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrack_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrack_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrack_PHI[i][jjj]=refit_tks[i].phi();
	}
      }
                   
      jjj++;
    }

    

    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByToken(RECOcollNameEEMM, CandidatesEEMM);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMass, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesEEMM,kk);
      //@//  cout << "Original 2e2mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_EEMM_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }


    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertex, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexX[jjj]=cand->position().x();
      KinFitVertexY[jjj]=cand->position().y();
      KinFitVertexZ[jjj]=cand->position().z();
      KinFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());

      //@//  
    /*  cout << "Kin Fit: " 
	   <<  KinFitVertexX[jjj] << " " 
	   <<  KinFitVertexY[jjj] << " " 
	   <<  KinFitVertexZ[jjj] << " " 
	   <<  KinFitVertexChi2r[jjj] << " " 
	   <<  KinFitVertexProb[jjj] << endl;*/
      jjj++;
    }
  }

  void fillConstraintVtx4mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexMMMM, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexXMMMM[jjj]=cand->position().x();
      StdFitVertexYMMMM[jjj]=cand->position().y();
      StdFitVertexZMMMM[jjj]=cand->position().z();
      StdFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      
      /* cout << "Std Fit MMMM: " 
	   <<  StdFitVertexXMMMM[jjj] << " " 
	   <<  StdFitVertexYMMMM[jjj] << " " 
	   <<  StdFitVertexZMMMM[jjj] << " " 
	   <<  StdFitVertexChi2rMMMM[jjj] << " " 
	   <<  StdFitVertexProbMMMM[jjj] << endl;*/

      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){
	  if (i > 3) break;
	  //@//  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrackMMMM_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrackMMMM_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrackMMMM_PHI[i][jjj]=refit_tks[i].phi();
	}
      }

      jjj++;

    }


    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByToken(RECOcollNameMMMM_, CandidatesMMMM);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMassMMMM, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesMMMM,kk);
      //@//  cout << "Original 4mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_MMMM_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }



    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexMMMM, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexXMMMM[jjj]=cand->position().x();
      KinFitVertexYMMMM[jjj]=cand->position().y();
      KinFitVertexZMMMM[jjj]=cand->position().z();
      KinFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
    

      /*cout << "Kin Fit MMMM: " 
	   <<  KinFitVertexXMMMM[jjj] << " " 
	   <<  KinFitVertexYMMMM[jjj] << " " 
	   <<  KinFitVertexZMMMM[jjj] << " " 
	   <<  KinFitVertexChi2rMMMM[jjj] << " " 
	   <<  KinFitVertexProbMMMM[jjj] << endl;*/
      jjj++;
    }
  }

  void fillConstraintVtx4e(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexEEEE, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexXEEEE[jjj]=cand->position().x();
      StdFitVertexYEEEE[jjj]=cand->position().y();
      StdFitVertexZEEEE[jjj]=cand->position().z();
      StdFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      /*  cout << "Std Fit EEEE: " 
	   <<  StdFitVertexXEEEE[jjj] << " " 
	   <<  StdFitVertexYEEEE[jjj] << " " 
	   <<  StdFitVertexZEEEE[jjj] << " " 
	   <<  StdFitVertexChi2rEEEE[jjj] << " " 
	   <<  StdFitVertexProbEEEE[jjj] << endl;*/

      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){
	  if (i > 3) break;
	  //@//  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrackEEEE_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrackEEEE_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrackEEEE_PHI[i][jjj]=refit_tks[i].phi();
	}
      }

      jjj++;      

    }



    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByToken(RECOcollNameEEEE, CandidatesEEEE);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMassEEEE, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesEEEE,kk);
      //@// cout << "Original 4e mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_EEEE_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }


    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexEEEE, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexXEEEE[jjj]=cand->position().x();
      KinFitVertexYEEEE[jjj]=cand->position().y();
      KinFitVertexZEEEE[jjj]=cand->position().z();
      KinFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      //@//   
      /* cout << "Kin Fit EEEE: " 
	   <<  KinFitVertexXEEEE[jjj] << " " 
	   <<  KinFitVertexYEEEE[jjj] << " " 
	   <<  KinFitVertexZEEEE[jjj] << " " 
	   <<  KinFitVertexChi2rEEEE[jjj] << " " 
	   <<  KinFitVertexProbEEEE[jjj] << endl;*/
      jjj++;
    }
  }


  void fillConstraintVtxDiLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxDiLep_;
    iEvent.getByLabel(StandardFitVertexDiLep, StandardFitVtxDiLep_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxDiLep_->begin(); cand!=StandardFitVtxDiLep_->end(); ++cand){
      StdFitVertexChi2rDiLep[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbDiLep[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      //@//
      /*  cout << "Std Fit DiLeptons: " 
	   <<  StdFitVertexChi2rDiLep[jjj] << " " 
	   <<  StdFitVertexProbDiLep[jjj] << endl;*/
      jjj++;
    }
  }


  void fillConstraintVtxTriLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxMMM_;
    iEvent.getByLabel(StandardFitVertexMMM, StandardFitVtxMMM_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMMM_->begin(); cand!=StandardFitVtxMMM_->end(); ++cand){
      if (jjj > 39) break;
      StdFitVertexChi2rMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
   
      /*   cout << "Std Fit MMM: " 
	   <<  StdFitVertexChi2rMMM[jjj] << " " 
	   <<  StdFitVertexProbMMM[jjj] << endl;*/
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMME_;
    iEvent.getByLabel(StandardFitVertexMME, StandardFitVtxMME_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMME_->begin(); cand!=StandardFitVtxMME_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rMME[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMME[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      /*  cout << "Std Fit MME: " 
	   <<  StdFitVertexChi2rMME[jjj] << " " 
	   <<  StdFitVertexProbMME[jjj] << endl;*/
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxEEE_;
    iEvent.getByLabel(StandardFitVertexEEE, StandardFitVtxEEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxEEE_->begin(); cand!=StandardFitVtxEEE_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      /*  cout << "Std Fit EEE: " 
	   <<  StdFitVertexChi2rEEE[jjj] << " " 
	   <<  StdFitVertexProbEEE[jjj] << endl;*/
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMEE_;
    iEvent.getByLabel(StandardFitVertexMEE, StandardFitVtxMEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMEE_->begin(); cand!=StandardFitVtxMEE_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rMEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      /*  cout << "Std Fit MEE: " 
	   <<  StdFitVertexChi2rMEE[jjj] << " " 
	   <<  StdFitVertexProbMEE[jjj] << endl;*/
      jjj++;
    }

    
  }

  
  		      
  void filljets(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    edm::Handle<pat::JetCollection> pfjets,pfjetsmva;

   iEvent.getByToken(jetsTag_, pfjets);
   iEvent.getByToken(jetsMVATag_, pfjetsmva);

    //@ To retrieve JEC Uncertainty
     edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
     iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);//slimmedJets are ak4PFJetsCHS
     JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
     std::unique_ptr<JetCorrectionUncertainty> jecUnc = std::unique_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(JetCorPar));

/*    
    if (fillMCTruth == 1){
      cout << "Using Jet energy collection for MC " << endl;

      iEvent.getByToken(jetsTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaMCfullDiscr_,puJetIdMVAMC);
    }
    else{
      iEvent.getByToken(jetsDataTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaDatafullDiscr_,puJetIdMVAData);      
      
    }
*/

    RECO_PFJET_N = pfjets->size();
    //@// cout << "Number of PFJets in the event= " << RECO_PFJET_N << endl;
    
    int index_jets = 0;
    
    for ( pat::JetCollection::const_iterator i=pfjets->begin(); i!=pfjets->end(); i++) {  
      if (index_jets>99) continue;
      
      edm::Ref<pat::JetCollection> pfjetref(pfjets,index_jets);      
      edm::Ref<pat::JetCollection> pfjetrefmva(pfjetsmva,index_jets);
      
      float mva = 0.;
      // int pupass = 1;
  
//      mva = i->userFloat("vtxPx");    
      mva = i->userFloat("pileupJetId:fullDiscriminant");

 
      ////////////////////////////////////

     //test Reham    

      bool passLoose_PUID = (i->userInt("pileupJetId:fullId") & (1 << 2));
      bool passMedium_PUID = (i->userInt("pileupJetId:fullId") & (1 << 1));
      bool passTight_PUID = (i->userInt("pileupJetId:fullId") & (1 << 0));
      
      
      // cout<<"passLoose_PUID = "<<passLoose_PUID<<endl;
      // cout<<"passMedium_PUID = "<<passMedium_PUID<<endl;
      // cout<<"passTight_PUID = "<<passTight_PUID<<endl;
      
      
      RECO_PFJET_CHARGE[index_jets] = i->charge();
      RECO_PFJET_ET[index_jets]     = i->et();
      RECO_PFJET_PT[index_jets]     = i->pt();
      RECO_PFJET_ETA[index_jets]    = i->eta();
      RECO_PFJET_PHI[index_jets]    = i->phi();
      //RECO_PFJET_PUID[index_jets]     = pupass;
      RECO_PFJET_PUID_loose[index_jets]     = passLoose_PUID;
      RECO_PFJET_PUID_medium[index_jets]     = passMedium_PUID;
      RECO_PFJET_PUID[index_jets]     = passTight_PUID;
      RECO_PFJET_PUID_MVA[index_jets] = mva;

      //@//
      /*    cout 
	<< "PF Jet with ET= " << RECO_PFJET_ET[index_jets]   
	<< " PT="   << RECO_PFJET_PT[index_jets]   
	<< " ETA="  << RECO_PFJET_ETA[index_jets]  
	<< " PHI="  << RECO_PFJET_PHI[index_jets]  
	<< " PUID=" << RECO_PFJET_PUID[index_jets] 
	<< " PUID_MVA=" << RECO_PFJET_PUID_MVA[index_jets]
       //qier test
        << " Uncorrected Pt=" << i->correctedP4("Uncorrected").Pt()
        << " L1FastJet Pt=" << i->correctedP4("L1FastJet").Pt()
        << " L2Relative Pt=" << i->correctedP4("L2Relative").Pt()
        << " L3Absolute Pt=" << i->correctedP4("L3Absolute").Pt()
	<< endl;*/

      //Branches for Uncertinity

      jecUnc->setJetEta( i->eta() );
      jecUnc->setJetPt( i->pt() ); // here you must use the CORRECTED jet pt
      double unc = jecUnc->getUncertainty(true);
      RECO_PFJET_PT_UncUp[index_jets] = i->pt() + unc;
      RECO_PFJET_PT_UncDn[index_jets] = i->pt() - unc;
      RECO_PFJET_AREA[index_jets] = i->jetArea();
      RECO_PFJET_CHARGED_HADRON_ENERGY[index_jets] = i->chargedHadronEnergy();
      RECO_PFJET_NEUTRAL_HADRON_ENERGY[index_jets] = i->neutralHadronEnergy();
      RECO_PFJET_PHOTON_ENERGY[index_jets] = i->photonEnergy();
      RECO_PFJET_ELECTRON_ENERGY[index_jets] = i->electronEnergy();
      RECO_PFJET_MUON_ENERGY[index_jets] = i->muonEnergy();
      RECO_PFJET_HF_HADRON_ENERGY[index_jets] = i->HFHadronEnergy();
      RECO_PFJET_HF_EM_ENERGY[index_jets] = i->HFEMEnergy();
      RECO_PFJET_CHARGED_EM_ENERGY[index_jets] = i->chargedEmEnergy();
      RECO_PFJET_CHARGED_MU_ENERGY[index_jets] = i->chargedMuEnergy();
      RECO_PFJET_NEUTRAL_EM_ENERGY[index_jets] = i->neutralEmEnergy();
      RECO_PFJET_CHARGED_HADRON_MULTIPLICITY[index_jets] = i->chargedHadronMultiplicity();
      RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY[index_jets] = i->neutralHadronMultiplicity();
      RECO_PFJET_PHOTON_MULTIPLICITY[index_jets] = i->photonMultiplicity();
      RECO_PFJET_ELECTRON_MULTIPLICITY[index_jets] = i->electronMultiplicity();
      RECO_PFJET_MUON_MULTIPLICITY[index_jets] = i->muonMultiplicity();
      RECO_PFJET_HF_HADRON_MULTIPLICTY[index_jets] = i->HFHadronMultiplicity();
      RECO_PFJET_HF_EM_MULTIPLICITY[index_jets] = i->HFEMMultiplicity();
      RECO_PFJET_CHARGED_MULTIPLICITY[index_jets] = i->chargedMultiplicity();
      RECO_PFJET_NEUTRAL_MULTIPLICITY[index_jets] = i->neutralMultiplicity();

      //Stores jet components (jet substructure)

      int Ncomponents = i->numberOfDaughters();
      RECO_PFJET_NCOMPONENTS[index_jets] = Ncomponents;
      std::vector<float> o_jc_pt, o_jc_eta, o_jc_phi, o_jc_energy, o_jc_charge, o_jc_mt, o_jc_vx, o_jc_vy, o_jc_vz, o_jc_pdgid, o_jc_vchi2;
      float sumpt = 0, sumpt2 = 0;
      for(int idau=0; idau<Ncomponents; ++idau){
        const reco::Candidate *jetComponent = i->daughter(idau);
        float jetComponentPt = jetComponent->pt();
        sumpt  += jetComponentPt;
        sumpt2 += (jetComponentPt*jetComponentPt);
        o_jc_pt      .push_back( jetComponentPt );
        o_jc_eta     .push_back( jetComponent->eta() );
        o_jc_phi     .push_back( jetComponent->phi() );
        o_jc_energy  .push_back( jetComponent->energy() );
        o_jc_charge  .push_back( jetComponent->charge() );
        o_jc_mt      .push_back( jetComponent->mt() );//transverse mass
        o_jc_vx      .push_back( jetComponent->vx() );
        o_jc_vy      .push_back( jetComponent->vy() );
        o_jc_vz      .push_back( jetComponent->vz() );
        o_jc_pdgid   .push_back( jetComponent->pdgId() );
        o_jc_vchi2   .push_back( jetComponent->vertexNormalizedChi2() );
      }
      RECO_PFJET_PTD[index_jets] = sqrt(sumpt2)/sumpt;

      //Sort the jet components by pt (to keep thing organized for MVA)

      for(int i1=0; i1 < Ncomponents; ++i1){
        float maxpt = -1, ci2 = -1;
        for(int i2=0; i2 < Ncomponents; ++i2){
          if(o_jc_pt.at(i2) > maxpt){
            maxpt = o_jc_pt[i2];
            ci2 = i2;
          }
        }
        RECO_PFJET_COMPONENT_PDGID[index_jets][i1] = o_jc_pdgid[ci2];
        RECO_PFJET_COMPONENT_PT[index_jets][i1] = o_jc_pt[ci2];
        RECO_PFJET_COMPONENT_ETA[index_jets][i1] = o_jc_eta[ci2];
        RECO_PFJET_COMPONENT_PHI[index_jets][i1] = o_jc_phi[ci2];
        RECO_PFJET_COMPONENT_E[index_jets][i1] = o_jc_energy[ci2];
        RECO_PFJET_COMPONENT_CHARGE[index_jets][i1] = o_jc_charge[ci2];
        RECO_PFJET_COMPONENT_TRANSVERSE_MASS[index_jets][i1] = o_jc_mt[ci2];
        RECO_PFJET_COMPONENT_XVERTEX[index_jets][i1] = o_jc_vx[ci2];
        RECO_PFJET_COMPONENT_YVERTEX[index_jets][i1] = o_jc_vy[ci2];
        RECO_PFJET_COMPONENT_ZVERTEX[index_jets][i1] = o_jc_vz[ci2];
        RECO_PFJET_COMPONENT_VERTEX_CHI2[index_jets][i1] = o_jc_vchi2[ci2];
        o_jc_pt[ci2] = -2;
      }
    
      //@// 
     /*    cout 
	<< "PF Jet with ET= " << RECO_PFJET_ET[index_jets]   
	<< " PT="   << RECO_PFJET_PT[index_jets]   
	<< " ETA="  << RECO_PFJET_ETA[index_jets]  
	<< " PHI="  << RECO_PFJET_PHI[index_jets]  
	<< " PUID=" << RECO_PFJET_PUID[index_jets] 
	<< " PUID_MVA=" << RECO_PFJET_PUID_MVA[index_jets]
       //qier test
        << " Uncorrected Pt=" << i->correctedP4("Uncorrected").Pt()
        << " L1FastJet Pt=" << i->correctedP4("L1FastJet").Pt()
        << " L2Relative Pt=" << i->correctedP4("L2Relative").Pt()
        << " L3Absolute Pt=" << i->correctedP4("L3Absolute").Pt()
	<< endl;*/
      
      index_jets++;

      
    } // for loop on PFJets jets


    /////////////////////////////////////

    //QGTagger REHAM
    
    edm::Handle<edm::ValueMap<float> > QGjetmap;
    iEvent.getByToken(jetQGMapTag_, QGjetmap);

    edm::Handle<edm::ValueMap<float> > QGjetmap_axis2;
    iEvent.getByToken(jetQGMapTag_axis2_, QGjetmap_axis2);

    edm::Handle<edm::ValueMap<float> > QGjetmap_ptd;
    iEvent.getByToken(jetQGMapTag_ptd_, QGjetmap_ptd);

    edm::Handle<edm::ValueMap<int> > QGjetmap_mult;
    iEvent.getByToken(jetQGMapTag_mult_, QGjetmap_mult);

    edm::Handle<pat::JetCollection> pfjetsQG;
    iEvent.getByToken(jetsQGTag_, pfjetsQG); 
    
    int N_QGjets = 0;
    
    for ( pat::JetCollection::const_iterator i=pfjetsQG->begin(); i!=pfjetsQG->end(); i++) {  
      if (N_QGjets>99) continue;
    
      edm::Ref<pat::JetCollection> pfjetref(pfjetsQG , N_QGjets);     

      float qgLikelihood = (*QGjetmap)[pfjetref];
      float qgaxis2 = (*QGjetmap_axis2)[pfjetref];
      float qgptd = (*QGjetmap_ptd)[pfjetref];
      int qgmult = (*QGjetmap_mult)[pfjetref];

      // cout <<"qgLikelihood = "<<qgLikelihood<<endl;

      RECO_PFJET_QG_Likelihood[N_QGjets] = qgLikelihood;
      RECO_PFJET_QG_axis2[N_QGjets] = qgaxis2;
      RECO_PFJET_QG_ptd[N_QGjets] = qgptd;
      RECO_PFJET_QG_mult[N_QGjets] = qgmult;

      // cout<<"qgLikelihood branch = "<< RECO_PFJET_QG_Likelihood[N_QGjets]<<endl;

      N_QGjets++;

    }

    /////////////////////////////////////////////

    //@
    //For Q/G study
    
    if( fillMCTruth ){
      //Prepare LHE particles list (only the partons in final state and check if they have DR=0.4 clean)

      //cout<<"MCCCCCCCCCCCCC onlyyyyyyyyyyy"<<endl;

      edm::Handle<LHEEventProduct> lheInfo;
      iEvent.getByToken(LHE_, lheInfo);
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheInfo->hepeup().PUP;
      int NlheParticles = lheParticles.size();
      int lhep_index = 0;
      bool lhep_clear = true;
      TLorentzVector lheParton1, lheParton2;
      for(int ip1=0; ip1<NlheParticles; ++ip1){
    	int lhep1_pdgid = abs(lheInfo->hepeup().IDUP.at(ip1));
    	int lhep1_state = lheInfo->hepeup().ISTUP.at(ip1);
        if(lhep1_state == 1 && ((lhep1_pdgid >= 1 && lhep1_pdgid <= 8) || lhep1_pdgid == 21)){
          lheParton1.SetPxPyPzE(lheParticles[ip1][0],lheParticles[ip1][1],lheParticles[ip1][2],lheParticles[ip1][3]);
          for(int ip2=0; ip2<NlheParticles; ++ip2){
    	      int lhep2_pdgid = abs(lheInfo->hepeup().IDUP.at(ip2));
    	      int lhep2_state = lheInfo->hepeup().ISTUP.at(ip2);
            if(ip2 != ip1 && lhep2_state == 1 && ((lhep2_pdgid >= 1 && lhep2_pdgid <= 8) || lhep2_pdgid == 21)){
              lheParton2.SetPxPyPzE(lheParticles[ip2][0],lheParticles[ip2][1],lheParticles[ip2][2],lheParticles[ip2][3]);
              //Check parton isolation (is there different partons inside PF Jet (Ak4) radius?)
              //If so, then the PFJet-Parton matching is of course not the best (we might want to avoid those jets?)
              //Does this take in account the case of different PFJets matching better to same parton?
              if( deltaR2(lheParton1.Eta(), lheParton1.Phi(), lheParton2.Eta(), lheParton2.Phi()) < 0.4 ) lhep_clear = false;
            }
          }
          LHE_PARTON_CLEAR[lhep_index] = lhep_clear;
          LHE_PARTON_PDGID[lhep_index] = lhep1_pdgid;
    	  LHE_PARTON_PT[lhep_index]    = lheParton1.Pt();
    	  LHE_PARTON_ETA[lhep_index]   = lheParton1.Eta();
    	  LHE_PARTON_PHI[lhep_index]   = lheParton1.Phi();
    	  LHE_PARTON_E[lhep_index]     = lheParton1.E();
    	  ++lhep_index;
        }
      }
      LHE_PARTON_N = lhep_index;
    }

    ////////////////////////////////////
    
    edm::Handle<double> rhoHandle;

    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_mu=*rhoHandle;
       cout << "RHO mu fastjet= " << RHO_mu << endl; 
    }
    else {
       cout << "Not valid RHO mu collection" << endl;
    }
    
    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_ele=*rhoHandle;
      cout << "RHO ele fastjet= " << RHO_ele << endl; 
    }
    else {
      cout << "Not valid RHO ele collection" << endl;
    } 
    
  }
  


  void fillBTagging(const edm::Event& iEvent){

    // trackCountingHighEffBJetTags
    edm::Handle<pat::JetCollection> bTagHandle;
    iEvent.getByToken(jetsTag_ ,bTagHandle);
    int l=0;
    for (pat::JetCollection::const_iterator btagIter=bTagHandle->begin(); btagIter!=bTagHandle->end();++btagIter) {
      if(l>=49) continue;
      double discrCSV1 = btagIter->bDiscriminator(tCHighEff_bTag_);
      //@//     
      /*	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighEffBJetTags = "<< discrCSV1 
	    << " and jet Pt = "<<btagIter->pt()<<endl;   */  
	tCHighEff_BTagJet_PT[l]=btagIter->pt();
	tCHighEff_BTagJet_ETA[l]=btagIter->eta();
	tCHighEff_BTagJet_PHI[l]=btagIter->phi();
	tCHighEff_BTagJet_DISCR[l]=discrCSV1;

    // trackCountingHighPurBJetTags
     double discrCSV2 = btagIter->bDiscriminator(tCHighPur_bTag_);
     //@//
     /*	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighPurBJetTags = "<< discrCSV2
	    << " and jet Pt = "<<btagIter->pt()<<endl; */     
	tCHighPur_BTagJet_PT[l]=btagIter->pt();
	tCHighPur_BTagJet_ETA[l]=btagIter->eta();
	tCHighPur_BTagJet_PHI[l]=btagIter->phi();
	tCHighPur_BTagJet_DISCR[l]=discrCSV2;

    // Deep CSV tagger 2017 Reham
    double discrCSV3_1 = btagIter->bDiscriminator(cSV_bTag1_);
    double discrCSV3_2 = btagIter->bDiscriminator(cSV_bTag2_);
    double discrCSV3 = discrCSV3_1 + discrCSV3_2 ;

    //@//
    /*	cout<<" Jet "<< l
	    <<" has b tag discriminator pfDeepCSVJetTags:probb  = "<< discrCSV3_1<<" and pfDeepCSVJetTags:probbb  = "<<discrCSV3_2<<" and total = "
            <<discrCSV3 << " and jet Pt = "<<btagIter->pt()<<endl; */     
	cSV_BTagJet_PT[l]=btagIter->pt();
	cSV_BTagJet_ETA[l]=btagIter->eta();
	cSV_BTagJet_PHI[l]=btagIter->phi();
	cSV_BTagJet_ET[l]=btagIter->et();
	cSV_BTagJet_DISCR[l]=discrCSV3;
      l++;
    }

  }
  
  
  void fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const {
    double dp = PFEnergyResolution().getEnergyResolutionEm(c.energy(), c.eta());
    // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
    //            [ Px/P  ]                
    //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
    //            [ Pz/P  ]                
    AlgebraicMatrix31 ptop3;
    ptop3(0,0) = c.px()/c.p();
    ptop3(1,0) = c.py()/c.p();
    ptop3(2,0) = c.pz()/c.p();
    AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
	bigCov(offset+i,offset+j) = mat(i,j);
      } } 
  }
  
  
  
 private:
 
  std::vector<std::string> module_to_search;
  std::string par_to_search;
  
  // ROOT definition
  TFile *theFile_ ;
  TTree *theTree_ ;
  
  bool useRECOformat;
  std::string inputfileName;

  // Input tags
  std::string decaychannel;
  std::string flaginst;
  std::vector<std::string> flagtags;
  std::string rootFileName;

  // PU
  bool fillPUinfo;
  int num_PU_vertices;
  int PU_BunchCrossing; // bunch crossing for the PU event added
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSrc_;     

  // MC weight
  float MC_weighting;
  edm::EDGetTokenT<GenEventInfoProduct> generator_;

  // HLT
  bool fillHLTinfo;
  edm::InputTag HLTInfoFired;
  std::string HLTAnalysisinst;
  std::vector<edm::InputTag> flagHLTnames; 

  edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent;     

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  //Reham Read MET filter decision from trigger results

  edm::EDGetTokenT<edm::TriggerResults> noiseFilterTag_; 


  edm::InputTag triggerMatchObjectEle;
  edm::EDGetTokenT<edm::Association<std::vector<pat::TriggerObjectStandAlone> > > triggerMatchObject;
  std::string triggerFilter,triggerEleFilter,triggerHLTcollection;
  
  
  // SkimEarlyData
  std::string SkimEarlyDataAnalysisinst;
  std::vector<edm::InputTag> flagSkimEarlyDatanames; 
 
  // MC truth
  bool fillMCTruth;
  edm::EDGetTokenT<edm::View<reco::Candidate> > MCcollName;
  // GenParticles;
  edm::EDGetTokenT<vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > fourgenleptons_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > digenZ_;

	
  // RECO
  bool useAdditionalRECO;
  bool use2011EA;
  std::vector<edm::InputTag> RECOcollNameBest2e2mu,RECOcollNameBest4e,RECOcollNameBest4mu,
    RECOcollNameBestRestFrame2e2mu,RECOcollNameBestRestFrame4mu,RECOcollNameBestRestFrame4e;
   edm::EDGetTokenT<edm::View<reco::Candidate> > zToMuMu, zToEE,
    zToMuMussmerge, zToEEssmerge, zToCrossLeptons,
    RECOcollNameDiLep_,RECOcollNameMMMM_,RECOcollNameEEEE,RECOcollNameEEMM,
    quadLeptons4Mu, quadLeptons2Mu2E, quadLeptons4E,
    triLeptonsMuMuMu, triLeptonsMuMuE, triLeptonsMuEE, triLeptonsEEE,
    quadLeptons3Mu1E, quadLeptons3E1Mu, quadLeptonsSSOSele,
    quadLeptonsSSOSmu, quadLeptonsSSOSelemu, quadLeptonsSSOSmuele,
    RECOcollNameLLLL;
    std::vector<edm::InputTag> RECOcollNameZ,RECOcollNameZss,
 //   RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,
    RECOcollNameLLLLss,RECOcollNameLLL,RECOcollNameLLLl,RECOcollNameLLLLssos;
  edm::InputTag RECOcollNameDiLep;//RECOcollNameLLLL,RECOcollNameDiLep;
  
  // electron and muon tags
  bool useBestCandidate;
  edm::InputTag BestCandidatesLeptonsTag_;
  edm::InputTag clusterCollectionTag_,gsftrackCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonPFTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronEgmTag_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_;

  edm::EDGetTokenT<edm::ValueMap<float> > muonCorrPtErrorMapTag_;
    
  edm::InputTag electronEgmTkMapTag_;
  edm::InputTag electronEgmEcalMapTag_;
  edm::InputTag electronEgmHcalMapTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > mvaElectronTag_;
  ///@@@/// TEST Electron ID  REHAM
  // edm::EDGetTokenT<edm::ValueMap<float> > /*//@// mvaTrigV0MapTag_,*/mvaNonTrigV0MapTag_;
  //edm::EDGetToken eleIDToken_; 
  edm::EDGetToken electronsMiniAODToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleIDToken_;
  // edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
/*
  edm::EDGetTokenT<edm::ValueMap<double> >
    electronPFIsoValueChargedAllTag_,
    electronPFIsoValueChargedTag_,
    electronPFIsoValueNeutralTag_,
    electronPFIsoValueGammaTag_,
    electronPFIsoValuePUTag_;   


  edm::InputTag eleRegressionEnergyErrorTag_,eleRegressionEnergyTag_;

  edm::EDGetTokenT<edm::ValueMap<double> >
    muonPFIsoValueChargedAllTag_,
    muonPFIsoValueChargedTag_,
    muonPFIsoValueNeutralTag_,
    muonPFIsoValueGammaTag_,
    muonPFIsoValuePUTag_;   


  edm::EDGetTokenT<edm::ValueMap<double> >
    photonPFIsoValueChargedAllTag_,
    photonPFIsoValueChargedTag_,
    photonPFIsoValueNeutralTag_,
    photonPFIsoValueGammaTag_,
    photonPFIsoValuePUTag_;   

*/


  // vertexing 3D 
  edm::InputTag electronTag_Vert,muonTag_Vert;
  edm::EDGetTokenT<edm::ValueMap<float> > electronMapTag_Vert,electronMapTag_VertKF,muonMapTag_Vert,muonMapTag_VertKF,
    electronMapTag_VertValue,electronMapTag_VertValueKF,muonMapTag_VertValue,muonMapTag_VertValueKF,
    electronMapTag_VertError,electronMapTag_VertErrorKF,muonMapTag_VertError,muonMapTag_VertErrorKF;

  edm::InputTag electronMapTag_VertGD,muonMapTag_VertGD;
  edm::InputTag electronMapTag_VertGDEEEE,muonMapTag_VertGDMMMM;
  edm::InputTag electronMapTag_VertStd,muonMapTag_VertStd;
  edm::InputTag electronMapTag_VertStdEEEE,muonMapTag_VertStdMMMM;
  edm::InputTag electronMapTag_VertKin,muonMapTag_VertKin;
  edm::InputTag electronMapTag_VertKinEEEE,muonMapTag_VertKinMMMM;

  // vertexing STIP / SLIP
  edm::EDGetTokenT<edm::ValueMap<float> >  electronSTIPMapTag_Vert,muonSTIPMapTag_Vert,
    electronSLIPMapTag_Vert,muonSLIPMapTag_Vert,
    electronSTIPMapTag_VertValue,muonSTIPMapTag_VertValue,
    electronSLIPMapTag_VertValue,muonSLIPMapTag_VertValue,
    electronSTIPMapTag_VertError,muonSTIPMapTag_VertError,
    electronSLIPMapTag_VertError,muonSLIPMapTag_VertError;
  
  // vertexing GD
  edm::InputTag ftsigma_Vert,ftsigmalag_Vert,
                ftsigma_VertMMMM,ftsigmalag_VertMMMM,
                ftsigma_VertEEEE,ftsigmalag_VertEEEE;
  edm::InputTag 
    gdX_Vert,gdY_Vert,gdZ_Vert,
    gdX_VertMMMM,gdY_VertMMMM,gdZ_VertMMMM,
    gdX_VertEEEE,gdY_VertEEEE,gdZ_VertEEEE;
   edm::InputTag 
    gdlagX_Vert,gdlagY_Vert,gdlagZ_Vert,gdlagProb_Vert,gdlagNdof_Vert,
     gdlagX_VertMMMM,gdlagY_VertMMMM,gdlagZ_VertMMMM,gdlagProb_VertMMMM,gdlagNdof_VertMMMM,
     gdlagX_VertEEEE,gdlagY_VertEEEE,gdlagZ_VertEEEE,gdlagProb_VertEEEE,gdlagNdof_VertEEEE;
  //

   // ConstraintVtx
   edm::InputTag 
     StandardFitVertex,StandardFitVertexMMMM,StandardFitVertexEEEE,
     KinematicFitVertex,KinematicFitVertexMMMM,KinematicFitVertexEEEE,
     RefittedMass,RefittedMassMMMM,RefittedMassEEEE;

   edm::InputTag 
     StandardFitVertexEEE,StandardFitVertexMMM,StandardFitVertexMEE,StandardFitVertexMME,StandardFitVertexDiLep;

  //electronID
  std::vector<edm::InputTag> eleIDTag_;
  edm::InputTag EleID_VeryLooseTag_ ;
  edm::InputTag EleID_LooseTag_  ;
  edm::InputTag EleID_MediumTag_  ;
  edm::InputTag EleID_TightTag_ ;

  edm::InputTag EleID_HZZVeryLooseTag_ ;
  edm::InputTag EleID_HZZLooseTag_  ;
  edm::InputTag EleID_HZZMediumTag_  ;
  edm::InputTag EleID_HZZTightTag_ ;
  
  edm::InputTag allelectronsColl, allmuonsColl;
  edm::InputTag isoVarTagElectronsCal,isoVarTag,isoVarTagElectronsTracker,isoVarTagElectronsX;
  edm::InputTag isoVarTagMuonsHCalIso,isoVarTagMuonsECalIso,
    isoVarTagMuonsCalIso,isoVarTagMuonsTracker,isoVarTagMuonsX;
  
  edm::InputTag theECALIsoDepositLabel;    //EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel;    //Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel;   //Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel; //Tracker Isolation deposit label 
  
  // Photon, Tracks, Jets, Vertices
//  edm::EDGetTokenT<vector<reco::Track> > tracksTag_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfTag_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonsTag_;
  vector<reco::Vertex> PV;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesTag_;
  edm::EDGetTokenT<double> rhojetsTag_;
  edm::EDGetTokenT<vector<pat::Jet> > jetsTag_,jetsMVATag_, jetsQGTag_;
  edm::EDGetTokenT<edm::ValueMap<float> > jetQGMapTag_;
  edm::EDGetTokenT<edm::ValueMap<float> > jetQGMapTag_axis2_;
  edm::EDGetTokenT<edm::ValueMap<float> > jetQGMapTag_ptd_;
  edm::EDGetTokenT<edm::ValueMap<int> > jetQGMapTag_mult_;

//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaMCfullDiscr_,PuJetMvaMCfullId_;
//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaDatafullDiscr_,PuJetMvaDatafullId_;

  // genJet
  edm::EDGetTokenT<vector<reco::GenJet> >  genjetTag_;
  
  // MET
  edm::InputTag trackermetTag_;
  edm::EDGetTokenT<vector<pat::MET> >  pfmetTag_;
  edm::EDGetTokenT<vector<reco::GenMET> > genmetTag_;
  edm::InputTag calometTag_,calometoptTag_,calometoptnohfTag_,calometoptnohfhoTag_;
  edm::InputTag calometopthoTag_,calometnohfTag_,calometnohfhoTag_,calomethoTag_;
  bool useAdditionalMET_;
  edm::InputTag htmetic5Tag_,htmetkt4Tag_,htmetkt6Tag_,htmetsc5Tag_,htmetsc7Tag_;
  edm::InputTag jescormetic5Tag_,jescormetkt4Tag_,jescormetkt6Tag_,jescormetsc5Tag_,jescormetsc7Tag_;
  edm::InputTag cormetMuTag_;
  //edm::EDGetTokenT<bool>ecalBadCalibFilterUpdate_token ;// new met filter to be used (unser test)
  // edm::EDGetTokenT<bool> BadChCandFilterToken_;
  // edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

  // Conversion
//  edm::EDGetTokenT<edm::ValueMap<float> > ConvMapDistTag_,ConvMapDcotTag_;

  // Matching
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodElectronMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myElectrons_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodMuonMCMatch_; 
  edm::EDGetTokenT<reco::CandidateCollection> myMuons_; 
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodGammaMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myGammas_;

  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoMuMuMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoEEMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToMMMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEEEMCMatch_;
  // Beam Spot
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;

  // bTagging
//  edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >  
  std::string tCHighEff_bTag_, tCHighPur_bTag_,cSV_bTag1_,cSV_bTag2_;

//  const std::vector<std::string> bDiscriminators_;

  // counters
  int irun, ievt,ils,nevt;
  float Avginstlumi;

  // HLT
  int RECO_nMuHLTMatch,RECO_nEleHLTMatch;
  ArrayVector<float> RECOMU_PT_MuHLTMatch,RECOMU_ETA_MuHLTMatch,RECOMU_PHI_MuHLTMatch;
  ArrayVector<float> RECOELE_PT_EleHLTMatch,RECOELE_ETA_EleHLTMatch,RECOELE_PHI_EleHLTMatch;
  
  char HLTPathsFired[20000];

 
  // tmp candidate collections
//  reco::CandidateCollection *leptonscands2e2mu_;//unused branch
//  reco::CandidateCollection *leptonscands2e2murf_;//unused branch
//  reco::CandidateCollection *leptonscands4mu_;//unused branch
//  reco::CandidateCollection *leptonscands4murf_;//unused branch
//  reco::CandidateCollection *leptonscands4e_;//unused branch
//  reco::CandidateCollection *leptonscands4erf_;//unused branch
  
  reco::CandidateCollection *leptonscands_Z0;
  reco::CandidateCollection *leptonscands_Z1;
  reco::CandidateCollection *leptonscands_Zss0;
  reco::CandidateCollection *leptonscands_Zss1;
  reco::CandidateCollection *leptonscands_Zcross;
  reco::CandidateCollection *leptonscands_DiLep;
  reco::CandidateCollection *leptonscands_MMMM;
  reco::CandidateCollection *leptonscands_EEEE;
  reco::CandidateCollection *leptonscands_EEMM;
//  reco::CandidateCollection *leptonscands_LLL0;//unused branch
//  reco::CandidateCollection *leptonscands_LLL1;//unused branch
//  reco::CandidateCollection *leptonscands_LLL2;//unused branch
//  reco::CandidateCollection *leptonscands_LLL3;//unused branch
//  reco::CandidateCollection *leptonscands_LLLLss0;//unused branch
//  reco::CandidateCollection *leptonscands_LLLLss1;//unused branch
//  reco::CandidateCollection *leptonscands_LLLLss2;//unused branch
//  reco::CandidateCollection *leptonscands_LLLl0;//unused branch
//  reco::CandidateCollection *leptonscands_LLLl1;//unused branch
//  reco::CandidateCollection *leptonscands_LLLL;//unused branch


  // MC info
  edm::ESHandle<ParticleDataTable>  pdt_;
  
  // MC truth
  ArrayVector<float> MC_E,MC_PT,MC_ETA,MC_THETA,MC_PHI,MC_MASS,MC_PDGID;

  ArrayVector<double> MCRF_cosTheta1_spin, MCRF_cosTheta2_spin, MCRF_cosThetaStar_spin, MCRF_Phi_spin, 
    MCRF_Phi1_spin, MCRF_Phi2_spin, MCRF_phi1RF_spin, MCRF_phi2RF_spin, MCRF_MELA;
  
  ArrayVector<float> MC_LEPT_PT,MC_LEPT_ETA,MC_LEPT_PHI,MC_LEPT_THETA,MC_LEPT_PDGID;
  ArrayVector<ArrayVector<float> > MC_Z_MASS,MC_Z_PT,MC_Z_ETA,MC_Z_PHI,MC_Z_THETA,MC_Z_PDGID;

  ArrayVector<ArrayVector<float> > MC_fourl_MASS,MC_fourl_PT,MC_fourl_PDGID;
  ArrayVector<ArrayVector<float> > MC_ZZ_MASS,MC_ZZ_PT,MC_ZZ_ETA,MC_ZZ_PHI,MC_ZZ_THETA,MC_ZZ_PDGID;

  // RECO collection
 
  
  // RECORF
    
  ArrayVector<double> RECORF_2e2mu_cosTheta1_spin, RECORF_2e2mu_cosTheta2_spin, RECORF_2e2mu_cosThetaStar_spin, RECORF_2e2mu_Phi_spin, 
    RECORF_2e2mu_Phi1_spin, RECORF_2e2mu_Phi2_spin, RECORF_2e2mu_phi1RF_spin, RECORF_2e2mu_phi2RF_spin,RECORF_2e2mu_MELA;
  
    
  ArrayVector<double> RECORF_4e_cosTheta1_spin, RECORF_4e_cosTheta2_spin, RECORF_4e_cosThetaStar_spin, RECORF_4e_Phi_spin,RECORF_4e_MELA, 
    RECORF_4e_Phi1_spin, RECORF_4e_Phi2_spin, RECORF_4e_phi1RF_spin, RECORF_4e_phi2RF_spin,RECORF_4mu_MELA;
  
    
  ArrayVector<double> RECORF_4mu_cosTheta1_spin, RECORF_4mu_cosTheta2_spin, RECORF_4mu_cosThetaStar_spin, RECORF_4mu_Phi_spin, 
    RECORF_4mu_Phi1_spin, RECORF_4mu_Phi2_spin, RECORF_4mu_phi1RF_spin, RECORF_4mu_phi2RF_spin;
  

  int leptonflavor;

  // RECO additional
  ArrayVector<float> 
    RECO_ZMM_MASS,
    RECO_ZEE_MASS,
    RECO_ZMMss_MASS,
    RECO_ZEEss_MASS,
    RECO_ZEM_MASS,
    RECO_DiLep_MASS;
  ArrayVector<ArrayVector<float> >
    RECO_ZMM_PT,RECO_ZMM_ETA,RECO_ZMM_PHI,
    RECO_ZEE_PT,RECO_ZEE_ETA,RECO_ZEE_PHI,
    RECO_ZMMss_PT,RECO_ZMMss_ETA,RECO_ZMMss_PHI,
    RECO_ZEEss_PT,RECO_ZEEss_ETA,RECO_ZEEss_PHI,
    RECO_ZEM_PT,RECO_ZEM_ETA,RECO_ZEM_PHI,
    RECO_DiLep_PT,RECO_DiLep_ETA,RECO_DiLep_PHI;

  ArrayVector<ArrayVector<float> > 
    RECO_EEMM_MASS,RECO_MMMM_MASS,RECO_EEEE_MASS,
    RECO_EEMM_PT,  RECO_MMMM_PT,  RECO_EEEE_PT,
    RECO_EEMM_ETA,RECO_MMMM_ETA,RECO_EEEE_ETA,
    RECO_EEMM_PHI,  RECO_MMMM_PHI,  RECO_EEEE_PHI,
    RECO_LLLL_MASS,RECO_LLLL_PT,RECO_LLLL_ETA,RECO_LLLL_PHI;

  ArrayVector<float>
    RECO_MMMM_MASS_REFIT,RECO_EEMM_MASS_REFIT,RECO_EEEE_MASS_REFIT;

  ArrayVector<float> RECO_LLL0_MASS,RECO_LLL1_MASS,RECO_LLL2_MASS,RECO_LLL3_MASS;
  ArrayVector<ArrayVector<float> > RECO_LLL0_PT,RECO_LLL1_PT,RECO_LLL2_PT,RECO_LLL3_PT;

  ArrayVector<float> RECO_LLLl0_MASS,RECO_LLLl1_MASS;
  ArrayVector<ArrayVector<float> > RECO_LLLl0_PT,RECO_LLLl1_PT;

  ArrayVector<float> 
    RECO_LLLL0ss_MASS,
    RECO_LLLL1ss_MASS,
    RECO_LLLL2ss_MASS;
  ArrayVector<ArrayVector<float> > 
    RECO_LLLL0ss_PT,
    RECO_LLLL1ss_PT,
    RECO_LLLL2ss_PT;
  // RECO electrons
  edm::ESHandle<CaloGeometry> theCaloGeom_;  
  ArrayVector<float> RECOELE_E,RECOELE_PT,RECOELE_PTError,RECOELE_P,RECOELE_ETA,RECOELE_THETA,RECOELE_PHI,RECOELE_MASS;
  ArrayVector<float> RECOELE_CHARGE,RECOELE_ID, RECOELE_PT_uncorr ;

  ArrayVector<int > RECOELE_isEcalDriven, RECOELE_isTrackerDriven;
  ArrayVector<float> 
    RECOELE_gsftrack_NPixHits,
    RECOELE_gsftrack_NStripHits,
    RECOELE_gsftrack_chi2, 
    RECOELE_gsftrack_dxyB, RECOELE_gsftrack_dxy, RECOELE_gsftrack_dxyError,
    RECOELE_gsftrack_dzB, RECOELE_gsftrack_dz,RECOELE_gsftrack_dzError;
  ArrayVector<int> RECOELE_gsftrack_losthits,RECOELE_gsftrack_validhits,RECOELE_gsftrack_expected_inner_hits; 
  ArrayVector<float>
    RECOELE_scl_E,RECOELE_scl_Et,RECOELE_scl_Eta,RECOELE_scl_Phi; 
  //float RECOELE_eSuperClusterOverP[100], RECOELE_eSeedClusterOverPout[100], RECOELE_deltaEtaSuperClusterTrackAtVtx[100], RECOELE_deltaPhiSuperClusterTrackAtVtx[100];
  ArrayVector<float> 
    RECOELE_ep, RECOELE_eSeedp, RECOELE_eSeedpout, RECOELE_eElepout,
    RECOELE_deltaEtaIn,RECOELE_deltaEtaSeed,RECOELE_deltaEtaEle,RECOELE_deltaPhiIn,
    RECOELE_deltaPhiSeed,RECOELE_deltaPhiEle,RECOELE_ecalEnergy;
  ArrayVector<int> RECOELE_isbarrel, RECOELE_isendcap, RECOELE_isEBetaGap, RECOELE_isEBphiGap, RECOELE_isEEdeeGap, RECOELE_isEEringGap, RECOELE_isGap; 
  ArrayVector<float> RECOELE_sigmaIetaIeta, RECOELE_sigmaEtaEta, RECOELE_e15, RECOELE_e25max, RECOELE_e55, RECOELE_he, RECOELE_r9;
  ArrayVector<float> RECOELE_mva, RECOELE_fbrem,RECOELE_fbrem_mean,RECOELE_fbrem_mode;
  ArrayVector<int> RECOELE_nbrems, RECOELE_Class;
  
  ArrayVector<float> RECOELE_EGMTRACKISO,RECOELE_EGMECALISO,RECOELE_EGMHCALISO,RECOELE_EGMX,
    RECOELE_IP,RECOELE_SIP,RECOELE_IPERROR,
    RECOELE_IP_KF,RECOELE_SIP_KF,RECOELE_IPERROR_KF,
    RECOELE_STIP,RECOELE_TIP,RECOELE_TIPERROR,
    RECOELE_SLIP,RECOELE_LIP,RECOELE_LIPERROR,
    RECOELE_SIP_GD, RECOELE_SIP_GDEEEE,
    RECOELE_SIP_Std, RECOELE_SIP_StdEEEE, 
    RECOELE_SIP_Kin, RECOELE_SIP_KinEEEE;

  ArrayVector<double> RECOELE_PFchAllPart,RECOELE_PFchHad,RECOELE_PFneuHad,RECOELE_PFphoton,
    RECOELE_PFPUchAllPart,RECOELE_PFX_dB,RECOELE_PFX_rho,RECOELE_PF_RingsIsoMVA;

  ArrayVector<double> RECOELE_regEnergy,RECOELE_regEnergyError;

  ArrayVector<float> RECOELE_TLE_ParentSC_X,RECOELE_TLE_ParentSC_Y,RECOELE_TLE_ParentSC_Z;
  
  ArrayVector<int> RECOELE_EEEE_MATCHED,RECOELE_EEMM_MATCHED,RECOELE_ZEE_MATCHED,RECOELE_ZssEE_MATCHED,RECOELE_ZEM_MATCHED,
    RECOELE_LLL0_MATCHED,RECOELE_LLL1_MATCHED,RECOELE_LLL2_MATCHED,RECOELE_LLL3_MATCHED,
    RECOELE_LLLLss0_MATCHED,RECOELE_LLLLss1_MATCHED,RECOELE_LLLLss2_MATCHED,
    RECOELE_LLLl0_MATCHED,RECOELE_LLLl1_MATCHED,RECOELE_LLLL_MATCHED;
  ArrayVector<double> RECOELE_mvaTrigV0,RECOELE_mvaNonTrigV0;
  
  ArrayVector<double> ele_sclRawE ;
  ArrayVector<double> ele_sclX, ele_sclY, ele_sclZ;
  ArrayVector<int> ele_seedSubdet1;
  ArrayVector<double> ele_seedDphi1, ele_seedDrz1;
  ArrayVector<int> ele_seedSubdet2;
  ArrayVector<double> ele_seedDphi2, ele_seedDrz2;
  ArrayVector<double> ele_eidVeryLoose, ele_eidLoose, ele_eidMedium, ele_eidTight ;
  ArrayVector<double> ele_eidHZZVeryLoose, ele_eidHZZLoose, ele_eidHZZMedium, ele_eidHZZTight ;
  ArrayVector<ArrayVector<ArrayVector<double> > > RECOELE_COV;

  //Reham electron systematic variables

  ArrayVector<float> RECOELE_ecalTrkEnergyErrPostCorr,RECOELE_energyScaleValue,RECOELE_energySigmaValue, RECOELE_energyScaleUp, RECOELE_energyScaleDown, RECOELE_energyScaleStatUp, RECOELE_energyScaleStatDown, RECOELE_energyScaleSystUp, RECOELE_energyScaleSystDown, RECOELE_energyScaleGainUp, RECOELE_energyScaleGainDown,RECOELE_energyScaleEtUp, RECOELE_energyScaleEtDown, RECOELE_energySigmaUp, RECOELE_energySigmaDown, RECOELE_energySigmaPhiUp, RECOELE_energySigmaPhiDown, RECOELE_energySigmaRhoUp, RECOELE_energySigmaRhoDown,RECOELE_ecalTrkEnergyPreCorr, RECOELE_ecalTrkEnergyErrPreCorr;

  // RECO muons
  ArrayVector<int > RECOMU_isPFMu,RECOMU_isGlobalMu,RECOMU_isStandAloneMu,RECOMU_isTrackerMu,RECOMU_isCaloMu,RECOMU_isTrackerHighPtMu,RECOMU_isME0Muon;
  ArrayVector<float> RECOMU_E,RECOMU_PT,RECOMU_P,RECOMU_ETA,RECOMU_THETA,RECOMU_PHI,RECOMU_MASS,RECOMU_CHARGE;
  ArrayVector<ArrayVector<ArrayVector<double> > > RECOMU_COV;
  ArrayVector<double> /*,RECOMU_Roch_calib_error,*/ RECOMU_PT_uncorr;

  ArrayVector<float> 
    RECOMU_TRACKISO,RECOMU_TRACKISO_SUMPT,RECOMU_ECALISO,RECOMU_HCALISO, RECOMU_X,   
    RECOMU_IP,RECOMU_SIP,RECOMU_IPERROR,
    RECOMU_IP_KF,RECOMU_SIP_KF,RECOMU_IPERROR_KF,
    RECOMU_STIP,RECOMU_TIP,RECOMU_TIPERROR,
    RECOMU_SLIP,RECOMU_LIP,RECOMU_LIPERROR,
    RECOMU_SIP_GD, RECOMU_SIP_GDMMMM,
    RECOMU_SIP_Std, RECOMU_SIP_StdMMMM, 
    RECOMU_SIP_Kin, RECOMU_SIP_KinMMMM;

  ArrayVector<double> 
    RECOMU_PFchHad,RECOMU_PFneuHad,RECOMU_PFphoton,RECOMU_PFsumPUPt,RECOMU_PFX_dB,RECOMU_PFPUchAllPart,RECOMU_PFchAllPart,RECOMU_PFX_rho,
    RECOMU_PFchHad42,RECOMU_PFneuHad42,RECOMU_PFphoton42,RECOMU_PFPUchAllPart42,RECOMU_PFchAllPart42,
    RECOMU_PF_RingsIsoMVA,RECOMU_PF_RingsIDMVA;

  ArrayVector<float>
    RECOMU_caloCompatibility,RECOMU_segmentCompatibility;
  ArrayVector<int > RECOMU_glbmuPromptTight;
 
  ArrayVector<int> RECOMU_MMMM_MATCHED,RECOMU_EEMM_MATCHED,
      RECOMU_ZMM_MATCHED,RECOMU_ZssMM_MATCHED,RECOMU_ZEM_MATCHED,
      RECOMU_LLL0_MATCHED,RECOMU_LLL1_MATCHED,RECOMU_LLL2_MATCHED,RECOMU_LLL3_MATCHED,
    RECOMU_LLLLss0_MATCHED,RECOMU_LLLLss1_MATCHED,RECOMU_LLLLss2_MATCHED,
    RECOMU_LLLl0_MATCHED,RECOMU_LLLl1_MATCHED,RECOMU_LLLL_MATCHED;
  ArrayVector<int> RECOMU_numberOfMatches,RECOMU_numberOfMatchedStations;
  ArrayVector<int> RECOMU_mubesttrkType;
  
  ArrayVector<float> 
    RECOMU_mutrkPT,RECOMU_mutrkPTError,
    RECOMU_mutrkDxy,RECOMU_mutrkDxyError,RECOMU_mutrkDxyB,
    RECOMU_mutrkDz,RECOMU_mutrkDzError,RECOMU_mutrkDzB,
    RECOMU_mutrkChi2PerNdof,
    RECOMU_mutrktrackerLayersWithMeasurement,
    RECOMU_muInnertrkDxy,RECOMU_muInnertrkDxyError,RECOMU_muInnertrkDxyB,
    RECOMU_muInnertrkDz,RECOMU_muInnertrkDzError,RECOMU_muInnertrkDzB,
    RECOMU_mubesttrkDxy,RECOMU_mubesttrkDxyError,RECOMU_mubesttrkDxyB,
    RECOMU_mubesttrkDz,RECOMU_mubesttrkDzError,RECOMU_mubesttrkDzB, RECOMU_mubesttrkPTError, RECOMU_Rochester_Error,
    RECOMU_muInnertrkChi2PerNdof,
    RECOMU_muInnertrktrackerLayersWithMeasurement,RECOMU_muInnertrkPT,RECOMU_muInnertrkPTError,
    RECOMU_muInnertrkCharge,RECOMU_muInnertrkNHits,RECOMU_muInnertrkNPixHits,RECOMU_muInnertrkNStripHits,
    RECOMU_mutrkCharge,RECOMU_mutrkNHits,RECOMU_mutrkNPixHits,RECOMU_mutrkNStripHits,RECOMU_mutrkNMuonHits;
  ArrayVector<int > RECOMU_trkmuArbitration,RECOMU_trkmu2DCompatibilityLoose,RECOMU_trkmu2DCompatibilityTight;
  ArrayVector<int > RECOMU_trkmuOneStationLoose,RECOMU_trkmuOneStationTight;
  ArrayVector<int > RECOMU_trkmuLastStationLoose,RECOMU_trkmuLastStationTight;
  ArrayVector<int > RECOMU_trkmuLastStationAngLoose,RECOMU_trkmuLastStationAngTight;
  ArrayVector<int > RECOMU_trkmuOneStationAngLoose,RECOMU_trkmuOneStationAngTight;
  ArrayVector<int > RECOMU_trkmuLastStationOptimizedLowPtLoose,RECOMU_trkmuLastStationOptimizedLowPtTight;
  
   // Photons
  ArrayVector<float> RECOPHOT_PT,RECOPHOT_ETA,RECOPHOT_PHI,RECOPHOT_THETA,RECOPHOT_TLE_ParentSC_X,RECOPHOT_TLE_ParentSC_Y,RECOPHOT_TLE_ParentSC_Z;
  ArrayVector<float> RECOPFPHOT_PT,RECOPFPHOT_PTError,RECOPFPHOT_ETA,RECOPFPHOT_PHI,RECOPFPHOT_THETA;
  ArrayVector<double> RECOPFPHOT_PFchAllPart,RECOPFPHOT_PFchHad,RECOPFPHOT_PFneuHad,RECOPFPHOT_PFphoton,
    RECOPFPHOT_PFPUchAllPart,RECOPFPHOT_PFX_rho, RECOPFPHOT_PT_uncorr;

  //Reham 
  
  ArrayVector<float> RECOPFPHOT_ecalEnergyErrPostCorr,RECOPFPHOT_energyScaleValue,RECOPFPHOT_energySigmaValue, RECOPFPHOT_energyScaleUp, RECOPFPHOT_energyScaleDown, RECOPFPHOT_energyScaleStatUp, RECOPFPHOT_energyScaleStatDown, RECOPFPHOT_energyScaleSystUp, RECOPFPHOT_energyScaleSystDown, RECOPFPHOT_energyScaleGainUp, RECOPFPHOT_energyScaleGainDown,RECOPFPHOT_energyScaleEtUp, RECOPFPHOT_energyScaleEtDown, RECOPFPHOT_energySigmaUp, RECOPFPHOT_energySigmaDown, RECOPFPHOT_energySigmaPhiUp, RECOPFPHOT_energySigmaPhiDown, RECOPFPHOT_energySigmaRhoUp, RECOPFPHOT_energySigmaRhoDown,RECOPFPHOT_ecalEnergyPreCorr,RECOPFPHOT_ecalEnergyErrPreCorr;
  

  // Vertexing
  ArrayVector<double> ftsigma,ftsigmalag,ftsigmaMMMM,ftsigmalagMMMM,ftsigmaEEEE,ftsigmalagEEEE;
  ArrayVector<double> gdX,gdY,gdZ,gdXMMMM,gdYMMMM,gdZMMMM,gdXEEEE,gdYEEEE,gdZEEEE;
  ArrayVector<double> gdlagX,gdlagY,gdlagZ,gdlagProb,gdlagNdof,
    gdlagXMMMM,gdlagYMMMM,gdlagZMMMM,gdlagProbMMMM,gdlagNdofMMMM,
    gdlagXEEEE,gdlagYEEEE,gdlagZEEEE,gdlagProbEEEE,gdlagNdofEEEE;

  //ConstraintFit 4l
  ArrayVector<double>  StdFitVertexX, StdFitVertexY, StdFitVertexZ, StdFitVertexChi2r, StdFitVertexProb;
  ArrayVector<double>  KinFitVertexX, KinFitVertexY, KinFitVertexZ, KinFitVertexChi2r, KinFitVertexProb;
  ArrayVector<double>  StdFitVertexXMMMM, StdFitVertexYMMMM, StdFitVertexZMMMM, StdFitVertexChi2rMMMM, StdFitVertexProbMMMM;
  ArrayVector<double>  KinFitVertexXMMMM, KinFitVertexYMMMM, KinFitVertexZMMMM, KinFitVertexChi2rMMMM, KinFitVertexProbMMMM;
  ArrayVector<double>  StdFitVertexXEEEE, StdFitVertexYEEEE, StdFitVertexZEEEE, StdFitVertexChi2rEEEE, StdFitVertexProbEEEE;
  ArrayVector<double>  KinFitVertexXEEEE, KinFitVertexYEEEE, KinFitVertexZEEEE, KinFitVertexChi2rEEEE, KinFitVertexProbEEEE;

  ArrayVector<ArrayVector<float> > StdFitVertexTrack_PT,StdFitVertexTrack_ETA,StdFitVertexTrack_PHI,
    StdFitVertexTrackMMMM_PT,StdFitVertexTrackMMMM_ETA,StdFitVertexTrackMMMM_PHI,
    StdFitVertexTrackEEEE_PT,StdFitVertexTrackEEEE_ETA,StdFitVertexTrackEEEE_PHI;

   //ConstraintFit 2l
  ArrayVector<double>  StdFitVertexChi2rDiLep, StdFitVertexProbDiLep;

  //ConstraintFit 3l
  ArrayVector<double>  StdFitVertexChi2rMMM, StdFitVertexProbMMM;
  ArrayVector<double>  StdFitVertexChi2rMME, StdFitVertexProbMME;
  ArrayVector<double>  StdFitVertexChi2rEEE, StdFitVertexProbEEE;
  ArrayVector<double>  StdFitVertexChi2rMEE, StdFitVertexProbMEE;
   
  
  //Muons Matching
  ArrayVector<int > RECOMU_MatchingMCTruth;
  ArrayVector<float> RECOMU_MatchingMCpT;
  ArrayVector<float> RECOMU_MatchingMCEta;
  ArrayVector<float> RECOMU_MatchingMCPhi;
  
  //Electrons:
  ArrayVector<int > RECOELE_MatchingMCTruth;
  ArrayVector<float> RECOELE_MatchingMCpT;
  ArrayVector<float> RECOELE_MatchingMCEta;
  ArrayVector<float> RECOELE_MatchingMCPhi;
  //Gamma:
  ArrayVector<int > RECOPHOT_MatchingMCTruth;
  ArrayVector<float> RECOPHOT_MatchingMCpT;
  ArrayVector<float> RECOPHOT_MatchingMCEta;
  ArrayVector<float> RECOPHOT_MatchingMCPhi;

  //zToMuMu:
  ArrayVector<int > RECOzMuMu_MatchingMCTruth;
  ArrayVector<float> RECOzMuMu_MatchingMCpT;
  ArrayVector<float> RECOzMuMu_MatchingMCmass;
  ArrayVector<float> RECOzMuMu_MatchingMCEta;
  ArrayVector<float> RECOzMuMu_MatchingMCPhi;

  //zToEE:
  ArrayVector<int > RECOzEE_MatchingMCTruth;
  ArrayVector<float> RECOzEE_MatchingMCpT;
  ArrayVector<float> RECOzEE_MatchingMCmass;
  ArrayVector<float> RECOzEE_MatchingMCEta;
  ArrayVector<float> RECOzEE_MatchingMCPhi;

  //HtoZtoEEEE:
  ArrayVector<int > RECOHzzEEEE_MatchingMCTruth;
  ArrayVector<float> RECOHzzEEEE_MatchingMCpT;
  ArrayVector<float> RECOHzzEEEE_MatchingMCmass;
  ArrayVector<float> RECOHzzEEEE_MatchingMCEta;
  ArrayVector<float> RECOHzzEEEE_MatchingMCPhi;

  //HtoZtoMMMM:
  ArrayVector<int > RECOHzzMMMM_MatchingMCTruth;
  ArrayVector<float> RECOHzzMMMM_MatchingMCpT;
  ArrayVector<float> RECOHzzMMMM_MatchingMCmass;
  ArrayVector<float> RECOHzzMMMM_MatchingMCEta;
  ArrayVector<float> RECOHzzMMMM_MatchingMCPhi;

  //HtoZtoEEMM:
  ArrayVector<int > RECOHzzEEMM_MatchingMCTruth;
  ArrayVector<float> RECOHzzEEMM_MatchingMCpT;
  ArrayVector<float> RECOHzzEEMM_MatchingMCmass;
  ArrayVector<float> RECOHzzEEMM_MatchingMCEta;
  ArrayVector<float> RECOHzzEEMM_MatchingMCPhi;
  
 

  // RECO counters
  int RECO_NMU, RECO_NELE, RECO_NTRACK, RECO_NPHOT, RECO_NPFPHOT,RECO_NJET, RECO_NVTX;
  ArrayVector<float> RECO_TRACK_PT, RECO_TRACK_ETA, RECO_TRACK_PHI,
    RECO_TRACK_CHI2,RECO_TRACK_CHI2RED,RECO_TRACK_CHI2PROB, 
    RECO_TRACK_DXY,RECO_TRACK_DXYERR, 
    RECO_TRACK_DZ,RECO_TRACK_DZERR;
  ArrayVector<int> RECO_TRACK_NHITS;
  
  // Primary Vertices
  ArrayVector<float> RECO_VERTEX_x, RECO_VERTEX_y, RECO_VERTEX_z,RECO_VERTEX_ndof,RECO_VERTEX_chi2,RECO_VERTEXPROB;
  ArrayVector<ArrayVector<float > > RECO_VERTEX_TRACK_PT;
  ArrayVector<int > RECO_VERTEX_isValid;
  ArrayVector<int> RECO_VERTEX_ntracks;
  
  // RECO JETS
  int RECO_PFJET_N;
  ArrayVector<int>  RECO_PFJET_CHARGE,RECO_PFJET_PUID, RECO_PFJET_PUID_loose, RECO_PFJET_PUID_medium;
  ArrayVector<float> RECO_PFJET_ET, RECO_PFJET_PT, RECO_PFJET_ETA, RECO_PFJET_PHI,RECO_PFJET_PUID_MVA,RECO_PFJET_QG_Likelihood,RECO_PFJET_QG_axis2,RECO_PFJET_QG_ptd,RECO_PFJET_QG_mult;
  double RHO,RHO_ele,RHO_mu;

 //@
  bool isData;
  edm::EDGetTokenT<edm::View<pat::Muon>> slimmedMuonsTag_;
  edm::EDGetTokenT<LHEEventProduct> LHE_;
  ArrayVector<int > LHE_PARTON_CLEAR;
  int LHE_PARTON_N;
  ArrayVector<int> LHE_PARTON_PDGID;
  ArrayVector<float> LHE_PARTON_PT, LHE_PARTON_ETA, LHE_PARTON_PHI, LHE_PARTON_E;
  ArrayVector<int> RECO_PFJET_CHARGED_HADRON_MULTIPLICITY, RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY, RECO_PFJET_PHOTON_MULTIPLICITY, RECO_PFJET_ELECTRON_MULTIPLICITY;
  ArrayVector<ArrayVector<int > > RECO_PFJET_COMPONENT_PDGID;
  ArrayVector<int> RECO_PFJET_MUON_MULTIPLICITY, RECO_PFJET_HF_HADRON_MULTIPLICTY, RECO_PFJET_HF_EM_MULTIPLICITY, RECO_PFJET_CHARGED_MULTIPLICITY, RECO_PFJET_NEUTRAL_MULTIPLICITY, RECO_PFJET_NCOMPONENTS;
  ArrayVector<float> RECO_PFJET_PT_UncUp, RECO_PFJET_PT_UncDn, RECO_PFJET_AREA, RECO_PFJET_CHARGED_HADRON_ENERGY, RECO_PFJET_NEUTRAL_HADRON_ENERGY, RECO_PFJET_PHOTON_ENERGY, RECO_PFJET_ELECTRON_ENERGY, RECO_PFJET_PTD;
  ArrayVector<float> RECO_PFJET_MUON_ENERGY, RECO_PFJET_HF_HADRON_ENERGY, RECO_PFJET_HF_EM_ENERGY, RECO_PFJET_CHARGED_EM_ENERGY, RECO_PFJET_CHARGED_MU_ENERGY;
  ArrayVector<float> RECO_PFJET_NEUTRAL_EM_ENERGY;
  ArrayVector<ArrayVector<float > > RECO_PFJET_COMPONENT_PT, RECO_PFJET_COMPONENT_ETA, RECO_PFJET_COMPONENT_PHI, RECO_PFJET_COMPONENT_E, RECO_PFJET_COMPONENT_TRANSVERSE_MASS, RECO_PFJET_COMPONENT_CHARGE;
  ArrayVector<ArrayVector<float> > RECO_PFJET_COMPONENT_XVERTEX, RECO_PFJET_COMPONENT_YVERTEX, RECO_PFJET_COMPONENT_ZVERTEX, RECO_PFJET_COMPONENT_VERTEX_CHI2;

  // GenJET
  ArrayVector<float> MC_GENJET_PT, MC_GENJET_ETA, MC_GENJET_PHI;

  // RECO MET
  float genmet;
  float calomet;
    //calometopt,calometoptnohf,calometoptnohfho,calometoptho,calometnohf,calometnohfho,calometho;       
  float pfmet,pfmet_x,pfmet_y,pfmet_phi,pfmet_theta,pfmet_uncorr,pfmet_x_uncorr,pfmet_y_uncorr,pfmet_phi_uncorr,pfmet_theta_uncorr;
  float pfmet_JetEnUp, pfmet_JetEnDn, pfmet_ElectronEnUp, pfmet_ElectronEnDn, pfmet_MuonEnUp, pfmet_MuonEnDn;
  float pfmet_JetResUp, pfmet_JetResDn, pfmet_UnclusteredEnUp, pfmet_UnclusteredEnDn, pfmet_PhotonEnUp, pfmet_PhotonEnDn,pfmet_TauEnUp,pfmet_TauEnDown;
  int PassecalBadCalibFilterUpdated, filterbadChCandidate, filterbadPFMuon;

    //htmetic5,htmetkt4,htmetkt6,htmetsc5,htmetsc7;        
    float tcmet;
    //jescormetic5,jescormetkt4,jescormetkt6,jescormetsc5,jescormetsc7;    
  float cormetmuons; 

  //Read MET filter decisions

    std::string GoodVtxNoiseFilter_Selector_;
    std::string GlobalSuperTightHalo2016NoiseFilter_Selector_;
    std::string HBHENoiseFilter_Selector_;
    std::string HBHENoiseIsoFilter_Selector_;
    std::string EcalDeadCellTriggerPrimitiveNoiseFilter_Selector_;
    std::string BadPFMuonFilter_Selector_;
    std::string BadChargedCandidateFilter_Selector_;
    std::string EEBadScNoiseFilter_Selector_;                         
    std::string EcalBadCalibFilter_Selector_; 
  
    int  passFilterGoodVtxNoise, passFilterGlobalSuperTightHalo2016NoiseFilter, passFilterHBHENoise, passFilterHBHENoiseIso, passFilterEcalDeadCellTriggerPrimitiveNoise, passFilterBadPFMuon, passFilterBadChargedCandidate, passFilterEEBadScNoise, passFilterEcalBadCalib; 

  // Beam Spot
  double BeamSpot_X,BeamSpot_Y,BeamSpot_Z;	
  BeamSpot bs;
  
 
  
  ArrayVector<float> tCHighEff_BTagJet_PT,
    tCHighPur_BTagJet_PT,
    cSV_BTagJet_PT;
  ArrayVector<float> tCHighEff_BTagJet_ETA,
    tCHighPur_BTagJet_ETA,
    cSV_BTagJet_ETA;
  ArrayVector<float> tCHighEff_BTagJet_PHI,
    tCHighPur_BTagJet_PHI,
    cSV_BTagJet_PHI;
  ArrayVector<float> tCHighEff_BTagJet_DISCR,
    tCHighPur_BTagJet_DISCR,
    cSV_BTagJet_DISCR;
  ArrayVector<float> cSV_BTagJet_ET;

  ArrayVector<float> ConvMapDist,ConvMapDcot;

  // MVA Ring Isolation
  //MuonMVAEstimator *fMuonIsoMVA;

  // Magnetic Field
  edm::ESHandle<MagneticField> magfield_;
  std::unique_ptr<TRandom3> rgen_;
};

#endif
