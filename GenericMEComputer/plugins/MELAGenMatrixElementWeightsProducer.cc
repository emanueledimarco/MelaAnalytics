#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "JHUGenMELA/MELA/interface/Mela.h"
#include "JHUGenMELA/MELA/interface/Mela.h"
#include "MelaAnalytics/EventContainer/interface/MELAEvent.h"
#include "MelaAnalytics/GenericMEComputer/interface/MELAOptionParser.h"
#include "MelaAnalytics/GenericMEComputer/interface/MELAHypothesis.h"
#include "MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h"
#include "CommonLHETools/LHEHandler/interface/LHEHandler.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>
#include <regex>

namespace {

class MELAGenMatrixElementWeightsProducer : public edm::one::EDProducer<edm::one::WatchRuns, edm::one::SharedResources> {
    public:
        MELAGenMatrixElementWeightsProducer( edm::ParameterSet const & iConfig )  ;
        ~MELAGenMatrixElementWeightsProducer() override {} // nothing to do

        void beginRun(const edm::Run& iRun, const edm::EventSetup&) override ;
        void endRun(const edm::Run&, const edm::EventSetup&) override ;
        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override ;

    protected:
        const std::vector<edm::InputTag> lheLabel_;
        const std::vector<edm::EDGetTokenT<LHEEventProduct>> lheTag_;
        const std::vector<edm::EDGetTokenT<LHERunInfoProduct>> lheRunTag_;

        const std::string name_;

        const int year_;
        const float sqrts_, mH_;
        MELAEvent::CandidateVVMode candVVmode_;
        int decayVVmode_;
        enum MelaSetup {
            None=0, Decay=1, VBF_NLO=2, VBF_LO=3, ZH_NLO=4, WH_NLO=5, ZH_LO=6, WH_LO=7
        } melaSetup_;
        bool normalize_;

        bool firstRun_;
        std::unique_ptr<LHEHandler> lheHandler_;
        std::unique_ptr<Mela> mela_;

        class MEBranch {
            protected:
                std::unique_ptr<MELAOptionParser> opt_;
                const std::string stropt_, name_;
                double value_;
                bool dummy_;

            public:
                MEBranch(const std::string & line) ; 
                const std::string & name() const { return name_; }
                const std::string & doc() const { return stropt_; }
                double value() const { 
                    return dummy_ ? -1 : value_; 
                }
                double normValue(double norm) const { 
                    return norm ? ( dummy_ ? 1 : value_/norm  ) : -1;
                }
               
                void reset() { value_ = -1; } 
                void calc(Mela &mela) {
                    if (dummy_) return;
                    MELAHypothesis hypot(&mela, opt_.get());
                    hypot.computeP();
                    value_ = hypot.getVal(MELAHypothesis::UseME);
                }

                void setProduction(TVar::Production prod) {
                    opt_->prod = prod;
                }
                void optimizeForNormValue(const MEBranch & other) ;

        };
        std::vector<MEBranch> meBranches_;

    private:
        bool runNLOVHApprox() ;
        bool runNLOVBFApprox();
        bool runBestLOAssociatedV();
        bool runBestLOAssociatedVBF();
        bool runDecay();
}; // class


MELAGenMatrixElementWeightsProducer::MELAGenMatrixElementWeightsProducer( edm::ParameterSet const & iConfig ) :
    lheLabel_(iConfig.getParameter<std::vector<edm::InputTag>>("lheInfo")),
    lheTag_(edm::vector_transform(lheLabel_, [this](const edm::InputTag & tag) { return mayConsume<LHEEventProduct>(tag); })),
    lheRunTag_(edm::vector_transform(lheLabel_, [this](const edm::InputTag & tag) { return mayConsume<LHERunInfoProduct, edm::InRun>(tag); })),
    name_(iConfig.getParameter<std::string>("name")),
    year_(iConfig.getParameter<int>("year")),
    sqrts_(iConfig.getParameter<double>("sqrts")),
    mH_(iConfig.getParameter<double>("mH")),
    normalize_(iConfig.getParameter<bool>("normalize")),
    firstRun_(true)
{
    LHEHandler::set_maxlines_print_header(year_==2016 ? 1000 : -1);
    this->usesResource("MELA");
    mela_.reset(new Mela(sqrts_, mH_, TVar::ERROR));
    //mela_->setVerbosity(TVar::DEBUG);

    const std::string mode = iConfig.getParameter<std::string>("mode");

    candVVmode_ = MELAEvent::UndecayedMode; decayVVmode_ = -1;
    std::string commonString = "Process:SelfDefine_spin0 ";
    if (mode == "ZH_NLO") { 
        commonString += "MatrixElement:JHUGen Production:Lep_ZH "; melaSetup_ = ZH_NLO;  
    } else if (mode ==  "WH_NLO") { 
        commonString += "MatrixElement:JHUGen Production:Lep_WH "; melaSetup_ = WH_NLO;  
    } else if (mode == "ZH_LO") { 
        commonString += "MatrixElement:JHUGen Production:Lep_ZH "; melaSetup_ = ZH_LO; 
    } else if (mode ==  "WH_LO") { 
        commonString += "MatrixElement:JHUGen Production:Lep_WH "; melaSetup_ = WH_LO; 
    } else if (mode ==  "VBF_NLO") { 
        commonString += "MatrixElement:JHUGen Production:JJVBF "; melaSetup_ = VBF_NLO;  
    } else if (mode ==  "VBF_LO") { 
        commonString += "MatrixElement:JHUGen Production:JJVBF "; melaSetup_ = VBF_LO; 
    } else if (mode ==  "Decay_ZZ4l") { 
        candVVmode_ = MELAEvent::ZZMode; decayVVmode_ = 0;
        commonString += "MatrixElement:JHUGen Production:ZZINDEPENDENT "; melaSetup_ = Decay; 
    } else if (mode ==  "Decay_WW2l2v") { 
        candVVmode_ = MELAEvent::WWMode; decayVVmode_ = 0;
        commonString += "MatrixElement:JHUGen Production:ZZINDEPENDENT "; melaSetup_ = Decay; 
    } else if (mode ==  "Decay_gammagamma") { 
        candVVmode_ = MELAEvent::GammaGammaMode; decayVVmode_ = 0;
        commonString += "MatrixElement:JHUGen Production:ZZINDEPENDENT "; melaSetup_ = Decay; 
    } else {
        throw cms::Exception("Configuration", "mode '"+mode+"' not yet implemented\n");
    }
    commonString += "hmass:125 hwidth:0.00407 ";

    lheHandler_.reset(new LHEHandler(candVVmode_, decayVVmode_, LHEHandler::doHiggsKinematics, year_, LHEHandler::keepDefaultPDF, LHEHandler::keepDefaultQCDOrder));

    for (const std::string & element : iConfig.getParameter<std::vector<std::string>>("matrixElements")) {
        meBranches_.emplace_back(commonString+element);
    }
    if (normalize_) {
        const MEBranch & normTerm = meBranches_.front();
        for (auto it = meBranches_.begin()+1, ed = meBranches_.end(); it != ed; ++it) {
            it->optimizeForNormValue(normTerm);
        }

    }

    produces<std::vector<float> >("values");
    produces<std::vector<std::string> >("labels");
}


void MELAGenMatrixElementWeightsProducer::beginRun(const edm::Run& iRun, const edm::EventSetup&) {

}

void MELAGenMatrixElementWeightsProducer::endRun(const edm::Run& iRun, const edm::EventSetup&) {
    if (firstRun_) {
        edm::Handle<LHERunInfoProduct> lheInfo;
        for (const auto & lheLabel: lheLabel_) {
            iRun.getByLabel(lheLabel, lheInfo);
            if (lheInfo.isValid()) {
                break;
            }
        }
        lheHandler_->setHeaderFromRunInfo(&lheInfo);
        firstRun_ = false;
    }     // nothing to do
}

void MELAGenMatrixElementWeightsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    std::unique_ptr<std::vector<float> > outValues(new std::vector<float>);
    std::unique_ptr<std::vector<std::string> > outLabels(new std::vector<std::string>);

    edm::Handle<LHEEventProduct> lheInfo;
    for (const auto & lheTag: lheTag_) {
        iEvent.getByToken(lheTag, lheInfo);
        if (lheInfo.isValid()) {
            break;
        }
    }

    lheHandler_->setHandle(&lheInfo);
    lheHandler_->extract();
    MELACandidate* cand = lheHandler_->getBestCandidate();
    if (!cand) { edm::LogWarning("Mela") << "Didn't find best candidate. "; return; }
    else mela_->setCurrentCandidate(cand);

    bool runOk = false;
    switch (melaSetup_) {
        case ZH_LO: runOk = runBestLOAssociatedV(); break;
        case WH_LO: runOk = runBestLOAssociatedV(); break;
        case ZH_NLO: runOk = runNLOVHApprox(); break;
        case WH_NLO: runOk = runNLOVHApprox(); break;
        case VBF_LO: runOk = runBestLOAssociatedVBF(); break;
        case VBF_NLO: runOk = runNLOVBFApprox(); break;
        case Decay: runOk = runDecay(); break;
        case None: 
            runOk = (cand != nullptr);
            for (auto & me : meBranches_) me.calc(*mela_);
            break;
        default: 
            break;
    }
    if (!runOk) { 
        edm::LogWarning("MelaSetup") << "Error setting up MELA candidate for setup " << melaSetup_; 
        for (auto & me : meBranches_) me.reset();
    }

    if (normalize_) {
        double norm = meBranches_.front().value();
        if (norm == 0) edm::LogWarning("MelaSetup") << "Warning: ZERO matrix element for the normalization hypothesis.";
        for (auto it = meBranches_.begin()+1, ed = meBranches_.end(); it != ed; ++it) {
          outLabels->push_back(it->name());
          outValues->push_back(it->normValue(norm));
          //std::cout << "Putting inthe event MELA NORM weight: " << it->name() << " with value  = " << it->normValue(norm) << std::endl;
        }
    } else {
        for (const auto & me : meBranches_) {
          outLabels->push_back(me.name());
          outValues->push_back(me.value());
          // std::cout << "Putting in the event MELA weight: " << me.name() << " with value  = " << me.value() << std::endl;
        }
    }
    
    iEvent.put(std::move(outLabels), "labels");
    iEvent.put(std::move(outValues), "values");
    mela_->resetInputEvent();
}


MELAGenMatrixElementWeightsProducer::MEBranch::MEBranch(const std::string & line) :
        opt_(new MELAOptionParser(line)), 
        stropt_(line), 
        name_(opt_->getName()), 
        dummy_(false) 
{
    if (opt_->prod == TVar::Lep_WH || opt_->prod == TVar::Had_WH) {
        // MELA uses ghz1, ghz2, ghz4, ghz2_prime2 instead of ghw1, ghw2, ghw4, ghw2_prime2
        SpinZeroCouplings & coupl_H = opt_->coupl_H;
        for (int j = 0; j <= 1; ++j) {
            coupl_H.Hzzcoupl[gHIGGS_VV_1][j] = coupl_H.Hwwcoupl[gHIGGS_VV_1][j];
            coupl_H.Hzzcoupl[gHIGGS_VV_2][j] = coupl_H.Hwwcoupl[gHIGGS_VV_2][j];
            coupl_H.Hzzcoupl[gHIGGS_VV_4][j] = coupl_H.Hwwcoupl[gHIGGS_VV_4][j];
            coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME2][j] = coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME2][j];
        }
    } 
}

void MELAGenMatrixElementWeightsProducer::MEBranch::optimizeForNormValue(const MELAGenMatrixElementWeightsProducer::MEBranch & other) {
    SpinZeroCouplings & me = opt_->coupl_H;
    const SpinZeroCouplings & ref = other.opt_->coupl_H;
    if (opt_->prod == TVar::Lep_WH || opt_->prod == TVar::Had_WH) {
        // MELA uses ghz1, ghz2, ghz4, ghz2_prime2 instead of ghw1, ghw2, ghw4, ghw2_prime2
        dummy_ = true;
        for (int j = 0; j <= 1; ++j) {
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_1][j] == ref.Hzzcoupl[gHIGGS_VV_1][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_2][j] == ref.Hzzcoupl[gHIGGS_VV_2][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_4][j] == ref.Hzzcoupl[gHIGGS_VV_4][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_1_PRIME2][j] == ref.Hzzcoupl[gHIGGS_VV_1_PRIME2][j]);
        }
    } else if (opt_->prod == TVar::Lep_ZH || opt_->prod == TVar::Had_ZH) {
        dummy_ = true;
        for (int j = 0; j <= 1; ++j) {
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_1][j] == ref.Hzzcoupl[gHIGGS_VV_1][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_2][j] == ref.Hzzcoupl[gHIGGS_VV_2][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_4][j] == ref.Hzzcoupl[gHIGGS_VV_4][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_ZA_2][j] == ref.Hzzcoupl[gHIGGS_ZA_2][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_ZA_4][j] == ref.Hzzcoupl[gHIGGS_ZA_4][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_VV_1_PRIME2][j] == ref.Hzzcoupl[gHIGGS_VV_1_PRIME2][j]);
            dummy_ &= (me.Hzzcoupl[gHIGGS_ZA_1_PRIME2][j] == ref.Hzzcoupl[gHIGGS_ZA_1_PRIME2][j]);
        }
    }
    //if (dummy_) edm::LogInfo("MEOpt") << "Optimized dummy ME branch " << name_ << " with couplings " << stropt_ << "\n";
}



// All these methods below to setup the initial state particles are taken from
//     https://indico.cern.ch/event/873834/contributions/3697609/attachments/1968519/3273954/ExampleMELAHelpers.cc
//
bool MELAGenMatrixElementWeightsProducer::runDecay() {
    MELACandidate* melaCand = mela_->getCurrentCandidate();
    if (!melaCand) return false;

    for (int imot=0; imot<melaCand->getNMothers(); imot++) {
        melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) {
        melaCand->getAssociatedJet(ijet)->id = 0;
    }

    for (auto & me : meBranches_) {
        me.calc(*mela_);
    }

    return true;
}



bool MELAGenMatrixElementWeightsProducer::runNLOVHApprox() {
    MELACandidate* melaCand = mela_->getCurrentCandidate();
    if (!melaCand) return false;

    TVar::Production candScheme = (melaSetup_ == ZH_NLO || melaSetup_ == ZH_LO) ? TVar::Had_ZH : TVar::Had_WH;

    MELACandidateRecaster recaster(candScheme);
    MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(melaCand, candScheme);
    if (bestAV) {
        MELACandidate* candModified = nullptr;
        recaster.copyCandidate(melaCand, candModified);
        recaster.deduceLOVHTopology(candModified);
        mela_->setCurrentCandidate(candModified);

        int V_dauid = bestAV->getDaughter(0)->id;
        bool hasALepV = (PDGHelpers::isALepton(V_dauid) || PDGHelpers::isANeutrino(V_dauid));
        TVar::Production prod = (melaSetup_ == ZH_NLO || melaSetup_ == ZH_LO) ? 
                            (hasALepV ? TVar::Lep_ZH : TVar::Had_ZH) :
                            (hasALepV ? TVar::Lep_WH : TVar::Had_WH);
        for (auto & me : meBranches_) {
            me.setProduction(prod);
            me.calc(*mela_);
        }

        return true;
    } else {
        return false; 
    }
}

bool MELAGenMatrixElementWeightsProducer::runBestLOAssociatedV() {
    MELACandidate* melaCand = mela_->getCurrentCandidate();
    if (!melaCand) return false;

    int    VabsPdgId = (melaSetup_ == ZH_LO ? 23 : 24);
    double Vmass     = (melaSetup_ == ZH_LO ? PDGHelpers::Zmass : PDGHelpers::Wmass);

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    for (int imot = 0; imot < melaCand->getNMothers(); imot++) {
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet = 0; ijet < melaCand->getNAssociatedJets(); ijet++) {
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELAParticle*> associatedVs; // Vector of Ws or Zs to loop over
    for (MELAParticle* Vtmp : melaCand->getAssociatedSortedVs()) {
      if (Vtmp && (std::abs(Vtmp->id) == VabsPdgId) && Vtmp->getNDaughters()>=1) {
        bool passSelection = true;
        for (MELAParticle* dauVtmp : Vtmp->getDaughters()) passSelection &= dauVtmp->passSelection;
        if (!passSelection) continue;

        associatedVs.push_back(Vtmp);
      }
    }

    // Give precedence to leptonic V decays
    bool hasALepV = false;
    for (MELAParticle* Vtmp : associatedVs) {
        const int& Vtmp_dauid = Vtmp->getDaughter(0)->id;
        if (PDGHelpers::isALepton(Vtmp_dauid) || PDGHelpers::isANeutrino(Vtmp_dauid)) {
            hasALepV=true;
            break;
        }
    }

    MELAParticle* bestVbyMass = nullptr;
    double bestVMassDiff = -1;
    for (MELAParticle* Vtmp:associatedVs) {
        // switch off all hadronic V's (we will then turn on the best matching one only)
        if (PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)) {
            for (MELAParticle* dauVtmp : Vtmp->getDaughters()) dauVtmp->setSelected(false);
        }
        if (!hasALepV && (!bestVbyMass || std::abs(Vtmp->m()-Vmass) < bestVMassDiff)) {
            bestVMassDiff = std::abs(Vtmp->m()-Vmass);
            bestVbyMass = Vtmp;
        }
    }
    if (bestVbyMass) {
        for (MELAParticle* dauVtmp : bestVbyMass->getDaughters()) dauVtmp->setSelected(true);
    }

    if (hasALepV || bestVbyMass != nullptr) {
        TVar::Production prod = (melaSetup_ == ZH_NLO || melaSetup_ == ZH_LO) ? 
                            (hasALepV ? TVar::Lep_ZH : TVar::Had_ZH) :
                            (hasALepV ? TVar::Lep_WH : TVar::Had_WH);
        for (auto & me : meBranches_) {
            me.setProduction(prod);
            me.calc(*mela_);
        }
        return true;
    } else {
        return false; 
    }
}

bool MELAGenMatrixElementWeightsProducer::runNLOVBFApprox() {
    MELACandidate* melaCand = mela_->getCurrentCandidate();
    if (!melaCand) return false;

    // Need one recaster for VBF
    MELACandidateRecaster recaster(TVar::JJVBF);
    MELACandidate* candModified = nullptr;
    recaster.copyCandidate(melaCand, candModified);
    recaster.reduceJJtoQuarks(candModified);
    mela_->setCurrentCandidate(candModified);

    for (auto & me : meBranches_) {
        me.calc(*mela_);
    }

    return true;
}

bool MELAGenMatrixElementWeightsProducer::runBestLOAssociatedVBF() {
    MELACandidate* melaCand = mela_->getCurrentCandidate();
    if (!melaCand) return false;

    for (int imot=0; imot<melaCand->getNMothers(); imot++) {
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) {
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    for (auto & me : meBranches_) {
        me.calc(*mela_);
    }

    return true;
}
 
} // namespace



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MELAGenMatrixElementWeightsProducer);

