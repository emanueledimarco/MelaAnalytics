#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <iostream>

namespace {

class HiggsDecayLevelVariablesProducer : public edm::stream::EDProducer<> {
    public:
        HiggsDecayLevelVariablesProducer( edm::ParameterSet const & iConfig )  ;
        ~HiggsDecayLevelVariablesProducer() override {} // nothing to do

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override ;

        void collectDaus(const reco::GenParticle & mom, std::vector<const reco::GenParticle *> & daus) const {
            //std::cout << "collect daughters from genParticle pdgId " << mom.pdgId() << ", status " << mom.status() << ", lastCopy " << mom.isLastCopy() << ", numberOfDaughters " << mom.numberOfDaughters() << ". Already gave " <<  daus.size() << " daughters in list" << std::endl;
            if (std::abs(mom.pdgId()) != 23 && std::abs(mom.pdgId()) != 24 && std::abs(mom.pdgId()) != 25)  {
                daus.push_back(&mom);
                //std::cout << " - adding this to the list of daughters. new list size = " << daus.size()  << std::endl;
            } else if (!mom.isLastCopy()) {
                for (unsigned int i = 0; i < mom.numberOfDaughters(); ++i) {
                    const auto & di = mom.daughterRef(i);
                    //std::cout << " - examining daughter " << i << ", pdgId " << (di.isAvailable() ? di->pdgId() : 0) << ", status " << (di.isAvailable() ? di->status() : 0) << std::endl;
                    if (di.isAvailable() && di->pdgId() == mom.pdgId()) {
                        //std::cout << "     - recursive call on this " << daus.size()  << std::endl;
                        collectDaus(*di, daus);
                    }
                }
            } else {
                for (unsigned int i = 0; i < mom.numberOfDaughters(); ++i) {
                    const auto & di = mom.daughterRef(i);
                    //std::cout << " - recursive call on daughter " << i << ", pdgId " << (di.isAvailable() ? di->pdgId() : 0) << ", status " << (di.isAvailable() ? di->status() : 0) << std::endl;
                    if (di.isAvailable()) collectDaus(*di, daus);
                }
            }
        }
    protected:
        const edm::EDGetTokenT<reco::GenParticleCollection> genTag_;
}; // class


HiggsDecayLevelVariablesProducer::HiggsDecayLevelVariablesProducer( edm::ParameterSet const & iConfig ) :
    genTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles")))
{
    produces<nanoaod::FlatTable>();
}


void HiggsDecayLevelVariablesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    auto out = std::make_unique<nanoaod::FlatTable>(1, "GenHiggsDecay", true);
    out->setDoc("Higgs Gen-level decay information");

    edm::Handle<reco::GenParticleCollection> genInfo;
    iEvent.getByToken(genTag_, genInfo);

    bool found = false;
    for (const auto & h : *genInfo) {
        if (h.pdgId() == 25 && h.isLastCopy()) {
            found = true;
            out->addColumnValue<int>("ndauRaw", h.numberOfDaughters(), "number of daughters (Raw)", nanoaod::FlatTable::IntColumn);
            std::vector<const reco::GenParticle *> nonVDau;
            collectDaus(h, nonVDau);
            out->addColumnValue<int>("ndau", int(nonVDau.size()), "number of daughters (non-V)", nanoaod::FlatTable::IntColumn);
            if (nonVDau.size() == 4) {
                std::vector<const reco::GenParticle *> dau;
                float mZ1 = 0, mZ2 = 0, minMllSFOS = 999., minMllAFOS = 999.;
                int nlep = 0, nel = 0, nmu = 0, nnu = 0; float ptmin = 9999., etamax = 0;
                for (unsigned int i = 0; i < 4; ++i) {
                    const reco::GenParticle & di = * nonVDau[i];
                    if (abs(di.pdgId()) == 11 || abs(di.pdgId()) == 13) {
                        nlep++; 
                        nel += (abs(di.pdgId()) == 11);
                        nmu += (abs(di.pdgId()) == 13);
                        ptmin = std::min<float>(ptmin, di.pt()); 
                        etamax = std::max<float>(etamax, std::abs(di.eta()));
                    } else if (abs(di.pdgId()) == 12 || abs(di.pdgId()) == 14 || abs(di.pdgId()) == 16) {
                        nnu++;
                        continue;
                    } else {
                        continue;
                    }
                    for (unsigned int j = i+1; j < 4; ++j) {
                        const reco::GenParticle & dj = * nonVDau[j];
                        if (abs(dj.pdgId()) != 11 && abs(dj.pdgId()) != 13) {
                            continue;
                        }
                        float mij = (di.p4()+dj.p4()).M();
                        if (di.pdgId() == - dj.pdgId()) { // SFOS
                            minMllSFOS = std::min(minMllSFOS, mij);
                            if (std::abs(mij-91.1876) < std::abs(mZ1-91.1876)) {
                                mZ1 = mij;
                                for (unsigned int k = 0; k < 3; ++k) {
                                    if (k == i || k == j) continue;
                                    unsigned int l = 0+1+2+3 - (i+j+k);
                                    const reco::GenParticle & dk = * nonVDau[k];
                                    const reco::GenParticle & dl = * nonVDau[l];
                                    mZ2 = (dk.p4()+dl.p4()).M();
                                }
                            }
                        } 
                        if (di.pdgId() * dj.pdgId() < 0) { // AFOS
                            minMllAFOS = std::min(minMllAFOS, mij);
                        }
                    }
                }
                out->addColumnValue<int>("nlep", nlep, "number of leptons (e/mu)", nanoaod::FlatTable::IntColumn);
                out->addColumnValue<int>("nel", nel, "number of electrons (e/mu)", nanoaod::FlatTable::IntColumn);
                out->addColumnValue<int>("nmu", nmu, "number of muons (e/mu)", nanoaod::FlatTable::IntColumn);
                out->addColumnValue<int>("nnu", nnu, "number of neutrinos (e/mu/tau)", nanoaod::FlatTable::IntColumn);
                out->addColumnValue<float>("ptmin", ptmin, "trailing lepton pT (e/mu)", nanoaod::FlatTable::FloatColumn);
                out->addColumnValue<float>("etamax", etamax, "max lepton abs(eta) (e/mu)", nanoaod::FlatTable::FloatColumn);
                out->addColumnValue<float>("mZ1", nlep == 4 ? mZ1 : -1., "best Z mass (ZZ4l only)", nanoaod::FlatTable::FloatColumn);
                out->addColumnValue<float>("mZ2", nlep == 4 ? mZ2 : -1., "alternate Z mass (ZZ4l only)", nanoaod::FlatTable::FloatColumn);
                out->addColumnValue<float>("minMllSFOS", minMllSFOS, "min m(ll), same flavour, opposite sign", nanoaod::FlatTable::FloatColumn);
                out->addColumnValue<float>("minMllAFOS", minMllAFOS, "min m(ll), any flavour, opposite sign", nanoaod::FlatTable::FloatColumn);
            }  // ndau == 4
            break;
        }
    }

    out->addColumnValue<int>("nhig", found ? 1 : 0, "number of higgses", nanoaod::FlatTable::IntColumn);
    iEvent.put(std::move(out));
}

} // namespace

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HiggsDecayLevelVariablesProducer);

