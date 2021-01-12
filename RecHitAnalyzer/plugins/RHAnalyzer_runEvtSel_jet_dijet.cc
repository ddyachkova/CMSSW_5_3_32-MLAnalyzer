#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

using std::vector;

vector<float> v_gen_pT_;
vector<float> v_jet_pT_;

vector<float> v_gen_m0_;
vector<float> v_jet_m0_;

vector<float> v_jetIdxs;
vector<float> v_genIdxs;

vector<float> v_dR_jet_W;
vector<float> v_dR_jet_b;
vector<float> v_dR_jet_genTop;


//Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("jet_pT",  &v_jet_pT_);
  tree->Branch("jet_m0",   &v_jet_m0_);

  tree->Branch("gen_pT",  &v_gen_pT_);
  tree->Branch("gen_m0",   &v_gen_m0_);

  tree->Branch("jetIdxs",  &v_jetIdxs);

  tree->Branch("dR(jet, W)",   &v_dR_jet_W);
  tree->Branch("dR(jet, b)",  &v_dR_jet_b);
  tree->Branch("dR(jet, genTop)",   &v_dR_jet_genTop);

} // branchesEvtSel_jet_dijet()

// Run jet selection _____________________________________________________//


const reco::Candidate* get_parent_with_stable_daughter( const reco::Candidate* iC ) {

//    if ( iC->daughter(0)->status() == 1 )
    if ( std::abs(iC->daughter(0)->pdgId()) != 24 )
  return iC;
    else
        return get_parent_with_stable_daughter( iC->daughter(0) );
}



bool RecHitAnalyzer::runEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(jetCollectionT_, jets);
  iEvent.getByLabel(genParticleCollectionT_, genParticles);


  v_jetIdxs.clear();
  v_genIdxs.clear();

  v_gen_pT_.clear();
  v_gen_m0_.clear();

  v_jet_pT_.clear();
  v_jet_m0_.clear();

  v_dR_jet_W.clear();
  v_dR_jet_b.clear();
  v_dR_jet_genTop.clear();

  // int nJet = 0;
  int i=0;
  float dR;
  float dR_sum;
  int ir=0;
  // main loop
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); iGen++ ) {
    int id = iGen->pdgId();
    if ( abs(id) != 6 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;

    // Loop over reconstructed jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(), iGen->phi() );
      ir += 1;
      dR_sum +=dR;
      if ( dR > 0.8 ) continue;
      v_jetIdxs.push_back(iJ);
      v_genIdxs.push_back(iGen);

      if (abs(iGen -> daughter(0) -> pdgId()) == 24) { 
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) ->phi());
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
      }
      else if (abs(iGen -> daughter(0) -> pdgId()) == 5) {
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) ->phi());
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
        }
      float dR_jet_genTop = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> eta(), iGen -> phi());


      v_dR_jet_genTop.push_back(dR_jet_genTop);

 
      break;
      } // reco jets
    i++;
    dR = reco::deltaR( iGen -> daughter(1)->eta(),iGen -> daughter(0)->phi(), iGen -> daughter(1)->eta(), iGen -> daughter(1)->phi() );
 }
  return true;
}

void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByLabel( edm::InputTag("genParticles") , genparticles);

  // Fill branches and histogras
  for(int thisJetIdx : vJetIdxs){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    v_jet_pT_.push_back( std::abs(thisJet->pt()) );
    v_jet_m0_.push_back( thisJet->mass() );

}
  
    for(int thisGenIdx : vgenIdxs){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    v_jet_pT_.push_back( std::abs(thisJet->pt()) );
    v_jet_m0_.push_back( thisJet->mass() );

}

      v_gen_pT_.push_back( std::abs(iGen->pt()) );
    v_gen_m0_.push_back(iGen->mass() );

}
 // fillEvtSel_jet_dijet()
