#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

using std::vector;

vector<int> v_jet_Idxs;
vector<int> v_gen_Idxs;

vector<float> v_gen_pT_;
vector<float> v_jet_pT_;

vector<float> v_gen_m0_;
vector<float> v_jet_m0_;

vector<float> v_jetIdxs;
//vector<float> v_genIdxs;

vector<float> v_dR_jet_W;
vector<float> v_dR_jet_b;
vector<float> v_dR_jet_genTop;
vector<float> v_dR_jet_W_b;

vector<float> v_jetdR;
vector<float> v_jetdR_W;
vector<float> v_jetdR_b;
vector<float> v_jet_W_b;



TH1D *h_gen_m;
TH1D *h_jet_m;

TH1D *h_dR_jet_W;
TH1D *h_dR_jet_b;
TH1D *h_dR_jet_genTop;



//Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("jetIdx",  &v_jet_Idxs);
  tree->Branch("genIdx",   &v_gen_Idxs);

  tree->Branch("jet_pT",  &v_jet_pT_);
  tree->Branch("jet_m0",   &v_jet_m0_);

  tree->Branch("gen_pT",  &v_gen_pT_);
  tree->Branch("gen_m0",   &v_gen_m0_);

  tree->Branch("dR_jet_W",   &v_dR_jet_W);
  tree->Branch("dR_jet_b",  &v_dR_jet_b);
  tree->Branch("dR_jet_genTop",   &v_dR_jet_genTop);
  tree->Branch("dR_jet_W_b",   &v_dR_jet_W_b);
  
  h_gen_m    = fs->make<TH1D>("gen_m"  , "M", 84,  50., 550);
  h_jet_m    = fs->make<TH1D>("jet_m"  , "M", 84,  0., 450);
  h_dR_jet_W    = fs->make<TH1D>("dR_jet_W"  , "#DR;#DR;n_{top}", 25,  0., 0.087*25);
  h_dR_jet_b    = fs->make<TH1D>("dR_jet_b"  , "#DR;#DR;n_{top}", 25,  0., 0.087*25);
  h_dR_jet_genTop    = fs->make<TH1D>("dR_jet_genTop"  , "#DR;#DR;n_{top}", 25,  0., 0.087*25);

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


  vJetIdxs.clear();
  vGenIdxs.clear();

  v_jetdR.clear();
  v_jetdR_W.clear();
  v_jetdR_b.clear();
  v_jet_W_b.clear();
  
  // int nJet = 0;
  int i=0;
  //int gen_ind=0;
  float dR;
  //float dR_sum;
  //int ir=0;
  // main loop
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); iGen++ ) {
    bool passedGenSel = false;
    int gen_ind = iGen - genParticles-> begin();
    int id = iGen->pdgId();
    if ( abs(id) != 6 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;

    // Loop over reconstructed jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(), iGen->phi() );
      if ( dR > 0.8 ) continue;
      passedGenSel = true;
      
      if (passedGenSel) { 
      vJetIdxs.push_back( iJ );
      vGenIdxs.push_back(gen_ind);
      v_jetdR.push_back( dR );
//      h_dR_jet_genTop -> Fill(dR);
//      h_jet_m -> Fill(iJet -> mass());
//      h_gen_m -> Fill(iGen -> mass());

      if (abs(iGen -> daughter(0) -> pdgId()) == 24) { 
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) ->phi());
//        h_dR_jet_W -> Fill(dR_jet_W);
//        h_dR_jet_b -> Fill(dR_jet_b);
        v_jetdR_W.push_back( dR_jet_W );
        v_jetdR_b.push_back( dR_jet_b );
      }
      else if (abs(iGen -> daughter(0) -> pdgId()) == 5) {
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) ->phi());
//        h_dR_jet_W -> Fill(dR_jet_W);
//        h_dR_jet_b -> Fill(dR_jet_b);
        v_jetdR_W.push_back( dR_jet_W );
        v_jetdR_b.push_back( dR_jet_b );
        }
    
   float dR_jet_W_b = reco::deltaR( iGen -> daughter(0) -> eta(), iGen -> daughter(0) ->phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) ->phi());
        v_jet_W_b.push_back(dR_jet_W_b);
      }
      //float dR_jet_genTop = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> eta(), iGen -> phi());

      //v_dR_jet_genTop.push_back(dR_jet_genTop);
 //     break;
      } // reco jets
    i++;
    //dR = reco::deltaR( iGen -> daughter(1)->eta(),iGen -> daughter(0)->phi(), iGen -> daughter(1)->eta(), iGen -> daughter(1)->phi() );
 }
  return true;
}

void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleCollectionT_, genParticles);

  v_dR_jet_genTop.clear();
  v_dR_jet_W.clear();
  v_dR_jet_b.clear();
  v_gen_pT_.clear();
  v_gen_m0_.clear();
  v_jet_pT_.clear();
  v_jet_m0_.clear();
  
  v_jet_Idxs.clear();
  v_gen_Idxs.clear();
  v_dR_jet_W_b.clear();

  for ( unsigned iG(0); iG != vGenIdxs.size(); ++iG ) {
       reco::GenParticleCollection::const_iterator thisGen = genParticles-> begin() + vGenIdxs[iG];
       v_gen_pT_.push_back( std::abs(thisGen->pt()) );
       v_gen_m0_.push_back( thisGen->mass() );
       v_gen_Idxs.push_back(vGenIdxs[iG]);
  }

  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {
    reco::PFJetRef thisJet( jets, vJetIdxs[iJ] );
    v_dR_jet_genTop.push_back( v_jetdR[iJ] );
    v_dR_jet_W.push_back( v_jetdR_W[iJ] );
    v_dR_jet_b.push_back( v_jetdR_b[iJ] );
    v_dR_jet_W_b.push_back(v_jet_W_b[iJ]);

    v_jet_pT_.push_back( std::abs(thisJet->pt()));
    v_jet_m0_.push_back( thisJet->mass() );
    v_jet_Idxs.push_back(vJetIdxs[iJ]);
    
  }
    /*
       float dR_jet_genTop = reco::deltaR (thisJet -> eta(), thisJet -> phi(), thisGen -> eta(), thisGen -> phi());
       v_dR_jet_genTop.push_back(dR_jet_genTop);
       
      if (abs(thisGen -> daughter(0) -> pdgId()) == 24) { 
        float dR_jet_W = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(0) -> eta(), thisGen -> daughter(0) -> phi());
        float dR_jet_b = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(1) -> eta(), thisGen -> daughter(1) ->phi());
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
      }
      else if (abs(thisGen -> daughter(0) -> pdgId()) == 5) {
        float dR_jet_W = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(1) -> eta(), thisGen -> daughter(1) -> phi());
        float dR_jet_b = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(0) -> eta(), thisGen -> daughter(0) ->phi());
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
        }
  }
*/
  
}
// fillEvtSel_jet_dijet()
