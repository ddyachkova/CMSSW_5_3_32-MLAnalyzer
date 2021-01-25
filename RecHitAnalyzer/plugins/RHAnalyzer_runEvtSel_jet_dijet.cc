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

TH1D *h_gen_m;
TH1D *h_jet_m;

TH1D *h_dR_jet_W;
TH1D *h_dR_jet_b;
TH1D *h_dR_jet_genTop;



//Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("jet_ind",  &v_jet_Idxs);
  tree->Branch("genpart_ind",   &v_gen_Idxs);

  tree->Branch("jet_pT",  &v_jet_pT_);
  tree->Branch("jet_m0",   &v_jet_m0_);

  tree->Branch("gen_pT",  &v_gen_pT_);
  tree->Branch("gen_m0",   &v_gen_m0_);

  tree->Branch("jetIdxs",  &v_jetIdxs);

  tree->Branch("dR_jet_W",   &v_dR_jet_W);
  tree->Branch("dR_jet_b",  &v_dR_jet_b);
  tree->Branch("dR_jet_genTop",   &v_dR_jet_genTop);

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

  v_jet_Idxs.clear();
  v_gen_Idxs.clear();
  
  v_gen_pT_.clear();
  v_gen_m0_.clear();

  v_jet_pT_.clear();
  v_jet_m0_.clear();

  v_dR_jet_W.clear();
  v_dR_jet_b.clear();
  v_dR_jet_genTop.clear();
  v_jetIdxs.clear(); 

  // int nJet = 0;
  int i=0;
  //int gen_ind=0;
  float dR;
  float dR_sum;
  int ir=0;
  // main loop
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); iGen++ ) {
    //std::cout << "Gen Part Ind " << gen_ind << std::endl;
    int gen_ind = iGen - genParticles-> begin();
    //std::cout << "gen index" << gen_ind << std::endl;
    int id = iGen->pdgId();
    if ( abs(id) != 6 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    //for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {
      //const reco::Candidate* topd = iGen->daughter(iD);
      //top_daughter ->Fill(abs(topd -> pdgId()));
      //if (abs(topd -> pdgId()) != 24) continue;
      //const reco::Candidate *w = get_parent_with_stable_daughter(topd);
      // W daughters loop
      //for (unsigned  int iw = 0; iw < w -> numberOfDaughters(); iw++ ){
      //  const reco::Candidate *w_daughter = w -> daughter(iw);
      //  w_daughters -> Fill(std::abs(w_daughter -> pdgId()));
      //}
    //} // Top daughter loop

    // Loop over reconstructed jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(), iGen->phi() );
      ir += 1;
      dR_sum +=dR;
      if ( dR < 0.8 ) {
      //std::cout << "Gen Sel eta " << iGen -> eta() << std::endl;
      //std::cout << "Jet Sel eta " << iJet -> eta() << std::endl;
      
      //vJetIdxs.push_back(gen_ind);
      
      vJetIdxs.push_back(iJ);
      vGenIdxs.push_back(gen_ind);
      h_jet_m -> Fill(iJet -> mass());
      h_gen_m -> Fill(iGen -> mass());
      std::cout << "Selected gen ind " << gen_ind << std::endl;
      std::cout << "Jet ind " << iJ << std::endl;

      if (abs(iGen -> daughter(0) -> pdgId()) == 24) { 
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) ->phi());
        //v_dR_jet_W.push_back(dR_jet_W);
        //v_dR_jet_b.push_back(dR_jet_b);
        h_dR_jet_W -> Fill(dR_jet_W);
        h_dR_jet_b -> Fill(dR_jet_b);

      }
      else if (abs(iGen -> daughter(0) -> pdgId()) == 5) {
        float dR_jet_W = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(1) -> eta(), iGen -> daughter(1) -> phi());
        float dR_jet_b = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> daughter(0) -> eta(), iGen -> daughter(0) ->phi());
        //v_dR_jet_W.push_back(dR_jet_W);
        //v_dR_jet_b.push_back(dR_jet_b);
        h_dR_jet_W -> Fill(dR_jet_W);
        h_dR_jet_b -> Fill(dR_jet_b);

        }
      //float dR_jet_genTop = reco::deltaR( iJet -> eta(),iJet -> phi(), iGen -> eta(), iGen -> phi());
      h_dR_jet_genTop -> Fill(dR);


      //v_dR_jet_genTop.push_back(dR_jet_genTop);
      }
 
 //     break;
      } // reco jets
    //gen_ind += 1;
    i++;
    //dR = reco::deltaR( iGen -> daughter(1)->eta(),iGen -> daughter(0)->phi(), iGen -> daughter(1)->eta(), iGen -> daughter(1)->phi() );
 }
  return true;
}

void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);
  
  std::cout << "vJetIdxs" << std::endl;
  int n1 = vJetIdxs.size(); 
  for (int i=0; i<n1; i++){
  std::cout << vJetIdxs[i] << std::endl; 
  }
  std::cout << "vGenIdxs" << std::endl;
  int n2 = vGenIdxs.size(); 
  for (int i=0; i<n2; i++){
  std::cout << vGenIdxs[i] << " "; 
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleCollectionT_, genParticles);
  //edm::Handle<std::vector<reco::GenParticle> > genparticles;
  //iEvent.getByLabel( edm::InputTag("genParticles") , genparticles);

  int s = vGenIdxs.size();
  for (int count=0; count < s; count ++){
       int GenId = vGenIdxs[count];
       int JetId = vJetIdxs[count];
       reco::PFJetRef thisJet( jets, JetId );
       reco::GenParticleCollection::const_iterator thisGen = genParticles-> begin() + GenId;
       
       v_gen_pT_.push_back( std::abs(thisGen->pt()) );
       v_gen_m0_.push_back( thisGen->mass() );
       v_jet_pT_.push_back( std::abs(thisJet->pt()));
       v_jet_m0_.push_back( thisJet->mass() );
       
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

  /*
  vector<float> jet_eta;
  vector<float> jet_phi;

  vector<float> gen_eta;
  vector<float> gen_phi;


  std::cout << "vJetIdxs" << std::endl;
  int n1 = vJetIdxs.size(); 
  for (int i=0; i<n1; i++){
  std::cout << vJetIdxs[i] << std::endl; 
  }
  std::cout << "vGenIdxs" << std::endl;
  int n2 = vGenIdxs.size(); 
  for (int i=0; i<n2; i++){
  std::cout << vGenIdxs[i] << " "; 
  }
  //std::cout << vGenIdxs[n2 - 1] << std::endl;
  // Fill branches and histogras


  for(int thisGenIdx : vGenIdxs){
    reco::GenParticleCollection::const_iterator thisGen = genparticles-> begin() + thisGenIdx;
    if ( debug ) std::cout << " >> Gen[" << thisGenIdx << "] Pt:" << thisGen->pt() << std::endl;
    v_gen_pT_.push_back( std::abs(thisGen->pt()) );
    v_gen_m0_.push_back( thisGen->mass() );
    v_gen_Idxs.push_back(thisGenIdx);
    std::cout << "gen eta" << thisGen -> eta() << std::endl;
    std::cout << "gen phi" << thisGen -> phi() << std::endl;
    gen_eta.push_back(thisGen -> eta());
    gen_phi.push_back(thisGen -> phi());

  }
  for(int thisJetIdx : vJetIdxs){
      reco::PFJetRef thisJet( jets, thisJetIdx );
      if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
      v_jet_pT_.push_back( std::abs(thisJet->pt()));
      v_jet_m0_.push_back( thisJet->mass() );
      //jet_eta.push_back(thisJet -> eta());
      //jet_phi.push_back(thisJet -> phi());
      v_jet_Idxs.push_back(thisJetIdx);
      std::cout << "jet eta" << thisJet -> eta() << std::endl;
      std::cout << "jet phi" << thisJet -> phi() << std::endl;
      jet_eta.push_back(thisJet -> eta());
      jet_phi.push_back(thisJet -> phi());
  }

  if (abs(thisGen -> daughter(0) -> pdgId()) == 24) { 
        float dR_jet_W = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(0) -> eta(), thisGen -> daughter(0) -> phi());
        float dR_jet_b = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(1) -> eta(), thisGen -> daughter(1) ->phi());
        std::cout << "dr_jet_W_1 " << dR_jet_W << std::endl;
        std::cout << "dr_jet_b_1 " << dR_jet_b << std::endl;
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
      }
      else if (abs(thisGen -> daughter(0) -> pdgId()) == 5) {
        float dR_jet_W = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(1) -> eta(), thisGen -> daughter(1) -> phi());
        float dR_jet_b = reco::deltaR( thisJet -> eta(), thisJet -> phi(), thisGen -> daughter(0) -> eta(), thisGen -> daughter(0) ->phi());
        std::cout << "dr_jet_W_2 " << dR_jet_W << std::endl;
        std::cout << "dr_jet_b_2 " << dR_jet_b << std::endl;
        v_dR_jet_W.push_back(dR_jet_W);
        v_dR_jet_b.push_back(dR_jet_b);
        }
/
  int n = jet_eta.size();
  for (int i=0; i<n; i++){
    float dR_jet_genTop = reco::deltaR( jet_eta[i], jet_phi[i], gen_eta[i], gen_phi[i]);
     v_dR_jet_genTop.push_back(dR_jet_genTop);
  
  }


}
*/
}
// fillEvtSel_jet_dijet()
