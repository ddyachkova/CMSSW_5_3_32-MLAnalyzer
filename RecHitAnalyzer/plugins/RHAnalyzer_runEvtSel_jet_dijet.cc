#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

using std::vector;

TH1D *h_dijet_jet_pT;
TH1D *h_dijet_jet_E;
TH1D *h_dijet_jet_eta;
TH1D *h_dijet_jet_m0;
TH1D *h_dijet_jet_nJet;
TH1D *w_daughters;
TH1D *top_daughter;


TH1D *reco_Jet_pT;
TH1D *reco_Jet_eta;
TH1D *reco_Jet_phi;
TH1D *reco_Jet_R;
TH1D *dRwb;
TH1D *reco_Jet_m;

TProfile2D  *meanGenLevelDeltaR;
TH2F  *top_genptvm_occupancy;


vector<float> vDijet_jet_pT_;
vector<float> vDijet_jet_m0_;
vector<float> vDijet_jet_eta_;



std::vector<std::vector<int> > seljet_genpart_collid;
std::vector<std::vector<int> > seljet_genpart_pdgid;
std::vector<std::vector<int> > seljet_genpart_charge;

std::vector<std::vector<float> > seljet_genpart_px;
std::vector<std::vector<float> > seljet_genpart_py;
std::vector<std::vector<float> > seljet_genpart_pz;
std::vector<std::vector<float> > seljet_genpart_energy;

std::vector<std::vector<int> > seljet_genpart_status;

std::vector<std::vector<int> > seljet_genpart_motherpdgid;
std::vector<std::vector<int> > seljet_genpart_dau1pdgid;
std::vector<std::vector<int> > seljet_genpart_dau2pdgid;

std::vector<float> seljet_px;
std::vector<float> seljet_py;
std::vector<float> seljet_pz;
std::vector<float> seljet_energy;

std::vector<std::vector<float> > seljet_pfcand_px;
std::vector<std::vector<float> > seljet_pfcand_py;
std::vector<std::vector<float> > seljet_pfcand_pz;
std::vector<std::vector<float> > seljet_pfcand_energy;
std::vector<std::vector<int> > seljet_pfcand_type;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  h_dijet_jet_pT    = fs->make<TH1D>("top_pT"  , "p_{T};p_{T};Particles", 350,  0., 1500.);
  h_dijet_jet_E     = fs->make<TH1D>("top_E"   , "E;E;Particles"        , 400,  0., 2000.);
  h_dijet_jet_eta   = fs->make<TH1D>("top_eta" , "#eta;#eta;Particles"  , 350, -5., 5.);
  h_dijet_jet_nJet  = fs->make<TH1D>("top_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_dijet_jet_m0    = fs->make<TH1D>("top_m0"  , "m_0;m_0;Events"         , 100,  70., 520.);
  w_daughters = fs->make<TH1D>("w_daughters"  , "w_daughters;w_daughters;Events"         , 30,  0., 30.);
  top_daughter = fs->make<TH1D>("top_daughters"  , "top_daughters;top_daughters;Events"         , 30,  0., 30.);

  reco_Jet_pT    = fs->make<TH1D>("reco_Jet_pT"  , "p_{T};p_{T};Particles", 350,  0., 1500.);
  reco_Jet_eta    = fs->make<TH1D>("reco_Jet_eta"  , "#eta; #eta;Particles", 100,  -5, 5.);
  reco_Jet_phi    = fs->make<TH1D>("reco_Jet_phi"  , "phi;phi;Particles", 100,  -5., 100.);
  reco_Jet_R    = fs->make<TH1D>("reco_Jet_R"  , "#DR;#DR;n_{top}", 25,  0., 0.087*25);
  dRwb    = fs->make<TH1D>("delta_R_w_b"  , "#DR;#DR;n_{top}", 25,  0., 0.087*25);
  reco_Jet_m    = fs->make<TH1D>("reco_Jet_m"  , "m;m;Events", 100, 70., 520.);

  meanGenLevelDeltaR = fs->make<TProfile2D>("meanGenLevelDeltaR", "Profile of mean Gen_Level Delta R",10, 70.,520.,10,0.,1500.);

top_genptvm_occupancy = fs->make<TH2F>("top_genptvm_occupancy", "Profile of mean Gen_Level Delta R",10, 70.,520.,10,0.,1500.);

  tree->Branch("jetPt",  &vDijet_jet_pT_);
  tree->Branch("jetM",   &vDijet_jet_m0_);
  tree->Branch("jetEta", &vDijet_jet_eta_);



  tree->Branch("seljet_genpart_collid", &seljet_genpart_collid);
  tree->Branch("seljet_genpart_pdgid", &seljet_genpart_pdgid);
  tree->Branch("seljet_genpart_charge", &seljet_genpart_charge);

  tree->Branch("seljet_genpart_px", &seljet_genpart_px);
  tree->Branch("seljet_genpart_py", &seljet_genpart_py);
  tree->Branch("seljet_genpart_pz", &seljet_genpart_pz);
  tree->Branch("seljet_genpart_energy", &seljet_genpart_energy);

  tree->Branch("seljet_genpart_status", &seljet_genpart_status);

  tree->Branch("seljet_genpart_motherpdgid", &seljet_genpart_motherpdgid);
  tree->Branch("seljet_genpart_dau1pdgid", &seljet_genpart_dau1pdgid);
  tree->Branch("seljet_genpart_dau2pdgid", &seljet_genpart_dau2pdgid);


  tree->Branch("seljet_px", &seljet_px);
  tree->Branch("seljet_py", &seljet_py);
  tree->Branch("seljet_pz", &seljet_pz);
  tree->Branch("seljet_energy", &seljet_energy);

  tree->Branch("seljet_pfcand_px", &seljet_pfcand_px);
  tree->Branch("seljet_pfcand_py", &seljet_pfcand_py);
  tree->Branch("seljet_pfcand_pz", &seljet_pfcand_pz);
  tree->Branch("seljet_pfcand_energy", &seljet_pfcand_energy);
  tree->Branch("seljet_pfcand_type", &seljet_pfcand_type);

} // branchesEvtSel_jet_dijet()

// Run jet selection _____________________________________________________//

const reco::Candidate* get_parent_with_stable_daughter( const reco::Candidate* iC ) {

//    if ( iC->daughter(0)->status() == 1 )
    if ( std::abs(iC->daughter(0)->pdgId()) != 24 )
  return iC;
    else
        return get_parent_with_stable_daughter( iC->daughter(0) );
}


int sum(vector <int> dist) {
    return std::accumulate(dist.begin(), dist.end(), 0);
}

int max_element(vector <int> dist) {
    int max = 0;
    int s = dist.size();
      for (int i = 0; i < s; i++) {
        int el = dist[i];
        if (max < el){max = el;}
              }
        return max;
}

vector <float> get_inverse_pdf(vector <int> dist) {
    vector <float> invpdf(dist.size());
      float summax = sum(dist) / max_element(dist);
      int s = dist.size();
      for (int i = 0; i < s; i++) {
              if (dist[i] != 0 ) {invpdf[i] = summax / dist[i];}
              else {invpdf[i] = 1;}
                }
          return invpdf;
}

float lookup_pt_invpdf(int pTgen, vector <int> pT_bins, vector <float> pT_invpdf) {
    int ipt = 0;
    int s1 = pT_bins.size();
    int s2 = pT_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
            ipt = ib;
            if (ib + 1 >  s2 - 1) { break; }
            if (pTgen <= pT_bins[ib]) { break; }
                        }
        return pT_invpdf[ipt];
}


float get_rand_el(vector <int> dist) {
    int randomIndex = rand() % dist.size();
      return dist[randomIndex];
}

bool RecHitAnalyzer::runEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(jetCollectionT_, jets);
  iEvent.getByLabel(genParticleCollectionT_, genParticles);

  vJetIdxs.clear();
  vDijet_jet_pT_.clear();
  vDijet_jet_m0_.clear();
  vDijet_jet_eta_.clear();

  // std::vector<TLorentzVector> had_tops,bdau,wdau;
  // int nJet = 0;
  int i=0;
  float dR;
  float dR_sum;
  int ir=0;
  vector <int> pT_bins = {0, 24, 32, 54, 30, 54, 56, 42, 36, 36, 36};
  vector <int> m_bins = {0, 0, 0, 38, 98, 108, 90, 66, 0, 0, 0};
  vector <float> m_invpdf = get_inverse_pdf(m_bins);
  vector <float> pT_invpdf = get_inverse_pdf(pT_bins);
  //std::cout << "pT" << sum(pT_bins) <<std::endl;
  //std::cout << "m" << sum(m_bins) <<std::endl;
  //for (int i=0; i<11; i++){
  //  std::cout << "pT" << pT_invpdf[i] <<std::endl; 
  //  std::cout << "m" << m_invpdf[i] << std::endl; 
 // }
  //vector <int> unbiased_pts(100);
  // main loop
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); iGen++ ) {
    double rand_sampler_pT = rand() / double(RAND_MAX);
    //std::cout << " rand pT" << rand_sampler_pT  <<std::endl; 
    float rand_sampler_m = rand() / double(RAND_MAX);
    int pT_gen = iGen -> pt();
    int m_gen = iGen -> mass();
    float pT_wgt = lookup_pt_invpdf(pT_gen, pT_bins, pT_invpdf);
    std::cout << " wgt pT" << pT_wgt  <<std::endl;
    float m_wgt = lookup_pt_invpdf(m_gen, m_bins, m_invpdf);
    if (rand_sampler_pT > pT_wgt) continue;
    if (rand_sampler_m > m_wgt) continue;
    int id = iGen->pdgId();
    if ( abs(id) != 6 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;

    // Top daughter loop
    for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {
      const reco::Candidate* topd = iGen->daughter(iD);
      top_daughter ->Fill(abs(topd -> pdgId()));
      if (abs(topd -> pdgId()) != 24) continue;
      const reco::Candidate *w = get_parent_with_stable_daughter(topd);
      // W daughters loop
      for (unsigned  int iw = 0; iw < w -> numberOfDaughters(); iw++ ){
        const reco::Candidate *w_daughter = w -> daughter(iw);
        w_daughters -> Fill(std::abs(w_daughter -> pdgId()));
      }
    } // Top daughter loop

    // Loop over reconstructed jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      reco::PFJetRef iJet( jets, iJ );
      //if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      //if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(), iGen->phi() );
      ir += 1;
      dR_sum +=dR;
      //std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl;
      
      
      //vDijet_jet_pT_.push_back( std::abs(p.pt()) );
      //vDijet_jet_m0_.push_back(p.mass() );
      //vDijet_jet_eta_.push_back(p.eta() );
      if ( dR > 0.8 ) continue;
      //std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId:" << std::abs(iGen -> pdgId()) << std::endl;
      
      
      reco_Jet_eta->Fill( iJet -> eta() );
      reco_Jet_phi->Fill(iJet -> phi());
      reco_Jet_R->Fill(dR);
      reco_Jet_m->Fill(iJet -> mass());
      reco_Jet_pT->Fill(iJet -> pt());
       break;
      } // reco jets
    i++;
    h_dijet_jet_pT->Fill( std::abs(iGen->pt()) );
    h_dijet_jet_E->Fill( iGen->energy() );
    h_dijet_jet_m0->Fill( iGen->mass() );
    h_dijet_jet_eta->Fill( iGen->eta() );
    vDijet_jet_pT_.push_back( std::abs(iGen->pt()) );
    vDijet_jet_m0_.push_back(iGen->mass() );
    vDijet_jet_eta_.push_back(iGen->eta() );
    dR = reco::deltaR( iGen -> daughter(1)->eta(),iGen -> daughter(0)->phi(), iGen -> daughter(1)->eta(), iGen -> daughter(1)->phi() );
    dRwb -> Fill(dR);
    meanGenLevelDeltaR -> Fill(iGen -> mass(), iGen-> pt(), dR);
    top_genptvm_occupancy -> Fill(iGen -> mass(), iGen-> pt(), 1.);
 }
  return true;
}

void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  seljet_genpart_collid.clear();
  seljet_genpart_pdgid.clear();
  seljet_genpart_charge.clear();

  seljet_genpart_px.clear();
  seljet_genpart_py.clear();
  seljet_genpart_pz.clear();
  seljet_genpart_energy.clear();

  seljet_pz.clear();
  seljet_energy.clear();

  seljet_pfcand_px.clear();
  seljet_pfcand_py.clear();
  seljet_pfcand_pz.clear();
  seljet_pfcand_energy.clear();
  seljet_pfcand_type.clear();


  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByLabel( edm::InputTag("genParticles") , genparticles);

  h_dijet_jet_nJet->Fill( vJetIdxs.size() );
  // Fill branches and histogras
  for(int thisJetIdx : vJetIdxs){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    h_dijet_jet_pT->Fill( std::abs(thisJet->pt()) );
    h_dijet_jet_E->Fill( thisJet->energy() );
    h_dijet_jet_m0->Fill( thisJet->mass() );
    h_dijet_jet_eta->Fill( thisJet->eta() );
    vDijet_jet_pT_.push_back( std::abs(thisJet->pt()) );
    vDijet_jet_m0_.push_back( thisJet->mass() );
    vDijet_jet_eta_.push_back( thisJet->eta() );


    h_dijet_jet_E->Fill( thisJet->energy() );
    h_dijet_jet_m0->Fill( thisJet->mass() );
    h_dijet_jet_eta->Fill( thisJet->eta() );
    vDijet_jet_pT_.push_back( std::abs(thisJet->pt()) );
    vDijet_jet_m0_.push_back( thisJet->mass() );
    vDijet_jet_eta_.push_back( thisJet->eta() );
/*

    seljet_px.push_back(thisJet->px());
    seljet_py.push_back(thisJet->py());
    seljet_pz.push_back(thisJet->pz());
    seljet_energy.push_back(thisJet->energy());

    std::vector<float> pfcand_px;
    std::vector<float> pfcand_py;
    std::vector<float> pfcand_pz;
    std::vector<float> pfcand_energy;
    std::vector<int> pfcand_type;

    for(auto & pfcand : thisJet->getPFConstituents())
    {

      pfcand_px.push_back(pfcand->px());
      pfcand_py.push_back(pfcand->py());
      pfcand_pz.push_back(pfcand->pz());
      pfcand_energy.push_back(pfcand->energy());
      pfcand_type.push_back((int)pfcand->particleId());

    }

    seljet_pfcand_px.push_back(pfcand_px);
    seljet_pfcand_py.push_back(pfcand_py);
    seljet_pfcand_pz.push_back(pfcand_pz);
    seljet_pfcand_energy.push_back(pfcand_energy);
    seljet_pfcand_type.push_back(pfcand_type);

    TLorentzVector TLVJet(thisJet->px(),thisJet->py(),thisJet->pz(),thisJet->energy());
    double cosTheta = TLVJet.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVJet.SetPx(0.0001);

    std::vector<int> genpart_collid;
    std::vector<int> genpart_pdgid;
    std::vector<int> genpart_charge;

    std::vector<float> genpart_px;
    std::vector<float> genpart_py;
    std::vector<float> genpart_pz;
    std::vector<float> genpart_energy;

    std::vector<int> genpart_status;

    std::vector<int> genpart_motherpdgid;
    std::vector<int> genpart_dau1pdgid;
    std::vector<int> genpart_dau2pdgid;

    std::vector<reco::GenParticle>::const_iterator genpartIterator      = (genparticles.product())->begin();
    std::vector<reco::GenParticle>::const_iterator genpartIteratorEnd   = (genparticles.product())->end();
    for ( ; genpartIterator != genpartIteratorEnd; genpartIterator++)
    {

      TLorentzVector TLVgenpart(genpartIterator->px(),genpartIterator->py(),genpartIterator->pz(),genpartIterator->energy());
      cosTheta = TLVgenpart.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVgenpart.SetPx(0.0001); 
      if (TLVJet.DeltaR(TLVgenpart)<0.8)
      {
        genpart_collid.push_back(genpartIterator->collisionId());
        genpart_pdgid.push_back(genpartIterator->pdgId());
        genpart_charge.push_back(genpartIterator->charge());

        genpart_px.push_back(genpartIterator->px());
        genpart_py.push_back(genpartIterator->py());
        genpart_pz.push_back(genpartIterator->pz());
        genpart_energy.push_back(genpartIterator->energy());

        genpart_status.push_back(genpartIterator->status());

        if (genpartIterator->numberOfMothers()>0)
        {
          genpart_motherpdgid.push_back(genpartIterator->mother(0)->pdgId());
        }
        genpart_energy.push_back(genpartIterator->energy());

        genpart_status.push_back(genpartIterator->status());

        if (genpartIterator->numberOfMothers()>0)
        {
          genpart_motherpdgid.push_back(genpartIterator->mother(0)->pdgId());
        }
        else
        {
          genpart_motherpdgid.push_back(-9999);
        }

        switch (genpartIterator->numberOfDaughters())
        {
          case 0:
            genpart_dau1pdgid.push_back(-9999);
            genpart_dau2pdgid.push_back(-9999);
          break;

          case 1:
            genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
            genpart_dau2pdgid.push_back(-9999);
    }//genpart loop
    seljet_genpart_collid.push_back(genpart_collid);
    seljet_genpart_pdgid.push_back(genpart_pdgid);
    seljet_genpart_charge.push_back(genpart_charge);

    seljet_genpart_px.push_back(genpart_px);
    seljet_genpart_py.push_back(genpart_py);
    seljet_genpart_pz.push_back(genpart_pz);
    seljet_genpart_energy.push_back(genpart_energy);

    seljet_genpart_status.push_back(genpart_status);

    seljet_genpart_motherpdgid.push_back(genpart_motherpdgid);
    seljet_genpart_dau1pdgid.push_back(genpart_dau1pdgid);
    seljet_genpart_dau2pdgid.push_back(genpart_dau2pdgid);

  }//jet loop
*/ 
}
}
 // fillEvtSel_jet_dijet()
