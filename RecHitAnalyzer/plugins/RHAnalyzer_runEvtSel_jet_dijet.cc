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


float sum(vector <float> dist) {
  return std::accumulate(dist.begin(), dist.end(), 0.0);
}

float max_element(vector <float> dist) {
    float max = 0;
    int s = dist.size();
      for (int i = 0; i < s; i++) {
        float el = dist[i];
        if (max < el){max = el;}}
        return max;
}

vector <float> get_inverse_pdf(vector <float> dist) {
      vector <float> invpdf(dist.size());
      float sum_hist = sum(dist);
      std::cout << "sum_hist " << sum_hist << std::endl;
      int s = dist.size();
      std::cout << "dist size " << s << std::endl;
      for (int i = 0; i < s; i++) {
              //if (dist[i] != 0 ) {invpdf[i] = sum_hist / dist[i];}
              invpdf[i] = sum_hist / dist[i];
              //else {invpdf[i] = 1;}
                }
      float max_invpdf = max_element(invpdf);
      std::cout << "max " << max_invpdf << std::endl;
      for (int i = 0; i < s; i++) {
              invpdf[i] = invpdf[i] / max_invpdf;}
      return invpdf;
}

float lookup_pt_invpdf(float pTgen, vector <float> pT_bins, vector <float> pT_invpdf) {
    int ipt = 0;
    int s1 = pT_bins.size();
    int s2 = pT_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
            ipt = ib;
            if (ib + 1 >  s2 - 1) { break; }
            if (pTgen <= pT_bins[ib]) { break; }}
    return pT_invpdf[ipt];
}


float get_rand_el(vector <int> dist) {
  int randomIndex = rand() % dist.size();
  return dist[randomIndex];
}


bool resampling(float val, vector<float> bins, vector<float> invpdf){
    double rand_sampler_pT = rand() / double(RAND_MAX);
    //std::cout << "rand" << rand_sampler_pT << std::endl;
    float wgt = lookup_pt_invpdf(val, bins, invpdf);
    //std::cout << "wgt" << wgt << std::endl;
    if (rand_sampler_pT > wgt) {return false;}
    else {return true;}
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
  //vector <int> pT_hist = {603,583,526,458,502,472,461,446,448,487};
  //vector <int> m_hist = {166,380,441,505,600,581,571,562,578,602};
  vector <float> m_hist = {0.0332932,0.0762134,0.0884477,0.101284,0.120337,0.116526,0.114521,0.112716,0.115925,0.12073};
  vector <float> pT_hist = {0.120939,0.116927,0.105495,0.0918572,0.100682,0.0946651,0.0924589,0.0894505,0.0898516,0.0976735};
  //vector <float> pT_bins = {400, 430,490,550,610,670,730,790,850,910,970, 1000};
  //vector <float> m_bins = {85, 105.75,147.25,188.75,230.25,271.75,313.25,354.75,396.25,437.75,479.25,500};
  vector <float> pT_bins = {400.,  460.,  520.,  580.,  640.,  700.,  760.,  820.,  880.,
          940., 1000.};
  vector <float> m_bins = {85. , 126.5, 168. , 209.5, 251. , 292.5, 334. , 375.5, 417. ,
         458.5, 500. };
  
  vector <float> m_invpdf = get_inverse_pdf(m_hist);
  int m_hist_size = m_hist.size();
  std::cout << "hist " << std::endl;
  for (int i=0; i< m_hist_size; i++) {std::cout << m_hist[i] << " ";}
  std::cout << std::endl;
  
  int m_size = m_invpdf.size();
  for (int i=0; i< m_size; i++) {std::cout << m_invpdf[i] << " ";}
  std::cout << std::endl;
  
  vector <float> pT_invpdf = get_inverse_pdf(pT_hist);
  int pT_size = pT_invpdf.size();
  for (int i=0; i< pT_size; i++) {std::cout << pT_invpdf[i] << " ";}
  std::cout << std::endl;
  // main loop
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); iGen++ ) {
    float pT_gen = iGen -> pt();
    //std::cout << "pT_gen" << pT_gen << "   ";
    float m_gen = iGen -> mass();
    if (resampling(pT_gen, pT_bins, pT_invpdf) != true) continue;
    if (resampling(m_gen, pT_bins, pT_invpdf) != true) continue;
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
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(), iGen->phi() );
      ir += 1;
      dR_sum +=dR;
      //std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl;
      //vDijet_jet_pT_.push_back( std::abs(p.pt()) );
      //vDijet_jet_m0_.push_back(p.mass() );
      //vDijet_jet_eta_.push_back(p.eta() );
      if ( dR > 0.8 ) continue;
      //std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId:" << std::abs(iGen -> pdgId()) << std::endl;
      top_genptvm_occupancy -> Fill(iGen -> mass(), iGen-> pt(), 1.);
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
    //top_genptvm_occupancy -> Fill(iGen -> mass(), iGen-> pt(), 1.);
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
