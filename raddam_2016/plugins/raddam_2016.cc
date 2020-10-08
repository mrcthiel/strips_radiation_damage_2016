// -*- C++ -*-
// 
// Package:    raddam_2016/raddam_2016
// Class:      raddam_2016
//
/**\class raddam_2016 raddam_2016.cc raddam_2016/raddam_2016/plugins/raddam_2016.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Mon, 03 Aug 2020 21:38:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

#include "RecoCTPPS/TotemRPLocal/interface/FastLineRecognition.h"

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "RecoCTPPS/TotemRPLocal/interface/TotemRPLocalTrackFitterAlgorithm.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TF1.h>


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class raddam_2016 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit raddam_2016(const edm::ParameterSet&);
		~raddam_2016();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::InputTag tagRecHit;
		edm::EDGetTokenT<CTPPSLocalTrackLiteCollection> tracksToken_;

		FastLineRecognition *lrcgn;

		edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;

		void recognizeAndSelect(TotemRPUVPattern::ProjectionType proj, double z0, double threshold,unsigned int planes_required,const edm::DetSetVector<TotemRPRecHit> &hits, edm::DetSet<TotemRPUVPattern> &patterns);

		TotemRPLocalTrackFitterAlgorithm fitter_;

		edm::EDGetTokenT<edm::DetSetVector<TotemRPRecHit > > detSetVectorTotemRPRecHitToken;

		TH2D * x_y_all_planes;
		TH2D * x_y_0;
		TH2D * x_y_1;
		TH2D * x_y_2;
		TH2D * x_y_3;
		TH2D * x_y_4;
		TH2D * x_y_5;
		TH2D * x_y_6;
		TH2D * x_y_7;
		TH2D * x_y_8;
		TH2D * x_y_9;

		TH2D * x_corr;
		TH2D * y_corr;
		TH1D * x_diff;
		TH1D * y_diff;

                double mean_02, sigma_02, mean_04, sigma_04, mean_06, sigma_06, mean_08, sigma_08, mean_24, sigma_24, mean_26, sigma_26, mean_28, sigma_28, mean_46, sigma_46, mean_48, sigma_48, mean_68, sigma_68, mean_13, sigma_13, mean_15, sigma_15, mean_17, sigma_17, mean_19, sigma_19, mean_35, sigma_35, mean_37, sigma_37, mean_39, sigma_39, mean_57, sigma_57, mean_59, sigma_59, mean_79, sigma_79;

};


using namespace std;
using namespace edm;


raddam_2016::raddam_2016(const edm::ParameterSet& iConfig):
	tagRecHit(iConfig.getParameter<edm::InputTag >("tagRecHit")),
	tracksToken_(consumes<CTPPSLocalTrackLiteCollection>(iConfig.getParameter<edm::InputTag>("tagLocalTrackLite"))),
	lrcgn(new FastLineRecognition(iConfig.getParameter<double>("clusterSize_a"), iConfig.getParameter<double>("clusterSize_b"))),
	fitter_(iConfig){
		usesResource("TFileService");
		edm::Service<TFileService> fs;
		detSetVectorTotemRPRecHitToken = consumes<edm::DetSetVector<TotemRPRecHit>>(tagRecHit);


		TFile *file2 = TFile::Open("/eos/project/c/ctpps/subsystems/Strips/StripsTracking/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root");

		TH2F * h_eff_jonathan;
		file2->GetObject("Strips/2016/2016C/h56_2016C_RP102_all_2D",h_eff_jonathan);

		x_y_all_planes = fs->make<TH2D>("x_y_all_planes", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_0 = fs->make<TH2D>("x_y_0", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_1 = fs->make<TH2D>("x_y_1", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_2 = fs->make<TH2D>("x_y_2", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_3 = fs->make<TH2D>("x_y_3", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_4 = fs->make<TH2D>("x_y_4", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_5 = fs->make<TH2D>("x_y_5", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_6 = fs->make<TH2D>("x_y_6", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_7 = fs->make<TH2D>("x_y_7", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_8 = fs->make<TH2D>("x_y_8", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());
		x_y_9 = fs->make<TH2D>("x_y_9", ";x   (mm);y   (mm);", h_eff_jonathan->GetNbinsX(),h_eff_jonathan->GetXaxis()->GetXbins()->GetArray(),h_eff_jonathan->GetNbinsY(),h_eff_jonathan->GetYaxis()->GetXbins()->GetArray());

		x_corr = fs->make<TH2D>("x_corr", ";x all;x -1;", 100, 0., 15., 100, 1., 15.);
		y_corr = fs->make<TH2D>("y_corr", ";y all;y -1;", 100, -5., +5., 100, -5., +5.);

		x_diff = fs->make<TH1D>("x_diff", ";(x_{all}-x_{all-1})/x_{all};", 100, -0.1, 0.1);
		y_diff = fs->make<TH1D>("y_diff", ";(y_{all}-y_{all-1})/y_{all};", 100, -0.1, 0.1);

                TFile *rootfile_ = TFile::Open("/afs/cern.ch/work/m/mthiel/private/PPS_efficiency_2019/04_02_2020/optimization/CMSSW_10_6_12/src/raddam_2016/raddam_2016/parameters.root");

                TH1D * diff_U_02;
                TH1D * diff_U_04;
                TH1D * diff_U_06;
                TH1D * diff_U_08;
                TH1D * diff_U_24;
                TH1D * diff_U_26;
                TH1D * diff_U_28;
                TH1D * diff_U_46;
                TH1D * diff_U_48;
                TH1D * diff_U_68;
                TH1D * diff_V_13;
                TH1D * diff_V_15;
                TH1D * diff_V_17;
                TH1D * diff_V_19;
                TH1D * diff_V_35;
                TH1D * diff_V_37;
                TH1D * diff_V_39;
                TH1D * diff_V_57;
                TH1D * diff_V_59;
                TH1D * diff_V_79;

	       	rootfile_->GetObject("demo/diffg_U_02",diff_U_02);
                rootfile_->GetObject("demo/diffg_U_04",diff_U_04);
                rootfile_->GetObject("demo/diffg_U_06",diff_U_06);
                rootfile_->GetObject("demo/diffg_U_08",diff_U_08);
                rootfile_->GetObject("demo/diffg_U_24",diff_U_24);
                rootfile_->GetObject("demo/diffg_U_26",diff_U_26);
                rootfile_->GetObject("demo/diffg_U_28",diff_U_28);
                rootfile_->GetObject("demo/diffg_U_46",diff_U_46);
                rootfile_->GetObject("demo/diffg_U_48",diff_U_48);
                rootfile_->GetObject("demo/diffg_U_68",diff_U_68);
                rootfile_->GetObject("demo/diffg_V_13",diff_V_13);
                rootfile_->GetObject("demo/diffg_V_15",diff_V_15);
                rootfile_->GetObject("demo/diffg_V_17",diff_V_17);
                rootfile_->GetObject("demo/diffg_V_19",diff_V_19);
                rootfile_->GetObject("demo/diffg_V_35",diff_V_35);
                rootfile_->GetObject("demo/diffg_V_37",diff_V_37);
                rootfile_->GetObject("demo/diffg_V_39",diff_V_39);
                rootfile_->GetObject("demo/diffg_V_57",diff_V_57);
                rootfile_->GetObject("demo/diffg_V_59",diff_V_59);
                rootfile_->GetObject("demo/diffg_V_79",diff_V_79);


              	diff_U_02->Fit("gaus");
                TF1 *g_02 = (TF1*)diff_U_02->GetListOfFunctions()->FindObject("gaus");
                mean_02 = g_02->GetParameter(1);
                sigma_02 = g_02->GetParameter(2);

                diff_U_04->Fit("gaus");
                TF1 *g_04 = (TF1*)diff_U_04->GetListOfFunctions()->FindObject("gaus");
                mean_04 = g_04->GetParameter(1);
                sigma_04 = g_04->GetParameter(2);

                diff_U_06->Fit("gaus");
                TF1 *g_06 = (TF1*)diff_U_06->GetListOfFunctions()->FindObject("gaus");
                mean_06 = g_06->GetParameter(1);
                sigma_06 = g_06->GetParameter(2);

                diff_U_08->Fit("gaus");
                TF1 *g_08 = (TF1*)diff_U_08->GetListOfFunctions()->FindObject("gaus");
                mean_08 = g_08->GetParameter(1);
                sigma_08 = g_08->GetParameter(2);

                diff_U_24->Fit("gaus");
                TF1 *g_24 = (TF1*)diff_U_24->GetListOfFunctions()->FindObject("gaus");
                mean_24 = g_24->GetParameter(1);
                sigma_24 = g_24->GetParameter(2);

                diff_U_26->Fit("gaus");
                TF1 *g_26 = (TF1*)diff_U_26->GetListOfFunctions()->FindObject("gaus");
                mean_26 = g_26->GetParameter(1);
                sigma_26 = g_26->GetParameter(2);

                diff_U_28->Fit("gaus");
                TF1 *g_28 = (TF1*)diff_U_28->GetListOfFunctions()->FindObject("gaus");
                mean_28 = g_28->GetParameter(1);
                sigma_28 = g_28->GetParameter(2);

                diff_U_46->Fit("gaus");
                TF1 *g_46 = (TF1*)diff_U_46->GetListOfFunctions()->FindObject("gaus");
                mean_46 = g_46->GetParameter(1);
                sigma_46 = g_46->GetParameter(2);

                diff_U_48->Fit("gaus");
                TF1 *g_48 = (TF1*)diff_U_48->GetListOfFunctions()->FindObject("gaus");
                mean_48 = g_48->GetParameter(1);
                sigma_48 = g_48->GetParameter(2);

                diff_U_68->Fit("gaus");
                TF1 *g_68 = (TF1*)diff_U_68->GetListOfFunctions()->FindObject("gaus");
                mean_68 = g_68->GetParameter(1);
                sigma_68 = g_68->GetParameter(2);

                diff_V_13->Fit("gaus");
                TF1 *g_13 = (TF1*)diff_V_13->GetListOfFunctions()->FindObject("gaus");
                mean_13 = g_13->GetParameter(1);
                sigma_13 = g_13->GetParameter(2);

                diff_V_15->Fit("gaus");
                TF1 *g_15 = (TF1*)diff_V_15->GetListOfFunctions()->FindObject("gaus");
                mean_15 = g_15->GetParameter(1);
                sigma_15 = g_15->GetParameter(2);

                diff_V_17->Fit("gaus");
                TF1 *g_17 = (TF1*)diff_V_17->GetListOfFunctions()->FindObject("gaus");
                mean_17 = g_17->GetParameter(1);
                sigma_17 = g_17->GetParameter(2);

                diff_V_19->Fit("gaus");
                TF1 *g_19 = (TF1*)diff_V_19->GetListOfFunctions()->FindObject("gaus");
                mean_19 = g_19->GetParameter(1);
                sigma_19 = g_19->GetParameter(2);

                diff_V_35->Fit("gaus");
                TF1 *g_35 = (TF1*)diff_V_35->GetListOfFunctions()->FindObject("gaus");
                mean_35 = g_35->GetParameter(1);
                sigma_35 = g_35->GetParameter(2);

                diff_V_37->Fit("gaus");
                TF1 *g_37 = (TF1*)diff_V_37->GetListOfFunctions()->FindObject("gaus");
                mean_37 = g_37->GetParameter(1);
                sigma_37 = g_37->GetParameter(2);

                diff_V_39->Fit("gaus");
                TF1 *g_39 = (TF1*)diff_V_39->GetListOfFunctions()->FindObject("gaus");
                mean_39 = g_39->GetParameter(1);
                sigma_39 = g_39->GetParameter(2);

                diff_V_57->Fit("gaus");
                TF1 *g_57 = (TF1*)diff_V_57->GetListOfFunctions()->FindObject("gaus");
                mean_57 = g_57->GetParameter(1);
                sigma_57 = g_57->GetParameter(2);

                diff_V_59->Fit("gaus");
                TF1 *g_59 = (TF1*)diff_V_59->GetListOfFunctions()->FindObject("gaus");
                mean_59 = g_59->GetParameter(1);
                sigma_59 = g_59->GetParameter(2);

                diff_V_79->Fit("gaus");
                TF1 *g_79 = (TF1*)diff_V_79->GetListOfFunctions()->FindObject("gaus");
                mean_79 = g_79->GetParameter(1);
                sigma_79 = g_79->GetParameter(2);

}


raddam_2016::~raddam_2016(){ 
	delete lrcgn;
}

void raddam_2016::recognizeAndSelect(TotemRPUVPattern::ProjectionType proj,double z0, double threshold_loc, unsigned int planes_required,const DetSetVector<TotemRPRecHit> &hits, DetSet<TotemRPUVPattern> &patterns){
	DetSet<TotemRPUVPattern> newPatterns;
	lrcgn->getPatterns(hits, z0, threshold_loc, newPatterns);

	for (auto &p : newPatterns)
	{
		p.setProjection(proj);

		p.setFittable(true);

		set<unsigned int> planes;
		for (const auto &ds : p.getHits())
			planes.insert(TotemRPDetId(ds.detId()).plane());

		if (planes.size() < planes_required)
			p.setFittable(false);

		if (fabs(p.getA()) > 10.0) //max_a_toFit
			p.setFittable(false);

		patterns.push_back(p);
	}
}


//
// member functions
//

// ------------ method called for each event  ------------
	void
raddam_2016::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	double n_sigma = 3;

	sigma_02 = n_sigma*sigma_02;
	sigma_04 = n_sigma*sigma_04;
	sigma_06 = n_sigma*sigma_06;
	sigma_08 = n_sigma*sigma_08;
	sigma_24 = n_sigma*sigma_24;
	sigma_26 = n_sigma*sigma_26;
	sigma_28 = n_sigma*sigma_28;
	sigma_46 = n_sigma*sigma_46;
	sigma_48 = n_sigma*sigma_48;
	sigma_68 = n_sigma*sigma_68;
	sigma_13 = n_sigma*sigma_13;
	sigma_15 = n_sigma*sigma_15;
	sigma_17 = n_sigma*sigma_17;
	sigma_19 = n_sigma*sigma_19;
	sigma_35 = n_sigma*sigma_35;
	sigma_37 = n_sigma*sigma_37;
	sigma_39 = n_sigma*sigma_39;
	sigma_57 = n_sigma*sigma_57;
	sigma_59 = n_sigma*sigma_59;
	sigma_79 = n_sigma*sigma_79;

	edm::Handle<CTPPSLocalTrackLiteCollection> hTracks;
	iEvent.getByToken(tracksToken_, hTracks);

	bool pass = false;
	double x_official = -999.;
	double y_official = -999.;
	double x_mine = -999.;
	double y_mine = -999.;

	std::vector<TH2D*> histos;
	histos.clear();
	histos.push_back(x_y_0);
	histos.push_back(x_y_1);
	histos.push_back(x_y_2);
	histos.push_back(x_y_3);
	histos.push_back(x_y_4);
	histos.push_back(x_y_5);
	histos.push_back(x_y_6);
	histos.push_back(x_y_7);
	histos.push_back(x_y_8);
	histos.push_back(x_y_9);



	for (unsigned int idx = 0; idx < hTracks->size(); ++idx){
		const auto& tr = hTracks->at(idx);
		const CTPPSDetId rpId(tr.getRPId());
		unsigned int rpId_n = rpId.getRPId();
		if(rpId_n==1997537280) x_official =  tr.getX();
		if(rpId_n==1997537280) y_official =  tr.getY();
		if(rpId_n==1997537280) pass = true;
	}

	if(pass){

		x_y_all_planes->Fill(x_official,y_official);

		edm::Handle<edm::DetSetVector<TotemRPRecHit>> input;
		iEvent.getByToken(detSetVectorTotemRPRecHitToken, input);
		ESHandle<CTPPSGeometry> geometry;
		iSetup.get<VeryForwardRealGeometryRecord>().get(geometry);
		if (geometryWatcher.check(iSetup)) lrcgn->resetGeometry(geometry.product());

		std::vector<double> plane_0;
		plane_0.clear();
		std::vector<double> plane_1;
		plane_1.clear();
		std::vector<double> plane_2;
		plane_2.clear();
		std::vector<double> plane_3;
		plane_3.clear();
		std::vector<double> plane_4;
		plane_4.clear();
		std::vector<double> plane_5;
		plane_5.clear();
		std::vector<double> plane_6;
		plane_6.clear();
		std::vector<double> plane_7;
		plane_7.clear();
		std::vector<double> plane_8;
		plane_8.clear();
		std::vector<double> plane_9;
		plane_9.clear();


		for (auto &ids : *input) {
			TotemRPDetId detId(ids.detId());
			if(!(detId.arm()==1 && detId.rp()==2)) continue;
			unsigned int plane = detId.plane();
			if(plane==0){
				for (auto &h : ids){
					plane_0.push_back(h.getPosition());
				}
			}
			if(plane==1){
				for (auto &h : ids){
					plane_1.push_back(h.getPosition());
				}
			}
			if(plane==2){
				for (auto &h : ids){
					plane_2.push_back(h.getPosition());
				}
			}
			if(plane==3){
				for (auto &h : ids){
					plane_3.push_back(h.getPosition());
				}
			}
			if(plane==4){
				for (auto &h : ids){
					plane_4.push_back(h.getPosition());
				}
			}
			if(plane==5){
				for (auto &h : ids){
					plane_5.push_back(h.getPosition());
				}
			}
			if(plane==6){
				for (auto &h : ids){
					plane_6.push_back(h.getPosition());
				}
			}
			if(plane==7){
				for (auto &h : ids){
					plane_7.push_back(h.getPosition());
				}
			}
			if(plane==8){
				for (auto &h : ids){
					plane_8.push_back(h.getPosition());
				}
			}
			if(plane==9){
				for (auto &h : ids){
					plane_9.push_back(h.getPosition());
				}
			}
		}


		for(int nn=0; nn<5; nn++){
			struct RPData
			{
				DetSetVector<TotemRPRecHit> hits_U, hits_V;
				map<uint8_t, uint16_t> planeOccupancy_U, planeOccupancy_V;
			};
			map<unsigned int, RPData> rpData;

			for (auto &ids : *input) {
				TotemRPDetId detId(ids.detId());
				if(!(detId.arm()==1 && detId.rp()==2)) continue;
				unsigned int plane = detId.plane();
				if(nn==0 && (plane==0 || plane==1)) continue;
				if(nn==1 && (plane==2 || plane==3)) continue;
				if(nn==2 && (plane==4 || plane==5)) continue;
				if(nn==3 && (plane==5 || plane==7)) continue;
				if(nn==4 && (plane==8 || plane==9)) continue;

				bool uDir = detId.isStripsCoordinateUDirection();
				unsigned int rpId = detId.getRPId();
				unsigned int arm = detId.arm();
				RPData &data = rpData[rpId];


				for (auto &h : ids){
					if (uDir){
						auto &ods = data.hits_U.find_or_insert(ids.detId());
						ods.push_back(h);
						data.planeOccupancy_U[plane]++;
					} else {
						auto &ods = data.hits_V.find_or_insert(ids.detId());
						ods.push_back(h);
						data.planeOccupancy_V[plane]++;
					}
				}
			}


			DetSetVector<TotemRPUVPattern> patternsVector;
			for (auto it : rpData){
				CTPPSDetId rpId(it.first);
				RPData &data = it.second;
				auto &uColl = data.planeOccupancy_U;
				auto &vColl = data.planeOccupancy_V;

				unsigned int uPlanes = 0, vPlanes = 0;
				for (auto pit : uColl)
					if (pit.second <= 5)
						uPlanes++;

				for (auto pit : vColl)
					if (pit.second <= 5)
						vPlanes++;

				if (uPlanes < 3/*minPlanesPerProjectionToSearch*/ || vPlanes < 3/*minPlanesPerProjectionToSearch*/) continue;

				DetSet<TotemRPUVPattern> &patterns = patternsVector.find_or_insert(rpId);
				double z0 = geometry->getRP(rpId)->translation().z();

				recognizeAndSelect(TotemRPUVPattern::projU, z0, 2.99, 3, data.hits_U, patterns);
				recognizeAndSelect(TotemRPUVPattern::projV, z0, 2.99, 3, data.hits_V, patterns);

			}


			if (geometryWatcher.check(iSetup)) fitter_.reset();

			for (const auto &rpv : patternsVector){
				if(rpv.size()==0) break;
				CTPPSDetId rpId(rpv.detId());
				unsigned int n_U=0, n_V=0;
				unsigned int idx_U=0, idx_V=0;
				for (unsigned int pi = 0; pi < rpv.size(); pi++){
					const TotemRPUVPattern &pattern = rpv[pi];
					switch (pattern.getProjection())
					{
						case TotemRPUVPattern::projU:
							n_U++;
							idx_U=pi;
							break;
						case TotemRPUVPattern::projV:
							n_V++;
							idx_V=pi;
							break;
						default:
							break;
					}   
				}


				if (!rpv[idx_U].getFittable() || !rpv[idx_V].getFittable()) continue;

				DetSetVector<TotemRPRecHit> hits;
				for (auto &ids : rpv[idx_U].getHits())
				{
					auto &ods = hits.find_or_insert(ids.detId());
					for (auto &h : ids)
						ods.push_back(h);
				}

				for (auto &ids : rpv[idx_V].getHits())
				{
					auto &ods = hits.find_or_insert(ids.detId());
					for (auto &h : ids)
						ods.push_back(h);
				}

				double z0 = geometry->getRP(rpId)->translation().z();

				TotemRPLocalTrack track;
				fitter_.fitTrack(hits, z0, *geometry, track);

				x_mine = track.getX0();
				y_mine = track.getY0();

				histos.at(nn)->Fill(x_mine,y_mine);
//                                histos.at(nn)->Fill(x_official,y_official);

				bool save_U = false;
                                bool save_V = false;

				if(nn==0){
				        for (unsigned int mm = 0; mm < plane_0.size(); ++mm){
						for (unsigned int jj = 0; jj < plane_2.size(); ++jj){
							if(abs((plane_0.at(mm) - plane_2.at(jj))-mean_02)<sigma_02) save_U = true;
						}
                                                for (unsigned int jj = 0; jj < plane_4.size(); ++jj){
                                                        if(abs((plane_0.at(mm) - plane_4.at(jj))-mean_04)<sigma_04) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_6.size(); ++jj){
                                                        if(abs((plane_0.at(mm) - plane_6.at(jj))-mean_06)<sigma_06) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_8.size(); ++jj){
                                                        if(abs((plane_0.at(mm) - plane_8.at(jj))-mean_08)<sigma_08) save_U = true;
						}
					}
                                        for (unsigned int mm = 0; mm < plane_1.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_3.size(); ++jj){
                                                        if(abs((plane_1.at(mm) - plane_3.at(jj))-mean_13)<sigma_13) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_5.size(); ++jj){
                                                        if(abs((plane_1.at(mm) - plane_5.at(jj))-mean_15)<sigma_15) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_7.size(); ++jj){
                                                        if(abs((plane_1.at(mm) - plane_7.at(jj))-mean_17)<sigma_17) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_9.size(); ++jj){
                                                        if(abs((plane_1.at(mm) - plane_9.at(jj))-mean_19)<sigma_19) save_V = true;
                                                }
                                        }
				}

                                if(nn==1){
                                        for (unsigned int mm = 0; mm < plane_2.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_0.size(); ++jj){
                                                        if(abs((plane_0.at(jj) - plane_2.at(mm))-mean_02)<sigma_02) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_4.size(); ++jj){
                                                        if(abs((plane_2.at(mm) - plane_4.at(jj))-mean_24)<sigma_24) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_6.size(); ++jj){
                                                        if(abs((plane_2.at(mm) - plane_6.at(jj))-mean_26)<sigma_26) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_8.size(); ++jj){
                                                        if(abs((plane_2.at(mm) - plane_8.at(jj))-mean_28)<sigma_28) save_U = true;
                                                }
                                        }
                                        for (unsigned int mm = 0; mm < plane_3.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_1.size(); ++jj){
                                                        if(abs((plane_1.at(jj) - plane_3.at(mm))-mean_13)<sigma_13) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_5.size(); ++jj){
                                                        if(abs((plane_3.at(mm) - plane_5.at(jj))-mean_35)<sigma_35) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_7.size(); ++jj){
                                                        if(abs((plane_3.at(mm) - plane_7.at(jj))-mean_37)<sigma_37) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_9.size(); ++jj){
                                                        if(abs((plane_3.at(mm) - plane_9.at(jj))-mean_39)<sigma_39) save_V = true;
                                                }
                                        }
                                }

                                if(nn==2){
                                        for (unsigned int mm = 0; mm < plane_4.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_0.size(); ++jj){
                                                        if(abs((plane_0.at(jj) - plane_4.at(mm))-mean_04)<sigma_04) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_2.size(); ++jj){
                                                        if(abs((plane_2.at(jj) - plane_4.at(mm))-mean_24)<sigma_24) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_6.size(); ++jj){
                                                        if(abs((plane_4.at(mm) - plane_6.at(jj))-mean_46)<sigma_46) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_8.size(); ++jj){
                                                        if(abs((plane_4.at(mm) - plane_8.at(jj))-mean_48)<sigma_48) save_U = true;
                                                }
                                        }
                                        for (unsigned int mm = 0; mm < plane_5.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_1.size(); ++jj){
                                                        if(abs((plane_1.at(jj) - plane_5.at(mm))-mean_15)<sigma_15) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_3.size(); ++jj){
                                                        if(abs((plane_3.at(jj) - plane_5.at(mm))-mean_35)<sigma_13) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_7.size(); ++jj){
                                                        if(abs((plane_5.at(mm) - plane_7.at(jj))-mean_57)<sigma_57) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_9.size(); ++jj){
                                                        if(abs((plane_5.at(mm) - plane_9.at(jj))-mean_59)<sigma_59) save_V = true;
                                                }
                                        }
                                }

                                if(nn==3){
                                        for (unsigned int mm = 0; mm < plane_6.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_0.size(); ++jj){
                                                        if(abs((plane_0.at(jj) - plane_6.at(mm))-mean_06)<sigma_06) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_2.size(); ++jj){
                                                        if(abs((plane_2.at(jj) - plane_6.at(mm))-mean_26)<sigma_26) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_4.size(); ++jj){
                                                        if(abs((plane_4.at(jj) - plane_6.at(mm))-mean_46)<sigma_46) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_8.size(); ++jj){
                                                        if(abs((plane_6.at(mm) - plane_8.at(jj))-mean_68)<sigma_68) save_U = true;
                                                }
                                        }
                                        for (unsigned int mm = 0; mm < plane_7.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_1.size(); ++jj){
                                                        if(abs((plane_1.at(jj) - plane_7.at(mm))-mean_17)<sigma_17) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_3.size(); ++jj){
                                                        if(abs((plane_3.at(jj) - plane_7.at(mm))-mean_37)<sigma_37) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_5.size(); ++jj){
                                                        if(abs((plane_5.at(jj) - plane_7.at(mm))-mean_57)<sigma_57) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_9.size(); ++jj){
                                                        if(abs((plane_7.at(mm) - plane_9.at(jj))-mean_79)<sigma_79) save_V = true;
                                                }
                                        }
                                }

                                if(nn==4){
                                        for (unsigned int mm = 0; mm < plane_8.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_0.size(); ++jj){
                                                        if(abs((plane_0.at(jj) - plane_8.at(mm))-mean_08)<sigma_08) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_2.size(); ++jj){
                                                        if(abs((plane_2.at(jj) - plane_8.at(mm))-mean_28)<sigma_28) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_4.size(); ++jj){
                                                        if(abs((plane_4.at(jj) - plane_8.at(mm))-mean_48)<sigma_48) save_U = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_6.size(); ++jj){
                                                        if(abs((plane_6.at(jj) - plane_8.at(mm))-mean_68)<sigma_68) save_U = true;
                                                }
                                        }
                                        for (unsigned int mm = 0; mm < plane_9.size(); ++mm){
                                                for (unsigned int jj = 0; jj < plane_1.size(); ++jj){
                                                        if(abs((plane_1.at(jj) - plane_9.at(mm))-mean_19)<sigma_19) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_3.size(); ++jj){
                                                        if(abs((plane_3.at(jj) - plane_9.at(mm))-mean_39)<sigma_39) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_5.size(); ++jj){
                                                        if(abs((plane_5.at(jj) - plane_9.at(mm))-mean_59)<sigma_59) save_V = true;
                                                }
                                                for (unsigned int jj = 0; jj < plane_7.size(); ++jj){
                                                        if(abs((plane_7.at(jj) - plane_9.at(mm))-mean_79)<sigma_79) save_V = true;
                                                }
                                        }
                                }



				if(save_V && save_U){
				  histos.at(nn+5)->Fill(x_mine,y_mine);
                                  //histos.at(nn+5)->Fill(x_official,y_official);
				}
			}
		}
	}


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
	void
raddam_2016::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void
raddam_2016::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
raddam_2016::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(raddam_2016);
