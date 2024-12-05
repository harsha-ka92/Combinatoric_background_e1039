#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TMultiGraph.h>
#include <ktracker/SRecEvent.h>
#include <ktracker/SRawEvent.h>
using namespace std;

TFile *f_file;
TTree *_mixed;
TTree *_sorted;
TTree *_ana;

void plots(const string list_run_ana)
{   
	gSystem->Exec("rm -rf plots");
    	gSystem->mkdir("plots", 1);
    	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
  	gStyle->SetPadTickY(1);
	
	vector<string> list_in;
  	ifstream ifs(list_run_ana);

  	if (!ifs.is_open()) {
        	std::cout<<"___file for run ID "<<list_run_ana<<" could not be opened"<<std::endl;
        	std::exit(1);
    	}
	
	const string dir_in = "scratch_all_good_fpga1_only/run_00";
	const string fname = "/out/output.root";
  	int run_id;
	int run_id_low=999999;
	int run_id_high=0;
  	int n_runs =0;

	while (ifs >> run_id){ 
		list_in.push_back(dir_in + std::to_string(run_id) + fname);
		if(run_id_low > run_id) run_id_low = run_id;
		if(run_id_high < run_id) run_id_high = run_id;
		n_runs++;
	}
	ifs.close();

	ostringstream oss;	
	oss<<"Dimuon Mass Spectra (For "<<n_runs<<" runs between : "<<run_id_low<<" - "<<run_id_high<<")";
	TH1D* m_sorted = new TH1D("m_sorted",oss.str().c_str(), 90 , 0, 9);
	TH1D* m_mixed = new TH1D("m_mixed","m_mixed", 90 , 0, 9);
	
	oss.str("");
        oss<<"Target-Like Dimuon Mass Spectra (For "<<n_runs<<" runs between : "<<run_id_low<<" - "<<run_id_high<<")";
	TH1D* m_sorted_target = new TH1D("m_sorted_target",oss.str().c_str(), 90 , 0, 9);
        TH1D* m_mixed_target = new TH1D("m_mixed_target","m_mixed_target", 90 , 0, 9);

	oss.str("");
        oss<<"D0 Occupancy (For "<<n_runs<<" runs between : "<<run_id_low<<" - "<<run_id_high<<")";
	TH1D* h_occuD1 = new TH1D("h_occuD1", oss.str().c_str(), 400, 0, 2000); 
	
	int bins = 6;
	int max_tm = 6;
	oss.str("");
        oss<<"Positive Track Multiplicity (For "<<n_runs<<" runs between : "<<run_id_low<<" - "<<run_id_high<<")";
        TH1D* h_tm_pl = new TH1D("h_tm_pl", oss.str().c_str(), bins, 0, max_tm);
	TH1D* h_tm_pm = new TH1D("h_tm_pm", oss.str().c_str(), bins, 0, max_tm);
	TH1D* h_tm_ph = new TH1D("h_tm_ph", oss.str().c_str(), bins, 0, max_tm);

	oss.str("");
        oss<<"Negative Track Multiplicity (For "<<n_runs<<" runs between : "<<run_id_low<<" - "<<run_id_high<<")";
        TH1D* h_tm_nl = new TH1D("h_tm_nl", oss.str().c_str(), bins, 0, max_tm);
        TH1D* h_tm_nm = new TH1D("h_tm_nm", oss.str().c_str(), bins, 0, max_tm);
        TH1D* h_tm_nh = new TH1D("h_tm_nh", oss.str().c_str(), bins, 0, max_tm);

	TCanvas* c1 = new TCanvas("c1", "");
	c1->SetGrid();

	for (auto it = list_in.begin(); it != list_in.end(); it++) {
		std::cout << "Input = " << *it << std::endl;
    		const string fn_in = *it;	
		
		f_file = new TFile(fn_in.c_str());
		if(!f_file){
                	cout << "Cannot get the DST tree.  Abort." << endl;
        		exit(1);
    		}
		
		_sorted = (TTree*) f_file->Get("save_sorted");
		_mixed = (TTree*) f_file->Get("save_mix");	
	    	_ana = (TTree*) f_file->Get("ana_tree");

		SRecEvent* sorted_event = new SRecEvent();
		SRecEvent* mixed_event = new SRecEvent();
		int occuD1=0;
		std::vector<SRecTrack> * pos_tracks =0;
		std::vector<SRecTrack> * neg_tracks =0;

		_sorted->SetBranchAddress("recEvent", &sorted_event);
		_ana->SetBranchAddress("occuD1", &occuD1);
		_ana->SetBranchAddress("pos_tracks", &pos_tracks);
  		_ana->SetBranchAddress("neg_tracks", &neg_tracks);
		_mixed->SetBranchAddress("recEvent", &mixed_event);
		
		int n_sorted = _sorted ->GetEntries();
		int n_mixed = _mixed ->GetEntries();
		int n_ana = _ana ->GetEntries();
		
		for (int k= 0; k < n_ana; k++){

                        _ana->GetEntry(k);
			
			h_occuD1->Fill(occuD1);

                        //positive track multiplicity
                        if (occuD1 <=150) h_tm_pl->Fill(pos_tracks->size());
                        else if (150 < occuD1 && occuD1 <=250) h_tm_pm->Fill(pos_tracks->size());
                        else if (occuD1 >250) h_tm_ph->Fill(pos_tracks->size());

                        //negative track multuolicity
                        if (occuD1 <=150) h_tm_nl->Fill(neg_tracks->size());
                        else if (150 < occuD1 && occuD1 <=250) h_tm_nm->Fill(neg_tracks->size());
                        else if (occuD1 >250) h_tm_nh->Fill(neg_tracks->size());
		}

		for (int i= 0; i < n_sorted; i++){
			
			_sorted->GetEntry(i);
			std::cout<<"entry_"<<i<<std::endl;	
			for (int n_dims=0; n_dims <sorted_event->getNDimuons(); n_dims++){
				
				SRecDimuon s_dim = sorted_event->getDimuon(n_dims);
				double mass = s_dim.get_mass();
				
				m_sorted->Fill(mass);
				
				SRecTrack* trk_pos = &(sorted_event->getTrack(s_dim.get_track_id_pos()));
				std::cout<<"checking Ncols()"<<std::endl;
				if (trk_pos->getStateVector(0).GetNcols()==0) continue;
				std::cout<<"Ncols for pos track is nonzero"<<std::endl;
				SRecTrack* trk_neg = &(sorted_event->getTrack(s_dim.get_track_id_neg()));
				if (trk_neg->getStateVector(0).GetNcols()==0) continue;

				//if (!trk_pos->isValid() || !trk_neg->isValid()) continue;
				if (trk_pos->get_pos_vtx().Z() < -690 || trk_neg->get_pos_vtx().Z() < -690) continue;

				int pos_chisq_t = trk_pos->getChisqTarget();
				int pos_chisq_d = trk_pos->getChisqDump();
				int pos_chisq_us = trk_pos->get_chsiq_upstream();
				
				int neg_chisq_t = trk_neg->getChisqTarget();
                                int neg_chisq_d = trk_neg->getChisqDump();
                                int neg_chisq_us = trk_neg->get_chsiq_upstream();	
				
				bool pos_ok = pos_chisq_t >=0 && pos_chisq_d >=0 && pos_chisq_us >=0 && 
					pos_chisq_t - pos_chisq_d <0 && pos_chisq_t - pos_chisq_us <0;
				bool neg_ok = neg_chisq_t >=0 && neg_chisq_d >=0 && neg_chisq_us >=0 && 
					neg_chisq_t - neg_chisq_d <0 && neg_chisq_t - neg_chisq_us <0;

				if (pos_ok && neg_ok){
					m_sorted_target->Fill(mass);
				}
			}
		}
		
		for (int j= 0; j < n_mixed; j++){
			
			_mixed->GetEntry(j);
			for (int n_dims=0; n_dims <mixed_event->getNDimuons(); n_dims++){
				SRecDimuon m_dim = mixed_event->getDimuon(n_dims);
				m_mixed->Fill(m_dim.get_mass());
				
				//getting target like dimuons
                                SRecTrack* m_trk_pos = &(sorted_event->getTrack(m_dim.get_track_id_pos()));
                                SRecTrack* m_trk_neg = &(sorted_event->getTrack(m_dim.get_track_id_neg()));

                                if (m_trk_pos->get_pos_vtx().Z() < -690 || m_trk_neg->get_pos_vtx().Z() < -690) continue;

                                int pos_chisq_t = m_trk_pos->getChisqTarget();
                                int pos_chisq_d = m_trk_pos->getChisqDump();
                                int pos_chisq_us = m_trk_pos->get_chsiq_upstream();

                                int neg_chisq_t = m_trk_neg->getChisqTarget();
                                int neg_chisq_d = m_trk_neg->getChisqDump();
                                int neg_chisq_us = m_trk_neg->get_chsiq_upstream();

                                bool pos_ok = pos_chisq_t >=0 && pos_chisq_d >=0 && pos_chisq_us >=0 &&
                                        pos_chisq_t - pos_chisq_d <0 && pos_chisq_t - pos_chisq_us <0;
                                bool neg_ok = neg_chisq_t >=0 && neg_chisq_d >=0 && neg_chisq_us >=0 &&
                                        neg_chisq_t - neg_chisq_d <0 && neg_chisq_t - neg_chisq_us <0;

                                if (pos_ok && neg_ok){
                                        m_mixed_target->Fill(m_dim.get_mass());
                                }
			}
		}
	}
	/////
	//plot All dimuon mass spectra
	/////
	c1->Clear();
	m_sorted->SetFillColorAlpha(kBlue-7, 0.35);
	m_sorted->GetYaxis()->SetRangeUser(0,m_sorted->GetMaximum()+50);
	m_sorted->GetYaxis()->SetTitle("Yield");
	m_sorted->GetXaxis()->SetTitle("Mass (GeV)");
	m_sorted->Draw("HIST");

	c1->SaveAs("plots/mass_sorted.png");

	TH1D *m_signal = new TH1D(*m_sorted); // TH1F or TH1D, same as h_hit_hitsperchannel1
	m_signal->SetNameTitle("m_signal", "Mass Spectrum for Signal Dimuons;Mass (GeV);Yield");
	if (!(m_signal->GetSumw2N() > 0)) m_signal->Sumw2(kTRUE); // ensure proper error propagation
	m_signal->Add(m_mixed, -1.0);
	m_signal->SetFillColorAlpha(kGreen-9, 0.85);
	m_signal->Draw("SAME HISTE");
	
	//m_mixed->SetMarkerStyle(kFullCircle);
        //m_mixed->SetMarkerColor(kRed);
        m_mixed->SetFillColorAlpha(kRed-9, 0.75);
        //m_mixed->GetYaxis()->SetRangeUser(0,60);
        m_mixed->Draw("SAME HIST");

	auto legend = new TLegend(0.55,0.6,0.85,0.8);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(m_sorted,"m_sorted (Measured)","f");
   	legend->AddEntry(m_mixed,"m_mixed (Combinatoric)","f");
   	legend->AddEntry(m_signal,"Measured - Combinatoric","f");
   	legend->Draw();

	c1->SaveAs("plots/dim_mass.png");
	
	////
	//plot D0 occupancy
	////
	//
	c1->Clear();
	h_occuD1->SetFillColorAlpha(kBlue-7, 0.35);
	h_occuD1->GetYaxis()->SetTitle("Yield");
        h_occuD1->GetXaxis()->SetTitle("D0 Occupancy");
	h_occuD1->Draw("HIST");
	c1->SaveAs("plots/occuD1.png");
	
	////
	//plot target-like dimuon mass spectra
	////
	c1->Clear();
        m_sorted_target->SetFillColorAlpha(kBlue-7, 0.35);
        m_sorted_target->GetYaxis()->SetRangeUser(0,m_sorted_target->GetMaximum()+50);
        m_sorted_target->GetYaxis()->SetTitle("Yield");
        m_sorted_target->GetXaxis()->SetTitle("Mass (GeV)");
        m_sorted_target->Draw("HIST");

        TH1D *m_signal_target = new TH1D(*m_sorted_target); // TH1F or TH1D, same as h_hit_hitsperchannel1
        m_signal_target->SetNameTitle("m_signal_target", "Mass Spectrum for Signal Dimuons;Mass (GeV);Yield");
        if (!(m_signal_target->GetSumw2N() > 0)) m_signal_target->Sumw2(kTRUE); // ensure proper error propagation
        m_signal_target->Add(m_mixed_target, -1.0);
        m_signal_target->SetFillColorAlpha(kGreen-9, 0.85);
        m_signal_target->Draw("SAME HISTE");

        //m_mixed->SetMarkerStyle(kFullCircle);
        //m_mixed->SetMarkerColor(kRed);
        m_mixed_target->SetFillColorAlpha(kRed-9, 0.75);
        //m_mixed->GetYaxis()->SetRangeUser(0,60);
        m_mixed_target->Draw("SAME HIST");

        auto legend_tar = new TLegend(0.55,0.6,0.85,0.8);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend_tar->AddEntry(m_sorted_target,"m_sorted (Measured)","f");
        legend_tar->AddEntry(m_mixed_target,"m_mixed (Combinatoric)","f");
        legend_tar->AddEntry(m_signal_target,"Measured - Combinatoric","f");
        legend_tar->Draw();

        c1->SaveAs("plots/targetdim_mass.png");

	////
	//plot track muliplicity
	////
	TCanvas* c2 = new TCanvas("c2", "");
        c2->SetGrid();
	c2->SetLogy();

	h_tm_pl->SetLineColor(kBlue-9);
	h_tm_pl->GetYaxis()->SetTitle("Yield");
	h_tm_pl->GetYaxis()->SetRangeUser(0.2,h_tm_ph->GetMaximum()+200000);
	h_tm_pl->GetXaxis()->SetTitle("Number of Tracks / Event");
	h_tm_pl->SetLineWidth(3);
	h_tm_pl->Draw();
	h_tm_pm->SetLineColor(kRed-9);
	h_tm_pm->SetLineWidth(3);
	h_tm_pm->Draw("SAME");
	h_tm_ph->SetLineColor(kGreen-9);
	h_tm_ph->SetLineWidth(3);
        h_tm_ph->Draw("SAME");

	auto legend_p = new TLegend(0.55,0.6,0.85,0.8);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend_p->AddEntry(h_tm_pl,"0 #leq D0 #leq 150","l");
        legend_p->AddEntry(h_tm_pm,"150 < D0 #leq 250","l");
        legend_p->AddEntry(h_tm_ph,"250 < D0","l");
        legend_p->Draw();

	c2->SaveAs("plots/pos_track_mul.png");

	c2->Clear();
        h_tm_nl->SetLineColor(kBlue-9);
        h_tm_nl->GetYaxis()->SetTitle("Yield");
	h_tm_nl->GetYaxis()->SetRangeUser(0.2,h_tm_nh->GetMaximum()+200000);
        h_tm_nl->GetXaxis()->SetTitle("Number of Tracks / Event");
	h_tm_nl->SetLineWidth(3);
        h_tm_nl->Draw();
        h_tm_nm->SetLineColor(kRed-9);
        h_tm_nm->SetLineWidth(3);
	h_tm_nm->Draw("SAME");
        h_tm_nh->SetLineColor(kGreen-9);
        h_tm_nh->SetLineWidth(3);
	h_tm_nh->Draw("SAME");

	auto legend_n = new TLegend(0.55,0.6,0.85,0.8);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend_n->AddEntry(h_tm_nl,"0 #leq D0 #leq 150","l");
        legend_n->AddEntry(h_tm_nm,"150 < D0 #leq 250","l");
        legend_n->AddEntry(h_tm_nh,"250 < D0","l");
        legend_n->Draw();
        c2->SaveAs("plots/neg_track_mul.png");
	
		
}
