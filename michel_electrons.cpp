#include <string>
#include <vector> 
#include <iostream> 
#include <fstream>


// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

// ANNIE plot style 
#include "rootlogon.hpp"


int TEST_RUN_DATA = -1;





void MakeRatioPlot(TH1D *data_hist, TH1D* mc_hist, std::string fig_name){

    double max_data = data_hist->GetMaximum();
    double max_mc = mc_hist->GetMaximum();
    double y_max = std::max(max_data, max_mc);

    TCanvas *c = new TCanvas("c", "MC/Data Comparison", 800,800);
    c->Divide(1,2);

    // Create two pads, one for histograms other for the ratio
    TPad *pad1 = new TPad("pad1", "Hist plot", 0, 0.32, 1, 1.0);
    pad1->SetBottomMargin(0.02); // No x-axis labels
    pad1->Draw();
    pad1->cd();


    mc_hist->SetMaximum(1.1*y_max);
    mc_hist->GetXaxis()->SetLabelSize(0);
    mc_hist->Draw("HIST");
    data_hist->Draw("E1 SAME");
    TLegend *l = new TLegend(0.73, 0.7, 0.85, 0.8);
    l->AddEntry(data_hist, "Data");
    l->AddEntry(mc_hist, "MC");
    PreliminarySide();
    l->Draw("SAME");

    c->cd();
    TPad *pad2 = new TPad("pad2", "Ratio plot", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    TH1D *ratio = (TH1D*)mc_hist->Clone("ratio");
    ratio->Divide(data_hist);
    ratio->SetStats(0);
    ratio->SetTitle("");
    ratio->SetLineColor(kBlack);

    ratio->GetYaxis()->SetTitle("MC / Data");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.12);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetLabelSize(0.10);

    ratio->SetMinimum(0.0);
    ratio->SetMaximum(2.0);
    ratio->Draw("HIST");

    // Draw a line at y=1
    TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1.0,ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->SaveAs(fig_name.c_str());

    // Recall that objects created with new don't go away after the scope ends!! 
    delete c; 
}




void read_data(std::string file_list, std::vector<float>& cluster_time, std::vector<float>& cluster_pe,
    std::vector<float>& cluster_qb,std::vector<std::vector<double>>& hitID, 
    std::vector<std::vector<double>>& hitPE, std::vector<std::vector<double>>& hiT){

    /*

    This function reads data n-tuples and fills
    vectors with relevant information 

    */

    std::cout << "Reading data n-tuples..." << std::endl;
    if(TEST_RUN_DATA > 0){
    std::cout << "Test run, iterating over " << TEST_RUN_DATA << " files" << std::endl;

    }

    std::ifstream listFile(file_list);
    if(!listFile){
    std::cerr << "Error: could not open: " << file_list << std::endl;
    }



    TChain *data_tree = new TChain("data");
    // Iterate over files 
    std::string filePath;
    int file_counter = 0;
    while(std::getline(listFile, filePath)){

    if(filePath.empty()) continue;
    if(file_counter == TEST_RUN_DATA){
    break;
    }
    data_tree->Add(filePath.c_str());
    file_counter+=1;
    }

    data_tree->SetBranchStatus("*", 0);

    data_tree->SetBranchStatus("run_number", 1);
    data_tree->SetBranchStatus("number_of_clusters",1);
    data_tree->SetBranchStatus("event_number", 1);
    data_tree->SetBranchStatus("cluster_Qb", 1);
    data_tree->SetBranchStatus("cluster_time", 1); 
    data_tree->SetBranchStatus("cluster_PE", 1);
    data_tree->SetBranchStatus("cluster_Number", 1);
    data_tree->SetBranchStatus("cluster_Hits", 1);
    data_tree->SetBranchStatus("isBrightest", 1);
    data_tree->SetBranchStatus("TankMRDCoinc", 1);
    data_tree->SetBranchStatus("MRD_activity", 1);
    data_tree->SetBranchStatus("hadExtended", 1);
    data_tree->SetBranchStatus("NoVeto", 1);
    data_tree->SetBranchStatus("hitZ", 1);
    data_tree->SetBranchStatus("hitPE", 1);
    data_tree->SetBranchStatus("hitT", 1);
    data_tree->SetBranchStatus("hitID", 1);


    float ct2;
    float cpe1;
    float cb1;
    int   ch1;
    Long64_t    cn1;
    Long64_t    noc1;
    Long64_t    evn1; 
    int    rn1; 
    int    b1;
    int    tmrd1;
    int    mrd_yes;
    int    t1;
    int    nv1;
    int    had_ext;
    std::vector<double>  *hz1=nullptr;
    std::vector<double>  *hpe1=nullptr;
    std::vector<double>  *ht1=nullptr;
    std::vector<double>  *hid1=nullptr;


    data_tree->SetBranchAddress("run_number", &rn1);
    data_tree->SetBranchAddress("number_of_clusters", &noc1);
    data_tree->SetBranchAddress("event_number", &evn1);
    data_tree->SetBranchAddress("cluster_Qb", &cb1);
    data_tree->SetBranchAddress("cluster_time", &ct2); 
    data_tree->SetBranchAddress("cluster_PE", &cpe1);
    data_tree->SetBranchAddress("cluster_Number", &cn1);
    data_tree->SetBranchAddress("cluster_Hits", &ch1);
    data_tree->SetBranchAddress("isBrightest", &b1);
    data_tree->SetBranchAddress("TankMRDCoinc", &tmrd1);
    data_tree->SetBranchAddress("MRD_activity", &mrd_yes);
    data_tree->SetBranchAddress("hadExtended", &had_ext);
    data_tree->SetBranchAddress("NoVeto", &nv1);
    data_tree->SetBranchAddress("hitZ", &hz1);
    data_tree->SetBranchAddress("hitPE", &hpe1);
    data_tree->SetBranchAddress("hitT", &ht1);
    data_tree->SetBranchAddress("hitID", &hid1);


    int nEvents = data_tree->GetEntries();
    std::cout << "Number of entries in data chain: " << data_tree->GetEntries() << std::endl;


    int selected_muons = 0;
    int selected_michels = 0; 

    for(int i=0; i < nEvents; i++){
        if(i%100000==0){std::cout << "Processing entry: " << i << std::endl;}
        data_tree->GetEntry(i);


        double *hz1_v1 = hz1->data();
        double *hpe1_v1 = hpe1->data();
        double *ht1_v1 = ht1->data();
        double *hid1_v1 = hid1->data();
        float a = 0;
        for(int j=0; j<ch1; j++){
        a+=hz1_v1[j]*hpe1_v1[j];
        }

        bool dirt_muon = (tmrd1 == 0) &&  // No tank-MRD coincidence
                        (nv1 == 0)   &&  // Vetoed (FMV yes)
                        (b1 == 1)    &&  // Brightest object
                        (noc1 > 1) && // More than one cluster
                        (had_ext != 0) && // Extended window to look for Michels
                        (ch1 >= 50) && // at least 50 pmt hits
                        (cpe1 >= 1000 && cpe1 <= 4000) && // 1000 < cluster PE < 4000 (why (?))
                        (cb1 > 0. && cb1 < 0.2)    && // Cluster charge balance < 0.2
                        (ct2 >= 200   && ct2 <= 1800) && // Inside spill time
                        (a>0.);   // Charge barycenter downstream 

        if(dirt_muon){

            float parent_time = ct2;
            Long64_t parent_evn = evn1;
            int parent_run = rn1;
            selected_muons+=1;


            // Loop again to find michels

            bool temp_michel = false;
            for(int k = i+1; k < nEvents; k++){
                data_tree->GetEntry(k);

                // Break if no longer same event or if we found a Michel candidate
                if(evn1 != parent_evn || rn1!= parent_run){break;}
            
                    
                    

                float adj_time = ct2 - parent_time;
                //std::cout << "Adjusted time: " << adj_time << std::endl;
                //std::cout << cb2 << " " << cpe2 << std::endl;
                bool michel_electron = (adj_time >= 200. && adj_time <= 5000) &&
                                        (cpe1 > 0. && cpe1 < 650.) &&
                                        (ch1 >= 20) &&
                                        (cb1 > 0. && cb1 < 0.2);
                if(michel_electron){


                    temp_michel=true;
                    selected_michels+=1;
                    cluster_time.push_back(adj_time);
                    cluster_pe.push_back(cpe1);
                    cluster_qb.push_back(cb1);
                    //hitID.push_back(*hid1);
                    //hitPE.push_back(*hpe1);
                    //hiT.push_back(*ht1);
                    break;
                }

            }
                
    
                
            //std::cout << "--------------------" << std::endl;
            //delete sub_tree;
        }
        else{continue;}

    }

    std::cout << "Number of selected dirt muons in data: " << selected_muons << std::endl;
    std::cout << "Number of selected Michels in data: " << selected_michels << std::endl;

}


void read_michel_mc(std::string file_name, std::vector<double> &cluster_pe_mc,
                    std::vector<double> &cluster_qb_mc, std::vector<std::vector<double>> &hitID_mc,
                    std::vector<std::vector<double>> &hitPE_mc, std::vector<double> &true_vx,
                    std::vector<double> &true_vy, std::vector<double> &true_vz, std::vector<double> &true_E){

    std::cout << "Processing Michel electron MC sample..." << std::endl;
    TFile *f = TFile::Open(file_name.c_str(), "READ"); 
    TTree *beam_cluster_tree = (TTree*) f->Get("phaseIITankClusterTree");
    TTree *truth_tree = (TTree*) f->Get("phaseIITriggerTree");

    std::cout << "Number of entries in tank tree: " << beam_cluster_tree->GetEntries() << std::endl;
    std::cout << "Number of entries in truth tree: " << truth_tree->GetEntries() << std::endl;

    double cpe_mc;
    double cb_mc;
    int    ev_mc; 
    UInt_t    ch_mc;
    std::vector<double>  *hpe_mc=nullptr;
    std::vector<double>  *hid_mc=nullptr;

    beam_cluster_tree->SetBranchAddress("eventNumber", &ev_mc);
    beam_cluster_tree->SetBranchAddress("clusterPE", &cpe_mc);
    beam_cluster_tree->SetBranchAddress("clusterChargeBalance", &cb_mc);
    beam_cluster_tree->SetBranchAddress("clusterHits", &ch_mc);
    beam_cluster_tree->SetBranchAddress("hitPE", &hpe_mc);
    beam_cluster_tree->SetBranchAddress("hitDetID", &hid_mc);

    double vx_mc;
    double vy_mc;
    double vz_mc;
    double E_mc;

    truth_tree->SetBranchAddress("trueVtxX", &vx_mc);
    truth_tree->SetBranchAddress("trueVtxX", &vy_mc);
    truth_tree->SetBranchAddress("trueVtxX", &vz_mc);
    truth_tree->SetBranchAddress("trueMuonEnergy", &E_mc);


    // Iterate and apply selection cuts
    int n_selected = 0; 
    double corrective_factor = 1.2;
    for(int i=0; i < beam_cluster_tree->GetEntries(); i++){
        beam_cluster_tree->GetEntry(i);
        truth_tree->GetEntry(ev_mc);

        double vz_plane = TMath::Sqrt(TMath::Power(vx_mc/100.,2) + TMath::Power(vz_mc/100.,2));
        double vy_plane = TMath::Abs(vy_mc/100.); 
        
        bool select_michel = (vz_plane < 1.1) && (vy_plane < 1.5) &&
                             (cb_mc < 0.18) && (cpe_mc < 650);

        if(select_michel){
                    //std::cout << cpe1 << " " << cb1 <<" " << ct2 << std::endl;
    

            double *hpe1_v1 = hpe_mc->data();
            float a = 0;
            for(int j=0; j<ch_mc; j++){
                //std::cout << hz1_v1[j]<< " " << hpe1_v1[j] << std::endl;
                a+=hpe1_v1[j]*corrective_factor;
            }
            //std::cout << ch1 << std::endl;

        




            cluster_pe_mc.push_back(a);
            true_vx.push_back(vx_mc);
            true_vy.push_back(vy_mc);
            true_vz.push_back(vz_mc);
            true_E.push_back(E_mc);
            n_selected+=1;
        }

    }
    std::cout << "Number of selected MC michel electrons: " << n_selected << std::endl;


}

void michel_electrons(){
    rootlogon();
    gStyle->SetOptStat(0);

    std::cout << "Running michel electrons selection..." << std::endl;


    std::vector<double> cluster_pe_mc;
    std::vector<double> cluster_qb_mc;
    std::vector<std::vector<double>> hitID_mc;
    std::vector<std::vector<double>> hitPE_mc;
    std::vector<double> true_vx;
    std::vector<double> true_vy;
    std::vector<double> true_vz;
    std::vector<double> true_E;


    std::string michel_mc_file = "/exp/annie/data/users/lmoralep/michels_full_volume_default.root";
    read_michel_mc(michel_mc_file, cluster_pe_mc, cluster_qb_mc,
                  hitID_mc, hitPE_mc, true_vx, true_vy, true_vz, true_E);


    TH1D *true_E_hist = new TH1D("True_E", " Michel electrons; Energy [MeV]; Number of entries",100,0., 60. );
    true_E_hist->SetLineColor(kRed);
    for(double E_val: true_E){
        true_E_hist->Fill(E_val);
    }

    TCanvas *c = new TCanvas("c", "trueE", 800,600);
    true_E_hist->Draw("HIST");
    SimulationSide();
    c->SaveAs("./figures/michel_electron_E_spectrum.pdf");
    delete c; 


    // Data arrays
    std::vector<float> cluster_time; 
    std::vector<float> cluster_pe;
    std::vector<float> cluster_qb;
    std::vector<std::vector<double>> hitID;
    std::vector<std::vector<double>> hitPE;
    std::vector<std::vector<double>> hitT;

    // Read data
    std::string data_list = "./samples/2023_beamdata_v1.txt";
    read_data(data_list, cluster_time, cluster_pe, cluster_qb,
        hitID, hitPE, hitT);

    TH1D *cpe_mc_hist = new TH1D("cpe_mc", " ; charge [p.e]; Number of entries", 50, 0., 650);
    cpe_mc_hist->SetLineColor(kRed);
    
    int mc_counter = 0.;
    std::cout << "Cluster pe size: " << cluster_pe.size() << std::endl;
    for(double cpe_mc_val: cluster_pe_mc){
        if(mc_counter == cluster_pe.size()){
            break;
        }
        cpe_mc_hist->Fill(cpe_mc_val);
        mc_counter+=1;
    }

    TH1D *cpe_data_hist = new TH1D("cpe_data", " ; charge[p.e]; Number of entries", 50, 0., 650.);
    for(float cpe_data_val: cluster_pe){
        cpe_data_hist->Fill(cpe_data_val);
    }

    MakeRatioPlot(cpe_data_hist, cpe_mc_hist, "./figures/cpe_michel_electrons.pdf");

    // Save to a new file (optional)
    TFile* outFile = new TFile("michel_electron_selection.root", "RECREATE");
    cpe_mc_hist->Write();
    cpe_data_hist->Write();
    outFile->Close();

    


}