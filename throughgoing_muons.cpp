#include <string>
#include <vector> 
#include <iostream> 
#include <fstream>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStopwatch.h"

// ANNIE plot style
#include "rootlogon.hpp"





int TEST_RUN = -1;
int TEST_RUN_MC = -1;

void read_data(std::string file_list, std::vector<float>& cluster_time, std::vector<float>& cluster_pe,
               std::vector<float>& cluster_qb,std::vector<std::vector<double>>& hitID, 
               std::vector<std::vector<double>>& hitPE, std::vector<std::vector<double>>& hiT,
               std::vector<int>& tru){

    /*
    
    This function reads data n-tuples and fills
    vectors with relevant information 
    
    */

    std::cout << "Reading data n-tuples..." << std::endl;
    if(TEST_RUN > 0){
        std::cout << "Test run, iterating over " << TEST_RUN << " files" << std::endl;

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
        if(file_counter == TEST_RUN){
            break;
        }
        data_tree->Add(filePath.c_str());
        file_counter+=1;
    }

    /*
    // We disable all branches
    data_tree->SetBranchStatus("*",0);

    // Enable only the following branches
    data_tree->SetBranchStatus("cluster_Qb", 1);
    data_tree->SetBranchStatus("cluster_time_BRF", 1); 
    data_tree->SetBranchStatus("cluster_PE", 1);
    data_tree->SetBranchStatus("cluster_Hits", 1);
    data_tree->SetBranchStatus("isBrightest", 1);
    data_tree->SetBranchStatus("TankMRDCoinc", 1);
    data_tree->SetBranchStatus("MRD_activity", 1);
    data_tree->SetBranchStatus("MRD_Track", 1);
    data_tree->SetBranchStatus("NoVeto", 1);
    data_tree->SetBranchStatus("hitZ", 1);
    data_tree->SetBranchStatus("hitPE", 1);
    data_tree->SetBranchStatus("hitT", 1);
    data_tree->SetBranchStatus("hitID", 1);
    data_tree->SetBranchStatus("MRDThrough", 1);
    */ 


    float ct2;
    float cpe1;
    float cb1;
    int    ch1;
    int    b1;
    int    tmrd1;
    int    mrd_yes;
    int    t1;
    int    nv1;
    std::vector<double>  *hz1=nullptr;
    std::vector<double>  *hpe1=nullptr;
    std::vector<double>  *ht1=nullptr;
    std::vector<double>  *hid1=nullptr;
    int    tg1;


    data_tree->SetBranchAddress("cluster_Qb", &cb1);
    data_tree->SetBranchAddress("cluster_time_BRF", &ct2); 
    data_tree->SetBranchAddress("cluster_PE", &cpe1);
    data_tree->SetBranchAddress("cluster_Hits", &ch1);
    data_tree->SetBranchAddress("isBrightest", &b1);
    data_tree->SetBranchAddress("TankMRDCoinc", &tmrd1);
    data_tree->SetBranchAddress("MRD_activity", &mrd_yes);
    data_tree->SetBranchAddress("MRD_Track", &t1);
    data_tree->SetBranchAddress("NoVeto", &nv1);
    data_tree->SetBranchAddress("hitZ", &hz1);
    data_tree->SetBranchAddress("hitPE", &hpe1);
    data_tree->SetBranchAddress("hitT", &ht1);
    data_tree->SetBranchAddress("hitID", &hid1);
    data_tree->SetBranchAddress("MRDThrough", &tg1);

    int nEvents = data_tree->GetEntries();
    std::cout << "Number of entries in data chain: " << data_tree->GetEntries() << std::endl;


    int selected_thru = 0; 

    for(int i=0; i < nEvents; i++){
        data_tree->GetEntry(i);
        //std::cout << cpe1 << " " << cb1 <<" " << ct2 << std::endl;
        
        double *hz1_v1 = hz1->data();
        double *hpe1_v1 = hpe1->data();
        double *ht1_v1 = ht1->data();
        double *hid1_v1 = hid1->data();
        float a = 0;
        for(int j=0; j<ch1; j++){
            //std::cout << hz1_v1[j]<< " " << hpe1_v1[j] << std::endl;
            a+=hz1_v1[j]*hpe1_v1[j];
        }
        //std::cout << ch1 << std::endl;


        bool througgoing = (t1 == 1)    &&  // 1 track
                           (tmrd1 == 1) &&  // Tank-MRD coincidence
                           (nv1 == 0)   &&  // Vetoed (FMV yes)
                           (b1 == 1)    &&  // Brightest object
                           (cpe1 >= 1000 && cpe1 <= 6000) && // 1000 < cluster PE < 6000 (why (?))
                           (cb1 > 0. && cb1 < 0.2)    && // Cluster charge balance < 0.2
                           (ct2 >= 200   && ct2 <= 1800) && // Inside spill time
                           (a>0.);   // Charge barycenter downstream 

        if(througgoing){
            //std::cout << cpe1 << std::endl;
            cluster_time.push_back(ct2);
            cluster_pe.push_back(cpe1);
            cluster_qb.push_back(cb1);
            hitID.push_back(*hid1);
            hitPE.push_back(*hpe1);
            hiT.push_back(*ht1);
            selected_thru+=1;

        }
   }
        
   std::cout << "Number of selected throughgoing muons in data: " << selected_thru << std::endl;

}

void read_MC(std::string file_list, std::vector<double>& cluster_time, std::vector<double>& cluster_pe,
    std::vector<double>& cluster_qb,std::vector<std::vector<double>>& hitID, 
    std::vector<std::vector<double>>& hitPE, std::vector<std::vector<double>>& hiT,
    std::vector<int>& tru){


    /*
    
    This function reads MC n-tuples and fills
    vectors with relevant information 
    
    */


    std::cout << "Reading MC n-tuples..." << std::endl;
    if(TEST_RUN_MC > 0){
        std::cout << "Test run, iterating over " << TEST_RUN_MC << " files" << std::endl;

    }
    
    std::ifstream listFile(file_list);
    if(!listFile){
        std::cerr << "Error: could not open: " << file_list << std::endl;
    }




    // Iterate over files 
    std::string filePath;
    int file_counter = 0;
    int selected_thru_mc = 0; 
    while(std::getline(listFile, filePath)){

        if(filePath.empty()) continue;
        if(file_counter == TEST_RUN_MC){
            break;
        }
        //std::cout << "Reading: " << filePath << std::endl;
        TFile *f = TFile::Open(filePath.c_str(), "READ");
        TTree *beam_cluster_tree = (TTree*) f->Get("phaseIITankClusterTree");


        // If tree is not present we skip to next file
        if(!beam_cluster_tree){
            //std::cout << "PhaseIITankClusterTree not present in: " << filePath << std::endl;
            //std::cout << "Skipping to next file..." << std::endl;
            continue;
        }

        // Disable all branches
        beam_cluster_tree->SetBranchStatus("*", 0);

        // Enable only the following branches: 

        beam_cluster_tree->SetBranchStatus("eventNumber", 1);
        beam_cluster_tree->SetBranchStatus("clusterPE", 1);
        beam_cluster_tree->SetBranchStatus("clusterChargeBalance", 1);
        beam_cluster_tree->SetBranchStatus("clusterTime", 1);
        beam_cluster_tree->SetBranchStatus("hitX", 1);
        beam_cluster_tree->SetBranchStatus("hitY", 1);
        beam_cluster_tree->SetBranchStatus("hitZ", 1);
        beam_cluster_tree->SetBranchStatus("hitT", 1);
        beam_cluster_tree->SetBranchStatus("hitPE", 1);
        beam_cluster_tree->SetBranchStatus("hitDetID", 1);
        beam_cluster_tree->SetBranchStatus("TankMRDCoinc", 1);
        beam_cluster_tree->SetBranchStatus("NoVeto",1);
            
        double ct_mc;
        double cpe_mc;
        double cb_mc;
        int    ev_mc; 
        int    ch_mc;
        int    b_mc;
        int    tmrd_mc;
        int    mrd_yes_mc;
        int    nv_mc;
        std::vector<double>  *hx_mc=nullptr;
        std::vector<double>  *hy_mc=nullptr;
        std::vector<double>  *hz_mc=nullptr;
        std::vector<double>  *hpe_mc=nullptr;
        std::vector<double>  *ht_mc=nullptr;
        std::vector<double>  *hid_mc=nullptr;

        beam_cluster_tree->SetBranchAddress("eventNumber", &ev_mc);
        beam_cluster_tree->SetBranchAddress("clusterPE", &cpe_mc);
        beam_cluster_tree->SetBranchAddress("clusterChargeBalance", &cb_mc);
        beam_cluster_tree->SetBranchAddress("clusterTime", &ct_mc);
        beam_cluster_tree->SetBranchAddress("hitX", &hx_mc);
        beam_cluster_tree->SetBranchAddress("hitY", &hy_mc);
        beam_cluster_tree->SetBranchAddress("hitZ", &hz_mc);
        beam_cluster_tree->SetBranchAddress("hitT", &ht_mc);
        beam_cluster_tree->SetBranchAddress("hitPE", &hpe_mc);
        beam_cluster_tree->SetBranchAddress("hitDetID", &hid_mc);
        beam_cluster_tree->SetBranchAddress("TankMRDCoinc", &tmrd_mc);
        beam_cluster_tree->SetBranchAddress("NoVeto", &nv_mc);

        TTree *truth_tree = (TTree*) f->Get("phaseIITriggerTree");

        // If tree is not present we skip to next file
        if(!truth_tree){
            //std::cout << "PhaseIITriggerTree not present in: " << filePath << std::endl;
            //std::cout << "Skipping to next file..." << std::endl;
            continue;
        }

        // Disable all branches:
        truth_tree->SetBranchStatus("*", 0);

        // Enable only these branches:
        truth_tree->SetBranchStatus("numMRDTracks", 1);
        truth_tree->SetBranchStatus("MRDThrough", 1);


        int t_mc; 
        std::vector<int> *tg_mc=nullptr; 

        truth_tree->SetBranchAddress("numMRDTracks", &t_mc);
        truth_tree->SetBranchAddress("MRDThrough", &tg_mc);

        int nEvents = beam_cluster_tree->GetEntries();
        int nTruthEvs = truth_tree->GetEntries();



        if(file_counter%500==0){
            std::cout << "Processing file number: " << file_counter << std::endl;
        }

        for(int i=0; i < nEvents; i++){
            beam_cluster_tree->GetEntry(i);
            truth_tree->GetEntry(ev_mc);



            double *hz1_v1 = hz_mc->data();
            double *hpe1_v1 = hpe_mc->data();
            double *ht1_v1 = ht_mc->data();
            double *hid1_v1 = hid_mc->data();

            //std::cout << hz_mc->size() << std::endl;
        
            float a = 0;
            for(int j=0; j< hz_mc->size() ; j++){
                a+=hz1_v1[j]*hpe1_v1[j];
            }

            //std::cout << cpe_mc << " " << cb_mc << " " << ct_mc << std::endl;

            bool througgoing_mc = (t_mc == 1)    &&  // 1 track
                               (tmrd_mc == 1) &&  // Tank-MRD coincidence
                               (nv_mc == 0)   &&  // Vetoed (FMV yes)
                               (cpe_mc >= 1000 && cpe_mc <= 6000) && // 1000 < cluster PE < 6000 (why (?))
                               (cb_mc > 0. && cb_mc < 0.2)    && // Cluster charge balance < 0.2
                               (a>0.);   // Charge barycenter downstream 

            if(througgoing_mc){
            //std::cout << cpe1 << std::endl;
            cluster_time.push_back(ct_mc);
            cluster_pe.push_back(cpe_mc);
            cluster_qb.push_back(cb_mc);
            hitID.push_back(*hid_mc);
            hitPE.push_back(*hpe_mc);
            hiT.push_back(*ht_mc);
            selected_thru_mc+=1;

            }
    

        }



        // We only count a file if it we processed it 
        file_counter+=1;


        // Close file and deletes owned objects
        f->Close();
        //delete beam_cluster_tree;
        //delete truth_tree;
        delete f;
    }

    std::cout << "Number of selected throughgoing muons in MC: " << selected_thru_mc << std::endl;

}


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



void throughgoing_muons(){

    TStopwatch timer;
    timer.Start();

    rootlogon();
    gStyle->SetOptStat(0);

    // Data arrays
    std::vector<float> cluster_time; 
    std::vector<float> cluster_pe;
    std::vector<float> cluster_qb;
    std::vector<std::vector<double>> hitID;
    std::vector<std::vector<double>> hitPE;
    std::vector<std::vector<double>> hitT;
    std::vector<int> thru;

    std::vector<double> cluster_time_mc; 
    std::vector<double> cluster_pe_mc;
    std::vector<double> cluster_qb_mc;
    std::vector<std::vector<double>> hitID_mc;
    std::vector<std::vector<double>> hitPE_mc;
    std::vector<std::vector<double>> hitT_mc;
    std::vector<int> thru_mc;

    std::string data_list = "./samples/2023_beamdata_v1.txt";
    //std::string mc_list = "./samples/mc_thru_muons.list";
    std::string mc_list = "./samples/mc_world.txt";

    // Read data
    read_data(data_list, cluster_time, cluster_pe, cluster_qb,
              hitID, hitPE, hitT, thru);


    // Read MC 
    read_MC(mc_list, cluster_time_mc, cluster_pe_mc, cluster_qb_mc,
                hitID_mc, hitPE_mc, hitT_mc, thru_mc);


    // Create histogram and fill it
    TH1D *cpe_hist = new TH1D("cpe", " ", 30, 1000., 6000.);
    TH1D *cpe_hist_mc = new TH1D("cpe_mc", " ", 30, 1000. ,6000.);
    

    int mc_max = cluster_pe_mc.size();
    int data_counter = 0;

    for(float cp_val: cluster_pe){
        //std::cout << cp_val << std::endl;
        if(data_counter > mc_max){break;}
        cpe_hist->Fill(cp_val);
        data_counter+=1;
    }
    
    for(double cp_val_mc: cluster_pe_mc){
        cpe_hist_mc->Fill(cp_val_mc);
    }

    cpe_hist_mc->SetXTitle("Cluster charge [p.e]");
    cpe_hist_mc->SetYTitle("Number of entries");
    cpe_hist->SetLineColor(kBlack);
    cpe_hist_mc->SetLineColor(kRed);
    cpe_hist->SetXTitle("Cluster charge [pe]");
    cpe_hist->SetYTitle("Number of entries");
    TCanvas *c1 = new TCanvas("c1", "canvas", 800, 600);
    TLegend *l1 = new TLegend(0.73, 0.7, 0.85, 0.8);
    l1->AddEntry(cpe_hist, "Data");
    l1->AddEntry(cpe_hist_mc, "MC");
    cpe_hist_mc->Draw("HIST");
    cpe_hist->Draw("E SAME");
    l1->Draw("SAME");
    PreliminarySide();
    c1->SaveAs("./figures/charge_pe_hist_mc_world.pdf");


    // Generate Data/MC comparison plot
    std::string cpe_fig_name = "./figures/charge_pe_hist_ratio_mc_world.pdf";
    MakeRatioPlot(cpe_hist, cpe_hist_mc, cpe_fig_name);



    // Create histogram and fill it
    TH1D *cb_hist = new TH1D("cb", " ", 30, 0.07, 0.2);
    TH1D *cb_hist_mc = new TH1D("cb_mc", " ", 30, 0.07, 0.2);

    int data_counter_cb = 0;
    for(float cb_val: cluster_qb){
        if(data_counter_cb > mc_max){break;}
        cb_hist->Fill(cb_val);
        data_counter_cb+=1;
    }

    for(double cb_val_mc: cluster_qb_mc){
        cb_hist_mc->Fill(cb_val_mc);
    }

    cb_hist->SetLineColor(kBlack);
    cb_hist->SetXTitle("Cluster charge balance");
    cb_hist->SetYTitle("Number of entries");


    cb_hist_mc->SetLineColor(kRed);
    cb_hist_mc->SetXTitle("Cluster charge balance");
    cb_hist_mc->SetYTitle("Number of entries");


    TCanvas *c2 = new TCanvas("c2", "canvas", 800, 600);
    TLegend *l2 = new TLegend(0.6, 0.7, 0.8, 0.8);
    cb_hist->Draw("E");
    PreliminarySide();
    c2->SaveAs("./figures/charge_balance_hist_mc_world.pdf");

    std::string cb_fig_name = "./figures/charge_cb_hist_ratio_mc_world.pdf";
    MakeRatioPlot(cb_hist, cb_hist_mc, cb_fig_name);

    timer.Stop();
    timer.Print();


}
