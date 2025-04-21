#include <string>
#include <vector> 
#include <iostream> 
#include <fstream>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "rootlogon.hpp"


void michel_electrons(std::string file_name){

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

    std::string data_list = "2023_beamdata_v1.txt";
    std::string mc_list = "mc_thru_muons.list";
    
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
    c1->SaveAs("./figures/charge_pe_hist.pdf");


    // Generate Data/MC comparison plot
    std::string cpe_fig_name = "./figures/charge_pe_hist_ratio.pdf";
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
    c2->SaveAs("./figures/charge_balance_hist.pdf");

    std::string cb_fig_name = "./figures/charge_cb_hist_ratio.pdf";
    MakeRatioPlot(cb_hist, cb_hist_mc, cb_fig_name);


}