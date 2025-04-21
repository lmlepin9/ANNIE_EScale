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
    for(int i=0; i < beam_cluster_tree->GetEntries(); i++){
        beam_cluster_tree->GetEntry(i);
        truth_tree->GetEntry(ev_mc);

        double vz_plane = TMath::Sqrt(TMath::Power(vx_mc/100.,2) + TMath::Power(vz_mc/100.,2));
        double vy_plane = TMath::Abs(vy_mc/100.); 
        
        bool select_michel = (vz_plane < 1.1) && (vy_plane < 1.5) &&
                             (cb_mc < 0.18) && (cpe_mc < 650);

        if(select_michel){
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


  

}