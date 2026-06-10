#!/usr/bin/env python3
import os, sys, math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# --- Define here the helpers ---
def safe_get(obj, attr_name, default=0):
    """Safely retrieves attributes from NanoAOD events or collections, 
       catching the NanoAODTools RuntimeError if the branch is missing."""
    try:
        return getattr(obj, attr_name)
    except RuntimeError:
        return default
        
def deltaR(obj1, obj2):
  deta = obj1.eta - obj2.eta
  dphi = ROOT.TVector2.Phi_mpi_pi(obj1.phi - obj2.phi)
  return math.hypot(deta, dphi)
  
def get_nu_p4(lep_vec, met_pt, met_phi):
    """Reconstructs the neutrino 4-vector using the W mass constraint."""
    MW = 80.379
    px_nu = met_pt * math.cos(met_phi)
    py_nu = met_pt * math.sin(met_phi)

    Lambda = (MW**2) / 2.0 + lep_vec.Px() * px_nu + lep_vec.Py() * py_nu
    A = lep_vec.Pt()**2
    B = -2.0 * Lambda * lep_vec.Pz()
    C = (lep_vec.E()**2) * (met_pt**2) - Lambda**2

    delta = B**2 - 4 * A * C

    # Solve the quadratic equation for pz
    if delta >= 0:
        pz1 = (-B + math.sqrt(delta)) / (2.0 * A)
        pz2 = (-B - math.sqrt(delta)) / (2.0 * A)
        # Take the solution with the smallest absolute value
        pz_nu = pz1 if abs(pz1) < abs(pz2) else pz2
    else:
        # If complex, take the real part
        pz_nu = -B / (2.0 * A)

    nu_vec = ROOT.TLorentzVector()
    e_nu = math.sqrt(met_pt**2 + pz_nu**2)
    nu_vec.SetPxPyPzE(px_nu, py_nu, pz_nu, e_nu)
    return nu_vec  

class AsymmetryModule(Module):
    def __init__(self, channel="mu", year=2026):
        self.channel = channel
        self.year = year
        self.rp_ids = {"45": [3, 23], "56": [103, 123]}
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        # MPI Summaries 
        self.out.branch("nano_NMPI05", "I") # nMPI at PV with pT > 0.5
        self.out.branch("nano_NMPI09", "I") # nMPI at PV with pT > 0.9   
        
        # Proton Summaries
        self.out.branch("nano_nProtons", "I")
        self.out.branch("nano_pps_arm", "I", lenVar="nano_nProtons")
        self.out.branch("nano_pps_rpid", "I", lenVar="nano_nProtons")
        self.out.branch("nano_pps_x", "F", lenVar="nano_nProtons")
        self.out.branch("nano_pps_y", "F", lenVar="nano_nProtons")
        
        # Jet Summaries
        self.out.branch("nano_nJets", "I")
        self.out.branch("nano_jet_pt", "F", lenVar="nano_nJets")
        self.out.branch("nano_jet_eta", "F", lenVar="nano_nJets")
        self.out.branch("nano_jet_phi", "F", lenVar="nano_nJets")
        self.out.branch("nano_Jet_ntrk05", "I", lenVar="nano_nJets")
        self.out.branch("nano_Jet_ntrk09", "I", lenVar="nano_nJets")
        self.out.branch("nano_mJets", "F")
        self.out.branch("nano_yJets", "F")
        
        # Lepton Summaries
        self.out.branch("nano_nLeptons", "I")
        self.out.branch("nano_lep_pt", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_eta", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_phi", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_charge", "I", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_ntrk05", "I", lenVar="nano_nLeptons") 
        self.out.branch("nano_lep_ntrk09", "I", lenVar="nano_nLeptons")
        self.out.branch("nano_w_mT", "F")
        self.out.branch("nano_w_pt", "F")
        self.out.branch("nano_w_phi", "F")
        self.out.branch("nano_w_y", "F")
        self.out.branch("nano_w_m", "F")
        self.out.branch("nano_mll", "F")
        self.out.branch("nano_yll", "F")
        self.out.branch("nano_ptll", "F")
        
        # event branches
        self.out.branch("nano_Mall", "F")
        self.out.branch("nano_Yall", "F")
        
        # cutflow histogram
        self.h_cutflow = ROOT.TH1D("Cutflow", "Event Cutflow", 3, 0, 3)
        self.h_cutflow.GetXaxis().SetBinLabel(1, "All Initial Events")
        self.h_cutflow.GetXaxis().SetBinLabel(2, "Pass HLT")
        self.h_cutflow.GetXaxis().SetBinLabel(3, "Pass Module Filter")
        
        self.h_cutflow.SetBinContent(1, inputTree.GetEntries())
        
        self.events_seen = 0
        self.events_passed = 0

    def analyze(self, event):
        
        # count events that survived the "-c CUT" in nano_postproc.py
        self.events_seen += 1
        
        # Load Collections
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        jets = Collection(event, "Jet")
        protons = Collection(event, "PPSLocalTrack")
               
        # Object Selections & IDs
        # ----------------------------------------------------------------------
        # Muons: Loose (for ZCR/veto) and Tight (for WCR)
        loose_mu = [m for m in muons if m.pt > 15 and abs(m.eta) < 2.5 and m.looseId]
        tight_mu = [m for m in loose_mu if m.pt > 15 and m.tightId and m.pfRelIso04_all < 0.15]
        soft_mu = [m for m in muons if m.pt > 3 and abs(m.eta) < 2.4 and m.looseId]
        
        # Electrons: Loose (cutBased=2) and Tight (cutBased=4)
        loose_el = [e for e in electrons if e.pt > 15 and abs(e.eta) < 2.5 and e.cutBased >= 2]
        tight_el = [e for e in loose_el if e.pt > 15  and abs(e.eta) < 2.5 and e.cutBased >= 4]

        # Combine and sort by pT
        loose_leps = sorted(loose_mu + loose_el, key=lambda x: x.pt, reverse=True)
        tight_leps = sorted(tight_mu + tight_el, key=lambda x: x.pt, reverse=True)

        # Object Overlap Removal (Jets vs Leptons)
        # ----------------------------------------------------------------------
        # Remove jets that fall within dR < 0.4 of any loose lepton
        raw_jets = [j for j in jets if j.pt > 25 and abs(j.eta) < 4.7]
        sel_jets = []
        for j in raw_jets:
            has_overlap = False
            for l in loose_leps:
                if deltaR(j, l) < 0.4:
                    has_overlap = True
                    break
            if not has_overlap:
                sel_jets.append(j)

        jet_sum = ROOT.TLorentzVector()
        for j in sel_jets: jet_sum += j.p4()

        # Event Overlap Removal & Signal Region Definitions
        # ----------------------------------------------------------------------
        # By using strict elif conditions based on lepton multiplicity, 
        # events are forced into ONE region exclusively.
        
        is_ZCR = False
        is_WCR = False
        is_mj = False
        is_dimuon_inclusive = False
        leptons_to_save = []

        # DY: Exactly 2 loose leptons, opposite sign, same flavor
        if len(loose_leps) == 2 and loose_leps[0].charge != loose_leps[1].charge and loose_leps[0].pdgId == -loose_leps[1].pdgId:
            is_ZCR = True
            leptons_to_save = loose_leps
            
        if len(soft_mu) >= 2 and soft_mu[0].charge != soft_mu[1].charge:
            is_dimuon_inclusive = True
            leptons_to_save = soft_mu   
            
        # W+jets: Exactly 1 tight lepton, AND exactly 1 loose lepton (vetoes events with a 2nd loose lepton)
        elif len(tight_leps) == 1 and len(loose_leps) == 1:
            is_WCR = True
            leptons_to_save = tight_leps
            
        # MJ Control Region: >= 2 isolated jets, strictly 0 loose leptons
        elif len(sel_jets) >= 2 and len(loose_leps) == 0:
            is_mj = True
            leptons_to_save = []

        # 5. Event Filtering based on Channel
        # ----------------------------------------------------------------------
        if self.channel == "mu":
            # Only keep Muon WCR or ZCR
            if not (is_ZCR or is_WCR): return False
            if abs(leptons_to_save[0].pdgId) != 13: return False

        elif self.channel == "dimuon":
            # Only keep soft muons in di muon channel
            if not is_dimuon_inclusive: return False

        elif self.channel == "el":
            # Only keep Electron WCR or ZCR
            if not (is_ZCR or is_WCR): return False
            if abs(leptons_to_save[0].pdgId) != 11: return False
            
        elif self.channel == "mj":
            if not is_mj: return False
            
        elif self.channel == "zb":
            # Save everything
            leptons_to_save = loose_leps
            
        else:
            return False

        # Calculations: Dileptons & W variables
        # ----------------------------------------------------------------------
        mll = yll = ptll = -1.0
        w_mT = w_pt = w_phi = w_y = w_m = -1.0
        Mall = Yall = -999.0
        
        v_all = ROOT.TLorentzVector()
        
        # 1. Add all jets to the global system
        for j in sel_jets:
            v_all += j.p4()
        
        if len(leptons_to_save) >= 2:
            dilep = leptons_to_save[0].p4() + leptons_to_save[1].p4()
            mll, yll, ptll = dilep.M(), dilep.Rapidity(), dilep.Pt()
            v_all += leptons_to_save[0].p4()
            v_all += leptons_to_save[1].p4()
            
        if len(leptons_to_save) >= 1:
            met_pt = event.PuppiMET_pt
            met_phi = event.PuppiMET_phi
            dphi = ROOT.TVector2.Phi_mpi_pi(leptons_to_save[0].phi - met_phi)
            w_mT = math.sqrt(2 * leptons_to_save[0].pt * met_pt * (1 - math.cos(dphi)))
            
            w_vec = ROOT.TVector2(leptons_to_save[0].pt * math.cos(leptons_to_save[0].phi) + met_pt * math.cos(met_phi),
                                  leptons_to_save[0].pt * math.sin(leptons_to_save[0].phi) + met_pt * math.sin(met_phi))
            w_pt, w_phi = w_vec.Mod(), w_vec.Phi()
            
            # Full reconstruction using W mass constraint
            lep_p4 = ROOT.TLorentzVector()
            lep_p4.SetPtEtaPhiM(leptons_to_save[0].pt, leptons_to_save[0].eta, leptons_to_save[0].phi, leptons_to_save[0].mass)            
            nu_p4 = get_nu_p4(lep_p4, met_pt, met_phi)
            w_p4 = lep_p4 + nu_p4
            w_pt, w_phi = w_p4.Pt(), w_p4.Phi()
            w_y, w_m = w_p4.Rapidity(), w_p4.M()
            
            v_all += lep_p4
            v_all += nu_p4
            
        elif is_mj or self.channel == "zb":
            for l in leptons_to_save:
                v_all += l.p4()
                
        # Calculate Global Variables
        if v_all.E() > 0:
            Mall = v_all.M()
            if v_all.E() > abs(v_all.Pz()): # Protection for rapidity calculation
                Yall = v_all.Rapidity()

        
        # Protons Logic (LocalTrack mapping)
        pps_arm = []
        pps_rpid = []
        pps_x = []
        pps_y = []
        for p in protons:
            # Map decRPId to arm (0 for 45, 1 for 56)
            arm = 0 if p.decRPId < 100 else 1
            pps_arm.append(arm)
            pps_rpid.append(p.decRPId)
            pps_x.append(p.x)
            pps_y.append(p.y)

        
        # TRACK MULTIPLICITIES & MPI LOGIC
        
        # Read PV0 tracks
        pv0_ntrk05 = safe_get(event, "PV_ntrk0p5", 0)
        pv0_ntrk09 = safe_get(event, "PV_ntrk0p9", 0)
        
        # Extract object tracks (Default to 0 if branch missing)
        jet_trk05 = [safe_get(j, "ntrk0p5", 0) for j in sel_jets]
        jet_trk09 = [safe_get(j, "ntrk0p9", 0) for j in sel_jets]
        
        lep_trk05 = [safe_get(l, "ntrk0p5", 0) for l in leptons_to_save]
        lep_trk09 = [safe_get(l, "ntrk0p9", 0) for l in leptons_to_save]
        
        # N_trk^MPI = N_trk^PV - Sum(N_trk^objects)
        nmpi05 = pv0_ntrk05 - sum(jet_trk05) - sum(lep_trk05)
        nmpi09 = pv0_ntrk09 - sum(jet_trk09) - sum(lep_trk09)
        
        # Safety catch: prevent negative MPI in case of overlapping object cones
        if nmpi05 < 0: nmpi05 = 0
        if nmpi09 < 0: nmpi09 = 0
        
        # 7. Fill Branches
        self.out.fillBranch("nano_NMPI05", nmpi05)
        self.out.fillBranch("nano_NMPI09", nmpi09)

        self.out.fillBranch("nano_nProtons", len(pps_arm))
        self.out.fillBranch("nano_pps_arm", pps_arm)
        self.out.fillBranch("nano_pps_rpid", pps_rpid)
        self.out.fillBranch("nano_pps_x", pps_x)
        self.out.fillBranch("nano_pps_y", pps_y)

        self.out.fillBranch("nano_nJets", len(sel_jets))
        self.out.fillBranch("nano_jet_pt", [j.pt for j in sel_jets])
        self.out.fillBranch("nano_jet_eta", [j.eta for j in sel_jets])
        self.out.fillBranch("nano_jet_phi", [j.phi for j in sel_jets])
        self.out.fillBranch("nano_Jet_ntrk05", jet_trk05)
        self.out.fillBranch("nano_Jet_ntrk09", jet_trk09)
        self.out.fillBranch("nano_mJets", jet_sum.M())
        self.out.fillBranch("nano_yJets", jet_sum.Rapidity())

        self.out.fillBranch("nano_nLeptons", len(leptons_to_save))
        self.out.fillBranch("nano_lep_pt", [l.pt for l in leptons_to_save])
        self.out.fillBranch("nano_lep_eta", [l.eta for l in leptons_to_save])
        self.out.fillBranch("nano_lep_phi", [l.phi for l in leptons_to_save])
        self.out.fillBranch("nano_lep_charge", [l.charge for l in leptons_to_save])
        self.out.fillBranch("nano_lep_ntrk05", lep_trk05)
        self.out.fillBranch("nano_lep_ntrk09", lep_trk09)
        self.out.fillBranch("nano_mll", mll)
        self.out.fillBranch("nano_yll", yll)
        self.out.fillBranch("nano_ptll", ptll)
        self.out.fillBranch("nano_w_mT", w_mT)
        self.out.fillBranch("nano_w_pt", w_pt)
        self.out.fillBranch("nano_w_phi", w_phi)
        self.out.fillBranch("nano_w_y", w_y)
        self.out.fillBranch("nano_w_m", w_m)
        
        self.out.fillBranch("nano_Mall", Mall)
        self.out.fillBranch("nano_Yall", Yall)
        
        # total accepted events
        self.events_passed += 1

        return True

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        # Fill the final bins and write to the file
        self.h_cutflow.SetBinContent(2, self.events_seen)
        self.h_cutflow.SetBinContent(3, self.events_passed)
        
        outputFile.cd()
        self.h_cutflow.Write()

asymmetry_mu  = lambda : AsymmetryModule(channel="mu")
asymmetry_el  = lambda : AsymmetryModule(channel="el")
asymmetry_mj  = lambda : AsymmetryModule(channel="mj")
asymmetry_zb  = lambda : AsymmetryModule(channel="zb")

